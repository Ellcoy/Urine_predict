############################################################
# Reproducibility analysis for manuscript results
############################################################

# ---- Packages:
library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(caret)
library(ranger)
library(future)
library(furrr)
library(ggplot2)
library(mltools)
library(pROC)
library(PRROC)
library(nnet)
library(reshape2)
library(fastshap)
library(shapviz)

# Used only for the reported fixed-feature benchmark
library(xgboost)
library(e1071)


# ============================================================
# 1. INPUT
# ============================================================

INPUT_DIR <- "SET_ME!"

Batch2_innoc <- read.csv(
  file.path(INPUT_DIR, "targets_innoc_batch2.csv"),
  check.names = FALSE
)
Pure_samples <- read.csv(
  file.path(INPUT_DIR, "targets_cutures.csv"),
  check.names = FALSE
)
Batch1_Patient <- read.csv(
  file.path(INPUT_DIR, "targets_Clinical_batch1.csv"),
  check.names = FALSE
)
Batch2_Patient <- read.csv(
  file.path(INPUT_DIR, "targets_Clinical_batch2.csv"),
  check.names = FALSE
)

CLASS_LEVELS <- c("Blk","Eco","Efa","Kpn","Pae","Pmi","Sag","Sau","Ssa")
ID_COLS <- c("ID","Class")

keep_targets <- function(x) x[x$Class %in% CLASS_LEVELS, , drop = FALSE]

B2I <- keep_targets(Batch2_innoc)
PURE <- keep_targets(Pure_samples)
CLIN_DEV <- keep_targets(Batch2_Patient)
CLIN_HELD <- keep_targets(Batch1_Patient)

feature_names <- function(x) setdiff(names(x), ID_COLS)

# Fixed feature space is defined from DEVELOPMENT data only.
common_features <- Reduce(intersect, list(
  feature_names(B2I),
  feature_names(PURE),
  feature_names(CLIN_DEV)
))

Big_train <- bind_rows(B2I, PURE, CLIN_DEV)[
  , c("ID", common_features, "Class"), drop = FALSE
]
Big_train[common_features] <- lapply(Big_train[common_features], \(x) {
  x[is.na(x)] <- 0
  x
})
Big_train$Class <- factor(Big_train$Class, levels = CLASS_LEVELS)

# Align held-out cohort after development feature space has been fixed.
external_data <- CLIN_HELD
for(f in setdiff(common_features, names(external_data))) external_data[[f]] <- 0
external_data <- external_data[, c("ID", common_features, "Class"), drop = FALSE]
external_data[common_features] <- lapply(external_data[common_features], \(x) {
  x[is.na(x)] <- 0
  x
})
external_data$Class <- factor(external_data$Class, levels = CLASS_LEVELS)

cat("\nDevelopment samples:", nrow(Big_train))
cat("\nDevelopment feature space:", length(common_features))
cat("\nHeld-out target-positive samples:", nrow(external_data), "\n")


# ============================================================
# 2. CORE HELPERS
# ============================================================

align_features <- function(df, features) {
  for(f in setdiff(features, names(df))) df[[f]] <- 0
  out <- df[, features, drop = FALSE]
  out[] <- lapply(out, \(x) {
    x[is.na(x)] <- 0
    as.numeric(x)
  })
  out
}

safe_mcc <- function(pred, truth) {
  mltools::mcc(
    preds = factor(pred, levels = CLASS_LEVELS),
    actuals = factor(truth, levels = CLASS_LEVELS)
  )
}

per_class_mcc <- function(pred, truth) {
  sapply(CLASS_LEVELS, \(cl) {
    mltools::mcc(
      preds = factor(pred == cl, levels = c(FALSE, TRUE)),
      actuals = factor(truth == cl, levels = c(FALSE, TRUE))
    )
  })
}

select_features_rf <- function(train_df, top_n) {
  fit <- ranger(
    dependent.variable.name = "Class",
    data = train_df[, c(common_features, "Class"), drop = FALSE],
    probability = TRUE,
    importance = "permutation",
    num.trees = 500,
    seed = 123
  )
  imp <- sort(fit$variable.importance, decreasing = TRUE)
  names(imp)[seq_len(min(top_n, length(imp)))]
}

run_ova_rf <- function(train_df, test_df, features, mtry_frac, min_node, trees) {
  Xtr <- align_features(train_df, features)
  Xte <- align_features(test_df, features)

  scores <- matrix(NA_real_, nrow(Xte), length(CLASS_LEVELS),
                   dimnames = list(NULL, CLASS_LEVELS))
  models <- vector("list", length(CLASS_LEVELS))
  names(models) <- CLASS_LEVELS

  for(cl in CLASS_LEVELS) {
    y <- factor(ifelse(train_df$Class == cl, "Yes", "No"),
                levels = c("No","Yes"))

    models[[cl]] <- ranger(
      dependent.variable.name = "y",
      data = data.frame(Xtr, y = y, check.names = FALSE),
      probability = TRUE,
      num.trees = as.integer(trees),
      min.node.size = as.integer(min_node),
      mtry = max(1, min(ncol(Xtr), round(length(features) * mtry_frac))),
      importance = "permutation",
      seed = 123
    )

    scores[, cl] <- predict(models[[cl]], data = Xte)$predictions[, "Yes"]
  }

  pred <- factor(
    colnames(scores)[max.col(scores, ties.method = "first")],
    levels = CLASS_LEVELS
  )

  list(pred_class = pred, prob_matrix = scores, models = models)
}

run_multinom <- function(train_df, test_df, features) {
  Xtr <- align_features(train_df, features)
  Xte <- align_features(test_df, features)

  fit <- nnet::multinom(
    Class ~ .,
    data = data.frame(Xtr, Class = train_df$Class, check.names = FALSE),
    trace = FALSE,
    MaxNWts = 100000,
    maxit = 500
  )

  factor(predict(fit, newdata = Xte), levels = CLASS_LEVELS)
}


# ============================================================
# 3. GROUPED CV STRUCTURE
# ============================================================

extract_urine_id <- function(x) {
  stringr::str_match(x, "_(U[0-9]+)(?:_|\\.)")[,2]
}

is_b2_inoc <- grepl(
  "20250209_EXB_FS_inj2-300SPD",
  Big_train$ID,
  fixed = TRUE
)

urine_id <- rep(NA_character_, nrow(Big_train))
urine_id[is_b2_inoc] <- extract_urine_id(Big_train$ID[is_b2_inoc])
if(any(is.na(urine_id[is_b2_inoc]))) stop("Could not extract a Batch 2 donor ID.")

Big_train$CV_Group <- ifelse(
  is_b2_inoc,
  paste0("BATCH2_URINE_", urine_id),
  paste0("UNIQUE_ROW_", seq_len(nrow(Big_train)))
)

make_grouped_folds <- function(data, k, seed, proportions) {
  proportions <- proportions / sum(proportions)
  groups <- unique(as.character(data$CV_Group))

  gc <- as.matrix(xtabs(
    ~ factor(CV_Group, levels = groups) +
      factor(Class, levels = CLASS_LEVELS),
    data = data
  ))
  rownames(gc) <- groups
  sizes <- rowSums(gc)
  target_class <- outer(proportions, colSums(gc))
  target_n <- proportions * sum(sizes)

  set.seed(seed)
  best <- NULL
  best_score <- Inf

  for(i in 1:1000) {
    a <- sample(seq_len(k), length(groups), replace = TRUE, prob = proportions)
    fc <- matrix(0, k, length(CLASS_LEVELS))
    fn <- numeric(k)
    fg <- numeric(k)

    for(f in seq_len(k)) {
      idx <- which(a == f)
      if(length(idx)) {
        fc[f,] <- colSums(gc[idx,,drop = FALSE])
        fn[f] <- sum(sizes[idx])
        fg[f] <- length(idx)
      }
    }

    score <-
      sum(((fc - target_class)^2) / pmax(target_class, 1)) +
      0.25 * sum(((fn - target_n)^2) / pmax(target_n, 1)) +
      1e6 * sum(fc == 0) +
      1e6 * sum(fg == 0)

    if(score < best_score) {
      best <- a
      best_score <- score
    }
  }

  map <- setNames(best, groups)
  fold_id <- unname(map[as.character(data$CV_Group)])

  for(f in seq_len(k)) {
    if(length(intersect(
      unique(data$CV_Group[fold_id == f]),
      unique(data$CV_Group[fold_id != f])
    ))) stop("Grouped-CV leakage detected.")
  }

  fold_id
}


# ============================================================
# 4. 20 x GROUPED NESTED CV
# ============================================================

TOP_N <- seq(50, 500, 50)
GRID <- expand.grid(
  top_n = TOP_N,
  mtry_frac = c(0.05, 0.10, 0.20),
  min_node = c(1, 5, 10),
  num_trees = c(300, 500)
)

one_outer_run <- function(data, seed) {

  outer_fold <- make_grouped_folds(
    data, 2, seed, c(0.70, 0.30)
  )
  outer_train <- data[outer_fold == 1,,drop = FALSE]
  outer_test  <- data[outer_fold == 2,,drop = FALSE]

  inner_fold <- make_grouped_folds(
    outer_train, 5, seed + 10000L, rep(0.2, 5)
  )

  inner <- lapply(1:5, \(f) {
    tr <- outer_train[inner_fold != f,,drop = FALSE]
    va <- outer_train[inner_fold == f,,drop = FALSE]

    grid_res <- bind_rows(lapply(seq_len(nrow(GRID)), \(g) {
      p <- GRID[g,]
      features <- select_features_rf(tr, p$top_n)
      fit <- run_ova_rf(
        tr, va, features,
        p$mtry_frac, p$min_node, p$num_trees
      )
      tibble(
        top_n = p$top_n,
        mtry_frac = p$mtry_frac,
        min_node = p$min_node,
        num_trees = p$num_trees,
        mcc = safe_mcc(fit$pred_class, va$Class)
      )
    }))

    best <- grid_res[which.max(grid_res$mcc),,drop = FALSE]
    list(
      features = select_features_rf(tr, best$top_n),
      best = best
    )
  })

  freq <- table(unlist(map(inner, "features")))
  consensus <- names(freq[freq >= 3])
  if(!length(consensus)) {
    consensus <- names(sort(freq, decreasing = TRUE))[seq_len(min(50,length(freq)))]
  }

  param_summary <- bind_rows(map(inner, "best")) %>%
    group_by(top_n, mtry_frac, min_node, num_trees) %>%
    summarise(mean_mcc = mean(mcc), freq = n(), .groups = "drop") %>%
    arrange(desc(freq), desc(mean_mcc))

  best <- param_summary[1,,drop = FALSE]

  ova <- run_ova_rf(
    outer_train, outer_test, consensus,
    best$mtry_frac, best$min_node, best$num_trees
  )

  multinom_pred <- run_multinom(
    outer_train, outer_test, consensus
  )

  list(
    outer_mcc = safe_mcc(ova$pred_class, outer_test$Class),
    baseline_mcc = safe_mcc(multinom_pred, outer_test$Class),
    class_mcc = per_class_mcc(ova$pred_class, outer_test$Class),
    features = consensus,
    best_params = best,
    pred = ova$pred_class,
    truth = outer_test$Class,
    probs = ova$prob_matrix
  )
}

plan(multisession, workers = max(1, parallel::detectCores() - 1))
set.seed(999)

results <- future_map(
  1:20,
  ~ one_outer_run(Big_train, .x),
  .progress = TRUE,
  .options = furrr_options(seed = TRUE)
)

plan(sequential)


# ============================================================
# 5. PRIMARY CV OUTPUTS
# ============================================================

ova_mcc <- map_dbl(results, "outer_mcc")
baseline_mcc <- map_dbl(results, "baseline_mcc")

primary_cv <- tibble(
  Model = c("OVA Random Forest","Multinomial baseline"),
  Mean_MCC = c(mean(ova_mcc), mean(baseline_mcc)),
  SD_MCC = c(sd(ova_mcc), sd(baseline_mcc)),
  Minimum = c(min(ova_mcc), min(baseline_mcc)),
  Maximum = c(max(ova_mcc), max(baseline_mcc))
)

class_mcc <- bind_rows(lapply(seq_along(results), \(i) {
  tibble(
    Iteration = i,
    Class = names(results[[i]]$class_mcc),
    MCC = as.numeric(results[[i]]$class_mcc)
  )
})) %>%
  group_by(Class) %>%
  summarise(Mean_MCC = mean(MCC), SD_MCC = sd(MCC), .groups = "drop")

cat("\n================ PRIMARY GROUPED CV ================\n")
print(primary_cv, n = Inf)
cat("\nClass-specific MCC:\n")
print(class_mcc, n = Inf)


# ============================================================
# 6. CV CLASS METRICS + 95% CIs
# ============================================================

class_metrics <- function(truth, pred, score, cl) {
  pos <- truth == cl
  pp <- pred == cl
  TP <- sum(pos & pp); FN <- sum(pos & !pp)
  TN <- sum(!pos & !pp); FP <- sum(!pos & pp)

  sens <- TP/(TP+FN)
  spec <- TN/(TN+FP)
  ba <- mean(c(sens,spec))

  auc <- as.numeric(pROC::auc(pROC::roc(
    as.numeric(pos), as.numeric(score),
    quiet = TRUE, direction = "<"
  )))

  tibble(Sensitivity = sens, Specificity = spec,
         Balanced_Accuracy = ba, ROC_AUC = auc)
}

outer_class_metrics <- bind_rows(lapply(seq_along(results), \(i) {
  r <- results[[i]]
  bind_rows(lapply(CLASS_LEVELS, \(cl) {
    class_metrics(r$truth, r$pred, r$probs[,cl], cl) %>%
      mutate(Iteration = i, Class = cl, .before = 1)
  }))
}))

mean_ci <- function(x) {
  n <- sum(is.finite(x))
  x <- x[is.finite(x)]
  m <- mean(x); s <- sd(x)
  e <- qt(.975, n - 1) * s/sqrt(n)
  c(mean = m, sd = s, lower = m-e, upper = m+e)
}

cv_class_table <- outer_class_metrics %>%
  group_by(Class) %>%
  summarise(
    N_Iterations = n(),
    Sensitivity_Mean = mean_ci(Sensitivity)["mean"],
    Sensitivity_SD = mean_ci(Sensitivity)["sd"],
    Sensitivity_CI_Lower = mean_ci(Sensitivity)["lower"],
    Sensitivity_CI_Upper = mean_ci(Sensitivity)["upper"],
    Specificity_Mean = mean_ci(Specificity)["mean"],
    Specificity_SD = mean_ci(Specificity)["sd"],
    Specificity_CI_Lower = mean_ci(Specificity)["lower"],
    Specificity_CI_Upper = mean_ci(Specificity)["upper"],
    Balanced_Accuracy_Mean = mean_ci(Balanced_Accuracy)["mean"],
    Balanced_Accuracy_SD = mean_ci(Balanced_Accuracy)["sd"],
    Balanced_Accuracy_CI_Lower = mean_ci(Balanced_Accuracy)["lower"],
    Balanced_Accuracy_CI_Upper = mean_ci(Balanced_Accuracy)["upper"],
    ROC_AUC_Mean = mean_ci(ROC_AUC)["mean"],
    ROC_AUC_SD = mean_ci(ROC_AUC)["sd"],
    ROC_AUC_CI_Lower = mean_ci(ROC_AUC)["lower"],
    ROC_AUC_CI_Upper = mean_ci(ROC_AUC)["upper"],
    .groups = "drop"
  )

cat("\n================ CV CLASS METRICS ================\n")
print(cv_class_table %>% mutate(across(where(is.numeric), ~round(.x,3))),
      n = Inf, width = Inf)


# ============================================================
# 7. FINAL FEATURES, PARAMETERS, HELD-OUT MODEL
# ============================================================

feature_freq <- table(unlist(map(results, "features")))
final_features <- names(feature_freq[feature_freq >= 10])

global_params <- bind_rows(map(results, "best_params")) %>%
  group_by(top_n, mtry_frac, min_node, num_trees) %>%
  summarise(mean_mcc = mean(mean_mcc), freq = n(), .groups = "drop") %>%
  arrange(desc(freq), desc(mean_mcc))

global_best <- global_params[1,,drop = FALSE]

final_model <- run_ova_rf(
  Big_train, external_data, final_features,
  global_best$mtry_frac,
  global_best$min_node,
  global_best$num_trees
)

heldout_truth <- external_data$Class
heldout_pred <- final_model$pred_class
heldout_mcc <- safe_mcc(heldout_pred, heldout_truth)
heldout_correct <- sum(heldout_pred == heldout_truth)

cat("\n================ FINAL MODEL ================\n")
cat("Final features:", length(final_features), "\n")
print(global_best)
cat("Held-out correct:", heldout_correct, "/", nrow(external_data),
    "=", round(100*heldout_correct/nrow(external_data),1), "%\n")
cat("Held-out MCC:", round(heldout_mcc,4), "\n")


# ============================================================
# 8. HELD-OUT CLASS METRICS + 95% CIs
# ============================================================

binom_ci <- function(x,n) {
  if(n == 0) return(c(est=NA,lo=NA,hi=NA))
  z <- binom.test(x,n)$conf.int
  c(est=x/n, lo=z[1], hi=z[2])
}

bootstrap_ba <- function(truth,pred,cl,B=10000) {
  pos <- which(truth == cl)
  neg <- which(truth != cl)
  if(!length(pos) || !length(neg)) return(c(NA,NA))

  vals <- replicate(B, {
    idx <- c(sample(pos,length(pos),TRUE), sample(neg,length(neg),TRUE))
    t <- truth[idx]; p <- pred[idx]
    sens <- sum(t==cl & p==cl)/sum(t==cl)
    spec <- sum(t!=cl & p!=cl)/sum(t!=cl)
    mean(c(sens,spec))
  })
  unname(quantile(vals,c(.025,.975),na.rm=TRUE))
}

set.seed(12345)

heldout_class_table <- bind_rows(lapply(CLASS_LEVELS, \(cl) {
  pos <- heldout_truth == cl
  pred_pos <- heldout_pred == cl
  npos <- sum(pos)

  if(npos == 0) {
    return(tibble(
      Class=cl, N_Positive=0,
      Sensitivity=NA, Sensitivity_CI_Lower=NA, Sensitivity_CI_Upper=NA,
      Specificity=NA, Specificity_CI_Lower=NA, Specificity_CI_Upper=NA,
      Balanced_Accuracy=NA, Balanced_Accuracy_CI_Lower=NA,
      Balanced_Accuracy_CI_Upper=NA,
      ROC_AUC=NA, ROC_AUC_CI_Lower=NA, ROC_AUC_CI_Upper=NA
    ))
  }

  TP <- sum(pos & pred_pos); FN <- sum(pos & !pred_pos)
  TN <- sum(!pos & !pred_pos); FP <- sum(!pos & pred_pos)

  s1 <- binom_ci(TP,TP+FN)
  s2 <- binom_ci(TN,TN+FP)
  ba_ci <- bootstrap_ba(heldout_truth,heldout_pred,cl)

  roc_obj <- pROC::roc(
    as.numeric(pos),
    final_model$prob_matrix[,cl],
    quiet=TRUE, direction="<"
  )
  auc <- as.numeric(pROC::auc(roc_obj))
  auc_ci <- suppressWarnings(as.numeric(
    pROC::ci.auc(roc_obj, conf.level=.95, method="delong")
  ))

  tibble(
    Class=cl, N_Positive=npos,
    Sensitivity=s1["est"], Sensitivity_CI_Lower=s1["lo"],
    Sensitivity_CI_Upper=s1["hi"],
    Specificity=s2["est"], Specificity_CI_Lower=s2["lo"],
    Specificity_CI_Upper=s2["hi"],
    Balanced_Accuracy=mean(c(s1["est"],s2["est"])),
    Balanced_Accuracy_CI_Lower=ba_ci[1],
    Balanced_Accuracy_CI_Upper=ba_ci[2],
    ROC_AUC=auc, ROC_AUC_CI_Lower=auc_ci[1], ROC_AUC_CI_Upper=auc_ci[3]
  )
}))

cat("\n================ HELD-OUT CLASS METRICS ================\n")
print(heldout_class_table %>% mutate(across(where(is.numeric), ~round(.x,3))),
      n = Inf, width = Inf)


# ============================================================
# 9. Model score values between Correct, incorrect and exploratory predictions
# ============================================================

score_stats <- function(scores) {
  max_score <- apply(scores,1,max)
  margin <- apply(scores,1,\(x) {
    z <- sort(x,decreasing=TRUE)
    z[1]-z[2]
  })
  tibble(Max_Score=max_score, Margin=margin)
}

cv_scores <- bind_rows(lapply(results, \(r) {
  score_stats(r$probs) %>%
    mutate(Outcome=ifelse(r$pred==r$truth,"Correct","Incorrect"))
}))

held_scores <- score_stats(final_model$prob_matrix) %>%
  mutate(Outcome=ifelse(heldout_pred==heldout_truth,"Correct","Incorrect"))

exploratory <- bind_rows(
  Batch2_innoc[!Batch2_innoc$Class %in% CLASS_LEVELS,,drop=FALSE],
  Batch2_Patient[!Batch2_Patient$Class %in% CLASS_LEVELS,,drop=FALSE],
  Batch1_Patient[!Batch1_Patient$Class %in% CLASS_LEVELS,,drop=FALSE]
)

exploratory_fit <- run_ova_rf(
  Big_train, exploratory, final_features,
  global_best$mtry_frac, global_best$min_node, global_best$num_trees
)

exploratory_scores <- score_stats(exploratory_fit$prob_matrix) %>%
  mutate(Outcome="Exploratory")

score_summary <- bind_rows(
  cv_scores %>% mutate(Dataset="Grouped outer-CV"),
  held_scores %>% mutate(Dataset="Held-out clinical"),
  exploratory_scores %>% mutate(Dataset="Non-target / polymicrobial")
) %>%
  group_by(Dataset,Outcome) %>%
  summarise(
    N=n(),
    Median_Max_Score=median(Max_Score),
    Median_Margin=median(Margin),
    .groups="drop"
  )

cat("\n================ SCORE / MARGIN VALUES ================\n")
print(score_summary,n=Inf)


# ============================================================
# 10. FIXED-278-FEATURE CLASSIFIER BENCHMARK
# ============================================================

X_all <- align_features(Big_train, final_features)
y_all <- factor(Big_train$Class, levels=CLASS_LEVELS)

metric3 <- function(truth,pred) {
  cm <- caret::confusionMatrix(
    factor(pred,levels=CLASS_LEVELS),
    factor(truth,levels=CLASS_LEVELS)
  )
  tibble(
    MCC=safe_mcc(pred,truth),
    Accuracy=as.numeric(cm$overall["Accuracy"]),
    Balanced_Accuracy=mean(cm$byClass[,"Balanced Accuracy"],na.rm=TRUE)
  )
}

fit_bench <- function(name,X,y) {
  if(name=="OVA_RF") {
    out <- list()
    for(cl in CLASS_LEVELS) {
      yy <- factor(ifelse(y==cl,"Yes","No"),levels=c("No","Yes"))
      out[[cl]] <- ranger(
        dependent.variable.name="y",
        data=data.frame(X,y=yy,check.names=FALSE),
        probability=TRUE,num.trees=500,importance="none",seed=123
      )
    }
    return(out)
  }

  if(name=="Multiclass_RF") return(ranger(
    Class~.,data=data.frame(X,Class=y,check.names=FALSE),
    probability=TRUE,num.trees=500,importance="none",seed=123
  ))

  if(name=="XGBoost") return(xgboost::xgb.train(
    params=list(
      objective="multi:softprob",num_class=length(CLASS_LEVELS),
      eval_metric="mlogloss",max_depth=6,eta=.05,
      subsample=.8,colsample_bytree=.8
    ),
    data=xgboost::xgb.DMatrix(as.matrix(X),label=as.integer(y)-1),
    nrounds=200,verbose=0
  ))

  if(name=="SVM") return(e1071::svm(
    x=X,y=y,kernel="radial",scale=TRUE,probability=TRUE
  ))

  if(name=="Naive_Bayes") return(e1071::naiveBayes(x=X,y=y))
}

pred_bench <- function(name,model,X) {
  if(name=="OVA_RF") {
    s <- sapply(CLASS_LEVELS,\(cl)
      predict(model[[cl]],data=X)$predictions[,"Yes"])
    return(factor(colnames(s)[max.col(s,ties.method="first")],
                  levels=CLASS_LEVELS))
  }

  if(name=="Multiclass_RF") {
    p <- predict(model,data=X)$predictions
    return(factor(colnames(p)[max.col(p,ties.method="first")],
                  levels=CLASS_LEVELS))
  }

  if(name=="XGBoost") {
    p <- predict(model,xgboost::xgb.DMatrix(as.matrix(X)))
    p <- matrix(p,ncol=length(CLASS_LEVELS),byrow=TRUE)
    return(factor(CLASS_LEVELS[max.col(p,ties.method="first")],
                  levels=CLASS_LEVELS))
  }

  if(name=="SVM")
    return(factor(predict(model,X,probability=TRUE),levels=CLASS_LEVELS))

  if(name=="Naive_Bayes")
    return(factor(predict(model,X,type="class"),levels=CLASS_LEVELS))
}

MODELS <- c("OVA_RF","Multiclass_RF","XGBoost","SVM","Naive_Bayes")

set.seed(2026)
folds <- caret::createFolds(y_all,k=10,returnTrain=FALSE)

benchmark_folds <- bind_rows(lapply(MODELS,\(m) {
  bind_rows(lapply(1:10,\(f) {
    va <- folds[[f]]
    tr <- setdiff(seq_len(nrow(X_all)),va)
    fit <- fit_bench(m,X_all[tr,,drop=FALSE],y_all[tr])
    pred <- pred_bench(m,fit,X_all[va,,drop=FALSE])
    metric3(y_all[va],pred) %>% mutate(Model=m,Fold=f,.before=1)
  }))
}))

benchmark_summary <- benchmark_folds %>%
  group_by(Model) %>%
  summarise(
    Mean_MCC=mean(MCC), SD_MCC=sd(MCC),
    Mean_Accuracy=mean(Accuracy), SD_Accuracy=sd(Accuracy),
    Mean_Balanced_Accuracy=mean(Balanced_Accuracy),
    SD_Balanced_Accuracy=sd(Balanced_Accuracy),
    .groups="drop"
  ) %>%
  arrange(desc(Mean_MCC))

cat("\n================ CLASSIFIER BENCHMARK ================\n")
print(benchmark_summary,n=Inf,width=Inf)


# Optional compact exports
dir.create("reproducibility_results",showWarnings=FALSE)
write.csv(primary_cv,"reproducibility_results/primary_grouped_cv.csv",row.names=FALSE)
write.csv(class_mcc,"reproducibility_results/class_specific_mcc.csv",row.names=FALSE)
write.csv(cv_class_table,"reproducibility_results/cv_class_metrics_with_CI.csv",row.names=FALSE)
write.csv(heldout_class_table,"reproducibility_results/heldout_class_metrics_with_CI.csv",row.names=FALSE)
write.csv(score_summary,"reproducibility_results/score_margin_summary.csv",row.names=FALSE)
write.csv(benchmark_summary,"reproducibility_results/classifier_benchmark.csv",row.names=FALSE)
write.csv(checks,"reproducibility_results/article_consistency_checks.csv",row.names=FALSE)
capture.output(sessionInfo(),file="reproducibility_results/sessionInfo.txt")
