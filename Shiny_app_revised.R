############################################################
# REVISED SHINY APP
# Direct mzML preprocessing + one-vs-all random-forest scoring
#
# Revised for manuscript repository:
# - Uses the current 278-feature model bundle
# - Uses the current filtered target table
# - Uses OVA "score" terminology only
# - No entropy
# - No probability/confidence terminology
# - No accuracy/rejection thresholds
# - Prediction = class with highest independent OVA score
# - Reports descriptive top1-top2 score margin
############################################################

############################################################
# 0. REQUIRED PACKAGES
############################################################

required_packages <- c(
  "shiny",
  "readr",
  "ranger",
  "dplyr",
  "tidyr",
  "tibble",
  "Spectra",
  "data.table"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_packages, collapse = ", "),
      ". Install them before running the app."
    )
  )
}

library(shiny)
library(readr)
library(ranger)
library(dplyr)
library(tidyr)
library(tibble)
library(Spectra)
library(data.table)

options(shiny.maxRequestSize = 20000 * 1024)


############################################################
# 1. FILE LOCATIONS
############################################################
#
# Expected directory structure:
#
#   app.R (or this file)
#   app_bundle/
#       app_model_bundle.rds
#       target_data_revised.csv
#
############################################################

MODEL_BUNDLE_PATH <- "app_model_bundle.rds"

TARGET_DATA_PATH <- "target_data_revised.csv"

if (!file.exists(MODEL_BUNDLE_PATH)) {
  stop(
    "Model bundle not found: ",
    MODEL_BUNDLE_PATH
  )
}

if (!file.exists(TARGET_DATA_PATH)) {
  stop(
    "Target data not found: ",
    TARGET_DATA_PATH
  )
}


############################################################
# 2. LOAD CURRENT MODEL BUNDLE
############################################################

app_model_bundle <- readRDS(
  MODEL_BUNDLE_PATH
)

models <- app_model_bundle$models
selected_features <- app_model_bundle$final_features
class_levels <- app_model_bundle$classes

if (is.null(models) ||
    is.null(selected_features) ||
    is.null(class_levels)) {
  stop(
    paste0(
      "app_model_bundle.rds must contain: ",
      "models, final_features, and classes."
    )
  )
}

if (length(selected_features) != 278) {
  stop(
    "Expected 278 final model features, found ",
    length(selected_features),
    "."
  )
}

if (!identical(names(models), class_levels)) {
  stop(
    "Model names and class levels are not identical."
  )
}


############################################################
# 3. LOAD CURRENT TARGET TABLE
############################################################

target_data <- fread(
  TARGET_DATA_PATH
)

required_target_columns <- c(
  "PeptideModifiedSequence",
  "PrecursorCharge",
  "PrecursorMz",
  "PrecursorId",
  "ProductMz",
  "FragmentIon",
  "RT_mean",
  "RT_sd"
)

missing_target_columns <- setdiff(
  required_target_columns,
  names(target_data)
)

if (length(missing_target_columns) > 0) {
  stop(
    paste0(
      "target_data_revised.csv is missing required columns: ",
      paste(missing_target_columns, collapse = ", ")
    )
  )
}

target_features <- unique(
  target_data$PrecursorId
)

missing_model_targets <- setdiff(
  selected_features,
  target_features
)

if (length(missing_model_targets) > 0) {
  stop(
    paste0(
      "Target table is missing ",
      length(missing_model_targets),
      " final model features."
    )
  )
}


############################################################
# 4. PREDICTION FUNCTION
############################################################
#
# The nine ranger models are independent binary OVA models.
# Their "Yes" outputs are therefore called OVA scores.
# They are NOT treated as calibrated multiclass probabilities.
#
############################################################

predict_single_sample <- function(
    sample_df,
    models,
    selected_features,
    class_levels
) {

  sample_df <- as.data.frame(
    sample_df,
    check.names = FALSE
  )

  missing_cols <- setdiff(
    selected_features,
    colnames(sample_df)
  )

  if (length(missing_cols) > 0) {
    sample_df[
      missing_cols
    ] <- 0
  }

  sample_df <- sample_df[
    ,
    selected_features,
    drop = FALSE
  ]

  sample_df[] <- lapply(
    sample_df,
    function(x) {
      x <- suppressWarnings(
        as.numeric(x)
      )
      x[
        is.na(x) |
          !is.finite(x)
      ] <- 0
      x
    }
  )

  score_vec <- vapply(
    class_levels,
    function(class_i) {

      pred <- predict(
        models[[class_i]],
        data = sample_df
      )$predictions

      if (!"Yes" %in% colnames(pred)) {
        stop(
          "Model for class ",
          class_i,
          " did not return a 'Yes' score."
        )
      }

      as.numeric(
        pred[
          1,
          "Yes"
        ]
      )
    },
    FUN.VALUE = numeric(1)
  )

  names(score_vec) <- class_levels

  score_order <- order(
    score_vec,
    decreasing = TRUE
  )

  pred_class <- names(
    score_vec
  )[
    score_order[1]
  ]

  pred_score <- score_vec[
    score_order[1]
  ]

  score_margin <- if (length(score_order) >= 2) {
    score_vec[
      score_order[1]
    ] -
      score_vec[
        score_order[2]
      ]
  } else {
    NA_real_
  }

  score_output <- as.list(
    round(
      score_vec,
      4
    )
  )

  names(score_output) <- paste0(
    "Score_",
    names(score_output)
  )

  tibble(
    Predicted_Class = pred_class,
    Predicted_Score = round(
      pred_score,
      4
    ),
    Score_Margin = round(
      score_margin,
      4
    ),
    !!!score_output
  )
}


############################################################
# 5. mzML -> ACQUISITION SUMMARY
############################################################

generate_summary_dt <- function(
    file_mzml
) {

  spectra_data <- Spectra(
    file_mzml
  )

  results_list <- lapply(
    seq_along(spectra_data),
    function(i) {

      spectrum <- spectra_data[i]

      mz_values <- as.numeric(
        unlist(
          mz(spectrum)
        )
      )

      intensity_values <- as.numeric(
        unlist(
          intensity(spectrum)
        )
      )

      data.table(
        acquisition = rep(
          i,
          length(mz_values)
        ),
        rtime = rep(
          rtime(spectrum),
          length(mz_values)
        ),
        mz = mz_values,
        intensity = intensity_values
      )
    }
  )

  rbindlist(
    results_list,
    use.names = TRUE,
    fill = TRUE
  )
}


############################################################
# 6. MATCH ACQUISITION PEAKS TO EXPECTED ISOTOPE m/z
############################################################
#
# This follows the revised development preprocessing:
# observed acquisition m/z is matched to expected ProductMz
# using +/-3 ppm.
#
############################################################

match_acquisition_to_targets <- function(
    acq_data,
    target_data
) {

  acq_data <- as.data.table(
    copy(acq_data)
  )

  targets <- as.data.table(
    copy(target_data)
  )

  targets[
    ,
    lower_isotope_mz :=
      ProductMz * (1 - 3e-6)
  ]

  targets[
    ,
    higher_isotope_mz :=
      ProductMz * (1 + 3e-6)
  ]

  matched <- targets[
    acq_data,
    on = .(
      lower_isotope_mz <= mz,
      higher_isotope_mz >= mz
    ),
    nomatch = 0,
    allow.cartesian = TRUE
  ]

  if (nrow(matched) == 0) {
    return(
      matched
    )
  }

  matched[
    ,
    c(
      "lower_isotope_mz",
      "higher_isotope_mz"
    ) := NULL
  ]

  matched
}


############################################################
# 7. REQUIRE ISOTOPE EVIDENCE
############################################################
#
# Preserve the current development-pipeline rule:
# >=3 matched isotope rows for the peptide/acquisition/
# precursor-charge/precursor-mz grouping.
#
############################################################

filter_isotope_presence <- function(
    matched_data
) {

  if (nrow(matched_data) == 0) {
    return(
      matched_data
    )
  }

  matched_data %>%
    as.data.frame() %>%
    group_by(
      PeptideModifiedSequence,
      acquisition,
      PrecursorCharge,
      PrecursorMz
    ) %>%
    mutate(
      isotope_match_count = n()
    ) %>%
    filter(
      isotope_match_count >= 3
    ) %>%
    ungroup()
}


############################################################
# 8. RETENTION-TIME FILTER
############################################################
#
# Preserve current fixed +/-5 second RT rule.
#
############################################################

filter_retention_time <- function(
    matched_data
) {

  if (nrow(matched_data) == 0) {
    return(
      matched_data
    )
  }

  dt <- as.data.table(
    copy(matched_data)
  )

  dt <- dt[
    !is.na(RT_mean) &
      !is.na(rtime)
  ]

  dt[
    ,
    RT_mean_sec :=
      RT_mean * 60
  ]

  dt[
    ,
    RT_sd_sec :=
      RT_sd * 60
  ]

  dt[
    rtime >= (
      RT_mean_sec - 5
    ) &
      rtime <= (
        RT_mean_sec + 5
      )
  ]
}


############################################################
# 9. EXTRACT MAXIMUM PRECURSOR INTENSITY
############################################################

extract_precursor_features <- function(
    filtered_data,
    selected_features
) {

  feature_values <- setNames(
    rep(
      0,
      length(selected_features)
    ),
    selected_features
  )

  if (nrow(filtered_data) == 0) {
    return(
      feature_values
    )
  }

  dt <- as.data.table(
    copy(filtered_data)
  )

  dt <- dt[
    FragmentIon == "precursor"
  ]

  if (nrow(dt) == 0) {
    return(
      feature_values
    )
  }

  summary_dt <- dt[
    ,
    .(
      MaxIntensity = {
        x <- intensity[
          is.finite(intensity)
        ]

        if (length(x) == 0) {
          0
        } else {
          max(
            x,
            na.rm = TRUE
          )
        }
      }
    ),
    by = PrecursorId
  ]

  summary_dt <- summary_dt[
    PrecursorId %in%
      selected_features
  ]

  if (nrow(summary_dt) > 0) {

    feature_values[
      summary_dt$PrecursorId
    ] <- summary_dt$MaxIntensity
  }

  feature_values
}


############################################################
# 10. PROCESS ONE mzML FILE
############################################################

process_one_mzml <- function(
    mzml_path,
    target_data,
    selected_features
) {

  acq_data <- generate_summary_dt(
    mzml_path
  )

  matched <- match_acquisition_to_targets(
    acq_data,
    target_data
  )

  isotope_filtered <- filter_isotope_presence(
    matched
  )

  rt_filtered <- filter_retention_time(
    isotope_filtered
  )

  feature_values <- extract_precursor_features(
    rt_filtered,
    selected_features
  )

  sample_df <- as.data.frame(
    as.list(
      feature_values
    ),
    check.names = FALSE
  )

  list(
    ID = basename(
      mzml_path
    ),
    sample_df = sample_df
  )
}


############################################################
# 11. PROCESS MULTIPLE mzML FILES
############################################################

process_mzml_files <- function(
    mzml_paths,
    target_data,
    selected_features
) {

  processed <- lapply(
    mzml_paths,
    function(path_i) {

      tryCatch(
        process_one_mzml(
          path_i,
          target_data,
          selected_features
        ),
        error = function(e) {
          structure(
            list(
              ID = basename(path_i),
              message = conditionMessage(e)
            ),
            class = "mzml_processing_error"
          )
        }
      )
    }
  )

  errors <- vapply(
    processed,
    inherits,
    logical(1),
    what = "mzml_processing_error"
  )

  if (any(errors)) {

    error_text <- vapply(
      processed[errors],
      function(x) {
        paste0(
          x$ID,
          ": ",
          x$message
        )
      },
      FUN.VALUE = character(1)
    )

    stop(
      paste(
        c(
          "mzML preprocessing failed:",
          error_text
        ),
        collapse = "\n"
      )
    )
  }

  feature_matrix <- bind_rows(
    lapply(
      processed,
      function(x) {
        x$sample_df
      }
    )
  )

  feature_matrix <- as.data.frame(
    feature_matrix,
    check.names = FALSE
  )

  feature_matrix <- feature_matrix[
    ,
    selected_features,
    drop = FALSE
  ]

  metadata <- tibble(
    ID = vapply(
      processed,
      function(x) {
        x$ID
      },
      FUN.VALUE = character(1)
    )
  )

  list(
    metadata = metadata,
    feature_matrix = feature_matrix
  )
}


############################################################
# 12. PREDICT A COMPLETE FEATURE MATRIX
############################################################

predict_feature_matrix <- function(
    feature_matrix,
    sample_ids = NULL
) {

  feature_matrix <- as.data.frame(
    feature_matrix,
    check.names = FALSE
  )

  if (is.null(sample_ids)) {
    sample_ids <- paste0(
      "Sample_",
      seq_len(
        nrow(feature_matrix)
      )
    )
  }

  result_list <- lapply(
    seq_len(
      nrow(feature_matrix)
    ),
    function(i) {

      pred <- predict_single_sample(
        feature_matrix[
          i,
          ,
          drop = FALSE
        ],
        models = models,
        selected_features = selected_features,
        class_levels = class_levels
      )

      bind_cols(
        tibble(
          Sample = sample_ids[i]
        ),
        pred
      )
    }
  )

  bind_rows(
    result_list
  )
}


############################################################
# 13. USER INTERFACE
############################################################

ui <- fluidPage(

  titlePanel(
    "Urine Pathogen Classification from High-Resolution MS1 Data"
  ),

  tags$p(
    style = "max-width: 1000px;",
    paste0(
      "The classifier reports one-vs-all (OVA) model scores. ",
      "These scores are outputs from independently fitted binary ",
      "classifiers and are not calibrated multiclass probabilities. ",
      "The predicted class is the class with the highest OVA score."
    )
  ),

  sidebarLayout(

    sidebarPanel(

      h4(
        "Option 1: mzML files"
      ),

      fileInput(
        "mzml_files",
        "Upload mzML file(s)",
        multiple = TRUE,
        accept = c(
          ".mzML",
          ".mzml"
        )
      ),

      actionButton(
        "preprocess_predict_btn",
        "Preprocess and Predict"
      ),

      tags$hr(),

      h4(
        "Option 2: model-ready CSV"
      ),

      fileInput(
        "sample_file",
        "Upload model-ready CSV",
        accept = ".csv"
      ),

      actionButton(
        "predict_csv_btn",
        "Predict from CSV"
      ),

      tags$hr(),

      downloadButton(
        "download_predictions",
        "Download Predictions"
      ),

      helpText(
        paste0(
          "For mzML input, target matching is performed automatically ",
          "using the bundled revised target table."
        )
      )
    ),

    mainPanel(

      textOutput(
        "status_text"
      ),

      br(),

      tableOutput(
        "prediction_table"
      )
    )
  )
)


############################################################
# 14. SERVER
############################################################

server <- function(
    input,
    output,
    session
) {

  predictions <- reactiveVal(
    NULL
  )

  status <- reactiveVal(
    ""
  )

  output$status_text <- renderText({
    status()
  })

  output$prediction_table <- renderTable({

    preds <- predictions()

    if (is.null(preds) ||
        nrow(preds) == 0) {

      return(
        data.frame(
          Message =
            "No predictions yet."
        )
      )
    }

    head(
      preds,
      20
    )
  })


  ##########################################################
  # DOWNLOAD RESULTS
  ##########################################################

  output$download_predictions <- downloadHandler(

    filename = function() {
      paste0(
        "OVA_score_predictions_",
        Sys.Date(),
        ".csv"
      )
    },

    content = function(file) {

      preds <- predictions()

      if (is.null(preds) ||
          nrow(preds) == 0) {

        writeLines(
          "No predictions available.",
          file
        )

      } else {

        write_csv(
          preds,
          file
        )
      }
    }
  )


  ##########################################################
  # mzML -> FEATURES -> PREDICTION
  ##########################################################

  observeEvent(
    input$preprocess_predict_btn,
    {

      req(
        input$mzml_files
      )

      predictions(
        NULL
      )

      status(
        "Preprocessing mzML file(s)..."
      )

      mzml_paths <- input$mzml_files$datapath

      withProgress(
        message =
          "Preprocessing and classification",
        value = 0,
        {

          incProgress(
            0.1,
            detail =
              "Reading mzML files"
          )

          processed <- tryCatch(

            process_mzml_files(
              mzml_paths = mzml_paths,
              target_data = target_data,
              selected_features = selected_features
            ),

            error = function(e) {
              showNotification(
                conditionMessage(e),
                type = "error",
                duration = NULL
              )
              NULL
            }
          )

          if (is.null(processed)) {
            status(
              "Preprocessing failed."
            )
            return()
          }

          incProgress(
            0.75,
            detail =
              "Calculating OVA scores"
          )

          pred_df <- predict_feature_matrix(
            processed$feature_matrix,
            sample_ids =
              input$mzml_files$name
          )

          predictions(
            pred_df
          )

          incProgress(
            1,
            detail =
              "Complete"
          )

          status(
            paste0(
              "Complete. ",
              nrow(pred_df),
              " sample(s) classified."
            )
          )
        }
      )
    }
  )


  ##########################################################
  # MODEL-READY CSV -> PREDICTION
  ##########################################################

  observeEvent(
    input$predict_csv_btn,
    {

      req(
        input$sample_file
      )

      predictions(
        NULL
      )

      status(
        "Reading model-ready data..."
      )

      sample_df <- read_csv(
        input$sample_file$datapath,
        show_col_types = FALSE
      )

      if (nrow(sample_df) == 0) {

        status(
          "No samples found in uploaded CSV."
        )

        return()
      }

      if ("ID" %in% colnames(sample_df)) {

        sample_ids <- as.character(
          sample_df$ID
        )

      } else if ("Sample" %in% colnames(sample_df)) {

        sample_ids <- as.character(
          sample_df$Sample
        )

      } else {

        sample_ids <- paste0(
          "Sample_",
          seq_len(
            nrow(sample_df)
          )
        )
      }

      status(
        "Calculating OVA scores..."
      )

      pred_df <- tryCatch(

        predict_feature_matrix(
          sample_df,
          sample_ids = sample_ids
        ),

        error = function(e) {

          showNotification(
            conditionMessage(e),
            type = "error",
            duration = NULL
          )

          NULL
        }
      )

      if (is.null(pred_df)) {

        status(
          "Prediction failed."
        )

        return()
      }

      predictions(
        pred_df
      )

      status(
        paste0(
          "Complete. ",
          nrow(pred_df),
          " sample(s) classified."
        )
      )
    }
  )
}


############################################################
# 15. RUN APP
############################################################

shinyApp(
  ui = ui,
  server = server
)
