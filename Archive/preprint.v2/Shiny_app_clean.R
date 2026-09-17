# ---- Ensure Required Packages ----
required_packages <- c(
  "shiny", "shinyjs", "readr", "ranger", "furrr", "future", 
  "dplyr", "tools", "stringr", "purrr", "here","Spectra","data.table","tidyverse"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, library, character.only = TRUE))


# ---- Load Model Object ----
load("models.RData")

models <- deployment_model$models
selected_features <- deployment_model$selected_features
class_levels <- deployment_model$class_levels
accuracy_thresholds <- deployment_model$accuracy_thresholds

options(shiny.maxRequestSize = 20000 * 1024)
# accuracy_thresholds <- c(Blk=0.25, Eco=0.25, Efa=0.25, Kpn=0.25, Pae=0.25, Sag=0.25, Sau=0.25, Ssa=0.25, Pmi=0.25)
# class_levels <- names(accuracy_thresholds)

get_expected_features <- function(model) {
  feats <- names(model$forest$xlevels)
  if (length(feats) == 0 && !is.null(model$data)) feats <- colnames(model$data)[-ncol(model$data)]
  if (length(feats) == 0) stop("Unable to retrieve feature names from model.")
  feats
}

predict_single_sample <- function(
    sample_df,
    models,
    selected_features,
    class_levels,
    accuracy_thresholds
){
  
  missing_cols <- setdiff(
    selected_features,
    colnames(sample_df)
  )
  
  if(length(missing_cols) > 0)
    sample_df[missing_cols] <- 0
  
  sample_df <- sample_df[
    ,
    selected_features,
    drop = FALSE
  ]
  
  sample_df[is.na(sample_df)] <- 0
  
  prob_vec <- sapply(
    class_levels,
    function(class_i){
      
      predict(
        models[[class_i]],
        data = sample_df
      )$predictions[, "Yes"]
    }
  )
  
  prob_vec <- setNames(
    as.numeric(prob_vec),
    class_levels
  )
  
  pred_class <- names(prob_vec)[which.max(prob_vec)]
  pred_prob <- max(prob_vec)
  
  threshold <- accuracy_thresholds[pred_class]
  
  is_accurate <- !is.na(threshold) &&
    pred_prob >= threshold
  
  notification <- if (is.na(threshold)) {
    
    paste0(
      "No threshold for class '",
      pred_class,
      "'."
    )
    
  } else if (is_accurate) {
    
    paste0(
      "Prediction is ABOVE threshold (",
      threshold,
      ")."
    )
    
  } else {
    
    paste0(
      "Prediction is BELOW threshold (",
      threshold,
      ")."
    )
    
  }
  
  tibble(
    Predicted_Class = pred_class,
    Predicted_Probability = round(pred_prob, 4),
    Above_Threshold = ifelse(
      is_accurate,
      "Yes",
      "No"
    ),
    Message = notification,
    Entropy = round(
      -sum(prob_vec * log(prob_vec + 1e-12)),
      4
    ),
    Margin = round(
      sort(prob_vec, decreasing = TRUE)[1] -
        sort(prob_vec, decreasing = TRUE)[2],
      4
    ),
    !!!as.list(round(prob_vec, 4))
  )
}

# ---- Preprocessing functions -----

generate_summary_dt <- function(file_MZML, output_dir = ".") {
  log_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file_MZML)), "_log.txt"))
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file_MZML)), "_summary.csv"))
  
  log_conn <- file(log_file, open = "wt")
  sink(log_conn, type = "output")
  sink(log_conn, type = "message", append = TRUE)
  
  tryCatch({
    spectra_data <- Spectra(file_MZML)
    results_list <- lapply(seq_along(spectra_data), function(i) {
      spectrum <- spectra_data[i]
      mz_values <- as.numeric(unlist(mz(spectrum)))
      intensity_values <- as.numeric(unlist(intensity(spectrum)))
      data.table(
        acquisition = rep(i, length(mz_values)),
        rtime = rep(rtime(spectrum), length(mz_values)),
        mz = mz_values,
        intensity = intensity_values,
        peakcount = length(mz_values)
      )
    })
    summary_dt <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
    
    if (all(sapply(results_list, function(dt) nrow(dt)) == sapply(results_list, function(dt) unique(dt$peakcount)))) {
      message("Data integrity check passed: All peak counts match.")
    } else {
      warning("Data integrity check failed: Mismatch in peak counts.")
    }
    if (length(unique(summary_dt$acquisition)) == length(spectra_data)) {
      message("All spectra present")
    } else {
      warning("Spectra are missing from the data")
    }
    
    fwrite(summary_dt, output_file)
  }, error = function(e) {
    cat("An error occurred:", e$message, "\n")
  }, finally = {
    sink(type = "message")
    sink()
    close(log_conn)
  })
  
  return(output_file)
}

raw_acquisition_to_targets <- function(acq_file_path, target_data_path, output_dir = ".") {
  acq_data <- fread(acq_file_path)
  target_data <- fread(target_data_path)
  names(target_data)[names(target_data) == '-3ppm'] <- 'lowerlimit'
  names(target_data)[names(target_data) == '+3ppm'] <- 'higherlimit'
  
  targets_in_acq <- full_join(acq_data, target_data, by = join_by(between(mz, lowerlimit, higherlimit)))
  filtered_targets_in_acq <- targets_in_acq[!rowSums(is.na(targets_in_acq)) == 19,]
  
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(acq_file_path)), "_targets_3ppm_isotope_8000.csv"))
  fwrite(filtered_targets_in_acq, output_file)
  
  return(output_file)
}

get_isotope_presence <- function(file_path, output_dir = ".") {
  filtered_targets_in_acq <- fread(file_path, data.table = FALSE)
  test <- filtered_targets_in_acq %>%
    group_by(PeptideModifiedSequence, acquisition, PrecursorCharge, PrecursorMz) %>%
    mutate(counter = n()) %>%
    filter(counter == 3)
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file_path)), "_3ppm_isotopes_8000_present.csv"))
  fwrite(test, output_file)
  return(output_file)
}

get_rtimes <- function(file_path, reduced_target_info, output_dir = ".") {
  data_test <- fread(file_path, data.table = FALSE)
  data_test <- data_test[!rowSums(is.na(data_test)) == 5,]
  testy <- full_join(data_test, reduced_target_info, by = c("PrecursorMz" = "mz", "PeptideModifiedSequence" = "Modified.Sequence"))
  testy$RT_mean_sec <- testy$RT_mean * 60
  testy$RT_sd_sec <- testy$RT_sd * 60
  testy_full_30 <- testy %>%
    filter(rtime >= (RT_mean_sec - 5) & rtime <= (RT_mean_sec + 5))
  output_file <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file_path)), "_3ppm_isotopes_8000_rtimes.csv"))
  fwrite(testy_full_30, output_file)
  return(output_file)
}


# Get model ready csv
process_max_intensity_summary <- function(input_files, output_dir, class_regex = ".*180K_([^_]+)_.*") {
  results <- list()
  
  for (file in input_files) {
    data <- fread(file)
    
    # Filter for only precursor ions
    data_precursor <- data[FragmentIon == "precursor"]
    
    # Get max intensity per Precursor.Id
    summary <- data_precursor[, .(MaxIntensity = max(intensity, na.rm = TRUE)), by = .(Precursor.Id)]
    
    # Extract file and class info
    filename <- basename(file)
    class_name <- ifelse(grepl(class_regex, filename),
                         sub(class_regex, "\\1", filename),
                         NA_character_)
    
    summary[, File := filename]
    summary[, Class := class_name]
    
    results[[filename]] <- summary
  }
  
  # Combine all results
  results_df <- rbindlist(results, fill = TRUE)
  setcolorder(results_df, c("File", "Class", "Precursor.Id", "MaxIntensity"))
  
  # Pivot to wide format: rows = Precursor.Id, columns = File, values = MaxIntensity
  precursor_matrix <- results_df %>%
    select(Precursor.Id, File, MaxIntensity) %>%
    pivot_wider(names_from = File, values_from = MaxIntensity, values_fill = list(MaxIntensity = 0))
  
  # Set Precursor.Id as rownames
  heatmap_data <- precursor_matrix %>%
    column_to_rownames("Precursor.Id")
  
  # Create file-to-class mapping for joining later
  file_class_mapping <- results_df %>%
    distinct(File, Class) %>%
    column_to_rownames("File")
  
  # Transpose the matrix: now rows = files, columns = Precursor.Ids
  heatmap_data_transposed <- as.data.frame(t(heatmap_data))
  
  # Add rownames as a column for filenames
  heatmap_data_transposed <- tibble::rownames_to_column(heatmap_data_transposed, var = "File")
  
  # Join with class info and rename columns explicitly
  heatmap_data_transposed <- heatmap_data_transposed %>%
    left_join(file_class_mapping %>% rownames_to_column("File"), by = "File")
  
  # Check if 'Class' column exists before rename
  if ("Class" %in% colnames(heatmap_data_transposed)) {
    heatmap_data_transposed <- dplyr::rename(heatmap_data_transposed,
                                             class = Class,
                                             ID = File)
  } else {
    # If Class missing, just rename File to ID
    heatmap_data_transposed <- dplyr::rename(heatmap_data_transposed, ID = File)
  }
  
  # Move class to last column
  heatmap_data_transposed <- heatmap_data_transposed %>%
    relocate(class, .after = last_col())
  
  # Write CSV in output_dir
  output_csv <- file.path(output_dir, "model_ready_data.csv")
  write.csv(heatmap_data_transposed, output_csv, row.names = FALSE)
  
  return(output_csv)
}


library(future.apply)

run_pipeline <- function(mzML_files, target_data_path, reduced_target_info_path, 
                         output_dir = "./output", cores = 4, class_regex = ".*180K_([^_]+)_.*") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Setup parallel plan cross-platform
  plan(multisession, workers = cores)  # works on Windows, Linux, macOS
  
  # Helper to safely run a function with error capture
  safe_run <- function(f, ...) {
    tryCatch(f(...), error = function(e) e)
  }
  
  # Step 1: Generate summary CSVs
  message("Step 1: Generating summary CSVs from mzML files...")
  summary_files <- future_lapply(mzML_files, function(f) safe_run(generate_summary_dt, f, output_dir = output_dir))
  
  error_indices_1 <- which(sapply(summary_files, inherits, "error"))
  if (length(error_indices_1) > 0) {
    message("Errors occurred in generate_summary_dt for the following files:")
    for (i in error_indices_1) {
      cat("File:", mzML_files[[i]], "\nError message:\n")
      print(summary_files[[i]])
      cat("\n")
    }
  } else {
    message("All summary CSVs generated successfully.")
  }
  
  summary_files_success <- summary_files[!sapply(summary_files, inherits, "error")]
  
  # Step 2: Map acquisitions to target data
  message("Step 2: Mapping acquisitions to target data...")
  target_csv_files <- future_lapply(summary_files_success, function(f) safe_run(raw_acquisition_to_targets, f, target_data_path = target_data_path, output_dir = output_dir))
  
  error_indices_2 <- which(sapply(target_csv_files, inherits, "error"))
  if (length(error_indices_2) > 0) {
    message("Errors occurred in raw_acquisition_to_targets for the following files:")
    for (i in error_indices_2) {
      cat("File:", summary_files_success[[i]], "\nError message:\n")
      print(target_csv_files[[i]])
      cat("\n")
    }
  } else {
    message("All acquisitions mapped to target data successfully.")
  }
  
  target_csv_files_success <- target_csv_files[!sapply(target_csv_files, inherits, "error")]
  
  # Step 3: Check isotope presence
  message("Step 3: Checking isotope presence...")
  isotope_presence_files <- future_lapply(target_csv_files_success, function(f) safe_run(get_isotope_presence, f, output_dir = output_dir))
  
  error_indices_3 <- which(sapply(isotope_presence_files, inherits, "error"))
  if (length(error_indices_3) > 0) {
    message("Errors occurred in get_isotope_presence for the following files:")
    for (i in error_indices_3) {
      cat("File:", target_csv_files_success[[i]], "\nError message:\n")
      print(isotope_presence_files[[i]])
      cat("\n")
    }
  } else {
    message("Isotope presence checking completed successfully.")
  }
  
  isotope_presence_files_success <- isotope_presence_files[!sapply(isotope_presence_files, inherits, "error")]
  
  # Step 4: Filter on retention time
  message("Step 4: Filtering on retention time...")
  reduced_target_info <- fread(reduced_target_info_path, data.table = FALSE)
  filtered_rt_files <- future_lapply(isotope_presence_files_success, function(f) safe_run(get_rtimes, f, reduced_target_info = reduced_target_info, output_dir = output_dir))
  
  error_indices_4 <- which(sapply(filtered_rt_files, inherits, "error"))
  if (length(error_indices_4) > 0) {
    message("Errors occurred in get_rtimes for the following files:")
    for (i in error_indices_4) {
      cat("File:", isotope_presence_files_success[[i]], "\nError message:\n")
      print(filtered_rt_files[[i]])
      cat("\n")
    }
  } else {
    message("All retention time filtering files generated successfully.")
  }
  
  # Step 5: Process max intensity summaries for final CSV output
  message("Step 5: Processing max intensity summaries for CSV output...")
  
  filtered_rt_files_success <- filtered_rt_files[!sapply(filtered_rt_files, inherits, "error")]
  
  if (length(filtered_rt_files_success) == 0) {
    warning("No filtered retention time files available to process for max intensity summary.")
    output_csv <- NULL
  } else {
    output_csv <- process_max_intensity_summary(filtered_rt_files_success, output_dir, class_regex)
    message("Max intensity summary CSV saved to: ", output_csv)
  }
  
  any_errors <- length(error_indices_1) > 0 || length(error_indices_2) > 0 ||
    length(error_indices_3) > 0 || length(error_indices_4) > 0
  
  if (!any_errors) {
    message("Pipeline completed successfully. Output files are in: ", normalizePath(output_dir))
  } else {
    message("Pipeline completed with some errors. Please check messages above.")
  }
  
  # Reset plan to sequential to avoid side effects
  plan(sequential)
  
  return(list(
    summary_files = summary_files,
    target_csv_files = target_csv_files,
    isotope_presence_files = isotope_presence_files,
    filtered_rt_files = filtered_rt_files,
    max_intensity_summary_csv = output_csv
  ))
}


# ---- UI ----
ui <- fluidPage(
  useShinyjs(),
  titlePanel("mzML Preprocessing + Multi-Sample Classifier"),
  sidebarLayout(
    sidebarPanel(
      h4("Preprocessing Section"),
      fileInput("mzml_files", "Upload mzML Files", multiple = TRUE, accept = ".mzML"),
      fileInput("target_csv", "Upload Target Data CSV", accept = ".csv"),
      fileInput("reduced_target_csv", "Upload Reduced Target Info CSV", accept = ".csv"),
      actionButton("preprocess_btn", "Run Preprocessing"),
      tags$hr(),
      h4("Prediction Section"),
      fileInput("sample_file", "Upload model_ready_data.csv", accept = ".csv"),
      actionButton("predict_btn", "Predict"),
      downloadButton("download_predictions", "Download Predictions"),
      helpText("Upload a CSV or run preprocessing first. Only the first batch of results is shown.")
    ),
    mainPanel(
      textOutput("status_text"),
      br(),
      tableOutput("prediction_table")
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  predictions <- reactiveVal(NULL)
  visible_batch <- reactiveVal(NULL)
  is_running <- reactiveVal(FALSE)
  preprocessed_file <- reactiveVal(NULL)
  
  output$status_text <- renderText({
    if (is_running()) {
      "Running... Please wait."
    } else if (!is.null(predictions()) && nrow(predictions()) > nrow(visible_batch())) {
      paste("Showing first batch of", nrow(visible_batch()), "samples. Download full results below.")
    } else {
      ""
    }
  })
  
  output$prediction_table <- renderTable({
    preds <- visible_batch()
    if (is.null(preds) || nrow(preds) == 0) {
      data.frame(Message = "No predictions yet. Upload data or run preprocessing.")
    } else {
      preds
    }
  })
  
  output$download_predictions <- downloadHandler(
    filename = function() {
      paste0("model_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      preds <- predictions()
      if (is.null(preds) || nrow(preds) == 0) {
        writeLines("No predictions available to download.", file)
      } else {
        readr::write_csv(preds, file)
      }
    }
  )
  
  observeEvent(input$predict_btn, {
    file_to_use <- if (!is.null(preprocessed_file())) preprocessed_file() else input$sample_file$datapath
    req(file_to_use)
    
    is_running(TRUE)
    on.exit(is_running(FALSE))
    
    # Update status text and flush UI immediately
    shinyjs::runjs("Shiny.setInputValue('flush', Math.random());")
    Sys.sleep(0.1)  # optional small pause to ensure flush
    
    all_results <- list()
    
    withProgress(message = "Running prediction...", value = 0, {
      
      # Now do the heavy work inside withProgress
      sample_df <- read_csv(file_to_use, show_col_types = FALSE)
      if (nrow(sample_df) < 1) {
        predictions(data.frame(Error = "No samples provided."))
        return()
      }
      
      n_samples <- nrow(sample_df)
      
      # Reset workers before prediction
      plan(sequential)
      
      if (n_samples < 20) {
        # Sequential prediction
        for (i in seq_len(n_samples)) {
          row_data <- sample_df[i, , drop = FALSE]
          pred <- predict_single_sample(
            row_data,
            models,
            selected_features,
            class_levels,
            accuracy_thresholds
          )
          all_results[[i]] <- cbind(Sample = i, pred)
          incProgress(i / n_samples)
        }
      } else {
        # Parallel prediction
        total_cores <- parallel::detectCores()
        safe_max_workers <- min(25, total_cores)
        used_cores <- max(1, floor(0.8 * safe_max_workers))
        
        plan(multisession, workers = used_cores)
        
        batch_size <- ceiling(n_samples / used_cores)
        sample_batches <- split(sample_df, ceiling(seq_len(n_samples) / batch_size))
        
        for (i in seq_along(sample_batches)) {
          this_batch <- sample_batches[[i]]
          batch_results <- future_map_dfr(
            1:nrow(this_batch),
            function(j) {
              row_data <- this_batch[j, , drop = FALSE]
              pred <- predict_single_sample(
                row_data,
                models,
                selected_features,
                class_levels,
                accuracy_thresholds
              )
              cbind(Sample = (i - 1) * batch_size + j, pred)
            },
            .options = furrr_options(packages = c("ranger"))
          )
          all_results[[i]] <- batch_results
          incProgress(i / length(sample_batches))
        }
      }
    })
    
    final_df <- bind_rows(all_results)
    predictions(final_df)
    output$status_text <- renderText("Prediction complete.")
    visible_batch(head(final_df, min(10, n_samples)))  # Show first 10 samples
  })
  
  observeEvent(input$sample_file, {
    preprocessed_file(NULL)        # Clear preprocessed file
    predictions(NULL)              # Clear full predictions
    visible_batch(NULL)            # Clear preview batch
    output$status_text <- renderText("")  # Clear status text if used
    shinyjs::runjs("Shiny.setInputValue('flush', Math.random());")
    Sys.sleep(0.1) 
  })
  
  
  
  # ---- Preprocessing Pipeline ----
  observeEvent(input$preprocess_btn, {
    req(input$mzml_files, input$target_csv, input$reduced_target_csv)
    
    is_running(TRUE)
    on.exit(is_running(FALSE))
    
    output$status_text <- renderText("Running preprocessing...")
    
    withProgress(message = "Running preprocessing...", value = 0.1, {
      mzml_paths <- sapply(input$mzml_files$datapath, normalizePath)
      target_path <- input$target_csv$datapath
      reduced_target_path <- input$reduced_target_csv$datapath
      output_dir <- file.path("preprocessed_results")
      dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
      
      cores <- max(1, floor(0.8 * parallel::detectCores()))
      
      results <- run_pipeline(
        mzML_files = mzml_paths,
        target_data_path = target_path,
        reduced_target_info_path = reduced_target_path,
        output_dir = output_dir,
        cores = cores,
        class_regex = ".*180K_([^_]+)_.*"
      )
      
      if (!is.null(results$max_intensity_summary_csv)) {
        preprocessed_file(results$max_intensity_summary_csv)
        sample_df <- read_csv(results$max_intensity_summary_csv, show_col_types = FALSE)
        
        output$status_text <- renderText("Preprocessing complete. Running classification...")
        predictions(NULL)
        visible_batch(NULL)
        
        n_samples <- nrow(sample_df)
        
        all_results <- list()
        if (n_samples < 20) {
          # Run sequential
          for (i in seq_len(n_samples)) {
            row_data <- sample_df[i, , drop = FALSE]
            pred <- predict_single_sample(
              row_data,
              models,
              selected_features,
              class_levels,
              accuracy_thresholds
            )
            all_results[[i]] <- cbind(Sample = i, pred)
            incProgress(i / n_samples)
          }
          final_df <- bind_rows(all_results)
        } else {
          # Run in parallel
          batch_size <- cores
          sample_batches <- split(sample_df, ceiling(seq_len(n_samples) / batch_size))
          plan(multisession, workers = batch_size)
          
          for (i in seq_along(sample_batches)) {
            this_batch <- sample_batches[[i]]
            batch_results <- future_map_dfr(
              1:nrow(this_batch),
              function(j) {
                row_data <- this_batch[j, , drop = FALSE]
                pred <- predict_single_sample(
                  row_data,
                  models,
                  selected_features,
                  class_levels,
                  accuracy_thresholds
                )
                cbind(Sample = (i - 1) * batch_size + j, pred)
              },
              .options = furrr_options(packages = c("ranger"))
            )
            all_results[[i]] <- batch_results
            incProgress(i / length(sample_batches))
          }
          final_df <- bind_rows(all_results)
        }
        
        predictions(final_df)
        visible_batch(head(final_df, 10))
        
        output$status_text <- renderText("Preprocessing and prediction complete.")
      } else {
        output$status_text <- renderText("Preprocessing failed. Check inputs.")
      }
    })
  })
}

shinyApp(ui, server)
