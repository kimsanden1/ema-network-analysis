# sensitivity-analysis.R
# author: Kim P.C. van der Sanden
# date: 12/08/2025

# Load required packages
library(dplyr)
library(mice)
library(imputeTS)
library(mlVAR)
library(qgraph)

# Load dataset
data <- readRDS("path/to/data.rds")

#-------------------------------------------------------------------------------
# Sensitivity analysis
#-------------------------------------------------------------------------------

#---- Helper function to filter participants based on compliance rate ----
filter_compliance <- function(data, threshold, var = "Fat") {
  compliant_ids <- data %>%
    group_by(subjno) %>%
    summarise(rate = sum(!is.na(.data[[var]])) / n()) %>%
    filter(rate >= threshold) %>%
    pull(subjno)
  return(compliant_ids)
}

# Define compliance thresholds (0%, 25%, 50%, 75%)
compliance_lists <- list(
  "00" = unique(data$subjno),
  "25" = filter_compliance(data, 0.25, "Fat"),
  "50" = filter_compliance(data, 0.50, "Fat"),
  "75" = filter_compliance(data, 0.75, "Fat")
)

vars <- c("Fat", "ConR", "TdoR", "PoAf", "NeAf", "Pain", "PA")

#--- IMPUTATION METHODS ---

#---- (A) No imputation reference datasets ----
# For each threshold, simply filter data
noimp <- list()
for (thr in names(compliance_lists)) {
  noimp[[paste0("noimp_", thr)]] <- data %>%
    filter(subjno %in% compliance_lists[[thr]])
}

# ---- (B) Kalman filtering imputation datasets ----
# Define the Kalman filtering function
kalman_filtering <- function(data_imp, participant_ids) {
  participant_ids <- as.character(participant_ids)
  data_klmn <- data %>%
    mutate(subjno = as.character(subjno)) %>%
    filter(subjno %in% participant_ids)
  optim_control_params = list(maxit=100)
  for (var in vars) {
    for (participant in participant_ids) {
      participant_rows <- data_klmn$subjno == participant
      time_series <- data_klmn[participant_rows, var]
      if (all(is.na(time_series))) {
        next 
      }
      min_val <- min(time_series, na.rm = TRUE)
      max_val <- max(time_series, na.rm = TRUE)
      if (min_val == max_val) {
        time_series[is.na(time_series)] <- min_val
      } else {
        normalized_series <- (time_series - min_val) / (max_val - min_val)
        imputed_series <- na_kalman(normalized_series, model = "StructTS",
                                    optim.control = optim_control_params)
        time_series <- imputed_series * (max_val - min_val) + min_val
      }
      data_klmn[participant_rows, var] <- time_series
    }
    data_klmn[[var]] <- pmax(pmin(data_klmn[[var]], 10), 0)
  }
  return(data_klmn)
}

# Apply Kalman filtering to all datasets
kalman_datasets <- list()
for (thr in names(compliance_lists)) {
  kalman_datasets[[paste0("klmn_", thr)]] <- kalman_filtering(data, compliance_lists[[thr]])
}

# ---- (C) Multilevel MICE imputation datasets ----
# Uses the mice package 
# Multilevel approach, treating subjects as a level-2 variable
run_mlMICE <- function(data_subset) {
  data_miceti <- data_subset %>%
    mutate(TimeIndex = (Daynumber - 1) * 5 + Bleepnumber) %>%
    select(subjno, TimeIndex, Bleepnumber, all_of(vars))
  
  subjno_lookup <- data_miceti %>% distinct(subjno) %>% mutate(subjno_factor = as.integer(as.factor(subjno)))
  
  data_miceti <- data_miceti %>%
    mutate(subjno_factor = as.integer(as.factor(subjno))) %>%
    select(subjno_factor, TimeIndex, Bleepnumber, all_of(vars))
  
  pred_matrix <- make.predictorMatrix(data_miceti)
  pred_matrix[,] <- 0
  pred_matrix[vars, vars] <- 1
  pred_matrix[vars, "subjno_factor"] <- -2
  pred_matrix[vars, "TimeIndex"] <- 1
  pred_matrix <- pred_matrix[vars, ]
  
  bl <- setNames(as.list(vars), vars)
  meth <- rep("2l.pan", length(vars))
  names(meth) <- vars
  
  impti <- mice(data_miceti,
                method = meth,
                predictorMatrix = pred_matrix,
                blocks = bl,
                m = 5,
                maxit = 10,
                seed = 123,
                printFlag = FALSE)
  
  data_micedti <- complete(impti, "long") %>%
    group_by(.id) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE))) %>%
    ungroup()
  
  data_micedti <- data_micedti %>%
    mutate(subjno_factor = ceiling(.id / 105),
           Bleepnumber = rep(1:5, times = 21, length.out = n()),
           Daynumber = rep(1:21, each = 5, length.out = n()),
           TimeIndex = (Daynumber - 1) * 5 + Bleepnumber) %>%
    left_join(subjno_lookup, by = "subjno_factor") %>%
    select(-subjno_factor, -.id) %>%
    relocate(subjno)
  
  for (v in vars) {
    data_micedti[[v]] <- pmax(pmin(data_micedti[[v]], 10), 0)
  }
  return(as.data.frame(data_micedti))
}

miceti_list <- list()
for (thr in names(compliance_lists)) {
  data_subset <- data %>% filter(subjno %in% compliance_lists[[thr]])
  miceti_list[[paste0("mice_", thr)]] <- run_mlMICE(data_subset)
}

# --- Plotting networks from the 12 imputed datasets ---
# Helper function to run mlVAR on a single dataset
run_mlvar_single <- function(dataset,
                             vars,
                             subjno = "subjno",
                             dayvar = "Daynumber",
                             beepvar = "Bleepnumber",
                             lags = 1,
                             estimator = "lmer",
                             temporal = "orthogonal",
                             plot = TRUE,
                             name = "Dataset") {
  missing_vars <- setdiff(vars, names(dataset))
  if (length(missing_vars) > 0) {
    stop("Error: Dataset is missing required variables: ", paste(missing_vars, collapse = ", "))
  }
  
  if (is.null(dataset)) {
    stop("Error: Dataset is NULL.")
  }
  
  # Run mlVAR
  result <- mlVAR(
    data = dataset,
    vars = vars,
    idvar = subjno,
    lags = lags,
    dayvar = dayvar,
    beepvar = beepvar,
    estimator = estimator,
    temporal = temporal,
    verbose = FALSE,
    scaleWithin = TRUE
  )
  
  summary_obj <- summary(result)
  
  return(list(
    result = result,
    temporal = result$results$Beta,
    contemporaneous = summary_obj$contemporaneous,
    between = summary_obj$between
  ))
}

# --- Streamlined mlVAR pipeline ---
# Define method and threshold lists
methods <- list(
  noimp = noimp,
  klmn = klmn,
  mice = miceti_list
)
method_names <- names(methods)
thresholds <- c("00", "25", "50", "75")
# Store results
networks_temp <- list()
networks_contemp <- list()
# Add a container
model_summaries <- list()
# Loop to run mlVAR on each dataset
Cnow <- Sys.time()
for (method in method_names) {
  networks_temp[[method]] <- list()
  networks_contemp[[method]] <- list()
  
  model_summaries[[method]] <- list()
  
  for (thr in thresholds) {
    data <- methods[[method]][[paste0(method, "_", thr)]]
    
    if (is.null(data)) {
      warning(paste("Missing dataset:", method, thr))
      next
    }
    
    if (!"ConR" %in% names(data)) data$ConR <- 10 - data$Con
    if (!"TdoR" %in% names(data)) data$TdoR <- 10 - data$Tdo
    
    model_out <- run_mlvar_single(
      dataset = data,
      vars = vars,
      plot = FALSE,
      name = paste(toupper(method), "-", thr, "%")
    )
    
    # Store network components
    model_summaries[[method]][[thr]] <- model_out
    networks_temp[[method]][[thr]] <- getNet(model_out$result, "temporal")
    networks_contemp[[method]][[thr]] <- getNet(model_out$result, "contemporaneous")
  }
}
str(model_summaries, max.level=2)
# Pick any available model as reference
ref_model <- model_summaries[["noimp"]][["00"]]$result

# Extract fixed layout from the contemporaneous network
layout_contemp <- qgraph(getNet(ref_model, "contemporaneous"), layout = "circle", DoNotPlot = TRUE)$layout
layout_temporal <- qgraph(getNet(ref_model, "temporal"), layout = "circle", DoNotPlot = TRUE)$layout


# ---- Plot all contemporaneous networks in a 3x4 grid ----
node_labels <- c(
  Fat       = "F1",
  ConR      = "F2",
  TdoR      = "F3",
  PoAf      = "A+",
  NeAf      = "A-",
  Pain      = "PA",
  PA        = "PH"
)

par(mfrow = c(3, 4), mar = c(2, 2, 3, 1))

for (method in method_names) {
  for (thr in thresholds) {
    model <- model_summaries[[method]][[thr]]$result
    
    if (is.null(model)) next
    
    plot_title <- paste(toupper(method), "-", thr, "%")
    
    cmat <- getNet(model, "contemporaneous")  # Extract matrix
    
    qgraph(cmat,
           layout = layout_contemp,
           labels = node_labels[colnames(cmat)],
           title = plot_title,
           title.cex = 1.2,
           theme = "classic",
           label.cex = 1.5,
           label.prop = 1,
           fade = FALSE,
           vsize = 12,
           shape = "circle",
           edge.width = 2,
           edge.labels = FALSE)
  }
}
grid_plot_contemp <- recordPlot()

# ---- Plot all temporal networks in a 3x4 grid ----
par(mfrow = c(3, 4), mar = c(2, 2, 3, 1))

for (method in method_names) {
  for (thr in thresholds) {
    model <- model_summaries[[method]][[thr]]$result
    
    if (is.null(model)) next
    
    plot_title <- paste(toupper(method), "-", thr, "%")
    
    tmat <- getNet(model, "temporal")  # Extract temporal network matrix
    
    qgraph(tmat,
           layout = layout_temporal,
           labels = node_labels[colnames(tmat)],
           title = plot_title,
           title.cex = 1.2,
           theme = "classic",
           label.cex = 1.5,
           label.prop = 1,
           fade = FALSE,
           edge.labels = FALSE,
           vsize = 12,
           shape = "circle",
           edge.width = 3)
  }
}
grid_plot_temp <- recordPlot()

file_out <- function(type) paste0("networks_", type, "_", Sys.Date(), ".png")

png(file_out("contemporaneous"), width = 1600, height = 1200, res = 200)
replayPlot(grid_plot_contemp)
dev.off()

png(file_out("temporal"), width = 1600, height = 1200, res = 200)
replayPlot(grid_plot_temp)
dev.off()

#---- Save plots as pdf
file_out_pdf <- function(type) paste0("networks_", type, "_", Sys.Date(), ".pdf")

# Save contemporaneous networks
pdf(file_out_pdf("contemporaneous"), width = 12, height = 9)
replayPlot(grid_plot_contemp)
dev.off()

# Save temporal networks
pdf(file_out_pdf("temporal"), width = 12, height = 9)
replayPlot(grid_plot_temp)
dev.off()


