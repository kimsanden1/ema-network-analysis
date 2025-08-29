# main-analysis.R
# author: Kim P.C. van der Sanden
# date: 12/08/2025

# Load required packages
library(tidyr)
library(dplyr)
library(mlVAR)
library(qgraph)
library(ggplot2)
library(mnet)

#-------------------------------------------------------------------------------
# Main analysis
#-------------------------------------------------------------------------------
# mlVAR modeling
# Using the imputed dataset with the following combination:
# Kalman filtering applied to participants with at least 75% compliance

# Nodes:
# Fat = Phyiscal fatigue
# ConR = Concentration problems (reverse coded)
# TdoR = Motivation problems (to do fun things, reverse-coded)
# PoAf = Positive affect (aggregate of I feel happy / relaxed / satisfied)
# NeAf_sqrt = sqrt transformed negative affect (frustrated, anxious, gloomy)
# Pain_sqrt = sqrt transformed pain
# PA = Physical activity (retrospective phrasing)

# Load imputed dataset
data <- readRDS("path/to/dataset.rds")

# Define variable sets for different network analyses:
# - vars_temp: Used to estimate the temporal network (excludes PA)
# - vars_full: Used for contemporaneous and between-subject networks
vars_temp <- c("Fat", "ConR", "TdoR", "PoAf", "NeAf_sqrt", "Pain_sqrt") 
vars_full <- c("Fat", "ConR", "TdoR", "PoAf", "NeAf_sqrt", "Pain_sqrt", "PA")

# Summarize mean and standard deviation of items
data %>%
  summarise(across(all_of(vars_full), list(mean = ~mean(.x, na.rm = TRUE),
                                           sd = ~sd(.x, na.rm = TRUE)))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "^([^_]+(?:_[^_]+)*)_(mean|sd)$"
  )

# --------- Network estimation with mlVAR for contemp and temp networks --------
# Estimate contemporaneous network
mlvar_con <- mlVAR(
  data = data,              # Use the imputed dataset
  vars = vars_full,         # Model all relevant vars
  idvar = "subjno",          # Participant ID col
  dayvar = "Daynumber",     # Day indicator
  beepvar = "Bleepnumber",  # Beep indicator within each day
  lags = 1,                 # Use one lag (ie previous beep)
  estimator = "lmer",       # Use linear mixed-effects regression
  temporal = "orthogonal",  # Ensure temporal effects are orthogonal
  verbose = FALSE,          # Suppress verbose output
  scaleWithin = TRUE        # Scale variables WITHIN each participant
)

# Estimate temporal network with lag 1
mlvar_temp <- mlVAR(
  data = data,
  vars = vars_temp,
  idvar = "subjno",
  dayvar = "Daynumber",
  beepvar = "Bleepnumber",
  lags = 1,
  estimator = "lmer",
  temporal = "orthogonal",
  verbose = FALSE,
  scaleWithin = TRUE
)

# Estimate temporal network with lag 2
mlvar_temp2 <- mlVAR(
  data = data,
  vars = vars_temp,
  idvar = "subjno",
  dayvar = "Daynumber",
  beepvar = "Bleepnumber",
  lags = 2,
  estimator = "lmer",
  temporal = "orthogonal",
  verbose = FALSE,
  scaleWithin = TRUE
)

# Extract matrices
cmat <- getNet(mlvar_con, "contemporaneous")
tmat <- getNet(mlvar_temp, "temporal")
tmat2 <-getNet(mlvar_temp2, "temporal")

# ---------------------------- Plotting the networks ---------------------------
# Define node labels and colors
node_labels <- c(
  "Fat" = "F1",
  "ConR" = "F2",
  "TdoR" = "F3",
  "PoAf" = "A+",
  "NeAf_sqrt" = "A-",
  "Pain_sqrt" = "PA",
  "PA" = "PH"
)

node_colors <- c(
  "Fat" = "#FFD700",
  "ConR" = "#FFD700",
  "TdoR" = "#FFD700",
  "PoAf" = "#1E90FF",
  "NeAf_sqrt" = "#1E90FF",
  "Pain_sqrt" = "#9370DB",
  "PA" = "#9370DB"
)

# Set up shared layout for plotting
layout_config_contemp <- qgraph(cmat, # just pick one
                                layout = "spring", DoNotPlot = TRUE)$layout
layout_config <- qgraph(tmat, # just pick one 
                        layout = "spring", DoNotPlot = TRUE)$layout
rownames(layout_config) <- colnames(tmat)

# Contemporaneous network plot
# png("figures/contemporaneous_network.png", width = 1250, height = 1200)
qgraph(cmat,
       layout = layout_config_contemp,
       labels = node_labels[colnames(cmat)],
       color = node_colors[colnames(cmat)],
       theme = "classic",
       vsize = 8,  
       label.cex = 1.5,  
       edge.labels = TRUE,
       edge.label.cex = 1.2,
       rescale = TRUE)
# dev.off()

# Temporal network plots
# lag-1 
# png("figures/temporal_lag1_network.png", width = 1250, height = 1200)
qgraph(tmat,
       layout = layout_config,
       labels = node_labels[colnames(tmat)],
       color = node_colors[colnames(tmat)],
       theme = "classic",
       vsize = 8,
       label.cex = 1.5,  
       edge.labels = TRUE,
       edge.width = 1.3,
       edge.label.cex = 1.2,
       esize = 12,   
       asize = 7.2,
       rescale = TRUE)
# dev.off()

# lag-2
# png("figures/temporal_lag2_network.png", width = 1250, height = 1200)
qgraph(tmat2,
       layout = layout_config,
       labels = node_labels[colnames(tmat2)],
       color = node_colors[colnames(tmat2)],
       theme = "classic",
       vsize = 8,  
       label.cex = 1.5,  
       edge.labels = TRUE,
       edge.width = 1.3,
       edge.label.cex = 1.2,
       esize = 12,   
       asize = 7.2,
       rescale = TRUE)
# dev.off()

# ----------------------------- Centrality indices ----------------------------- 
# Get qgraph objects for centrality
contemp_qgraph <- qgraph(cmat, layout = "spring", DoNotPlot = TRUE)
temporal_qgraph1 <- qgraph(tmat, layout = "spring", DoNotPlot = TRUE)
temporal_qgraph2 <- qgraph(tmat2, layout = "spring", DoNotPlot = TRUE)

# Contemporaneous strength plot
cp <- centralityPlot(contemp_qgraph,
                     scale = "raw0",
                     include = "Strength",
                     labels = node_labels[colnames(cmat)])
plotdata_c <- cp$data

node_order_c <- c("PA","A-","A+","F3","F2","F1", "PH") # Correct order of nodes
plotdata_c$node <- factor(plotdata_c$node, levels = node_order_c)

cplot <- ggplot(plotdata_c, aes(x = value, y = node)) +
  geom_path(group = 1, linewidth = 0.8, colour = "gray20") +
  geom_point(size = 2.5, colour = "gray20") +
  facet_grid(~measure, scales = "free") +
  scale_y_discrete(limits = node_order_c) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        strip.text = element_text(size = 16, face="bold")
  ) +
  xlim(0, NA) +
  labs(x = NULL, y = NULL)

# ggsave("figures/contemporaneous_centrality_plot.png", cplot, width = 4, height = 8, dpi = 100)

# Temporal InStrength + OutStrength
tp1 <- centralityPlot(temporal_qgraph1,
                      scale = "raw0",
                      include = c("InStrength", "OutStrength"),
                      labels = node_labels[colnames(tmat)])$data
tp2 <- centralityPlot(temporal_qgraph2,
                      scale = "raw0",
                      include = c("InStrength", "OutStrength"),
                      labels = node_labels[colnames(tmat2)])$data

# Tag and combine dataframes
tp1$model <- "lag1"
tp2$model <- "lag2"
plotdata_t <- rbind(tp1, tp2)

# Correct node order
node_order_t <- c("PA","A-","A+","F3","F2","F1") # PA, NeAf, PoAf, TdoR, ConR, Fat
plotdata_t$node <- factor(plotdata_t$node, levels = node_order_t)

tplot_combo <- ggplot(plotdata_t, aes(x = value, y = node, group = model, 
                                      colour = model, linetype = model, shape = model)) +
  geom_path(linewidth = 0.8) +
  geom_point(size = 3) +
  facet_grid(~measure, scales = "free") +
  scale_y_discrete(limits = node_order_t) +
  scale_colour_manual(
    name = "Temporal model",
    values = c(lag1 = "black", lag2 = "grey50")
  ) +
  scale_linetype_manual(
    name = "Temporal model",
    values = c(lag1 = "solid", lag2 = "dashed")
  ) +
  scale_shape_manual(
    name = "Temporal model",
    values = c(lag1 = 16, lag2 = 17)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 18),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold")
  ) +
  labs(x = NULL, y = NULL) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = function(x) {
      mx <- max(x, na.rm = TRUE)
      if (mx <= 0.25) {
        seq(0, 0.20, by = 0.05)
      } else {
        seq(0, 0.3, by = 0.1)
      }
    }
  ) +
  coord_cartesian(clip = "off")
# ggsave("figures/inout_strength_lag_comparison.png", tplot_combo, width = 7.5, height = 8, dpi = 100)

# ------------------ Stability analysis (case-drop bootstrap) ------------------
# Assess stability of centrality indices, ie, in- and out-strength, strength
# in the networks estimated with mlVAR using a case-drop bootstrap
# Some context on the procedure:
# 1. Resample 80% of pts for each bootstrap iteration
# 2. Estimate a new mlVAR model for each subsample
# 3. Extract centrality indices
# 4. Compute spearman correlations between bootstrap and full-model centrality
# 5. Histogram for visualization

# Rename columns for clarity in plots
data <- data %>% 
  rename(
    F1 = Fat,
    F2 = ConR,
    F3 = TdoR,
    Ap = PoAf,
    An = NeAf_sqrt,
    PA = Pain_sqrt,
    PH = PA
  )

# Define column order
new_order <- c(colnames(data)[1:3], "F1", "F2", "F3", "Ap", "An", "PA", "PH", 
               colnames(data)[!colnames(data) %in% c("F1", "F2", "F3", "Ap", "An", "PA", "PH", colnames(data)[1:3])])

# Reorder the dataset
data <- data[, new_order]
data$subjno <- as.integer(factor(data$subjno))

# --- Temporal network stability analysis ---
# Full model estimation
full_model <- mlVAR(
  data = data,
  vars = colnames(data)[4:9],
  idvar = "subjno",
  dayvar = "Daynumber",
  beepvar = "Bleepnumber",
  lags = 1,
  estimator = "lmer",
  temporal = "orthogonal",
  verbose = FALSE,
  scaleWithin = TRUE
)

# Get full model temporal network and centrality
full_model_net <- getWmat(plot(full_model, "temporal", DoNotPlot = TRUE))
full_centrality <- centrality_auto(full_model_net)

# Storage matrices for bootstrap results
n_nodes <- length(colnames(data)[4:9])
t_in_strength <- matrix(ncol = n_nodes)
t_out_strength <- matrix(ncol = n_nodes)

# Bootstrapping loop
set.seed(123)
n_iterations <- 100
all_ids <- unique(data$subjno)
sample_size <- floor(length(all_ids) * 0.8)
counter <- 0
Sys.time()
while (counter < n_iterations) {
  sampled_ids <- sample(all_ids, sample_size)
  sub_data <- data[data$subjno %in% sampled_ids, ]
  
  tryCatch({
    mod <- mlVAR(
      data = sub_data,
      vars = colnames(sub_data)[4:9],
      idvar = "subjno",
      dayvar = "Daynumber",
      beepvar = "Bleepnumber",
      lags = 1,
      estimator = "lmer",
      temporal = "orthogonal",
      verbose = FALSE,
      scaleWithin = TRUE
    )
    
    temp_net <- getWmat(plot(mod, "temporal", DoNotPlot = TRUE))
    centrality <- centrality_auto(temp_net)
    
    t_in_strength <- rbind(t_in_strength, centrality$node.centrality$InStrength)
    t_out_strength <- rbind(t_out_strength, centrality$node.centrality$OutStrength)
    
    counter <- counter + 1
  }, error = function(e) {
    message("skipped iteration due to error: ", e$message)
  })
}

# Compute spearman correlations with full model
cor_in <- apply(t_in_strength[-1, ], 1, function(x) cor(x, full_centrality$node.centrality$InStrength, method = "spearman"))
cor_out <- apply(t_out_strength[-1, ], 1, function(x) cor(x, full_centrality$node.centrality$OutStrength, method = "spearman"))


# --- Contemporaneous network stability analysis ---

# Full model estimation
full_cmodel <- mlVAR(
  data = data,
  vars = colnames(data)[4:10],
  idvar = "subjno",
  dayvar = "Daynumber",
  beepvar = "Bleepnumber",
  lags = 1,
  estimator = "lmer",
  temporal = "orthogonal",
  verbose = FALSE,
  scaleWithin = TRUE
)

# Get full model contemporaneous network and centrality
full_cmodel_net <- getWmat(plot(full_cmodel, "contemporaneous", DoNotPlot = TRUE))
full_c_centrality <- centrality_auto(full_cmodel_net)

# Storage matrix for bootstrap results
n_nodes <- length(colnames(data)[4:10])
c_strength <- matrix(ncol = n_nodes)

# Bootstrapping loop
set.seed(123)
counter <- 0
n_iterations <- 100
all_ids <- unique(data$subjno)
sample_size <- floor(length(all_ids) * 0.8)

while (counter < n_iterations) {
  sampled_ids <- sample(all_ids, sample_size)
  sub_data <- data[data$subjno %in% sampled_ids, ]
  
  tryCatch({
    cmod <- mlVAR(
      data = sub_data,
      vars = colnames(sub_data)[4:10],
      idvar = "subjno",
      dayvar = "Daynumber",
      beepvar = "Bleepnumber",
      lags = 1,
      estimator = "lmer",
      temporal = "orthogonal",
      verbose = FALSE,
      scaleWithin = TRUE
    )
    
    contemp_net <- getWmat(plot(cmod, "contemporaneous", DoNotPlot = TRUE))
    centrality <- centrality_auto(contemp_net)
    
    c_strength <- rbind(c_strength, centrality$node.centrality$Strength)
    
    counter <- counter + 1
  }, error = function(e) {
    message("Skipped iteration due to error: ", e$message)
  })
}
Sys.time()

# Compute Spearman correlations with full model
cor_str <- apply(c_strength[-1, ], 1, function(x) cor(x, full_c_centrality$node.centrality$Strength, method = "spearman"))

# --- Plotting stability results ---

# Create a tidy dataframe
df <- tibble(
  InStrength = cor_in,
  OutStrength = cor_out,
  ContemporaneousStrength = cor_str
) %>%
  pivot_longer(cols = everything(), names_to = "CentralityType", values_to = "Correlation")

# Label mapping for facets
centrality_labels <- c(
  InStrength = "Temporal in-strength",
  OutStrength = "Temporal out-strength",
  ContemporaneousStrength = "Contemporaneous strength"
)

# Histogram plot
histogram_plot <- ggplot(df, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.025, fill = "grey90", color = "black", boundary = 0) +
  facet_wrap(~ CentralityType, nrow = 1,
             labeller = labeller(CentralityType = centrality_labels)) +
  theme_bw(base_size = 10) +
  labs(
    x = "Spearman correlation with full model centrality",
    y = "Frequency"
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(strip.text = element_text(size = 10))

# Save plot
ggsave("figures/centrality_stability_histograms.png", histogram_plot, width = 12, height = 4)

