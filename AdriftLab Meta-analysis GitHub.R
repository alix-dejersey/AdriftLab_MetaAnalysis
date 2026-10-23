# ==============================================================================
# Title:  AdriftLab Meta-analysis GitHub V2
# Author: Alix de Jersey
# Date:   2026-06-07
# Paper:  Low toxicity thresholds for plastic exposure in Shearwaters
#
# Description:
#   Models the relationship between ingested plastic exposure (plastic mass in
#   the proventriculus/stomach) and sublethal biological responses across
#   tissues, organs, and biological fluids in seabirds.
#
#   For each response variable, the workflow compares:
#     1. A linear model using a power-transformed plastic mass predictor and
#        a Box-Cox transformed, Z-scored biological response.
#     2. A four-parameter sigmoidal dose-response (4PL) model on the Z-scored
#        biological response.
#
#   Model selection uses Jacobian-corrected AIC (projecting both models'
#   log-likelihoods back to the raw, untransformed response scale so they are
#   directly comparable) and checks whether the choice differs from the
#   traditional Adjusted R-squared rule.
# ==============================================================================


# ==============================================================================
# 0. Global Options, Packages, Directories, and Helper Functions
# ==============================================================================

# Prevent R from printing very small or very large numbers in scientific notation
options(scipen = 999)

# Load required packages (installing any that are missing)
packages <- c(
  "ggplot2", "viridisLite", "sjPlot", "MASS", "dplyr",
  "drc", "boot", "metafor", "forcats", "viridis", "ggfortify"
)

invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# Set working directory and output directory (defaults to current directory ".")
# setwd("path/to/working/directory")
output_dir <- "."

# Helper function to reverse the power transformation on plastic mass
# (plastic_mass_trans = plastic_mass^predictor_power) back to grams.
inverse_power_transform <- function(x, power) {
  x^(1 / power)
}


# ==============================================================================
# 1. Data Import and Variable Setup
# ==============================================================================

# Load raw dataset
raw_data <- read.csv("pdgfrb.csv", check.names = FALSE)
# metadata <- read.csv("path/to/metadata.csv", check.names = FALSE)

# Set labels and analysis parameters for this biological endpoint:
#   - response_variable: column name in raw_data
#   - response_name:     label used in plots and output CSV files
#   - year:              source study year
#   - predictor_power:   power transformation on plastic mass (chosen via diagnostics)
response_variable <- "response_variable"
response_name     <- "PDGFRB_2023"
year              <- 2023
predictor_power   <- 0.25

# Extract required columns and standardize the response column name
data_frame <- raw_data[, c("sample_id", response_variable, "plastic_mass")]
colnames(data_frame)[colnames(data_frame) == response_variable] <- "response_variable"


# ==============================================================================
# 2. Data Cleaning and Predictor Transformation
# ==============================================================================

set.seed(333)

# Ensure numeric types and remove rows with missing or non-finite values
data_frame <- data_frame[, c("sample_id", "response_variable", "plastic_mass")]
data_frame$response_variable <- as.numeric(as.character(data_frame$response_variable))
data_frame$plastic_mass      <- as.numeric(as.character(data_frame$plastic_mass))
data_frame <- data_frame[
  is.finite(data_frame$response_variable) & is.finite(data_frame$plastic_mass), 
]
rownames(data_frame) <- NULL

# Power-transform plastic mass for the linear model (raw plastic mass is kept for 4PL)
data_frame$plastic_mass_trans <- data_frame$plastic_mass^predictor_power


# ==============================================================================
# Step 1: Linear Model with Box-Cox Transformed Response
# ==============================================================================
# Fit a preliminary linear model to determine the direction of the effect and
# estimate the optimal Box-Cox lambda for the response variable.

lm0 <- lm(response_variable ~ plastic_mass_trans, data = data_frame)

# Align response direction so increasing biological effect is always positive
transform_variable <- ifelse(coef(lm0)["plastic_mass_trans"] < 0, -1, 1)

# Estimate optimal Box-Cox lambda across a 33-point grid from -5 to 5 (includes 0)
bc     <- MASS::boxcox(lm0, lambda = seq(-5, 5, length.out = 33), plotit = FALSE)
lambda <- bc$x[which.max(bc$y)]

# Apply Box-Cox power transformation (using log10 if lambda is effectively 0)
response_variable_bc_unscaled <- if (abs(lambda) < 1e-8) {
  log10(data_frame$response_variable)
} else {
  data_frame$response_variable^lambda
}

# Z-score scale the transformed response and align direction
data_frame$response_variable_bc <- as.numeric(scale(response_variable_bc_unscaled)) * transform_variable

# Fit the final Box-Cox linear model
lm_bc <- lm(response_variable_bc ~ plastic_mass_trans, data = data_frame)

# Diagnostic plots for the linear model
plot(response_variable_bc ~ plastic_mass_trans, data = data_frame, pch = 16, col = "#66000999")
autoplot(lm_bc)


# ==============================================================================
# Step 2: Four-Parameter Logistic (4PL) Sigmoidal Model
# ==============================================================================
# Fit a 4PL dose-response curve on the Z-scored raw response against raw plastic mass.

data_frame$response_variable_raw_scaled <- as.numeric(scale(data_frame$response_variable)) * transform_variable

four_pl <- drm(
  response_variable_raw_scaled ~ plastic_mass,
  data = data_frame,
  fct  = L.4(names = c("slope", "lower", "upper", "EC50"))
)

summary_4pl <- summary(four_pl)

# Basic diagnostic plots for the 4PL model
plot(four_pl, xlim = c(0, 100))
plot(response_variable_raw_scaled ~ plastic_mass_trans, data = data_frame)


# ==============================================================================
# Step 3: Model Fit Statistics & Jacobian Likelihood / AIC Correction
# ==============================================================================
# Calculate Adjusted R-squared, F-test p-values, and Jacobian-corrected
# log-likelihoods and AIC/AICc values projected back onto the raw y-scale.

# 3A. Adjusted R-squared for both models
adj_r2_lm <- summary(lm_bc)$adj.r.squared

yhat       <- predict(four_pl, newdata = data_frame)
n          <- length(data_frame$response_variable_raw_scaled)
RSS        <- sum((data_frame$response_variable_raw_scaled - yhat)^2)
TSS        <- sum((data_frame$response_variable_raw_scaled - mean(data_frame$response_variable_raw_scaled))^2)
adj_r2_4pl <- 1 - (RSS / (n - 4)) / (TSS / (n - 1))

# 3B. Jacobian Log-Likelihood & AIC Correction to the Original Raw Response Scale
y_raw          <- data_frame$response_variable
sd_bc_unscaled <- sd(response_variable_bc_unscaled, na.rm = TRUE)
sd_raw         <- sd(y_raw, na.rm = TRUE)

# Log-Jacobian for lm_bc (accounts for Box-Cox transformation u = y^lambda or log10(y) + scale())
if (abs(lambda) < 1e-8) {
  # For u = log10(y) and z = (u - mean(u)) / sd_u: |dz/dy| = 1 / (sd_u * y * ln(10))
  log_jac_lm <- -n * log(sd_bc_unscaled) - n * log(log(10)) - sum(log(y_raw))
} else {
  # For u = y^lambda and z = (u - mean(u)) / sd_u: |dz/dy| = |lambda| * y^(lambda - 1) / sd_u
  log_jac_lm <- -n * log(sd_bc_unscaled) + n * log(abs(lambda)) + (lambda - 1) * sum(log(y_raw))
}

# Log-Jacobian for four_pl (accounts for scale() on raw y: |dz/dy| = 1 / sd_raw)
log_jac_4pl <- -n * log(sd_raw)

# Project model log-likelihoods onto the raw untransformed response scale
logLik_lm_raw  <- as.numeric(logLik(lm_bc)) + log_jac_lm
logLik_4pl_raw <- as.numeric(logLik(four_pl)) + log_jac_4pl

# Parameter counts (k):
#   lm_bc:   intercept + slope + sigma^2 (df = 3) + 1 for estimated Box-Cox lambda = 4
#   four_pl: 4 curve parameters + sigma^2 (df = 5)
k_lm  <- attr(logLik(lm_bc), "df") + 1
k_4pl <- attr(logLik(four_pl), "df")

# Raw-scale AIC and small-sample corrected AICc
AIC_lm_raw  <- -2 * logLik_lm_raw + 2 * k_lm
AIC_4pl_raw <- -2 * logLik_4pl_raw + 2 * k_4pl

AICc_lm_raw  <- AIC_lm_raw + (2 * k_lm * (k_lm + 1)) / max(1, (n - k_lm - 1))
AICc_4pl_raw <- AIC_4pl_raw + (2 * k_4pl * (k_4pl + 1)) / max(1, (n - k_4pl - 1))

# Delta AIC: Positive = Linear beats 4PL; Negative = 4PL beats Linear
Delta_AIC_4pl_vs_lm <- AIC_4pl_raw - AIC_lm_raw

# 3C. Model significance p-values
Pvalue_Linear <- anova(lm_bc)$"Pr(>F)"[1]

m_null   <- lm(data_frame$response_variable_raw_scaled ~ 1)
RSS_null <- sum(residuals(m_null)^2)
RSS_dr   <- sum(residuals(four_pl)^2)
df_null  <- df.residual(m_null)
df_dr    <- df.residual(four_pl)

F_stat         <- ((RSS_null - RSS_dr) / (df_null - df_dr)) / (RSS_dr / df_dr)
Pvalue_Sigmoid <- 1 - pf(F_stat, df_null - df_dr, df_dr)


# ==============================================================================
# Step 4A: 4PL Prediction Grid, Bootstrap CI, Minimum Threshold, & Bootstrap EC50
# ==============================================================================
# Resample the dataset with replacement to generate 95% confidence intervals
# around the 4PL curve and EC50. The minimum effect threshold is the lowest
# plastic mass where the 95% bootstrap CI no longer overlaps the baseline mean.

grid_n_4pl  <- 2000
boot_number <- 2000

pm_grid <- seq(
  min(data_frame$plastic_mass, na.rm = TRUE),
  max(data_frame$plastic_mass, na.rm = TRUE),
  length.out = grid_n_4pl
)

boot_pred  <- matrix(NA_real_, nrow = boot_number, ncol = length(pm_grid))
boot_EC50  <- rep(NA_real_, boot_number)
x_max_root <- max(1, max(data_frame$plastic_mass, na.rm = TRUE) * 100)

for (i in seq_len(boot_number)) {
  idx_boot <- sample(seq_len(nrow(data_frame)), replace = TRUE)
  d        <- data_frame[idx_boot, ]

  fit_i <- try(
    suppressWarnings(
      drm(
        response_variable_raw_scaled ~ plastic_mass,
        data = d,
        fct  = L.4(names = c("slope", "lower", "upper", "EC50"))
      )
    ),
    silent = TRUE
  )

  if (inherits(fit_i, "try-error")) next

  pred_i <- try(
    suppressWarnings(
      predict(fit_i, newdata = data.frame(plastic_mass = pm_grid))
    ),
    silent = TRUE
  )

  if (inherits(pred_i, "try-error")) next

  boot_pred[i, ] <- as.numeric(pred_i)

  pars_i    <- coef(fit_i)
  asym1_i   <- unname(pars_i["lower:(Intercept)"])
  asym2_i   <- unname(pars_i["upper:(Intercept)"])
  plateau_i <- max(asym1_i, asym2_i)

  response0_i <- try(
    suppressWarnings(
      predict(fit_i, newdata = data.frame(plastic_mass = 0))
    ),
    silent = TRUE
  )

  if (inherits(response0_i, "try-error")) next

  response0_i <- as.numeric(response0_i)
  span_i      <- plateau_i - response0_i

  if (!is.finite(span_i) || abs(span_i) < 1e-12) next

  target_y_50_i <- response0_i + 0.50 * span_i

  f50_i <- function(x) {
    as.numeric(predict(fit_i, newdata = data.frame(plastic_mass = x))) - target_y_50_i
  }

  if (f50_i(0) >= 0) {
    boot_EC50[i] <- 0
  } else if (f50_i(x_max_root) >= 0) {
    boot_EC50[i] <- tryCatch(
      uniroot(f50_i, interval = c(0, x_max_root))$root,
      error = function(e) NA_real_
    )
  }
}

# Retain converged bootstrap iterations and compute pointwise 95% CIs
valid_rows <- rowSums(is.finite(boot_pred)) > 0
boot_mat   <- boot_pred[valid_rows, , drop = FALSE]
pred_fit   <- as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = pm_grid)))

if (nrow(boot_mat) > 0) {
  ci_low  <- apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE)
  ci_high <- apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE)
} else {
  ci_low  <- rep(NA_real_, length(pm_grid))
  ci_high <- rep(NA_real_, length(pm_grid))
}

pred_df <- data.frame(
  plastic_mass = pm_grid,
  fit          = pred_fit,
  lwr_boot     = ci_low,
  upr_boot     = ci_high
)

# Define baseline response (mean at 0 g plastic, or at minimum observed mass if no 0 g birds exist)
baseline <- mean(data_frame$response_variable_raw_scaled[data_frame$plastic_mass == 0], na.rm = TRUE)
if (!is.finite(baseline)) {
  min_mass <- min(data_frame$plastic_mass, na.rm = TRUE)
  baseline <- mean(data_frame$response_variable_raw_scaled[data_frame$plastic_mass == min_mass], na.rm = TRUE)
}

# Identify the minimum plastic mass where the 4PL 95% CI excludes the baseline
outside        <- (pred_df$upr_boot < baseline) | (pred_df$lwr_boot > baseline)
thr_idx        <- which(outside)[1]
threshold_mass <- if (length(thr_idx)) pred_df$plastic_mass[thr_idx] else NA_real_

# Build 4PL ggplot
four_pl_plot <- ggplot() +
  geom_point(
    data  = data_frame,
    aes(x = plastic_mass, y = response_variable_raw_scaled),
    alpha = 0.7,
    size  = 2
  ) +
  geom_ribbon(
    data  = pred_df,
    aes(x = plastic_mass, ymin = lwr_boot, ymax = upr_boot),
    alpha = 0.18
  ) +
  geom_line(
    data      = pred_df,
    aes(x = plastic_mass, y = fit),
    linewidth = 1.1
  ) +
  geom_hline(yintercept = baseline, linetype = "dashed") +
  { if (!is.na(threshold_mass)) geom_vline(xintercept = threshold_mass, linetype = "dotted") else NULL } +
  labs(
    x        = "Plastic mass (g)",
    y        = paste(response_name, "Z score"),
    title    = paste(response_name, "4PL fit with bootstrapped 95% CI vs real data"),
    subtitle = if (!is.na(threshold_mass)) {
      paste0("Threshold (CI excludes zero-plastic mean) \u2248 ", round(threshold_mass, 3), " g")
    } else {
      "No threshold: CI never excludes the zero-plastic mean within the observed range"
    }
  ) +
  scale_x_continuous(trans = "pseudo_log") +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# ==============================================================================
# Step 4B: Linear Prediction Grid and Minimum Threshold
# ==============================================================================
# Estimate the minimum threshold from the Box-Cox linear model as the first
# point where the 95% confidence band departs from the baseline response.

grid_n_lm <- 2000

pm_grid_linear <- seq(
  0,
  max(data_frame$plastic_mass_trans, na.rm = TRUE),
  length.out = grid_n_lm
)

predict_linear <- data.frame(
  predict(
    lm_bc,
    newdata  = data.frame(plastic_mass_trans = pm_grid_linear),
    interval = "confidence",
    level    = 0.95
  )
)

pred_df_linear <- data.frame(
  plastic_mass_trans = pm_grid_linear,
  fit                = predict_linear$fit,
  lwr                = predict_linear$lwr,
  upr                = predict_linear$upr
)

# Define baseline on the Box-Cox transformed scale
baseline_bc <- mean(data_frame$response_variable_bc[data_frame$plastic_mass == 0], na.rm = TRUE)
if (!is.finite(baseline_bc)) {
  min_mass    <- min(data_frame$plastic_mass, na.rm = TRUE)
  baseline_bc <- mean(data_frame$response_variable_bc[data_frame$plastic_mass == min_mass], na.rm = TRUE)
}

# Find first grid index where the 95% CI departs from baseline_bc in the slope direction
slope_bc <- coef(lm_bc)["plastic_mass_trans"]

if (is.finite(slope_bc) && slope_bc > 0) {
  idx <- which(pred_df_linear$lwr > baseline_bc)[1]
} else if (is.finite(slope_bc) && slope_bc < 0) {
  idx <- which(pred_df_linear$upr < baseline_bc)[1]
} else {
  idx <- integer(0)
}

threshold_trans  <- if (length(idx) && !is.na(idx)) pred_df_linear$plastic_mass_trans[idx] else NA_real_
threshold_mass_g <- if (is.finite(threshold_trans)) inverse_power_transform(threshold_trans, predictor_power) else NA_real_

cat("Threshold (transformed scale):", threshold_trans, "\n")
cat("Threshold (grams):", threshold_mass_g, "\n")

# Build Linear Model ggplot
lm_plot <- ggplot() +
  geom_point(
    data  = data_frame,
    aes(x = plastic_mass_trans, y = response_variable_bc),
    alpha = 0.7,
    size  = 2
  ) +
  geom_ribbon(
    data  = pred_df_linear,
    aes(x = plastic_mass_trans, ymin = lwr, ymax = upr),
    alpha = 0.18
  ) +
  geom_line(
    data      = pred_df_linear,
    aes(x = plastic_mass_trans, y = fit),
    linewidth = 1.05
  ) +
  geom_hline(yintercept = baseline_bc, linetype = "dashed") +
  { if (is.finite(threshold_trans)) geom_vline(xintercept = threshold_trans, linetype = "dotted") else NULL } +
  labs(
    x        = bquote("Plastic mass"^.(predictor_power)),
    y        = paste(response_name, "response (Box-Cox scale)"),
    title    = paste(response_name, "linear fit (95% CI)"),
    subtitle = if (is.finite(threshold_mass_g)) {
      paste0("First mass where 95% CI departs from zero-plastic mean: ~", round(threshold_mass_g, 3), " g")
    } else {
      "No threshold found within the observed range"
    }
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.5),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


# ==============================================================================
# Step 4C: Effective Concentrations (EC50, EC90, EC95, EC99 / ECMAX)
# ==============================================================================
# Estimate plastic mass where the fitted 4PL response reaches 50%, 90%, 95%,
# and 99% (ECMAX) of the modelled span from 0 g to the upper plateau.

pars  <- coef(four_pl)
asym1 <- unname(pars["lower:(Intercept)"])
asym2 <- unname(pars["upper:(Intercept)"])

Plateau_response_from0 <- max(asym1, asym2)
Response_0g            <- as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = 0)))
Plateau_span_from0     <- Plateau_response_from0 - Response_0g

target_y_50 <- Response_0g + 0.50 * Plateau_span_from0
target_y_90 <- Response_0g + 0.90 * Plateau_span_from0
target_y_95 <- Response_0g + 0.95 * Plateau_span_from0
target_y_99 <- Response_0g + 0.99 * Plateau_span_from0

f50 <- function(x) as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = x))) - target_y_50
f90 <- function(x) as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = x))) - target_y_90
f95 <- function(x) as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = x))) - target_y_95
f99 <- function(x) as.numeric(predict(four_pl, newdata = data.frame(plastic_mass = x))) - target_y_99

solve_ec <- function(fn, upper_bound) {
  if (fn(0) >= 0) return(0)
  if (fn(upper_bound) < 0) return(NA_real_)
  uniroot(fn, interval = c(0, upper_bound))$root
}

if (!is.finite(Plateau_span_from0) || abs(Plateau_span_from0) < 1e-12) {
  EC50 <- NA_real_
  EC90 <- NA_real_
  EC95 <- NA_real_
  EC99 <- NA_real_
} else {
  EC50 <- solve_ec(f50, x_max_root)
  EC90 <- solve_ec(f90, x_max_root)
  EC95 <- solve_ec(f95, x_max_root)
  EC99 <- solve_ec(f99, x_max_root)
}

ECMAX <- EC99

# Summarise bootstrapped EC50 uncertainty
boot_EC50_valid <- boot_EC50[is.finite(boot_EC50)]
EC50_SE         <- if (length(boot_EC50_valid) > 1) sd(boot_EC50_valid) else NA_real_
EC50_95CI_Lower <- if (length(boot_EC50_valid) > 0) as.numeric(quantile(boot_EC50_valid, 0.025)) else NA_real_
EC50_95CI_Upper <- if (length(boot_EC50_valid) > 0) as.numeric(quantile(boot_EC50_valid, 0.975)) else NA_real_

cat("EC50 from 0 g:", EC50, "\n")
cat("EC90 from 0 g:", EC90, "\n")
cat("EC95 from 0 g:", EC95, "\n")
cat("ECMAX plateau response:", EC99, "\n")
cat("EC50 bootstrap SE:", EC50_SE, "\n")
cat("EC50 bootstrap 95% CI:", EC50_95CI_Lower, "to", EC50_95CI_Upper, "\n")

# Extract slope and correlation statistics from the linear model
ct      <- summary(lm_bc)$coefficients
b       <- abs(ct["plastic_mass_trans", "Estimate"])
se      <- ct["plastic_mass_trans", "Std. Error"]
r       <- sqrt(summary(lm_bc)$r.squared)
r_error <- 1 / sqrt(summary(lm_bc)$df[2] + summary(lm_bc)$df[1])


# ==============================================================================
# Step 5: Model Selection (Jacobian AIC vs Adjusted R-Squared)
# ==============================================================================
# Selects the winning model using Jacobian-adjusted AIC on the raw scale and
# compares that decision directly against the original Adjusted R-squared rule.

max_plastic_mass <- max(data_frame$plastic_mass, na.rm = TRUE)

# Check that the 4PL model reaches EC50 (upper 95% CI) and ECMAX within the observed range
sigmoid_valid <- isTRUE(
  !is.na(EC50_95CI_Upper) && !is.na(ECMAX) &&
  EC50_95CI_Upper <= max_plastic_mass &&
  ECMAX <= max_plastic_mass
)

# 5A. Selection via Jacobian-Corrected AIC (Lower AIC = Better Fit)
model_choice_aic <- if (Pvalue_Linear < 0.05 && Pvalue_Sigmoid > 0.05) {
  "Linear"
} else if (Pvalue_Linear > 0.05 && Pvalue_Sigmoid < 0.05) {
  if (sigmoid_valid) "Sigmoid" else "No significant relationship detected"
} else if (Pvalue_Linear < 0.05 && Pvalue_Sigmoid < 0.05) {
  if (!sigmoid_valid) {
    "Linear"
  } else if (AIC_lm_raw <= AIC_4pl_raw) {
    "Linear"
  } else {
    "Sigmoid"
  }
} else {
  "No significant relationship detected"
}

# 5B. Selection via Adjusted R-squared (Higher Adjusted R2 = Better Fit)
model_choice_r2 <- if (Pvalue_Linear < 0.05 && Pvalue_Sigmoid > 0.05) {
  "Linear"
} else if (Pvalue_Linear > 0.05 && Pvalue_Sigmoid < 0.05) {
  if (sigmoid_valid) "Sigmoid" else "No significant relationship detected"
} else if (Pvalue_Linear < 0.05 && Pvalue_Sigmoid < 0.05) {
  if (!sigmoid_valid) {
    "Linear"
  } else if (adj_r2_lm >= adj_r2_4pl) {
    "Linear"
  } else {
    "Sigmoid"
  }
} else {
  "No significant relationship detected"
}

# Convert short choices to full output labels
map_model_label <- function(choice) {
  if (choice == "Linear")  return("Linear_BoxCox")
  if (choice == "Sigmoid") return("Four parameter Sigmoidal")
  return("No significant relationship detected")
}

model_aic      <- map_model_label(model_choice_aic)
model_r2       <- map_model_label(model_choice_r2)
Methods_Differ <- (model_aic != model_r2)

# Print comparison summary to console
cat("\n============================================================\n")
cat("MODEL COMPARISON SUMMARY FOR:", response_name, "\n")
cat("============================================================\n")
cat(sprintf(
  "Linear (Box-Cox, lambda = %.2f): Adj R2 = %.4f | Raw-scale AIC = %.2f (AICc = %.2f)\n",
  lambda, adj_r2_lm, AIC_lm_raw, AICc_lm_raw
))
cat(sprintf(
  "4PL Sigmoidal (Raw scale)      : Adj R2 = %.4f | Raw-scale AIC = %.2f (AICc = %.2f)\n",
  adj_r2_4pl, AIC_4pl_raw, AICc_4pl_raw
))
cat("------------------------------------------------------------\n")
cat("Selected by Jacobian AIC :", model_aic, "\n")
cat("Selected by Adjusted R2  :", model_r2, "\n")
cat("Do the methods differ?   :", ifelse(Methods_Differ, "YES - Model choice changed!", "NO - Both agree."), "\n")
cat("============================================================\n\n")

# Assign primary model and minimum threshold based on Jacobian AIC selection
model <- model_aic
Minimum_threshold <- if (model_choice_aic == "Linear") {
  threshold_mass_g
} else if (model_choice_aic == "Sigmoid") {
  threshold_mass
} else {
  NA_real_
}


# ==============================================================================
# Step 6: Compile, Export, and Plot Final Results
# ==============================================================================

result_dataframe <- data.frame(
  Response_variable      = response_name,
  Year                   = year,
  Model_Chosen           = model,
  Model_Chosen_AdjR2     = model_r2,
  Methods_Differ         = Methods_Differ,
  Lambda_BoxCox          = lambda,
  AIC_Linear_RawScale    = AIC_lm_raw,
  AIC_Sigmoid_RawScale   = AIC_4pl_raw,
  AICc_Linear_RawScale   = AICc_lm_raw,
  AICc_Sigmoid_RawScale  = AICc_4pl_raw,
  Delta_AIC_4PL_minus_LM = Delta_AIC_4pl_vs_lm,
  Rsquared_Linear        = adj_r2_lm,
  Pvalue_Linear          = Pvalue_Linear,
  Rsquared_Sigmoid       = adj_r2_4pl,
  Pvalue_Sigmoid         = Pvalue_Sigmoid,
  N                      = summary(lm_bc)$df[2] + summary(lm_bc)$df[1],
  R                      = r,
  R_error                = r_error,
  Slope_Z                = b,
  Slope_SE_Z             = se,
  Slope_CI95_Lower_Z     = b - 1.96 * se,
  Slope_CI95_Upper_Z     = b + 1.96 * se,
  EC50                   = EC50,
  EC50_SE                = EC50_SE,
  EC50_95CI_Lower        = EC50_95CI_Lower,
  EC50_95CI_Upper        = EC50_95CI_Upper,
  Minimum_threshold      = Minimum_threshold,
  EC95                   = EC95,
  EC90                   = EC90,
  ECMAX                  = ECMAX
)

# Print final results table
print(result_dataframe)

# Save analysis objects to global environment under endpoint-specific names
assign(response_variable, data_frame)
assign(paste0(response_variable, "_4pl_plot"), four_pl_plot)
assign(paste0(response_variable, "_lm_plot"), lm_plot)
assign(paste0("results_", response_variable), result_dataframe)

# Export results table to CSV
file_path <- file.path(output_dir, paste0(response_name, ".csv"))
write.csv(result_dataframe, file_path, row.names = FALSE)

# Render final plots
print(lm_plot)
print(four_pl_plot)