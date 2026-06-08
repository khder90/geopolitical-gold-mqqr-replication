## ============================================================
## Replication script
## Geopolitical Risk and Gold Return Predictability across Quantile States
## Quantile-on-Quantile Regression with Block Bootstrap and Scenario Forecasts
## ============================================================

## ------------------------------------------------------------
## 0. Project setup
## ------------------------------------------------------------

rm(list = ls())

project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir_input <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")

main_dir  <- file.path(output_dir, "01_Manuscript_Tables_Figures")
app_dir   <- file.path(output_dir, "02_Appendix_Tables_Figures")
model_dir <- file.path(output_dir, "03_Model_Objects")
audit_dir <- file.path(output_dir, "04_Data_Audit")

for (d in c(output_dir, main_dir, app_dir, model_dir, audit_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

required_pkgs <- c(
  "readxl", "openxlsx", "dplyr", "tidyr", "stringr", "lubridate", "zoo",
  "ggplot2", "scales", "patchwork", "ragg", "plot3D",
  "quantreg", "MASS", "mqqr", "QuantileOnQuantile"
)

missing_pkgs <- required_pkgs[!required_pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE)
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

set.seed(2026)

## ------------------------------------------------------------
## 1. General helper functions
## ------------------------------------------------------------

parse_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "M", "-")
  x <- ifelse(stringr::str_detect(x, "^\\d{4}-\\d{1,2}$"), paste0(x, "-01"), x)
  as.Date(lubridate::parse_date_time(x, orders = c("ymd", "ym", "dmy", "mdy")))
}

safe_min_ym <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  as.character(min(x, na.rm = TRUE))
}

safe_max_ym <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_character_)
  as.character(max(x, na.rm = TRUE))
}

theme_pub <- function(base_size = 16) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size - 1),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 2, color = "grey15"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 2),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.30, color = "grey86"),
      plot.caption = ggplot2::element_text(size = base_size - 3, color = "grey30")
    )
}

save_pub_png <- function(plot_obj, filename, width = 14, height = 8, dpi = 600) {
  ragg::agg_png(filename = filename, width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot_obj)
  dev.off()
}

matrix_to_long_safe <- function(mat, value_name) {
  mat <- as.matrix(mat)
  if (is.null(rownames(mat))) rownames(mat) <- seq(0.05, 0.95, by = 0.05)
  if (is.null(colnames(mat))) colnames(mat) <- seq(0.05, 0.95, by = 0.05)
  out <- as.data.frame(as.table(mat))
  colnames(out) <- c("theta", "tau", value_name)
  out$theta <- as.numeric(as.character(out$theta))
  out$tau <- as.numeric(as.character(out$tau))
  out
}

plot_heatmap_pub <- function(mat, title, fill_title, diverging = TRUE, base_size = 16) {
  d <- matrix_to_long_safe(mat, "Value")
  p <- ggplot2::ggplot(d, ggplot2::aes(x = tau, y = theta, fill = Value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_x_continuous(breaks = seq(0.05, 0.95, 0.10)) +
    ggplot2::scale_y_continuous(breaks = seq(0.05, 0.95, 0.10)) +
    ggplot2::labs(title = title, x = "Lagged GPR quantile", y = "Gold return quantile", fill = fill_title) +
    theme_pub(base_size)
  if (diverging) {
    p + ggplot2::scale_fill_gradient2(low = "#1F4E79", mid = "white", high = "#B22222", midpoint = 0,
                                      labels = scales::label_number(accuracy = 0.001))
  } else {
    p + ggplot2::scale_fill_viridis_c(option = "viridis", labels = scales::label_number(accuracy = 0.001))
  }
}

plot_3d_surface_pub <- function(M, main_title, z_title, file_name, zlim = NULL,
                                color_type = c("diverging", "sequential"),
                                theta_angle = 40, phi_angle = 25, width = 11,
                                height = 8, dpi = 600) {
  color_type <- match.arg(color_type)
  M <- as.matrix(M)
  x <- as.numeric(colnames(M))
  y <- as.numeric(rownames(M))
  z <- M
  if (is.null(zlim)) zlim <- range(z, na.rm = TRUE)
  col_pal <- if (color_type == "diverging") grDevices::hcl.colors(120, "Blue-Red 3") else grDevices::hcl.colors(120, "Viridis")
  ragg::agg_png(file.path(main_dir, file_name), width = width, height = height, units = "in", res = dpi, background = "white")
  par(mar = c(3.8, 3.8, 4.2, 5.2), cex.axis = 1.15, cex.lab = 1.25, cex.main = 1.45, font.main = 2)
  plot3D::persp3D(
    x = x, y = y, z = z, theta = theta_angle, phi = phi_angle,
    expand = 0.65, col = col_pal, border = "grey35", facets = TRUE,
    ticktype = "detailed", bty = "b2", shade = 0.25,
    xlab = "Lagged GPR quantile", ylab = "Gold return quantile", zlab = z_title,
    main = main_title, zlim = zlim,
    colkey = list(side = 4, length = 0.55, width = 0.85, cex.axis = 1.0, clab = z_title)
  )
  dev.off()
}

## ------------------------------------------------------------
## 2. Load and prepare data
## ------------------------------------------------------------

xlsx_file <- file.path(data_dir_input, "rd.xlsx")
csv_file  <- file.path(data_dir_input, "rd.csv")

if (file.exists(xlsx_file)) {
  df_raw <- readxl::read_excel(xlsx_file, sheet = "data")
} else if (file.exists(csv_file)) {
  df_raw <- read.csv(csv_file, check.names = FALSE)
} else {
  stop("No input data found. Place data/rd.xlsx with sheet 'data' or data/rd.csv in the data folder.")
}

df <- df_raw %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, "\\s+", " ")) %>%
  dplyr::rename_with(~ stringr::str_trim(.x))

required_vars <- c("date", "Gold Price", "Gold Return", "GPR", "Dollar Index", "VIX")
missing_vars <- setdiff(required_vars, names(df))
if (length(missing_vars) > 0) stop(paste("Missing required columns:", paste(missing_vars, collapse = ", ")))

df <- df %>%
  dplyr::mutate(
    date = parse_date_safe(date),
    ym = zoo::as.yearmon(date),
    Gold_Price  = as.numeric(`Gold Price`),
    Gold_Return = as.numeric(`Gold Return`),
    GPR = as.numeric(GPR),
    DXY = as.numeric(`Dollar Index`),
    VIX = as.numeric(VIX)
  ) %>%
  dplyr::arrange(ym) %>%
  dplyr::select(date, ym, Gold_Price, Gold_Return, GPR, DXY, VIX)

df_lagged <- df %>%
  dplyr::arrange(ym) %>%
  dplyr::mutate(
    GPR_l1 = dplyr::lag(GPR, 1),
    DXY_l1 = dplyr::lag(DXY, 1),
    VIX_l1 = dplyr::lag(VIX, 1),
    Gold_Return_l1 = dplyr::lag(Gold_Return, 1)
  ) %>%
  tidyr::drop_na(Gold_Return, GPR_l1, DXY_l1, VIX_l1)

est_end_ym <- zoo::as.yearmon("2026-01")
df_est <- df_lagged %>% dplyr::filter(ym <= est_end_ym)
df_available <- df_lagged

sample_audit <- data.frame(
  Item = c("Full raw sample start", "Full raw sample end", "Lagged sample start", "Lagged sample end",
           "Estimation sample start", "Estimation sample end", "Number of estimation observations"),
  Value = as.character(c(safe_min_ym(df$ym), safe_max_ym(df$ym), safe_min_ym(df_lagged$ym), safe_max_ym(df_lagged$ym),
                         safe_min_ym(df_est$ym), safe_max_ym(df_est$ym), nrow(df_est))),
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(list(Full_Data = df, Lagged_Data = df_lagged, Estimation_Sample = df_est, Sample_Audit = sample_audit),
                     file = file.path(audit_dir, "Data_Audit_and_Samples.xlsx"), overwrite = TRUE)

## ------------------------------------------------------------
## 3. Descriptive statistics and preliminary figures
## ------------------------------------------------------------

desc_vars <- df_est %>% dplyr::select(Gold_Price, Gold_Return, GPR, DXY, VIX)

desc_stats <- data.frame(
  Variable = names(desc_vars),
  N = sapply(desc_vars, function(x) sum(!is.na(x))),
  Mean = sapply(desc_vars, mean, na.rm = TRUE),
  SD = sapply(desc_vars, sd, na.rm = TRUE),
  Median = sapply(desc_vars, median, na.rm = TRUE),
  Min = sapply(desc_vars, min, na.rm = TRUE),
  Max = sapply(desc_vars, max, na.rm = TRUE),
  Skewness = sapply(desc_vars, function(x) mean((x - mean(x, na.rm = TRUE))^3, na.rm = TRUE) / sd(x, na.rm = TRUE)^3),
  Kurtosis = sapply(desc_vars, function(x) mean((x - mean(x, na.rm = TRUE))^4, na.rm = TRUE) / sd(x, na.rm = TRUE)^4 - 3),
  row.names = NULL
)

openxlsx::write.xlsx(desc_stats, file.path(main_dir, "Table_Descriptive_Statistics.xlsx"), overwrite = TRUE)

plot_df <- df_est %>%
  dplyr::select(date, Gold_Return, GPR, DXY, VIX) %>%
  tidyr::pivot_longer(-date, names_to = "Variable", values_to = "Value")

p_time <- ggplot2::ggplot(plot_df, ggplot2::aes(x = date, y = Value)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  ggplot2::labs(title = "Time Evolution of Study Variables", x = "Month", y = "Value") +
  theme_pub(16)
save_pub_png(p_time, file.path(main_dir, "Figure_1_Time_Evolution_600dpi.png"), width = 13, height = 11)

p_dist <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Value)) +
  ggplot2::geom_histogram(bins = 35, color = "white") +
  ggplot2::facet_wrap(~ Variable, scales = "free", ncol = 2) +
  ggplot2::labs(title = "Empirical Distributions of Study Variables", x = "Value", y = "Frequency") +
  theme_pub(16)
save_pub_png(p_dist, file.path(main_dir, "Figure_2_Empirical_Distribution_600dpi.png"), width = 13, height = 9)

## ------------------------------------------------------------
## 4. Baseline lagged quantile regression
## ------------------------------------------------------------

qr_taus <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
qr_results <- lapply(qr_taus, function(tau) {
  fit <- quantreg::rq(Gold_Return ~ GPR_l1 + DXY_l1 + VIX_l1, tau = tau, data = df_est, method = "fn")
  s <- summary(fit, se = "nid")
  out <- as.data.frame(s$coefficients)
  out$Variable <- rownames(out)
  names(out)[1:4] <- c("Coefficient", "Std_Error", "t_value", "p_value")
  out$Quantile <- tau
  out[, c("Quantile", "Variable", "Coefficient", "Std_Error", "t_value", "p_value")]
}) %>% dplyr::bind_rows()

openxlsx::write.xlsx(qr_results, file.path(main_dir, "Table_Baseline_Lagged_Quantile_Regression.xlsx"), overwrite = TRUE)

## ------------------------------------------------------------
## 5. Main mQQR model
## ------------------------------------------------------------

yq <- seq(0.05, 0.95, by = 0.05)
xq <- seq(0.05, 0.95, by = 0.05)

fit_mqqr_gpr <- mqqr::mqq_regression(
  y = df_est$Gold_Return,
  x = df_est$GPR_l1,
  moderators = list(DXY_l1 = df_est$DXY_l1, VIX_l1 = df_est$VIX_l1),
  y_quantiles = yq,
  x_quantiles = xq,
  bandwidth = 0.05,
  include_lag = FALSE,
  interactions = FALSE,
  se = "bootstrap",
  n_boot = 499,
  cdf_based_kernel = TRUE,
  x_name = "GPR_lag1",
  y_name = "Gold_Return",
  verbose = TRUE,
  seed = 2026
)

saveRDS(fit_mqqr_gpr, file.path(model_dir, "Main_mQQR_GPR_lag1_Model.rds"))

M_beta <- mqqr::mqq_to_matrix(fit_mqqr_gpr, value = "beta1")
M_se   <- mqqr::mqq_to_matrix(fit_mqqr_gpr, value = "se")
M_t    <- mqqr::mqq_to_matrix(fit_mqqr_gpr, value = "t_value")
M_p    <- mqqr::mqq_to_matrix(fit_mqqr_gpr, value = "p_value")
M_r2   <- mqqr::mqq_to_matrix(fit_mqqr_gpr, value = "r_squared")

res_main <- matrix_to_long_safe(M_beta, "beta1") %>%
  dplyr::left_join(matrix_to_long_safe(M_se, "se"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(M_t, "t_value"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(M_p, "p_value"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(M_r2, "r_squared"), by = c("theta", "tau")) %>%
  dplyr::mutate(abs_t = abs(t_value), Significant_5pct = p_value < 0.05, Significant_10pct = p_value < 0.10)

main_stats <- data.frame(
  Statistic = c("Mean beta1", "Median beta1", "Min beta1", "Max beta1", "SD beta1", "Mean R-squared", "Significant cells p<0.05", "Total cells"),
  Value = c(mean(res_main$beta1, na.rm = TRUE), median(res_main$beta1, na.rm = TRUE), min(res_main$beta1, na.rm = TRUE),
            max(res_main$beta1, na.rm = TRUE), sd(res_main$beta1, na.rm = TRUE), mean(res_main$r_squared, na.rm = TRUE),
            sum(res_main$Significant_5pct, na.rm = TRUE), nrow(res_main))
)

gpr_sd <- sd(df_est$GPR_l1, na.rm = TRUE)
res_main <- res_main %>%
  dplyr::mutate(One_SD_GPR_Effect_Percent = beta1 * gpr_sd,
                One_SD_GPR_Effect_Basis_Points = One_SD_GPR_Effect_Percent * 100)

main_surface_summary <- res_main %>%
  dplyr::summarise(Number_of_cells = dplyr::n(), Significant_share_5pct = mean(Significant_5pct, na.rm = TRUE),
                   Significant_share_10pct = mean(Significant_10pct, na.rm = TRUE), Mean_beta = mean(beta1, na.rm = TRUE),
                   Median_beta = median(beta1, na.rm = TRUE), Mean_R_squared = mean(r_squared, na.rm = TRUE),
                   Median_R_squared = median(r_squared, na.rm = TRUE), Max_abs_t = max(abs_t, na.rm = TRUE), Min_p_value = min(p_value, na.rm = TRUE))

top_main_cells <- res_main %>% dplyr::arrange(p_value, dplyr::desc(abs_t)) %>% dplyr::slice_head(n = 30)

openxlsx::write.xlsx(main_stats, file.path(main_dir, "Table_mQQR_Main_Model_Statistics.xlsx"), overwrite = TRUE)
openxlsx::write.xlsx(list(Surface_Summary = main_surface_summary, Top_Cells = top_main_cells, Full_Cell_Results = res_main),
                     file.path(main_dir, "Table_mQQR_Top_Cells_and_Economic_Magnitude.xlsx"), overwrite = TRUE)

p_beta <- plot_heatmap_pub(M_beta, "A. Main mQQR Coefficient Surface", "Coefficient", TRUE)
p_p    <- plot_heatmap_pub(M_p, "B. p-value Surface", "p-value", FALSE)
p_r2   <- plot_heatmap_pub(M_r2, "C. Pseudo R-squared Surface", "Pseudo R²", FALSE)
save_pub_png(p_beta | p_p | p_r2, file.path(main_dir, "Figure_3_mQQR_Main_Surface_Heatmaps_600dpi.png"), width = 21, height = 7.5)

plot_3d_surface_pub(M_beta, "3D Surface of Lagged GPR Effect on Gold Returns", "GPR effect", "Figure_3D_1_mQQR_Coefficient_Surface_600dpi.png", color_type = "diverging")
plot_3d_surface_pub(M_t, "3D mQQR t-statistic Surface", "t-statistic", "Figure_3D_2_mQQR_t_Statistic_Surface_600dpi.png", color_type = "diverging")
plot_3d_surface_pub(M_p, "3D mQQR p-value Surface", "p-value", "Figure_3D_3_mQQR_p_Value_Surface_600dpi.png", zlim = c(0, 1), color_type = "sequential")
plot_3d_surface_pub(M_r2, "3D mQQR Pseudo R-squared Surface", "Pseudo R²", "Figure_3D_4_mQQR_Pseudo_R2_Surface_600dpi.png", color_type = "sequential")

## ------------------------------------------------------------
## 6. Interactive mQQR model
## ------------------------------------------------------------

df_int <- df_est %>% tidyr::drop_na(Gold_Return, Gold_Return_l1, GPR_l1, DXY_l1, VIX_l1)

fit_mqqr_interactive <- mqqr::mqq_regression(
  y = df_int$Gold_Return,
  x = df_int$GPR_l1,
  moderators = list(Gold_Return_l1 = df_int$Gold_Return_l1, DXY_l1 = df_int$DXY_l1, VIX_l1 = df_int$VIX_l1),
  y_quantiles = yq,
  x_quantiles = xq,
  bandwidth = 0.05,
  include_lag = FALSE,
  interactions = TRUE,
  se = "bootstrap",
  n_boot = 499,
  cdf_based_kernel = TRUE,
  x_name = "GPR_lag1",
  y_name = "Gold_Return",
  verbose = TRUE,
  seed = 2026
)

saveRDS(fit_mqqr_interactive, file.path(model_dir, "Main_Interactive_mQQR_Model.rds"))

MI_beta <- mqqr::mqq_to_matrix(fit_mqqr_interactive, value = "beta1")
MI_se   <- mqqr::mqq_to_matrix(fit_mqqr_interactive, value = "se")
MI_t    <- mqqr::mqq_to_matrix(fit_mqqr_interactive, value = "t_value")
MI_p    <- mqqr::mqq_to_matrix(fit_mqqr_interactive, value = "p_value")
MI_r2   <- mqqr::mqq_to_matrix(fit_mqqr_interactive, value = "r_squared")

res_int <- matrix_to_long_safe(MI_beta, "beta_gpr") %>%
  dplyr::left_join(matrix_to_long_safe(MI_se, "se"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(MI_t, "t_value"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(MI_p, "p_value"), by = c("theta", "tau")) %>%
  dplyr::left_join(matrix_to_long_safe(MI_r2, "r_squared"), by = c("theta", "tau")) %>%
  dplyr::mutate(abs_t = abs(t_value), Significant_5pct = p_value < 0.05, Significant_10pct = p_value < 0.10)

int_summary <- res_int %>%
  dplyr::summarise(Number_of_cells = dplyr::n(), Significant_share_5pct = mean(Significant_5pct, na.rm = TRUE),
                   Significant_share_10pct = mean(Significant_10pct, na.rm = TRUE), Mean_beta_GPR = mean(beta_gpr, na.rm = TRUE),
                   Median_beta_GPR = median(beta_gpr, na.rm = TRUE), Mean_R_squared = mean(r_squared, na.rm = TRUE),
                   Median_R_squared = median(r_squared, na.rm = TRUE), Max_abs_t = max(abs_t, na.rm = TRUE), Min_p_value = min(p_value, na.rm = TRUE))

top_int_cells <- res_int %>% dplyr::arrange(p_value, dplyr::desc(abs_t)) %>% dplyr::slice_head(n = 30)

openxlsx::write.xlsx(list(Interactive_Model_Summary = int_summary, Top_30_Significant_Cells = top_int_cells, Full_Cell_Results = res_int),
                     file.path(main_dir, "Table_Main_Supporting_Interactive_mQQR_Results.xlsx"), overwrite = TRUE)

p_i_beta <- plot_heatmap_pub(MI_beta, "A. Interactive mQQR GPR Effect Surface", "GPR effect", TRUE)
p_i_p    <- plot_heatmap_pub(MI_p, "B. Interactive mQQR p-value Surface", "p-value", FALSE)
p_i_r2   <- plot_heatmap_pub(MI_r2, "C. Interactive mQQR Pseudo R-squared Surface", "Pseudo R²", FALSE)
save_pub_png(p_i_beta | p_i_p | p_i_r2, file.path(main_dir, "Figure_Main_Supporting_mQQR_Core_Surfaces_600dpi.png"), width = 21, height = 7.5)

extract_optional_surface <- function(fit, value, moderator_name) {
  tryCatch(mqqr::mqq_to_matrix(fit, value = value, moderator = moderator_name), error = function(e) NULL)
}

M_alpha_dxy <- extract_optional_surface(fit_mqqr_interactive, "alpha", "DXY_l1")
M_alpha_vix <- extract_optional_surface(fit_mqqr_interactive, "alpha", "VIX_l1")
if (!is.null(M_alpha_dxy) && !is.null(M_alpha_vix)) {
  p_ad <- plot_heatmap_pub(M_alpha_dxy, "A. Interaction Surface: GPR x DXY", "Interaction", TRUE)
  p_av <- plot_heatmap_pub(M_alpha_vix, "B. Interaction Surface: GPR x VIX", "Interaction", TRUE)
  save_pub_png(p_ad | p_av, file.path(main_dir, "Figure_Main_Supporting_mQQR_Interaction_Surfaces_600dpi.png"), width = 15, height = 7.5)
}

## ------------------------------------------------------------
## 7. Real-time and scenario forecasts
## ------------------------------------------------------------

emp_cdf_value <- function(x_hist, x_new) mean(x_hist <= x_new, na.rm = TRUE)
gaussian_kernel_weights <- function(x_hist, x_new, bandwidth = 0.05) {
  Fx <- ecdf(x_hist)
  u_hist <- Fx(x_hist)
  u_new <- emp_cdf_value(x_hist, x_new)
  w <- dnorm((u_hist - u_new) / bandwidth)
  if (sum(w, na.rm = TRUE) <= 0 || all(is.na(w))) w <- rep(1, length(x_hist))
  w
}
predict_local_mqqr <- function(train_data, new_data, theta, bandwidth = 0.05) {
  w <- gaussian_kernel_weights(train_data$GPR_l1, new_data$GPR_l1, bandwidth)
  fit <- quantreg::rq(Gold_Return ~ GPR_l1 + DXY_l1 + VIX_l1, tau = theta, data = train_data, weights = w, method = "fn")
  data.frame(theta = theta, qhat = as.numeric(predict(fit, newdata = new_data)), tau_hat = emp_cdf_value(train_data$GPR_l1, new_data$GPR_l1))
}

theta_forecast <- c(0.05, 0.10, 0.50, 0.90, 0.95)
last_obs <- df_est %>% dplyr::arrange(date) %>% dplyr::slice_tail(n = 1)
last_predictors <- data.frame(GPR_l1 = as.numeric(last_obs$GPR), DXY_l1 = as.numeric(last_obs$DXY), VIX_l1 = as.numeric(last_obs$VIX))
if (any(is.na(last_predictors))) last_predictors <- data.frame(GPR_l1 = last_obs$GPR_l1, DXY_l1 = last_obs$DXY_l1, VIX_l1 = last_obs$VIX_l1)

forecast_2026M02 <- lapply(theta_forecast, function(th) predict_local_mqqr(df_est, last_predictors, th, 0.05)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(date = as.Date("2026-02-01"), Scenario = "Observed predictors in 2026M01",
                GPR_l1 = last_predictors$GPR_l1, DXY_l1 = last_predictors$DXY_l1, VIX_l1 = last_predictors$VIX_l1)

scenario_probs <- c(Calm = 0.25, Baseline = 0.50, Tension = 0.75, Severe = 0.90)
scenario_levels <- data.frame(
  Scenario = names(scenario_probs), Probability = as.numeric(scenario_probs),
  GPR_l1 = as.numeric(stats::quantile(df_est$GPR_l1, probs = scenario_probs, na.rm = TRUE)),
  DXY_l1 = as.numeric(stats::quantile(df_est$DXY_l1, probs = scenario_probs, na.rm = TRUE)),
  VIX_l1 = as.numeric(stats::quantile(df_est$VIX_l1, probs = scenario_probs, na.rm = TRUE))
)

future_dates <- seq(as.Date("2026-03-01"), as.Date("2027-12-01"), by = "month")
scenario_rows <- list()
for (j in seq_along(future_dates)) {
  d <- future_dates[j]
  for (i in seq_len(nrow(scenario_levels))) {
    new_i <- data.frame(GPR_l1 = scenario_levels$GPR_l1[i], DXY_l1 = scenario_levels$DXY_l1[i], VIX_l1 = scenario_levels$VIX_l1[i])
    for (th in theta_forecast) {
      scenario_rows[[length(scenario_rows) + 1]] <- predict_local_mqqr(df_est, new_i, th, 0.05) %>%
        dplyr::mutate(date = as.Date(d), Scenario = scenario_levels$Scenario[i], GPR_l1 = new_i$GPR_l1, DXY_l1 = new_i$DXY_l1, VIX_l1 = new_i$VIX_l1)
    }
  }
}
scenario_forecasts <- dplyr::bind_rows(scenario_rows) %>% dplyr::mutate(date = as.Date(date))
all_forecasts <- dplyr::bind_rows(forecast_2026M02, scenario_forecasts) %>% dplyr::mutate(date = as.Date(date))
forecast_wide <- all_forecasts %>%
  dplyr::mutate(theta_label = paste0("Q", sprintf("%02d", round(theta * 100)))) %>%
  dplyr::select(date, Scenario, theta_label, qhat) %>%
  tidyr::pivot_wider(names_from = theta_label, values_from = qhat) %>%
  dplyr::arrange(date, Scenario)

forecast_design <- data.frame(
  Item = c("Available sample", "Estimation sample", "Predictive specification", "One-step-ahead forecast", "Scenario forecast horizon",
           "Scenario construction", "Forecast quantiles", "Ex-post accuracy evaluation"),
  Description = c("Monthly data ending at the last available observation", "Lagged model sample ending at the last available observation",
                  "Gold_Return_t is modeled using GPR_t-1, DXY_t-1, and VIX_t-1", "One-step-ahead forecast uses the latest observed predictors",
                  "2026M03 to 2027M12", "Calm, Baseline, Tension, and Severe scenarios use the 25th, 50th, 75th, and 90th percentiles of historical predictors",
                  "Q05, Q10, Q50, Q90, and Q95", "Not reported unless realized observations after the estimation end are available")
)

openxlsx::write.xlsx(list(Forecast_Design = forecast_design, Scenario_Levels = scenario_levels, Forecasts_Long = all_forecasts, Forecasts_Wide = forecast_wide),
                     file.path(main_dir, "Table_Real_Time_and_Scenario_Forecasts_2026M02_2027M12.xlsx"), overwrite = TRUE)

p_real <- forecast_wide %>% dplyr::filter(Scenario == "Observed predictors in 2026M01") %>%
  ggplot2::ggplot(ggplot2::aes(x = date)) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Q05, ymax = Q95), width = 12, color = "#1F4E79", linewidth = 1.2) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = Q10, ymax = Q90), width = 7, color = "#2E6F40", linewidth = 1.4) +
  ggplot2::geom_point(ggplot2::aes(y = Q50), color = "#B22222", size = 5) +
  ggplot2::labs(title = "One-step-ahead Forecast of Gold Returns", x = "Month", y = "Gold return, percent") + theme_pub(16)
save_pub_png(p_real, file.path(main_dir, "Figure_5_One_Step_Ahead_Forecast_2026M02_600dpi.png"), width = 11, height = 7)

p_scenario <- forecast_wide %>% dplyr::filter(Scenario != "Observed predictors in 2026M01") %>%
  ggplot2::ggplot(ggplot2::aes(x = date)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Q05, ymax = Q95, fill = Scenario), alpha = 0.18) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = Q10, ymax = Q90, fill = Scenario), alpha = 0.30) +
  ggplot2::geom_line(ggplot2::aes(y = Q50, color = Scenario), linewidth = 1.25) +
  ggplot2::facet_wrap(~ Scenario, ncol = 2) +
  ggplot2::labs(title = "Scenario-based Quantile Forecasts of Gold Returns", x = "Month", y = "Gold return, percent") +
  theme_pub(16) + ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
save_pub_png(p_scenario, file.path(main_dir, "Figure_6_Scenario_Forecasts_2026M03_2027M12_600dpi.png"), width = 16, height = 10)

## ------------------------------------------------------------
## 8. Appendix robustness: bivariate lagged QQR
## ------------------------------------------------------------

biv_dir <- file.path(app_dir, "Robustness_Bivariate_Lagged_QQR")
dir.create(biv_dir, showWarnings = FALSE, recursive = TRUE)

biv_results <- list()
biv_stats <- list()
for (x_name in c("GPR_l1", "DXY_l1", "VIX_l1")) {
  df_b <- df_est %>% dplyr::select(Gold_Return, dplyr::all_of(x_name)) %>% tidyr::drop_na()
  fit_b <- QuantileOnQuantile::qq_regression(y = df_b$Gold_Return, x = df_b[[x_name]], y_quantiles = yq, x_quantiles = xq,
                                             min_obs = 15, se_method = "boot", verbose = TRUE)
  saveRDS(fit_b, file.path(model_dir, paste0("Appendix_Bivariate_QQR_", x_name, ".rds")))
  M_coef <- QuantileOnQuantile::qq_to_matrix(fit_b, type = "coefficient")
  M_rsq  <- QuantileOnQuantile::qq_to_matrix(fit_b, type = "rsquared")
  M_pv   <- QuantileOnQuantile::qq_to_matrix(fit_b, type = "pvalue")
  biv_results[[x_name]] <- matrix_to_long_safe(M_coef, "coefficient") %>%
    dplyr::left_join(matrix_to_long_safe(M_rsq, "r_squared"), by = c("theta", "tau")) %>%
    dplyr::left_join(matrix_to_long_safe(M_pv, "p_value"), by = c("theta", "tau")) %>%
    dplyr::mutate(Predictor = x_name)
  biv_stats[[x_name]] <- data.frame(Predictor = x_name, Mean_Coefficient = mean(as.numeric(M_coef), na.rm = TRUE),
                                    Median_Coefficient = median(as.numeric(M_coef), na.rm = TRUE), Min_Coefficient = min(as.numeric(M_coef), na.rm = TRUE),
                                    Max_Coefficient = max(as.numeric(M_coef), na.rm = TRUE), SD_Coefficient = sd(as.numeric(M_coef), na.rm = TRUE),
                                    Mean_R2 = mean(as.numeric(M_rsq), na.rm = TRUE), Median_R2 = median(as.numeric(M_rsq), na.rm = TRUE),
                                    Max_R2 = max(as.numeric(M_rsq), na.rm = TRUE), Significant_Cells_5pct = sum(as.numeric(M_pv) < 0.05, na.rm = TRUE),
                                    Total_Cells = length(as.numeric(M_coef)))
}
openxlsx::write.xlsx(list(Bivariate_QQR_Statistics = dplyr::bind_rows(biv_stats), Bivariate_QQR_Cell_Results = dplyr::bind_rows(biv_results)),
                     file.path(biv_dir, "Appendix_Bivariate_Lagged_QQR_All_Results.xlsx"), overwrite = TRUE)

## ------------------------------------------------------------
## 9. Appendix robustness: bandwidth sensitivity
## ------------------------------------------------------------

bw_dir <- file.path(app_dir, "Robustness_mQQR_Bandwidth_Sensitivity")
dir.create(bw_dir, showWarnings = FALSE, recursive = TRUE)
bw_values <- c(0.05, 0.10, 0.15)
bw_summary_list <- list()

for (bw in bw_values) {
  fit_bw <- mqqr::mqq_regression(
    y = df_est$Gold_Return, x = df_est$GPR_l1,
    moderators = list(DXY_l1 = df_est$DXY_l1, VIX_l1 = df_est$VIX_l1),
    y_quantiles = yq, x_quantiles = xq, bandwidth = bw, include_lag = FALSE,
    interactions = FALSE, se = "bootstrap", n_boot = 199, cdf_based_kernel = TRUE,
    x_name = "GPR_lag1", y_name = "Gold_Return", verbose = TRUE, seed = 2026
  )
  saveRDS(fit_bw, file.path(model_dir, paste0("Robustness_mQQR_Bandwidth_", bw, ".rds")))
  B_beta <- mqqr::mqq_to_matrix(fit_bw, value = "beta1")
  B_se   <- mqqr::mqq_to_matrix(fit_bw, value = "se")
  B_t    <- mqqr::mqq_to_matrix(fit_bw, value = "t_value")
  B_p    <- mqqr::mqq_to_matrix(fit_bw, value = "p_value")
  B_r2   <- mqqr::mqq_to_matrix(fit_bw, value = "r_squared")
  res_bw <- matrix_to_long_safe(B_beta, "beta1") %>%
    dplyr::left_join(matrix_to_long_safe(B_se, "se"), by = c("theta", "tau")) %>%
    dplyr::left_join(matrix_to_long_safe(B_t, "t_value"), by = c("theta", "tau")) %>%
    dplyr::left_join(matrix_to_long_safe(B_p, "p_value"), by = c("theta", "tau")) %>%
    dplyr::left_join(matrix_to_long_safe(B_r2, "r_squared"), by = c("theta", "tau")) %>%
    dplyr::mutate(Bandwidth = bw)
  summary_bw <- res_bw %>% dplyr::summarise(Bandwidth = bw, Number_of_cells = dplyr::n(), Mean_Beta = mean(beta1, na.rm = TRUE),
                                            Median_Beta = median(beta1, na.rm = TRUE), Significant_Share_5pct = mean(p_value < 0.05, na.rm = TRUE),
                                            Significant_Share_10pct = mean(p_value < 0.10, na.rm = TRUE), Mean_R2 = mean(r_squared, na.rm = TRUE),
                                            Median_R2 = median(r_squared, na.rm = TRUE), Max_abs_t = max(abs(t_value), na.rm = TRUE), Min_p_value = min(p_value, na.rm = TRUE))
  bw_summary_list[[as.character(bw)]] <- summary_bw
  openxlsx::write.xlsx(list(Cell_Level_Results = res_bw, Summary = summary_bw), file.path(bw_dir, paste0("Bandwidth_", bw, "_Cell_Results.xlsx")), overwrite = TRUE)
}
openxlsx::write.xlsx(list(Bandwidth_Sensitivity_Summary = dplyr::bind_rows(bw_summary_list)),
                     file.path(bw_dir, "Appendix_mQQR_Bandwidth_Sensitivity_Summary.xlsx"), overwrite = TRUE)

## ------------------------------------------------------------
## 10. Appendix robustness: moving block bootstrap sensitivity
## ------------------------------------------------------------

mbb_dir <- file.path(app_dir, "Robustness_MBB_Block_Length_Sensitivity")
dir.create(mbb_dir, showWarnings = FALSE, recursive = TRUE)

mbb_index <- function(n, L) {
  starts <- sample.int(n, size = ceiling(n / L), replace = TRUE)
  idx <- unlist(lapply(starts, function(s) { j <- s:(s + L - 1); ((j - 1) %% n) + 1 }))
  idx[1:n]
}

run_mbb_symmetry_test <- function(data_in, taus, L = 12, B = 499, seed = 2026) {
  set.seed(seed)
  form <- Gold_Return ~ GPR_l1 + DXY_l1 + VIX_l1
  fit_list <- lapply(taus, function(tau) quantreg::rq(form, tau = tau, data = data_in, method = "fn"))
  coef_mat <- do.call(cbind, lapply(fit_list, coef))
  bhat <- as.vector(coef_mat)
  p <- nrow(coef_mat); n_tau <- length(taus)
  boot_mat <- matrix(NA_real_, nrow = B, ncol = length(bhat))
  for (b in seq_len(B)) {
    idx <- mbb_index(nrow(data_in), L)
    db <- data_in[idx, , drop = FALSE]
    boot_mat[b, ] <- tryCatch({
      fl <- lapply(taus, function(tau) quantreg::rq(form, tau = tau, data = db, method = "fn"))
      as.vector(do.call(cbind, lapply(fl, coef)))
    }, error = function(e) rep(NA_real_, length(bhat)))
  }
  boot_mat <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
  Vb <- stats::cov(boot_mat)
  pos <- function(k, j, p) (j - 1) * p + k
  sym_pairs <- c(0.10, 0.20, 0.30, 0.40)
  R_list <- list(); lab_list <- list()
  for (tL in sym_pairs) {
    tH <- 1 - tL; jL <- match(tL, taus); jH <- match(tH, taus)
    for (k in 1:p) {
      r <- rep(0, p * n_tau); r[pos(k, jL, p)] <- 1; r[pos(k, jH, p)] <- -1
      R_list[[length(R_list) + 1]] <- r
      lab_list[[length(lab_list) + 1]] <- data.frame(Block_Length = L, tau_low = tL, tau_high = tH, variable = rownames(coef_mat)[k])
    }
  }
  R <- do.call(rbind, R_list); lab <- do.call(rbind, lab_list)
  rhat <- as.vector(R %*% bhat); Vr <- R %*% Vb %*% t(R)
  inv_Vr <- tryCatch(solve(Vr), error = function(e) MASS::ginv(Vr))
  W <- as.numeric(t(rhat) %*% inv_Vr %*% rhat); df_w <- nrow(R); p_w <- 1 - stats::pchisq(W, df = df_w)
  se_r <- sqrt(diag(Vr)); z2 <- (rhat / se_r)^2; p_each <- 1 - stats::pchisq(z2, df = 1)
  list(summary = data.frame(Block_Length = L, Bootstrap_Requested = B, Bootstrap_Successful = nrow(boot_mat), Wald_Statistic = W, df = df_w, p_value = p_w),
       detail = cbind(lab, Restriction_Value = rhat, Std_Error = se_r, ChiSq_1df = z2, p_value = p_each))
}

taus_mbb <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)
mbb_6  <- run_mbb_symmetry_test(df_est, taus_mbb, L = 6,  B = 499, seed = 2026)
mbb_12 <- run_mbb_symmetry_test(df_est, taus_mbb, L = 12, B = 499, seed = 2026)
mbb_18 <- run_mbb_symmetry_test(df_est, taus_mbb, L = 18, B = 499, seed = 2026)
mbb_summary <- dplyr::bind_rows(mbb_6$summary, mbb_12$summary, mbb_18$summary)
mbb_detail <- dplyr::bind_rows(mbb_6$detail, mbb_12$detail, mbb_18$detail)
openxlsx::write.xlsx(list(MBB_Block_Length_Summary = mbb_summary, MBB_Restriction_Details = mbb_detail),
                     file.path(mbb_dir, "Appendix_MBB_Block_Length_Sensitivity.xlsx"), overwrite = TRUE)

## ------------------------------------------------------------
## 11. Output file index
## ------------------------------------------------------------

list_files_recursive <- function(path) {
  if (!dir.exists(path)) return(data.frame(File = character(0), Folder = character(0)))
  data.frame(File = list.files(path, recursive = TRUE, full.names = FALSE), Folder = path, stringsAsFactors = FALSE)
}
openxlsx::write.xlsx(list(Manuscript_Files = list_files_recursive(main_dir), Appendix_Files = list_files_recursive(app_dir),
                          Model_Objects = list_files_recursive(model_dir), Data_Audit = list_files_recursive(audit_dir)),
                     file.path(output_dir, "Output_File_Index.xlsx"), overwrite = TRUE)

cat("\nReplication completed successfully. Outputs saved in:\n", output_dir, "\n")
