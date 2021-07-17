# devtools::install_github("James-Thorson-NOAA/VAST")
library(VAST)
library(sdmTMB)

# load data set
# see `?load_example` for list of stocks with example data
# that are installed automatically with `FishStatsUtils`.


# Make settings (turning off bias.correct to save time for example)
settings <- make_settings(
  n_x = 100,
  Region = example$Region,
  purpose = "index2",
  strata.limits = example$strata.limits,
  bias.correct = FALSE
)

# Run model
fit <- fit_model(
  settings = settings,
  Lat_i = example$sampling_data[, "Lat"],
  Lon_i = example$sampling_data[, "Lon"],
  t_i = example$sampling_data[, "Year"],
  c_i = rep(0, nrow(example$sampling_data)),
  b_i = example$sampling_data[, "Catch_KG"],
  a_i = example$sampling_data[, "AreaSwept_km2"],
  v_i = example$sampling_data[, "Vessel"]
)

# Plot results
plot(fit)

example <- load_example(data_set = "BC_pacific_cod")

# Make settings (turning off bias.correct to save time for example)
settings <- make_settings(
  n_x = 100,
  Region = example$Region,
  purpose = "index2",
  fine_scale = TRUE,
  # RhoConfig = c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0),
  strata.limits = example$strata.limits,
  bias.correct = FALSE,
  knot_method="samples"
)

# Run model
fit <- fit_model(
  settings = settings,
  Lat_i = example$sampling_data[, "Lat"],
  Lon_i = example$sampling_data[, "Lon"],
  t_i = example$sampling_data[, "Year"],
  c_i = rep(0, nrow(example$sampling_data)),
  b_i = example$sampling_data[, "Catch_KG"],
  a_i = example$sampling_data[, "AreaSwept_km2"],
  v_i = example$sampling_data[, "Vessel"]
)

# Plot results
plot(fit)
