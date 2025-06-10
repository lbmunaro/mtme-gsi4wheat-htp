rm(list = ls())

library(asreml)
library(tidyverse)

# Handle SLURM array task ID from commandArgs
args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])

k_vals <- 1:5
workspace_vals <- c("16gb", "24gb", "72gb", "80gb", "96gb")

k <- k_vals[task_id]
workspace <- workspace_vals[task_id]
model_name <- paste0("NFA", k, "_Conv")

message("Fitting ", model_name, " with workspace = ", workspace)

# Load data
ILYT_Pheno_HTP.Dsplit <- readRDS("Data/ILYT_Pheno_HTP.Dsplit.rds")
Ginv <- readRDS("Data/Ginv.rds")

data <- ILYT_Pheno_HTP.Dsplit$Conv |>
  droplevels()

str(data)

# Fit model
fit <- tryCatch({
  # Begin logging to file
  log_file <- paste0("Data/log_", model_name, ".txt")
  sink(file = log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  
  model <- asreml(
    fixed = Pheno_z ~ TraitEnv,
    random = as.formula(substitute(
      ~ rr(TraitEnv, k_val):vm(Gkeep, Ginv_obj) +
        diag(TraitEnv):vm(Gkeep, Ginv_obj) +
        diag(TraitEnv):ide(Gkeep) +
        diag(TraitEnv):Block,
      list(k_val = k, Ginv_obj = Ginv)
    )),
    residual = ~ dsum(~ ar1(Col):ar1(Row) | TraitEnv),
    sparse = ~ TraitEnv:Gdrop,
    data = data,
    na.action = na.method(x = "include"),
    maxit = 13,
    workspace = workspace
  )
  
  count <- 0
  while (!isTRUE(model$converge) && count < 100) {
    count <- count + 1
    model <- update(model)
    message("↻ Update ", count, " (Convergence = ", model$converge, ")")
  }
  
  if (isTRUE(model$converge)) {
    message("✅ ", model_name, " converged after ", count, " update(s).")
  } else {
    message("⚠️ ", model_name, " did not converge after ", count, " update(s).")
  }
  
  save(model, file = paste0("Data/", model_name, ".RData"))
}, error = function(e) {
  warning("❌ Model failed for ", model_name, " → ", e$message)
})

