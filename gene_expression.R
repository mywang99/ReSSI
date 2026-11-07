# Run this script from the repository directory.
shhh <- suppressPackageStartupMessages
packages <- c("Rcpp", "RcppArmadillo", "nnet", "rpart", "caret", "xgboost", "foreach")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if(length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
# Keep this reduced rme package separate from other installed versions.
dir.create(".Rlibrary", showWarnings = FALSE)
.libPaths(c(file.path(getwd(), ".Rlibrary"), .libPaths()))
if(!dir.exists(".Rlibrary/rme"))
  install.packages("rme_1.0.tar.gz", repos = NULL, type = "source", lib = ".Rlibrary")
for(pkg in c("rme", "nnet", "rpart", "caret", "xgboost", "foreach"))
  shhh(library(pkg, character.only = TRUE))
source("ReSSI_function.R")

# Experiment function ----
gene_singlerun <- function(seed, Y, Y_pred, Y_pred_c, theta, method_train, method_opt, cv = T, n, eps, eps_all, intercept = F, method_rme, model = "logistic", alpha = 0.1, const = 1, iter_max = 20, tol = 1e-02, is_echo = F, is_plot = F, inference = "f_delta", fold = 5, h = 10){
  
  ## random split data into labeled data and unlabeled data ----
  id_labeled <- sample(c(1:length(Y)), n, replace = F)
  id_unlabeled <- setdiff(c(1:length(Y)), id_labeled)

  Y1 <- Y[id_labeled]
  Y1_hat <- Y_pred[id_labeled]
  
  Y2_hat0 <- Y_pred[id_unlabeled]
  Y2_hat <- Y2_hat0
  
  N <- length(id_unlabeled)

  if(floor(N * eps) > 0){
    b_index <- c(1: floor(N * eps))
    id_clean <- c((floor(N * eps) + 1):N)
    Y_c <- sample(Y_pred_c, size = length(b_index), replace = T)
    Y2_hat[b_index] <- Y_c
  }else{
    b_index <- c()
    id_clean <- c(1:N)
  }
  # result ----
  
  ## classical method ----
  res_classical <- classical(Y1 = Y1, X1 = NULL, theta = theta, model = model, intercept = intercept, alpha = alpha)
  res_classical
  
  X1 <- matrix(NA_real_, nrow = length(Y1), ncol = 1)
  X2 <- X2_c <- matrix(NA_real_, nrow = length(Y2_hat), ncol = 1)
  ## RPPI ----
  method_name <- c("PPI", "PPI++", "PDC", "RePPI")
  
  output <- function(seed, cv = T, X1, Y1, Y1_hat, X2, X2_c, Y2_hat, theta, method, model, intercept, eps_all, eps, alpha, const, b_index, method_train, method_rme, method_opt, iter_max, res_classical = NULL, theta_hat_PPI = NULL, theta_hat_RPPI = NULL){
    res_list <- list()
    res_PPI <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model, intercept = intercept, robust = F, th = 0.01, is_plot = is_plot, theta_hat_PPI = theta_hat_PPI, res_classical = res_classical,
                   b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps, method_opt = method_opt, iter_max = iter_max, tol = 1e-01, is_echo = is_echo)
    
    res_RPPI <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model,intercept = intercept, robust = T, th = 0.01, is_plot = is_plot, theta_hat_PPI = theta_hat_RPPI, res_classical = res_classical,
                    b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps, method_opt = method_opt, iter_max = iter_max, tol = 1e-01, is_echo = is_echo)

    res_list[[paste0("PPI_", method)]] <- res_PPI

    res_list[[paste0("RPPI_", method)]] <- res_RPPI
    eps_hat <- eps
    ## Choose eps by CV ----
    if(cv == T){
        res_cv_eps <- CV_eps(h = h, seed = seed, X1 = X1, Y1 = Y1, Y1_hat = Y1_hat, X2 = X2, X2_c = X2_c, Y2_hat = Y2_hat, theta = theta, method = method, robust = T, model = model, intercept = intercept, method_opt = method_opt, is_echo = is_echo,
                           eps_all = eps_all, alpha = alpha, const = 1, b_index = b_index, method_train = method_train, method_rme = method_rme, tol = 1e-01, iter_max = iter_max, fold = fold)
        eps_hat <- res_cv_eps$eps_hat

      res_RPPI_cv <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model, intercept = intercept, robust = T, th = 0.01, is_plot = is_plot, res_classical = res_classical, is_echo = is_echo,
                         b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps_hat, method_opt = method_opt, iter_max = iter_max, tol = 1e-01)

      res_list[[paste0("RPPI_cv_", method)]] <- res_RPPI_cv
    }
    return(list(res_list = res_list, eps_hat = eps_hat))
  }
  
  res_list <- list()
  eps_hat <- c()
  theta_hat_PPI <- theta_hat_RPPI <- NULL
  for(method in method_name){
    print(method)
    res <- output(seed = seed, cv = cv, X1 = X1, Y1 = Y1, Y1_hat = Y1_hat, X2 = X2, X2_c = X2_c, Y2_hat = Y2_hat, theta = theta, method = method, model = model, intercept = intercept, eps_all = eps_all, eps = eps, alpha = alpha, const = const, b_index = b_index, method_train = method_train, 
                  method_rme = method_rme, method_opt = method_opt, iter_max = iter_max, theta_hat_PPI = theta_hat_PPI, theta_hat_RPPI = theta_hat_RPPI, res_classical = res_classical)
    if(method == "PPI"){
      theta_hat_PPI <- res$res_list$PPI_PPI$theta_hat
      theta_hat_RPPI <- res$res_list$RPPI_PPI$theta_hat
    }
    res_list <- append(res_list, res$res_list)
    eps_hat <- c(eps_hat, res$eps_hat)
  }
  res_list <- append(res_list, list(classical = res_classical))
  
  # Summary ----
  theta_true <- if(model != "mean" && intercept) theta[-1] else theta
  method_name_list <- names(res_list)
  n_res <- length(res_list)
  est_error <- coverage <- ellip_vol <- width_CI <- numeric(n_res)
  for(i in seq_along(res_list)) {
    res <- res_list[[i]]
    theta_hat <- if(model != "mean" && intercept) res$theta_hat[-1] else res$theta_hat
    est_error[i] <- sqrt(sum((theta_hat - theta_true)^2))
    coverage[i] <- as.numeric(res$coverage)
    ellip_vol[i] <- res$ellip_vol
    width_CI[i] <- mean(res$CI[, 2] - res$CI[, 1])
  }
  out <- data.frame(method = method_name_list, est_error, coverage, ellip_vol, width_CI)
  names(eps_hat) <- method_name
  print(out)
  return(list(theta = theta, res_list = res_list, eps_hat = eps_hat, summary = out))
}

# Parameter setting ----
rep_times <- 5                    # Set to 500 for the full experiment.
seed_master <- 20261107
inference <- "f_delta"            # "f_delta" or "delta"
n <- 1000
alpha <- 0.1
const <- 1
method_rme <- "Filtering"         # "Filtering", "geom", or "trimming_mean"
h <- 10                          # MAD threshold for trimming_mean
method_opt <- "One-step"
method_train <- "fix"
iter_max <- 20
cv <- TRUE
fold <- 5
eps_all <- c(seq(0, 0.2, 0.025), 1)
eps <- 0.05
model <- "mean"
intercept <- FALSE
attack <- "high"                  # "high" or "middle"
level_l <- 7
level_h <- 13

# Data processing ----
data <- read.csv("data/Native_allchunks_test_tpu_model.csv")
data_Y <- data[, 3]
data_pred <- data[, 4]
id_low <- which(data_pred <= level_l)
id_c <- if(attack == "middle") which(data_pred > level_l & data_pred <= level_h) else which(data_pred > level_h)
Y <- data_Y[id_low]
Y_pred <- data_pred[id_low]
Y_pred_c <- data_pred[id_c]
theta <- mean(Y)
folder_name <- file.path("results", paste0("gene_", attack, "_eps", eps, "_", inference))

# Repeated experiments ----
set.seed(seed_master)
seed_all <- sample.int(.Machine$integer.max, rep_times)
dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
para <- mget(c("h", "rep_times", "seed_master", "inference", "n", "alpha", "const", "method_rme", "method_opt", "method_train", "iter_max", "cv", "fold", "eps_all", "eps", "model", "intercept", "attack", "level_l", "level_h"))
para_file <- file.path(folder_name, "parameters.rds")
if(file.exists(para_file) && !identical(readRDS(para_file), para))
  stop("Parameters changed: use a new folder_name to avoid mixing results.")
saveRDS(para, para_file)
saveRDS(seed_all, file.path(folder_name, "seeds.rds"))
res_all <- vector("list", rep_times)
for(r in seq_len(rep_times)) {
  seed <- seed_all[r]
  result_file <- file.path(folder_name, sprintf("rep_%04d.rds", r))
  if(file.exists(result_file)) {
    saved <- readRDS(result_file)
    stopifnot(saved$seed == seed)
    res_single <- saved$res_single
  } else {
    set.seed(seed)
    res_single <- gene_singlerun(seed = seed, Y = Y, Y_pred = Y_pred, Y_pred_c = Y_pred_c,
    theta = theta, method_train = method_train, method_opt = method_opt,
    cv = cv, n = n, eps = eps, eps_all = eps_all, intercept = intercept,
    method_rme = method_rme, model = model, alpha = alpha, const = const,
    iter_max = iter_max, h = h, inference = inference, fold = fold)
    saveRDS(list(seed = seed, res_single = res_single), result_file)
  }
  res_all[[r]] <- transform(res_single$summary, repetition = r)
  message("Completed ", r, "/", rep_times)
}

# Summary ----
results <- do.call(rbind, res_all)
metrics <- c("est_error", "coverage", "ellip_vol", "width_CI")
data_res <- aggregate(results[metrics], list(method = results$method), mean)
data_se <- aggregate(results[metrics], list(method = results$method),
                     function(x) sd(x) / sqrt(length(x)))
print(data_res)
write.csv(results, file.path(folder_name, "replicates.csv"), row.names = FALSE)
write.csv(data_res, file.path(folder_name, "summary.csv"), row.names = FALSE)
write.csv(data_se, file.path(folder_name, "standard_errors.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(folder_name, "sessionInfo.txt"))
