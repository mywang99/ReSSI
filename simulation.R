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
RPPI_singlerun <- function(seed, n, N, p, d, b, eps, eps_all, attack, method_rme, method_train, intercept = T, rho = 0.2, cv = F, model = "logistic", num_round = 100, alpha = 0.1, const = 1, iter_max = 20, is_echo = F, is_plot = F, inference = "f_delta", fold = 5, method_opt = "One-step", h = 10){
  # generate data ----
  ## labeled data ----
  data_labeled <- data_generator(n = n, p = p, d = d, model = model, rho = rho)
  X_labeled <- data_labeled$X
  Y_labeled <- data_labeled$Y
  
  ## unlabeled data ----
  data_unlabeled <- data_generator(n = N, p = p, d = d, model = model, rho = rho)
  X2 <- data_unlabeled$X
  Y2 <- data_unlabeled$Y
  ## target parameter ----
  theta <- target_est(p = p, d = d, model = model, rho = rho)
  
  ## contamination ----
  X2_c <- X2
  
  if(attack == "fix"){
    print("Fixed Attack")
    v <- rep(1, p)
    b_index <- sample(1:N, size = floor(N * eps))
    for(k in b_index){
      X2_c[k, ] <- t(t(X2_c[k, ]) + b * v)
    }
  }else if(attack == "conditional"){
    print("Conditional Attack")
    X2_col1 <- X2[, 1]
    
    threshold <- sort(X2_col1, decreasing = TRUE)[floor(N*eps)]
    b_index <- which(X2_col1 >= threshold)
    
    v <- rep(1, p)
    for (k in b_index) {
      X2_c[k, ] <- t(t(X2[k, ]) + b * X2[k, 1] * v)
    }
  }else if(attack == "adversial"){
    print("Adversial Attack")
    M <- matrix(rnorm(p * p, mean = 0, sd = 0.5), nrow = p, ncol = p)
    A <- diag(p) + b * M
    b_index <- sample(1:N, size = floor(N * eps))
    for (k in b_index) {
      X2_c[k, ] <- A %*% X2[k, ]
    }
  }else if(attack == "eigen"){
    pca <- prcomp(X2, scale = FALSE)
    b_index <- sample(1:N, size = floor(N * eps))
    for (k in b_index) { 
      v <- pca$rotation[, 1]  
      X2_c[k, ] <- X2[k, ] + b * v
    }
  }
  id_clean <- setdiff(1:N, b_index)
  data_predict <- predict_func(data_labeled = data_labeled, X2 = X2, X2_c = X2_c, method_train = method_train, model = model)
  X1 <- data_predict$X1
  Y1 <- data_predict$Y1
  Y1_hat <- data_predict$Y1_hat
  Y2_hat <- data_predict$Y2_hat
  Y2_hat0 <- data_predict$Y2_hat0
  
  
  # result ----
  ## classical method ----
  res_classical <- classical(Y1 = Y_labeled,  X1 = X_labeled, theta = theta, alpha = alpha, model = model, intercept = intercept)
  res_classical
  
  ## RPPI ----
  method_name <- c("PPI", "PPI++", "PDC", "RePPI")
  
  output <- function(seed, cv = T, X1, Y1, Y1_hat, X2, X2_c, Y2_hat, theta, method, model, eps_all, eps, alpha, const, b_index, method_train, method_rme, method_opt, iter_max,  intercept, is_echo = F, is_plot = F, res_classical = NULL, theta_hat_PPI = NULL, theta_hat_RPPI = NULL){
    res_list <- list()
    res_PPI <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model, robust = F, th = 0.01, is_plot = is_plot, theta_hat_PPI = theta_hat_PPI, res_classical = res_classical, intercept = intercept,
                   b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps, method_opt = method_opt, iter_max = iter_max, tol = 1e-03, is_echo = is_echo)
    
    res_RPPI <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model, robust = T, th = 0.01, is_plot = is_plot, theta_hat_PPI = theta_hat_RPPI, res_classical = res_classical, intercept = intercept,
                    b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps, method_opt = method_opt, iter_max = iter_max, tol = 1e-03, is_echo = is_echo)

    res_list[[paste0("PPI_", method)]] <- res_PPI

    res_list[[paste0("RPPI_", method)]] <- res_RPPI
    eps_hat <- eps
    ## Choose eps by CV ----
    if(cv == T){
        res_cv_eps <- CV_eps(h = h, seed = seed, X1 = X1, Y1 = Y1, Y1_hat = Y1_hat, X2 = X2, X2_c = X2_c, Y2_hat = Y2_hat, theta = theta, method = method, robust = T, model = model, method_opt = method_opt,  intercept = intercept,
                           eps_all = eps_all, alpha = alpha, const = 1, b_index = b_index, method_train = method_train, method_rme = method_rme, tol = 1e-02, iter_max = iter_max, fold = fold, is_echo = is_echo)
        eps_hat <- res_cv_eps$eps_hat

      res_RPPI_cv <- PPI(inference = inference, h = h, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, theta = theta, method = method, model = model, robust = T, th = 0.01, is_plot = is_plot, res_classical = res_classical, intercept = intercept,
                         b_index = b_index, method_rme = method_rme, alpha = alpha, const = const, eps = eps_hat, method_opt = method_opt, iter_max = iter_max, tol = 1e-03, is_echo = is_echo)

      res_list[[paste0("RPPI_cv_", method)]] <- res_RPPI_cv
    }
    return(list(res_list = res_list, eps_hat = eps_hat))
  }
  
  res_list <- list()
  eps_hat <- c()
  theta_hat_PPI <- theta_hat_RPPI <- NULL
  for(method in method_name){
    print(method)
    res <- output(seed = seed, cv = cv, X1 = X1, Y1 = Y1, Y1_hat = Y1_hat, X2 = X2, X2_c = X2_c, Y2_hat = Y2_hat, theta = theta, method = method, model = model, eps_all = eps_all, eps = eps, alpha = alpha, const = const, b_index = b_index, method_train = method_train, intercept = intercept,
                  method_rme = method_rme, method_opt = method_opt, iter_max = iter_max, theta_hat_PPI = theta_hat_PPI, theta_hat_RPPI = theta_hat_RPPI, res_classical = res_classical, is_echo = F)
    if(method == "PPI"){
      theta_hat_PPI <- res$res_list$PPI_PPI$theta_hat
      theta_hat_RPPI <- res$res_list$RPPI_PPI$theta_hat
    }
    res_list <- append(res_list, res$res_list)
    eps_hat <- c(eps_hat, res$eps_hat)
  }
  res_list <- append(res_list, list(classical = res_classical))
  method_name_list <- names(res_list)
  print(method_name_list)
  
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
N <- 10000
p <- 3
model <- "linear"                 # "mean", "linear", or "logistic"
d <- if(model == "mean") 3 else 1
b <- 20
attack <- "fix"
rho <- 0.2
intercept <- TRUE
folder_name <- file.path("results", paste0("simulation_", model, "_eps", eps, "_", inference))

# Repeated experiments ----
set.seed(seed_master)
seed_all <- sample.int(.Machine$integer.max, rep_times)
dir.create(folder_name, recursive = TRUE, showWarnings = FALSE)
para <- mget(c("h", "rep_times", "seed_master", "inference", "n", "alpha", "const", "method_rme", "method_opt", "method_train", "iter_max", "cv", "fold", "eps_all", "eps", "model", "intercept", "N", "p", "d", "b", "attack", "rho"))
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
    res_single <- RPPI_singlerun(seed = seed, n = n, N = N, p = p, d = d, b = b,
    eps = eps, eps_all = eps_all, attack = attack, method_rme = method_rme,
    method_train = method_train, intercept = intercept, rho = rho, cv = cv,
    model = model, alpha = alpha, const = const, iter_max = iter_max,
    h = h, inference = inference, fold = fold, method_opt = method_opt)
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
