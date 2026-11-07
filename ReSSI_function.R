# Generate simulation data for mean estimation, linear regression, or logistic regression; n is the sample size, p the feature dimension, and d the response dimension.
data_generator <- function(n, p, d, model, rho = 0.2) {
    if (model == "mean") {
        beta1 <- rep(1/p, p)
        beta2 <- c()
        for (i in 1:p) {
            beta2[i] <- (i - 1)/p
        }
        beta3 <- c(rep(1/2, floor(p/2)), rep(0, p - floor(p/2)))
        if (d == 1) {
            beta = beta1
        }
        else if (d == 2) {
            beta = cbind(beta1, beta2)
        }
        else if (d == 3) {
            beta = cbind(beta1, beta2, beta3)
        }
        X <- matrix(rnorm(n * p), n, p)
        Sig_0 <- diag(p)
        for (ii in 1:p) {
            for (jj in 1:p) {
                Sig_0[ii, jj] <- rho^(abs(ii - jj))
            }
        }
        eig = eigen(Sig_0)
        Sig_0.5 = eig$vectors %*% (t(eig$vectors) * eig$values^(0.5))
        mu_X = rep(0, p)
        X <- t(t(X %*% Sig_0.5) + mu_X)
        epsilon <- matrix(rnorm(n * d), n, d)
        mu <- rep(1, d)
        Y <- t(t(X %*% beta + epsilon) + mu)
    }
    if (model == "linear") {
        X <- matrix(rnorm(n * p), n, p)
        Sig_0 <- diag(p)
        for (ii in 1:p) {
            for (jj in 1:p) {
                Sig_0[ii, jj] <- rho^(abs(ii - jj))
            }
        }
        eig = eigen(Sig_0)
        Sig_0.5 = eig$vectors %*% (t(eig$vectors) * eig$values^(0.5))
        mu_X = rep(0, p)
        X <- t(t(X %*% Sig_0.5) + mu_X)
        epsilon <- matrix(rnorm(n * d), n, d)
        beta1 <- rep(1, p)
        beta2 <- rep(0.3, p)
        if (p > 3) {
            beta1 <- beta1/sqrt(p)
            beta2 <- beta2/sqrt(p)
        }
        Y <- X %*% beta1 + epsilon + X^2 %*% beta2
    }
    else if (model == "logistic") {
        Sig_0 <- diag(p)
        for (ii in 1:p) {
            for (jj in 1:p) {
                Sig_0[ii, jj] <- rho^(abs(ii - jj))
            }
        }
        eig = eigen(Sig_0)
        Sig_0.5 = eig$vectors %*% (t(eig$vectors) * eig$values^(0.5))
        mu_X = rep(1, p)
        X1 <- matrix(rnorm(n * p), n, p)
        X1 <- t(t(X1 %*% Sig_0.5) + mu_X)
        X2 <- matrix(rnorm(n * p), n, p)
        X2 <- t(t(X2 %*% Sig_0.5) - mu_X)
        X_partition_label <- rbinom(n, 1, 0.5)
        X <- X1 * X_partition_label + X2 * (1 - X_partition_label)
        beta0 <- 4
        beta1 <- rep(0.5, p)
        beta2 <- rep(-1, p)
        linear_combined <- beta0 + X %*% beta1 + X^2 %*% beta2
        mu_Y <- 1/(1 + exp(-linear_combined))
        Y <- rbinom(n, 1, mu_Y)
    }
    else if (model == "logistic1") {
        X <- matrix(rnorm(n * p), n, p)
        Sig_0 <- diag(p)
        rho <- 0.5
        for (ii in 1:p) {
            for (jj in 1:p) {
                Sig_0[ii, jj] <- rho^(abs(ii - jj))
            }
        }
        eig = eigen(Sig_0)
        Sig_0.5 = eig$vectors %*% (t(eig$vectors) * eig$values^(0.5))
        mu_X <- rep(0, p)
        X <- t(t(X %*% Sig_0.5) + mu_X)
        epsilon <- matrix(rnorm(n * d), n, d)
        beta1 <- rep(1, p)
        beta2 <- rep(1, p)
        mu_Y <- 1/(1 + exp(-X %*% beta1 - X^2 %*% beta2 - epsilon))
        Y <- rbinom(n, 1, mu_Y)
    }
    return(list(X = X, Y = Y))
}

# Estimate the scalar PPI++ weight lambda from observed and imputed gradients and the Hessian.
lam_est <- function(g1, g1_f, g2_f, Hessian) {
    N <- dim(g2_f)[1]
    n <- dim(g1)[1]
    Sig1 <- cov(g1, g1_f) + cov(g1_f, g1)
    Sig_f <- var(g1_f) * (n + 1)/n
    Hessian_inv <- solve(Hessian)
    lam <- sum(diag(Hessian_inv %*% Sig1 %*% Hessian_inv))/sum(diag(Hessian_inv %*% Sig_f %*% Hessian_inv))/(2 * 
        (1 + n/N))
    return(lam)
}

# Compute the volume of a chi-square Wald ellipsoid with covariance Sigma and confidence level 1 - alpha.
ellipsoid_vol <- function(Sigma, alpha = 0.05) {
    p <- ncol(Sigma)
    c <- qchisq(1 - alpha, df = p)
    vol <- (pi * c)^(p/2)/gamma(p/2 + 1) * sqrt(det(Sigma))
    return(vol)
}

# Compute the labeled-only baseline: sample mean, least squares, or logistic regression.
# With CI = TRUE, also return covariance, joint coverage, region volume, and marginal intervals; theta is used to assess coverage.
classical <- function(X1, Y1, theta, model, alpha, intercept = T, CI = TRUE) {
    if (is.vector(Y1) == T) {
        Y1 <- matrix(Y1, ncol = 1)
    }
    n = dim(Y1)[1]
    if (is.null(X1)) {
        X1 <- matrix(1, n, 1)
    }
    if (intercept == T) {
        X1_1 <- cbind(rep(1, n), X1)
    }
    else {
        X1_1 <- X1
    }
    if (model == "mean") {
        theta_hat <- apply(Y1, 2, mean)
    }
    else if (model == "linear") {
        theta_hat <- solve(t(X1_1) %*% X1_1) %*% t(X1_1) %*% Y1
    }
    else if (model == "logistic") {
        if (intercept == T) {
            theta_hat <- glm(Y1 ~ X1, family = binomial())$coefficients
        }
        else {
            theta_hat <- glm(Y1 ~ X1 - 1, family = binomial())$coefficients
        }
    }
    theta_hat_all <- theta_hat
    if (CI == T) {
        res <- gradient(X = X1_1, Y = Y1, theta = theta_hat, model = model)
        grad <- res$grad_mat
        Hessian <- res$Hessian
        Sig <- var(grad) * (n - 1)/n
        if (model != "mean" && intercept == T) {
            Hessian <- Hessian[-1, -1, drop = FALSE]
            Sig <- Sig[-1, -1, drop = FALSE]
            theta_hat = theta_hat[-1]
            theta = theta[-1]
        }
        V <- solve(Hessian) %*% Sig %*% solve(Hessian)
        p <- nrow(V)
        difference <- theta_hat - theta
        coverage <- as.numeric(n * t(difference) %*% solve(V) %*% difference <= qchisq(1 - alpha, p))
        ellip_vol <- ellipsoid_vol(Sigma = V/n, alpha = alpha)
        halfwidth <- qnorm(1 - alpha/2) * sqrt(diag(V)/n)
        CI_l <- theta_hat - halfwidth
        CI_r <- theta_hat + halfwidth
        CI <- cbind(CI_l, CI_r)
        return(list(theta_hat = theta_hat_all, coverage = coverage, Sig = V, ellip_vol = ellip_vol, CI = CI))
    }
    else {
        return(list(theta_hat = theta_hat_all))
    }
}

# Return per-sample gradients grad_mat, their mean grad, and the empirical Hessian.
# For mean estimation, the implementation uses grad_mat = Y, grad = colMeans(Y), and an identity Hessian.
gradient <- function(X, Y, theta, model) {
    n <- dim(X)[1]
    p <- length(theta)
    if (model == "mean") {
        if (is.vector(Y)) {
            Y <- matrix(Y, ncol = 1)
        }
        grad <- apply(Y, 2, mean)
        grad_mat <- Y
        Hessian <- diag(p)
    }
    else if (model == "linear") {
        grad_mat <- X * as.vector(X %*% theta - Y)
        grad <- apply(grad_mat, 2, mean)
        Hessian <- t(X) %*% X/n
    }
    else if (model == "logistic") {
        Y_hat <- 1/(1 + exp(-X %*% theta))
        grad <- t(X) %*% c(Y_hat - Y)/nrow(X)
        grad_mat <- cross_entropy_loss_gradient_mat(X, Y, Y_hat)
        Hessian <- (t(X) %*% (c(Y_hat * (1 - Y_hat)) * X))/nrow(X)
    }
    return(list(grad = grad, grad_mat = grad_mat, Hessian = Hessian))
}

# Fit a response model using features and existing predictions to construct RePPI imputed gradients.
# Return a prediction function scaled by 1 / (1 + r), where r is the labeled-to-unlabeled sample size ratio.
grad_fit <- function(X, Y, Y_hat, model, theta, r, intercept = T, fit_method = "mlp") {
    if (sum(is.na(X)) > 0) {
        covariates <- Y_hat
    }
    else {
        covariates <- cbind(X, Y_hat)
    }
    if (fit_method == "mlp") {
        fit_model <- nnet(Y ~ ., data = as.data.frame(covariates), size = 32, skip = TRUE, maxit = 1000)
    }
    else if (fit_method == "mlpclass") {
        fit_model <- nnet(Y ~ ., data = as.data.frame(covariates), size = 32, skip = TRUE, maxit = 1000, 
            linout = FALSE)
    }
    else if (fit_method == "linreg") {
        fit_model <- lm(Y ~ ., data = as.data.frame(covariates))
    }
    else if (fit_method == "logistic") {
        fit_model <- glm(Y ~ ., data = as.data.frame(covariates), family = binomial())
    }
    else if (fit_method == "xgboost") {
        dtrain <- xgb.DMatrix(data = as.matrix(covariates), label = Y)
        param <- list(max_depth = 3, eta = 0.3, silent = 1, objective = "reg:squarederror")
        fit_model <- xgb.train(param, dtrain, nrounds = 10)
    }
    else if (fit_method == "tree") {
        fit_model <- rpart(Y ~ ., data = as.data.frame(covariates), method = "anova")
    }
    else {
        stop(paste("Fit method", fit_method, "not recognized."))
    }
    # Predict responses on new samples, convert them to imputed gradients, and apply the scaling factor.
    f <- function(X_new, Y_hat_new) {
        if (sum(is.na(X)) > 0) {
            covariates_new <- as.data.frame(Y_hat_new)
        }
        else {
            covariates_new <- as.data.frame(cbind(X_new, Y_hat_new))
        }
        names(covariates_new) <- names(fit_model$model)[-1]
        if (fit_method == "linreg") {
            pred <- predict(fit_model, newdata = covariates_new)
            pred <- as.numeric(pred)
        }
        else if (fit_method == "logistic") {
            pred <- predict(fit_model, newdata = as.data.frame(covariates_new), type = "response")
        }
        else if (fit_method == "xgboost") {
            dtest <- xgb.DMatrix(data = as.matrix(covariates_new))
            pred <- predict(fit_model, dtest)
        }
        else {
            pred <- predict(fit_model, newdata = as.data.frame(covariates_new))
        }
        if (model == "mean") {
            grad <- pred
        }
        else {
            pred <- matrix(pred, ncol = 1)
            if (intercept == T) {
                X_new_1 <- cbind(rep(1, nrow(X_new)), X_new)
            }
            else {
                X_new_1 <- X_new
            }
            grad <- gradient(X_new_1, pred, theta, model)$grad_mat
        }
        return(1/(1 + r) * grad)
    }
    return(f)
}

# Construct or adjust imputed gradients for PPI, PPI++, PDC, or RePPI.
# With robust = TRUE, aggregate unlabeled gradients using the selected robust estimator; return adjusted gradients and reusable parameters para_list.
impute_grad <- function(model, fit_method, g1, g1_f, g2_f, X1, X2_c, Y1, Y1_hat, Y2_hat, theta, robust, 
    b_index, method, method_rme, eps, is_echo = F, intercept = T, theta_hat_PPI = NULL, para_list = NULL, 
    is_plot = F, const = 1, method_opt = "One-step", h = 2) {
    n <- dim(g1)[1]
    N <- dim(g2_f)[1]
    d <- dim(Y1)[2]
    if (intercept == T) {
        X1_1 <- cbind(rep(1, n), X1)
        X2_c_1 <- cbind(rep(1, N), X2_c)
    }
    else {
        X1_1 <- X1
        X2_c_1 <- X2_c
    }
    theta_hat <- theta
    if (method == "PDC") {
        if (is.null(para_list) == TRUE) {
            M <- cov(g1, g1_f) %*% solve(var(g1_f) * (n - 1)/n)
            M <- N * M/(N + n)
            para_list <- list(M = M)
        }
        M <- para_list$M
        g1_f <- g1_f %*% t(M)
        g2_f <- g2_f %*% t(M)
    }
    else if (method == "PPI++") {
        if (is.null(para_list) == TRUE) {
            if (model != "mean") {
                if (is.null(theta_hat_PPI) == T) {
                  theta_hat_PPI <- PPI(X1 = X1, X2 = X2_c, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, 
                    model = model, intercept = intercept, theta = theta, CI = F, method = "PPI", b_index = b_index, 
                    robust = robust, method_rme = method_rme, method_opt = method_opt, alpha = 0.05, 
                    const = 1, eps = eps, is_plot = F, h = h)$theta_hat
                }
                res_labeled_PPI <- gradient(X = X1_1, Y = Y1, theta = theta_hat_PPI, model = model)
                g1_PPI <- res_labeled_PPI$grad_mat
                Hessian <- res_labeled_PPI$Hessian
                g2_f_PPI <- gradient(X = X2_c_1, Y = Y2_hat, theta = theta_hat_PPI, model = model)$grad_mat
                g1_f_PPI <- gradient(X = X1_1, Y = Y1_hat, theta = theta_hat_PPI, model = model)$grad_mat
                lam <- lam_est(g1 = g1_PPI, g1_f = g1_f_PPI, g2_f = g2_f_PPI, Hessian = Hessian)
            }
            else {
                Hessian <- gradient(X = X1_1, Y = Y1, theta = theta, model = model)$Hessian
                lam <- lam_est(g1 = g1, g1_f = g1_f, g2_f = g2_f, Hessian = Hessian)
            }
            para_list <- list(lam = lam)
        }
        lam <- para_list$lam
        g1_f <- lam * g1_f
        g2_f <- lam * g2_f
    }
    else if (method == "RePPI") {
        fit_method <- "linreg"
        if (is.null(para_list) == TRUE) {
            r <- n/N
            if (model == "mean" && d > 1) {
                grad_model <- list()
                for (j in 1:d) {
                  grad_model[[j]] <- grad_fit(X = X1, Y = Y1[, j], Y_hat = Y1_hat[, j], model = model, 
                    theta = theta, r = r, fit_method = fit_method, intercept = intercept)
                  g1_f[, j] <- grad_model[[j]](X_new = X1, Y_hat_new = Y1_hat[, j])
                  g2_f[, j] <- grad_model[[j]](X_new = X2_c, Y_hat_new = Y2_hat[, j])
                }
            }
            else {
                grad_model <- grad_fit(X = X1, Y = Y1, Y_hat = Y1_hat, model = model, theta = theta, 
                  r = r, fit_method = fit_method, intercept = intercept)
                g1_f <- grad_model(X_new = X1, Y_hat_new = Y1_hat)
                g2_f <- grad_model(X_new = X2_c, Y_hat_new = Y2_hat)
            }
            M <- cov(g1, g1_f) %*% solve(var(g1_f) * (n - 1)/n)
            g1_f <- g1_f %*% t(M)
            g2_f <- g2_f %*% t(M)
            para_list <- list(g1_f = g1_f, g2_f = g2_f)
        }
        g1_f <- para_list$g1_f
        g2_f <- para_list$g2_f
    }
    if (robust == T) {
        if (eps == 0) {
            return(list(g1_f = g1_f, g2_f = g2_f, para_list = para_list))
        }
        else {
            g2_f_sd <- apply(g2_f, 2, sd)
            g2_f_sd <- max(g2_f_sd)
            g2_f_mean <- apply(g2_f, 2, mean)
            g2_f_scale <- t((t(g2_f) - g2_f_mean)/g2_f_sd)
            Sig_hat <- var(t((t(g1_f) - g2_f_mean)/g2_f_sd)) * (n - 1)/n
            if (method_rme == "Filtering") {
                w0 <- rep(1/N, N)
                res_gd <- rme_solver_refine(X = g2_f_scale, is_echo = F, w0 = w0, ep = const * eps, loop_max = 300, 
                  lr = 0.2, dc_rate = 0.1, iter_power = 300, Sig_0 = Sig_hat, tol = 1e-06, tol_pow = 1e-08)
                w_new <- res_gd[[2]]
                w <- w_new
                g2_hat <- t(t(w) %*% g2_f_scale) * g2_f_sd + g2_f_mean
            }
            else if (method_rme == "geom") {
                w <- rep(1/N, N)
                g2_hat <- geom_solver_C(X = g2_f_scale, tol = 1e-10, is_echo = FALSE)[[1]] * g2_f_sd + g2_f_mean
            }
            else if (method_rme == "trimming_mean") {
                w <- rep(1/N, N)
                g2_hat <- trimming_mean(g2_f_scale, h = h) * g2_f_sd + g2_f_mean
            }
            else {
                stop("Unknown robust method: use Filtering, geom, or trimming_mean")
            }
            id_det <- which(w < 0.01/N)
            id_fd <- setdiff(id_det, b_index)
            id_nd <- setdiff(b_index, id_det)
            if (is_echo == T) {
                cat(sprintf("Number of detection: %d false discovery: %d not detected: %d\n", length(id_det), 
                  length(id_fd), length(id_nd)))
            }
            if (is_plot == T) {
                N_clean_plot <- floor(N * eps * 3)
                par(mfrow = c(5, 1))
                cat(sprintf("w_max: %f\n ", max(w[1:floor(N * eps)]) * N))
                plot(x = c(1:N_clean_plot), y = w[1:N_clean_plot] * N) + points(x = b_index, y = w[b_index] * 
                  N, col = "red")
                for (j in 1:length(theta)) {
                  plot(x = c(1:N_clean_plot), y = g2_f[1:N_clean_plot, j]) + points(x = b_index, y = g2_f[b_index, 
                    j], col = "red")
                }
                id_clean1 <- intersect(id_clean, c(1:N_clean_plot))
                id_fd1 <- intersect(id_fd, c(1:N_clean_plot))
                plot(x = id_clean1, y = apply(t(g2_f[id_clean1, ]) - g2_f_mean, 2, function(x) sum(x^2))) + 
                  points(x = id_fd1, y = apply(t(matrix(g2_f[id_fd1, ], ncol = length(theta), byrow = T)) - 
                    g2_f_mean, 2, function(x) sum(x^2)), col = "red")
            }
            g2_f <- matrix(g2_hat, nrow = N, ncol = length(g2_hat), byrow = T)
        }
    }
    return(list(g1_f = g1_f, g2_f = g2_f, para_list = para_list))
}

# Compute a PPI-type estimate, using robust aggregation of unlabeled imputed gradients when robust = TRUE.
# Support one-step and iterative updates; CI = TRUE returns inference results, with f_delta or delta selected by inference.
# eps = 1 denotes the supervised rule; h is the fixed MAD threshold for trimming_mean.
PPI <- function(X1, X2_c, Y1, Y1_hat, Y2_hat, theta, alpha, b_index, eps, method, model, intercept = T, 
    theta_hat_PPI = NULL, res_classical = NULL, const = 1, robust = T, CI = T, method_rme = "Filtering", 
    th = 0.01, method_opt = "Newton", iter_max = 1000, tol = 1e-05, is_plot = F, is_echo = F, inference = "f_delta", 
    h = 2) {
    inference <- match.arg(inference, c("f_delta", "delta"))
    if (is.vector(Y1)) {
        Y1 <- matrix(Y1, nrow = length(Y1), ncol = 1)
    }
    if (is.vector(Y1_hat)) {
        Y1_hat <- matrix(Y1_hat, nrow = length(Y1_hat), ncol = 1)
    }
    if (is.vector(Y2_hat)) {
        Y2_hat <- matrix(Y2_hat, nrow = length(Y2_hat), ncol = 1)
    }
    n = dim(Y1_hat)[1]
    N = dim(Y2_hat)[1]
    if (is.null(X1)) {
        X1 <- matrix(NA, n, 1)
    }
    if (is.null(X2_c)) {
        X2_c <- matrix(NA, N, 1)
    }
    if (intercept == T) {
        X1_1 <- cbind(rep(1, n), X1)
        X2_c_1 <- cbind(rep(1, N), X2_c)
    }
    else {
        X1_1 <- X1
        X2_c_1 <- X2_c
    }
    if (is.null(res_classical) == T) {
        theta0 <- classical(X1 = X1, Y1 = Y1, model = model, intercept = intercept, CI = F)$theta_hat
    }
    else {
        theta0 <- res_classical$theta_hat
    }
    if (is_echo == T) {
        cat(sprintf("Iter %d: objective value: %3f\n", 0, sqrt(sum((theta0 - theta)^2))))
    }
    counter <- 0
    para_list <- NULL
    res_labeled <- gradient(X = X1_1, Y = Y1, theta = theta0, model = model)
    g1 <- res_labeled$grad_mat
    Hessian <- res_labeled$Hessian
    g2_f <- gradient(X = X2_c_1, Y = Y2_hat, theta = theta0, model = model)$grad_mat
    g1_f <- gradient(X = X1_1, Y = Y1_hat, theta = theta0, model = model)$grad_mat
    if (robust == T && eps == 1) {
        if (eps == 1) {
            if (is.null(res_classical) == T) {
                res <- classical(X1 = X1, Y1 = Y1, theta = theta, alpha = alpha, model = model, intercept = intercept)
            }
            else {
                res <- res_classical
            }
            return(res)
        }
    }
    if (model == "mean") {
        res_impute <- impute_grad(para_list = para_list, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, 
            Y2_hat = Y2_hat, model = model, intercept = intercept, b_index = b_index, is_echo = is_echo, 
            theta = theta0, theta_hat_PPI = theta_hat_PPI, g1 = g1, g1_f = g1_f, g2_f = g2_f, method = method, 
            robust = robust, method_rme = method_rme, eps = eps, is_plot = is_plot, const = const, method_opt = method_opt, 
            h = h)
        g1_f <- res_impute$g1_f
        g2_f <- res_impute$g2_f
        theta_hat <- apply(g1, 2, mean) + apply(g2_f, 2, mean) - apply(g1_f, 2, mean)
    }
    else if (model == "linear" && (robust == FALSE || eps == 0)) {
        res_impute <- impute_grad(para_list = para_list, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, 
            Y2_hat = Y2_hat, model = model, intercept = intercept, b_index = b_index, is_echo = is_echo, 
            theta = theta0, theta_hat_PPI = theta_hat_PPI, g1 = g1, g1_f = g1_f, g2_f = g2_f, method = method, 
            robust = robust, method_rme = method_rme, eps = eps, is_plot = is_plot, const = const, method_opt = method_opt, 
            h = h)
        g1_f <- res_impute$g1_f
        g2_f <- res_impute$g2_f
        if (method == "PPI") {
            Sig_inv <- solve(t(X2_c_1) %*% X2_c_1/N)
            theta_hat <- Sig_inv %*% (t(X2_c_1) %*% Y2_hat/N - t(X1_1) %*% (Y1_hat - Y1)/n)
        }
        else {
            Sig_inv <- solve(t(X1_1) %*% X1_1/n)
            theta_hat <- Sig_inv %*% (t(X1_1) %*% Y1/n + apply(g1_f, 2, mean) - apply(g2_f, 2, mean))
        }
    }
    else {
        for (t in 1:iter_max) {
            if (t > 1) {
                res_labeled <- gradient(X = X1_1, Y = Y1, theta = theta0, model = model)
                g1 <- res_labeled$grad_mat
                Hessian <- res_labeled$Hessian
                g2_f <- gradient(X = X2_c_1, Y = Y2_hat, theta = theta0, model = model)$grad_mat
                g1_f <- gradient(X = X1_1, Y = Y1_hat, theta = theta0, model = model)$grad_mat
            }
            res_impute <- impute_grad(para_list = para_list, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, 
                Y2_hat = Y2_hat, model = model, intercept = intercept, b_index = b_index, is_echo = is_echo, 
                theta = theta0, theta_hat_PPI = theta_hat_PPI, g1 = g1, g1_f = g1_f, g2_f = g2_f, method = method, 
                robust = robust, method_rme = method_rme, eps = eps, is_plot = is_plot, const = const, 
                method_opt = method_opt, h = h)
            grad <- apply(res_impute$g2_f, 2, mean) - apply(res_impute$g1_f, 2, mean) + apply(g1, 2, 
                mean)
            para_list <- res_impute$para_list
            if (method_opt == "Newton") {
                theta1 <- theta0 - solve(Hessian, grad)
            }
            else if (method_opt == "GD") {
                theta1 <- theta0 - lr * grad
            }
            else if (method_opt == "One-step") {
                theta1 <- theta0 - solve(Hessian, grad)
            }
            update_size <- sqrt(sum((theta1 - theta0)^2)/length(theta0))
            if (is_echo == T) {
                cat(sprintf("Iter %d: objective value: %3f update deviation: %3f\n", t, sqrt(sum((theta1 - 
                  theta)^2)), update_size))
            }
            if (method_opt == "One-step") 
                break
            if (update_size > 10 * sqrt(sum(theta0^2))) {
                theta1 <- theta0
                print("diverging")
                break
            }
            if (t == 1) {
                update_size0 <- update_size
            }
            if (update_size > update_size0) {
                counter <- counter + 1
                if (counter > 3) {
                  theta1 <- theta0
                  print("diverging")
                  break
                }
            }
            if (update_size < tol) {
                print("Converged")
                break
            }
            update_size0 <- update_size
            theta0 <- theta1
        }
        theta_hat <- theta1
    }
    g1_f <- gradient(X = X1_1, Y = Y1_hat, theta = theta_hat, model = model)$grad_mat
    res_impute <- impute_grad(para_list = para_list, X1 = X1, X2_c = X2_c, Y1 = Y1, Y1_hat = Y1_hat, 
        Y2_hat = Y2_hat, model = model, intercept = intercept, b_index = b_index, is_echo = F, theta = theta_hat, 
        theta_hat_PPI = theta_hat_PPI, g1 = g1, g1_f = g1_f, g2_f = g2_f, method = method, robust = F, 
        method_rme = method_rme, eps = eps, const = const, method_opt = method_opt, h = h)
    g1_f <- res_impute$g1_f
    res_labeled <- gradient(X = X1_1, Y = Y1, theta = theta_hat, model = model)
    g1 <- res_labeled$grad_mat
    Hessian <- res_labeled$Hessian
    Sig_f <- cov(g1_f) * (n - 1)/n
    Sig_delta <- cov(g1_f - g1) * (n - 1)/n
    V <- Sig_delta + if (inference == "f_delta") 
        (n/N) * Sig_f
    else 0
    if (!CI) 
        return(list(theta_hat = theta_hat))
    theta_hat_all <- theta_hat
    if (model != "mean" && intercept) {
        Hessian <- Hessian[-1, -1, drop = FALSE]
        V <- V[-1, -1, drop = FALSE]
        theta_hat <- theta_hat[-1]
        theta <- theta[-1]
    }
    H_inv <- solve(as.matrix(Hessian))
    Sigma_hat <- H_inv %*% V %*% H_inv
    p <- nrow(Sigma_hat)
    difference <- theta_hat - theta
    coverage <- as.numeric(n * t(difference) %*% solve(Sigma_hat) %*% difference <= qchisq(1 - alpha, 
        p))
    halfwidth <- qnorm(1 - alpha/2) * sqrt(diag(Sigma_hat)/n)
    list(theta_hat = theta_hat_all, coverage = coverage, Sig = Sigma_hat, ellip_vol = ellipsoid_vol(Sigma_hat/n, 
        alpha), CI = cbind(theta_hat - halfwidth, theta_hat + halfwidth))
}

# Train and return an XGBoost model for regression or binary classification.
# num_round controls boosting rounds
fit_xgboost <- function(X1, Y1, model = "linear", eta = 0.1, num_round = 100) {
    if (model == "linear") {
        param <- list(booster = "gbtree", max_depth = 7, eta = eta, objective = "reg:squarederror", eval_metric = c("rmse"))
    }
    else if (model == "logistic") {
        print("logistic")
        param <- list(booster = "gbtree", max_depth = 7, eta = eta, objective = "binary:logistic", eval_metric = c("logloss"))
    }
    dtrain <- xgb.DMatrix(data = X1, label = Y1)
    tree <- xgb.train(param, dtrain, nrounds = num_round, verbose = 0)
    return(tree)
}

# Select a robust-rule parameter from eps_all using held-out labeled validation loss.
# Return eps_hat, mean validation losses error, and diagnostic errors obj.
# The candidate 1 denotes the supervised rule.
CV_eps <- function(seed, X1, Y1, X2, X2_c, Y1_hat, Y2_hat, theta, b_index, method, model, robust, method_train, 
    method_rme, method_opt, alpha, const, eps_all, intercept = T, fold = 5, tol = 0.01, iter_max = 20, 
    is_echo = F, h = 2) {
    if (is.vector(Y1)) {
        Y1 <- matrix(Y1, nrow = length(Y1), ncol = 1)
    }
    if (is.vector(Y1_hat)) {
        Y1_hat <- matrix(Y1_hat, nrow = length(Y1_hat), ncol = 1)
    }
    if (is.vector(Y2_hat)) {
        Y2_hat <- matrix(Y2_hat, nrow = length(Y2_hat), ncol = 1)
    }
    n <- nrow(Y1)
    p <- ncol(X1)
    local_size <- n%/%fold
    sizes <- rep(local_size, fold)
    sizes[1] <- sizes[1] + n - sum(sizes)
    set.seed(seed)
    group <- sample(rep(1:fold, sizes), n, replace = F)
    error_eps <- obj_eps <- rep(0, length(eps_all))
    obj <- matrix(0, nrow = length(eps_all), ncol = fold)
    error <- matrix(0, nrow = length(eps_all), ncol = fold)
    for (i in 1:fold) {
        id_test <- which(group == i)
        id_train <- setdiff(c(1:n), id_test)
        X_train <- X1[id_train, , drop = FALSE]
        Y_train <- matrix(Y1[id_train, ], nrow = length(id_train), ncol = ncol(Y1))
        Y_train_hat <- matrix(Y1_hat[id_train, ], nrow = length(id_train), ncol = ncol(Y1))
        if (method_train == "crossfitting") {
            res_crossfitting <- cross_prediction(X1 = X_train, Y1 = Y_train, X2 = X2, X2_c = X2_c, model = model, 
                Kfold = 5, num_round = 100)
            Y2_hat <- res_crossfitting$Y2_hat
            Y_train_hat <- res_crossfitting$Y1_hat
        }
        X_test <- X1[id_test, , drop = FALSE]
        Y_test <- matrix(Y1[id_test, ], nrow = length(id_test), ncol = ncol(Y1))
        if (intercept == T) {
            X_test_1 <- cbind(1, X_test)
        }
        else {
            X_test_1 <- X_test
        }
        for (j in 1:length(eps_all)) {
            res_cv <- PPI(X1 = X_train, X2_c = X2_c, Y1 = Y_train, Y1_hat = Y_train_hat, Y2_hat = Y2_hat, 
                theta = theta, model = model, intercept = intercept, method = method, robust = robust, 
                CI = FALSE, res_classical = NULL, is_echo = is_echo, b_index = b_index, method_rme = method_rme, 
                alpha = alpha, const = const, eps = eps_all[j], method_opt = method_opt, iter_max = iter_max, 
                tol = tol, h = h)
            theta_hat <- res_cv$theta_hat
            if (model == "mean") {
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- mean(sweep(Y_test, 2, as.numeric(theta_hat), "-")^2)
            }
            else if (model == "linear") {
                Y_test_hat <- X_test_1 %*% theta_hat
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- mean((Y_test - Y_test_hat)^2)
            }
            else if (model == "logistic") {
                Y_test_hat_prob <- 1/(1 + exp(-(X_test_1 %*% theta_hat)))
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- -mean(c(Y_test) * log(apply(Y_test_hat_prob, 1, function(x) max(1e-09, 
                  x))) + c(1 - Y_test) * log(apply(1 - Y_test_hat_prob, 1, function(x) max(1e-09, x))))
            }
        }
    }
    obj_eps <- apply(obj, 1, mean)
    error_eps <- apply(error, 1, mean)
    eps_hat <- eps_all[which.min(error_eps)]
    if (is_echo == T) {
        out <- matrix(c(error_eps, obj_eps), nrow = 2, ncol = length(eps_all), byrow = T)
        colnames(out) <- eps_all
        rownames(out) <- c("error", "obj")
        print(out)
        print(paste("eps_hat = ", eps_hat))
    }
    return(list(eps_hat = eps_hat, error = error_eps, obj = obj_eps))
}

# Provide target parameters for simulation evaluation, using analytic values for mean and linear models.
# Approximate the logistic population target by fitting 100000 additional samples
target_est <- function(p, d, model, rho = 0.2) {
    N3 <- 1e+05
    data3 <- data_generator(n = N3, p = p, d = d, model = model, rho = rho)
    X3 <- data3$X
    Y3 <- data3$Y
    if (model == "mean") {
        theta <- rep(1, d)
    }
    else if (model == "linear") {
        theta <- c(0.3 * p, rep(1, p))
        if (p > 3) {
            theta <- theta/sqrt(p)
        }
    }
    else if (model == "logistic" | model == "logistic1") {
        model_large <- glm(Y3 ~ X3, family = quasibinomial(link = "logit"))
        theta <- model_large$coefficients
    }
    return(theta)
}

# Generate out-of-fold labeled predictions using K-fold cross-fitting and average unlabeled predictions across folds.
# Return Y1_hat, contaminated-feature predictions Y2_hat, and clean-feature predictions Y2_hat0; logistic predictions are thresholded at 0.5.
cross_prediction <- function(X1, Y1, X2, X2_c, model, Kfold = 5, num_round = 100) {
    n <- nrow(X1)
    p <- ncol(X1)
    N <- nrow(X2_c)
    if (is.vector(Y1) == 1) {
        d <- 1
        Y1 <- matrix(Y1, ncol = 1)
    }
    else {
        d <- ncol(Y1)
    }
    local_size <- n%/%Kfold
    sizes <- rep(local_size, Kfold)
    sizes[1] <- sizes[1] + n - sum(sizes)
    group <- sample(rep(1:Kfold, sizes), n, replace = F)
    Y1_hat <- matrix(Y1, nrow = n, ncol = d)
    Y2_hat <- matrix(0, nrow = N, ncol = d)
    Y2_hat0 <- matrix(0, nrow = N, ncol = d)
    for (i in 1:Kfold) {
        id_test <- which(group == i)
        id_train <- setdiff(c(1:n), id_test)
        X_train_i <- X1[id_train, , drop = FALSE]
        Y_train_i <- Y1[id_train, ]
        X_test_i <- X1[id_test, , drop = FALSE]
        if (model == "mean") {
            X_train_i_1 <- cbind(rep(1, length(id_train)), X_train_i)
            mu_beta_hat <- solve(t(X_train_i_1) %*% X_train_i_1) %*% t(X_train_i_1) %*% Y_train_i
            Y_test_i_hat <- cbind(rep(1, length(id_test)), X_test_i) %*% mu_beta_hat
            Y1_hat[id_test, ] <- Y_test_i_hat
            Y2_hat_i <- cbind(rep(1, N), X2_c) %*% mu_beta_hat
            Y2_hat <- Y2_hat + matrix(Y2_hat_i, ncol = d)/Kfold
            Y2_hat0_i <- cbind(rep(1, N), X2) %*% mu_beta_hat
            Y2_hat0 <- Y2_hat0 + matrix(Y2_hat0_i, ncol = d)/Kfold
        }
        else if (model == "linear") {
            tree <- fit_xgboost(X1 = X_train_i, Y1 = Y_train_i, model = model, num_round = 100)
            Y_test_i_hat <- predict(tree, newdata = xgb.DMatrix(X_test_i))
            Y1_hat[id_test, ] <- Y_test_i_hat
            Y2_hat_i <- predict(tree, newdata = xgb.DMatrix(X2_c))
            Y2_hat <- Y2_hat + matrix(Y2_hat_i, ncol = d)/Kfold
            Y2_hat0_i <- predict(tree, newdata = xgb.DMatrix(X2))
            Y2_hat0 <- Y2_hat0 + matrix(Y2_hat0_i, ncol = d)/Kfold
        }
        else if (model == "logistic") {
            tree <- fit_xgboost(X1 = X_train_i, Y1 = Y_train_i, model = model, num_round = 100)
            Y_test_i_hat_prob <- predict(tree, newdata = xgb.DMatrix(X_test_i))
            Y1_hat[id_test, ] <- Y_test_i_hat_prob
            Y2_hat_prob <- predict(tree, newdata = xgb.DMatrix(X2_c))
            Y2_hat <- Y2_hat + matrix(Y2_hat_prob, ncol = d)/Kfold
            Y2_hat0_prob <- predict(tree, newdata = xgb.DMatrix(X2))
            Y2_hat0 <- Y2_hat0 + matrix(Y2_hat0_prob, ncol = d)/Kfold
        }
    }
    if (model == "logistic") {
        Y1_hat <- matrix(ifelse(Y1_hat > 0.5, 1, 0), ncol = d)
        Y2_hat <- matrix(ifelse(Y2_hat > 0.5, 1, 0), ncol = d)
        Y2_hat0 <- matrix(ifelse(Y2_hat0 > 0.5, 1, 0), ncol = d)
    }
    return(list(Y1_hat = Y1_hat, Y2_hat = Y2_hat, Y2_hat0 = Y2_hat0))
}

# Generate predictions using fixed rules, sample splitting, or cross-fitting., as specified by method_train.
# Return the estimation samples X1 and Y1, plus predictions for labeled, contaminated unlabeled, and clean unlabeled data.
predict_func <- function(data_labeled, X2, X2_c, method_train, model, shift = -0.5, sig_f = 0.1, num_round = 100) {
    n <- nrow(data_labeled$X)
    p <- ncol(data_labeled$X)
    N <- nrow(X2)
    if (is.vector(data_labeled$Y) == 1) {
        d <- 1
        data_labeled$Y <- matrix(data_labeled$Y, ncol = 1)
    }
    else {
        d <- ncol(data_labeled$Y)
    }
    if (method_train == "splitting") {
        n0 <- floor(n * 0.3)
        X1 <- data_labeled$X[(n0 + 1):n, ]
        Y1 <- data_labeled$Y[(n0 + 1):n, ]
        X0 <- data_labeled$X[1:n0, ]
        Y0 <- data_labeled$Y[1:n0, ]
        if (model == "mean") {
            X0_1 <- cbind(rep(1, n0), X0)
            mu_beta_hat <- solve(t(X0_1) %*% X0_1) %*% t(X0_1) %*% Y0
            Y1_hat <- cbind(rep(1, n), X1) %*% mu_beta_hat
            Y2_hat <- cbind(rep(1, N), X2_c) %*% mu_beta_hat
            Y2_hat0 <- cbind(rep(1, N), X2) %*% mu_beta_hat
        }
        else if (model == "linear") {
            tree <- fit_xgboost(X1 = X0, Y1 = Y0, model = model, num_round = 100)
            Y1_hat <- predict(tree, newdata = xgb.DMatrix(X1))
            Y1_hat <- matrix(Y1_hat, ncol = d)
            Y2_hat <- predict(tree, newdata = xgb.DMatrix(X2_c))
            Y2_hat <- matrix(Y2_hat, ncol = d)
            Y2_hat0 <- predict(tree, newdata = xgb.DMatrix(X2))
            Y2_hat0 <- matrix(Y2_hat0, ncol = d)
        }
        else if (model == "logistic") {
            tree <- fit_xgboost(X1 = X0, Y1 = Y0, model = model, num_round = 100)
            Y1_hat_prob <- predict(tree, newdata = xgb.DMatrix(X1))
            Y1_hat <- matrix(ifelse(Y1_hat_prob > 0.5, 1, 0), ncol = d)
            Y2_hat_prob <- predict(tree, newdata = xgb.DMatrix(X2_c))
            Y2_hat <- matrix(ifelse(Y2_hat_prob > 0.5, 1, 0), ncol = d)
            Y2_hat0_prob <- predict(tree, newdata = xgb.DMatrix(X2))
            Y2_hat0 <- matrix(ifelse(Y2_hat0_prob > 0.5, 1, 0), ncol = d)
        }
    }
    else if (method_train == "crossfitting") {
        X1 <- data_labeled$X
        Y1 <- data_labeled$Y
        res_crossfitting <- cross_prediction(X1 = X1, Y1 = Y1, X2 = X2, X2_c = X2_c, model = model, Kfold = 5, 
            num_round = 100)
        Y1_hat <- res_crossfitting$Y1_hat
        Y2_hat <- res_crossfitting$Y2_hat
        Y2_hat0 <- res_crossfitting$Y2_hat0
    }
    else if (method_train == "fix") {
        X1 <- data_labeled$X
        Y1 <- data_labeled$Y
        if (model == "mean") {
            beta1 <- rep(1/p, p)
            beta2 <- c()
            for (i in 1:p) {
                beta2[i] <- (i - 1)/p
            }
            beta3 <- c(rep(1/2, floor(p/2)), rep(0, p - floor(p/2)))
            if (d == 1) {
                beta = beta1
            }
            else if (d == 2) {
                beta = cbind(beta1, beta2)
            }
            else if (d == 3) {
                beta = cbind(beta1, beta2, beta3)
            }
            # Fixed mean predictor: add biased random prediction error to the specified mean structure.
            pred_func <- function(x) {
                nrow <- dim(x)[1]
                epsilon <- matrix(rnorm(nrow * d, mean = -0.5, sd = 0.1), nrow, d)
                mu <- rep(1, d)
                Y <- t(t(x %*% beta + epsilon) + mu)
                return(Y)
            }
            Y1_hat <- pred_func(X1)
            Y2_hat <- pred_func(X2_c)
            Y2_hat0 <- pred_func(X2)
        }
        else if (model == "linear") {
            beta1 <- rep(1, p)
            beta2 <- rep(0.3, p)
            if (p > 3) {
                beta1 <- beta1/sqrt(p)
                beta2 <- beta2/sqrt(p)
            }
            # Fixed linear predictor: retain the linear and quadratic terms and add prediction error.
            pred_func <- function(x) {
                nrow <- dim(x)[1]
                epsilon <- rnorm(nrow, mean = shift, sd = sig_f)
                return(x %*% beta1 + epsilon + x^2 %*% beta2)
            }
            Y1_hat <- pred_func(X1)
            Y2_hat <- pred_func(X2_c)
            Y2_hat0 <- pred_func(X2)
        }
        else if (model == "logistic") {
            beta0 <- 4
            beta1 <- rep(0.5, p)
            beta2 <- rep(-1, p)
            # Fixed logistic predictor: omit the linear term and return logits for conversion to probabilities.
            pred_func <- function(x) {
                beta0 + x^2 %*% beta2
            }
            Y1_hat_prob <- 1/(1 + exp(-pred_func(X1)))
            Y1_hat <- matrix(ifelse(Y1_hat_prob > 0.5, 1, 0), ncol = d)
            Y2_hat_prob <- 1/(1 + exp(-pred_func(X2_c)))
            Y2_hat <- matrix(ifelse(Y2_hat_prob > 0.5, 1, 0), ncol = d)
            Y2_hat0_prob <- 1/(1 + exp(-pred_func(X2)))
            Y2_hat0 <- matrix(ifelse(Y2_hat0_prob > 0.5, 1, 0), ncol = d)
        }
    }
    return(list(X1 = X1, Y1 = Y1, Y1_hat = Y1_hat, Y2_hat = Y2_hat, Y2_hat0 = Y2_hat0))
}

# Compute mean logistic cross-entropy loss, its gradient, and the empirical Hessian.
# Return crit as the Euclidean gradient norm for convergence diagnostics.
gradient_logi <- function(X, theta, Y) {
    Y_hat <- 1/(1 + exp(-X %*% theta))
    loss <- -mean(c(Y) * log(apply(Y_hat, 1, function(x) max(1e-09, x))) + c(1 - Y) * log(apply(1 - Y_hat, 
        1, function(x) max(1e-09, x))))
    gradient <- t(X) %*% c(Y_hat - Y)/nrow(X)
    Hessian <- t(X) %*% diag(c(Y_hat * (1 - Y_hat))) %*% X/nrow(X)
    crit <- sqrt(sum(gradient^2))
    return(list(gradient = gradient, Hessian = Hessian, loss = loss, crit = crit))
}

# Compute logistic probabilities and per-sample gradients from X and theta.
# Return the gradient matrix, predicted probabilities Y_hat, and responses Y.
gradient_logi_mat <- function(X, theta, Y) {
    Y_hat <- 1/(1 + exp(-X %*% theta))
    gradient <- cross_entropy_loss_gradient_mat(X, Y, Y_hat)
    return(list(gradient = gradient, Y_hat = Y_hat, Y = Y))
}

# Compute per-sample cross-entropy gradients as X * (Y_hat - Y).
# For multiple response columns, compute each gradient block separately and concatenate the blocks by column.
cross_entropy_loss_gradient_mat <- function(X, Y, Y_hat) {
    res <- Y_hat - Y
    if (is.vector(res)) 
        res <- matrix(res, ncol = 1)
    col_y <- ncol(res)
    foreach(i = seq_len(col_y), .combine = cbind, .multicombine = TRUE) %do% X * res[, i]
}

# Remove entries outside each column's median +/- h * MAD interval, then average the remaining entries, where h is a fixed threshold
trimming_mean <- function(X, h) {
    p <- ncol(X)
    for (j in 1:p) {
        med_j <- median(X[, j])
        mad_j <- mad(X[, j])
        lower <- med_j - h * mad_j
        upper <- med_j + h * mad_j
        X[, j][X[, j] < lower | X[, j] > upper] <- NA
    }
    return(colMeans(X, na.rm = TRUE))
}

# Select a robust estimator using held-out labeled validation loss, with eps and h fixed.
# Return method_rme_choosen, validation losses, and diagnostic errors
CV_Meta <- function(seed, X1, Y1, X2, X2_c, Y1_hat, Y2_hat, theta, b_index, method, model, robust, method_train, 
    method_opt, alpha, const, method_rme_all, eps, h, intercept = T, fold = 5, tol = 0.01, iter_max = 20, 
    is_echo = F) {
    if (is.vector(Y1)) {
        Y1 <- matrix(Y1, nrow = length(Y1), ncol = 1)
    }
    if (is.vector(Y1_hat)) {
        Y1_hat <- matrix(Y1_hat, nrow = length(Y1_hat), ncol = 1)
    }
    if (is.vector(Y2_hat)) {
        Y2_hat <- matrix(Y2_hat, nrow = length(Y2_hat), ncol = 1)
    }
    n <- nrow(Y1)
    p <- ncol(X1)
    local_size <- n%/%fold
    sizes <- rep(local_size, fold)
    sizes[1] <- sizes[1] + n - sum(sizes)
    set.seed(seed)
    group <- sample(rep(1:fold, sizes), n, replace = F)
    error_method <- obj_method <- rep(0, length(method_rme_all))
    obj <- matrix(0, nrow = length(method_rme_all), ncol = fold)
    error <- matrix(0, nrow = length(method_rme_all), ncol = fold)
    for (i in 1:fold) {
        print(paste("Fold:", i))
        id_test <- which(group == i)
        id_train <- setdiff(c(1:n), id_test)
        X_train <- X1[id_train, , drop = FALSE]
        Y_train <- matrix(Y1[id_train, ], nrow = length(id_train), ncol = ncol(Y1))
        Y_train_hat <- matrix(Y1_hat[id_train, ], nrow = length(id_train), ncol = ncol(Y1))
        if (method_train == "crossfitting") {
            res_crossfitting <- cross_prediction(X1 = X_train, Y1 = Y_train, X2 = X2, X2_c = X2_c, model = model, 
                Kfold = 5, num_round = 100)
            Y2_hat <- res_crossfitting$Y2_hat
            Y_train_hat <- res_crossfitting$Y1_hat
        }
        X_test <- X1[id_test, , drop = FALSE]
        Y_test <- matrix(Y1[id_test, ], nrow = length(id_test), ncol = ncol(Y1))
        if (intercept == T) {
            X_test_1 <- cbind(1, X_test)
        }
        else {
            X_test_1 <- X_test
        }
        for (j in 1:length(method_rme_all)) {
            res_cv <- PPI(X1 = X_train, X2_c = X2_c, Y1 = Y_train, Y1_hat = Y_train_hat, Y2_hat = Y2_hat, 
                theta = theta, model = model, intercept = intercept, method = method, robust = robust, 
                CI = FALSE, res_classical = NULL, is_echo = is_echo, h = h, b_index = b_index, method_rme = method_rme_all[j], 
                alpha = alpha, const = const, eps = eps, method_opt = method_opt, iter_max = iter_max, 
                tol = tol)
            theta_hat <- res_cv$theta_hat
            if (model == "mean") {
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- mean(sweep(Y_test, 2, as.numeric(theta_hat), "-")^2)
            }
            else if (model == "linear") {
                Y_test_hat <- X_test_1 %*% theta_hat
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- mean((Y_test - Y_test_hat)^2)
            }
            else if (model == "logistic") {
                Y_test_hat_prob <- 1/(1 + exp(-(X_test_1 %*% theta_hat)))
                obj[j, i] <- sqrt(sum((theta_hat - theta)^2))
                error[j, i] <- -mean(c(Y_test) * log(apply(Y_test_hat_prob, 1, function(x) max(1e-09, 
                  x))) + c(1 - Y_test) * log(apply(1 - Y_test_hat_prob, 1, function(x) max(1e-09, x))))
            }
        }
    }
    obj_method <- apply(obj, 1, mean)
    error_method <- apply(error, 1, mean)
    method_rme_choosen <- method_rme_all[which.min(error_method)]
    if (is_echo == T) {
        out <- matrix(c(error_method, obj_method), nrow = 2, ncol = length(method_rme_all), byrow = T)
        colnames(out) <- method_rme_all
        rownames(out) <- c("error", "obj")
        print(out)
        print(paste("Method for RME choosen: ", method_rme_choosen))
    }
    return(list(method_rme_choosen = method_rme_choosen, error_method = error_method, obj_method = obj_method))
}

