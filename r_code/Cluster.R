library(bcf)
library(dplyr)

set.seed(42)
nburn <- 2000
nsim <- 2000
n <- 250

mu_types <- c("linear", "nonlinear")
tau_types <- c("homogeneous", "heterogeneous")

results_list <- list()
idx <- 1

for (mu_type in mu_types) {
  for (tau_type in tau_types) {
    cat(sprintf("Running: mu = %s | tau = %s\n", mu_type, tau_type))
    
    tryCatch({
      X1 <- rnorm(n)
      X2 <- rnorm(n)
      X3 <- rnorm(n)
      X4 <- rbinom(n, 1, 0.5)
      X <- data.frame(X1, X2, X3, X4)
      
      # tau
      tau <- 1 + 2 * X$X2 * X$X4
      if (tau_type == "homogeneous") tau <- rep(3, n)
      
      # mu
      g_X4 <- ifelse(X$X4 == 1, 2, -1)
      mu <- if (mu_type == "linear") {
        1 + g_X4 + X$X1 * X$X3
      } else {
        -6 + g_X4 + 6 * abs(X$X3 - 1)
      }
      
      s_mu <- sd(mu)
      u <- runif(n)
      pi <- pnorm(0.8 * ((3 * mu / s_mu) - 0.5 * X$X1) + 0.05 + u / 10)
      pi <- pmax(0, pmin(1, pi))
      Z <- rbinom(n, 1, pi)
      
      eps <- rnorm(n)
      Y <- mu + tau * Z + eps
      
      X_mat <- as.matrix(X)
      fit <- bcf(y = Y, z = Z,
                 x_control = X_mat,
                 x_moderate = X_mat,
                 pihat = pi,
                 nburn = nburn,
                 nsim = nsim)
      
      tau_hat <- colMeans(fit$tau)
      ate_hat <- mean(tau_hat)
      
      rmse_ate <- sqrt(mean((ate_hat - tau)^2))
      rmse_cate <- sqrt(mean((tau_hat - tau)^2))
      
      results_list[[idx]] <- data.frame(
        mu_type = mu_type,
        tau_type = tau_type,
        rmse_ate = rmse_ate,
        rmse_cate = rmse_cate
      )
      idx <- idx + 1
    }, error = function(e) {
      cat("⚠️ Error in fitting BCF for this scenario\n")
    })
  }
}

results <- do.call(rbind, results_list)
print(results)