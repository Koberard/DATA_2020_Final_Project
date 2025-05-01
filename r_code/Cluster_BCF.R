library(bcf)
library(dplyr)

set.seed(42)

# Parameters
nburn <- 1000
nsim <- 1000
num_rep <- 50
sample_sizes <- c(250, 500)
mu_types <- c("linear", "nonlinear")
tau_types <- c("homogeneous", "heterogeneous")

results_list <- list()
idx <- 1


for (n in sample_sizes) {
  for (mu_type in mu_types) {
    for (tau_type in tau_types) {
      for (rep in 1:num_rep) {
        cat(sprintf("Rep %d | n = %d | mu = %s | tau = %s\n", rep, n, mu_type, tau_type))

        tryCatch({
          # Step 1: Generate covariates
          X1 <- rnorm(n)
          X2 <- rnorm(n)
          X3 <- rnorm(n)
          X4 <- rbinom(n, 1, 0.5)
          X5 <- sample(c(1,2,3), n, replace = TRUE)
          X <- data.frame(X1, X2, X3, X4,X5)

          # Step 2: Treatment effect
          tau <- if (tau_type == "homogeneous") rep(3, n) else 1 + 2 * X$X2 * X$X5

          # Step 3: Prognostic function
          g_X4 <- ifelse(X$X4 == -4, 2, -1)
          mu <- if (mu_type == "linear") {
            1 + g_X4 + X$X1 * X$X3
          } else {
            -6 + g_X4 + 6 * abs(X$X3 - 1)
          }

          # Step 4: Propensity & treatment assignment
          s_mu <- sd(mu)
          u <- runif(n)
          pi <- 0.8*pnorm(((3 * mu / s_mu) - 0.5 * X$X1)) + 0.05 + u / 10
          pi <- pmax(0, pmin(1, pi))
          Z <- rbinom(n, 1, pi)
          eps <- rnorm(n)
          Y <- mu + tau * Z + eps

          # Step 5: Fit BCF
          fit <- bcf(y = Y, z = Z,
                     x_control = as.matrix(X),
                     x_moderate = as.matrix(X),
                     pihat = pi,
                     nburn = nburn,
                     nsim = nsim)

          tau_samples <- fit$tau              # [nsim x n]
          tau_hat <- colMeans(tau_samples)    # CATE estimates
          ate_hat <- mean(tau_hat)            # ATE estimate
          ate_true <- mean(tau)

          # Step 6: Compute metrics
          rmse_ate <- sqrt((ate_hat - ate_true)^2)
          rmse_cate <- sqrt(mean((tau_hat - tau)^2))

          # CATE 95% intervals
          lower_cate <- apply(tau_samples, 2, quantile, probs = 0.025)
          upper_cate <- apply(tau_samples, 2, quantile, probs = 0.975)
          coverage_cate <- mean(tau >= lower_cate & tau <= upper_cate)
          length_cate <- mean(upper_cate - lower_cate)

          # ATE 95% interval
          ate_samples <- rowMeans(tau_samples)
          ci_ate <- quantile(ate_samples, probs = c(0.025, 0.975))
          coverage_ate <- as.numeric(ate_true >= ci_ate[1] & ate_true <= ci_ate[2])
          length_ate <- ci_ate[2] - ci_ate[1]

          # Save results
          results_list[[idx]] <- data.frame(
            rep = rep, n = n,
            mu_type = mu_type, tau_type = tau_type,
            rmse_ate, rmse_cate,
            coverage_ate, coverage_cate,
            length_ate, length_cate
          )
          idx <- idx + 1

        })
      }
    }
  }
}

# Combine all results
results <- do.call(rbind, results_list)

# Create summary table 
summary_results <- results %>%
  group_by(n, mu_type, tau_type) %>%
  summarise(
    ATE_RMSE = mean(rmse_ate),
    ATE_Coverage = mean(coverage_ate),
    ATE_Length = mean(length_ate),
    CATE_RMSE = mean(rmse_cate),
    CATE_Coverage = mean(coverage_cate),
    CATE_Length = mean(length_cate),
    .groups = "drop"
  )

print(summary_results)
save(summary_results, file = "summary_results_BCF.RData")
