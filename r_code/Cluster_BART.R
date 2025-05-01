library(dbarts)  # for BART modeling
library(dplyr)

set.seed(42)

# Simulation Parameters
nburn    <- 1000    # Burn-in iterations for BART
nsim     <- 2000    # Posterior draws after burn-in
num_rep  <- 50
sample_sizes <- c(250, 500)
mu_types     <- c("linear", "nonlinear")
tau_types    <- c("homogeneous", "heterogeneous")

results_list <- list()
idx <- 1

# Loop over simulation scenarios
for (n in sample_sizes) {
  for (mu_type in mu_types) {
    for (tau_type in tau_types) {
      for (rep in 1:num_rep) {
        tryCatch({
          # --- Step 1: Generate Covariates ---
          X1 <- rnorm(n)
          X2 <- rnorm(n)
          X3 <- rnorm(n)
          X4 <- rbinom(n, 1, 0.5)
          X5 <- sample(c(1, 2, 3), n, replace = TRUE)
          X <- data.frame(X1, X2, X3, X4, X5)
          
          # --- Step 2: Specify the Treatment Effect ---
          tau <- if (tau_type == "homogeneous") rep(3, n) else 1 + 2 * X$X2 * X$X5
          
          # --- Step 3: Define the Prognostic Function ---
          # Adjusted logic: X4 is binary (0,1)
          g_X4 <- ifelse(X$X4 == 1, 2, -1)
          mu <- if (mu_type == "linear") {
            1 + g_X4 + X$X1 * X$X3
          } else {
            -6 + g_X4 + 6 * abs(X$X3 - 1)
          }
          
          # --- Step 4: Calculate Propensity Scores and Treatment Assignment ---
          s_mu <- sd(mu)
          u <- runif(n)
          pi <- 0.8 * pnorm((3 * mu / s_mu) - 0.5 * X$X1) + 0.05 + u/10
          pi <- pmax(0, pmin(1, pi))
          Z <- rbinom(n, 1, pi)
          if (length(unique(Z)) < 2) stop("Treatment assignment has no variation")
          
          eps <- rnorm(n)
          Y <- mu + tau * Z + eps
          
          # Prepare the data set
          dat <- data.frame(Y, Z, X)
          
          # --- Step 5: Fit the BART Model ---
          # Set up the training data. The model is analogous to Y ~ Z * (X1 + X2 + X3 + X4 + X5)
          x.train <- as.matrix(dat[, c("Z", "X1", "X2", "X3", "X4", "X5")])
          y.train <- dat$Y
          
          # Fit BART using dbarts.
          # We use:
          #   - keeptrees = TRUE so that the fitted trees are saved (this is necessary for later prediction)
          #   - ntree = 200 trees as a robust default
          #   - nskip = nburn and ndpost = nsim to control the MCMC iterations
          bart_fit <- bart(x.train = x.train, 
                           y.train = y.train,
                           keeptrees = TRUE,
                           nskip = nburn,
                           ndpost = nsim,
                           ntree = 200,
                           verbose = FALSE)
          
          # --- Step 6: Compute Potential Outcomes Using Posterior Predictions ---
          # Create test data for both counterfactual scenarios.
          newdata_control <- data.frame(Z = 0, X1 = X$X1, X2 = X$X2, X3 = X$X3, X4 = X$X4, X5 = X$X5)
          newdata_treatment <- data.frame(Z = 1, X1 = X$X1, X2 = X$X2, X3 = X$X3, X4 = X$X4, X5 = X$X5)
          
          # Convert the test data to matrices
          x_control <- as.matrix(newdata_control)
          x_treatment <- as.matrix(newdata_treatment)
          
          # Generate posterior predictions.
          # The predict function will return a matrix with dimensions:
          #   (number of posterior draws) x (number of observations)
          pred_control <- predict(bart_fit, newdata = x_control)
          pred_treatment <- predict(bart_fit, newdata = x_treatment)
          
          # Compute treatment effect samples for each observation across posterior draws.
          # (Assuming the output matrices are: nsim x n)
          tau_samples <- pred_treatment - pred_control
          
          # Estimate the individual treatment effects (CATEs) as the posterior mean.
          tau_hat <- colMeans(tau_samples)
          ate_hat <- mean(tau_hat)
          ate_true <- mean(tau)
          
          # --- Step 7: Compute Evaluation Metrics ---
          rmse_ate <- sqrt((ate_hat - ate_true)^2)
          rmse_cate <- sqrt(mean((tau_hat - tau)^2))
          
          # 95% credible intervals for CATEs
          lower_cate <- apply(tau_samples, 2, quantile, probs = 0.025)
          upper_cate <- apply(tau_samples, 2, quantile, probs = 0.975)
          coverage_cate <- mean(tau >= lower_cate & tau <= upper_cate)
          length_cate <- mean(upper_cate - lower_cate)
          
          # 95% credible interval for ATE
          ate_samples <- colMeans(tau_samples)  # same as tau_hat
          ci_ate <- quantile(ate_samples, probs = c(0.025, 0.975))
          coverage_ate <- as.numeric(ate_true >= ci_ate[1] & ate_true <= ci_ate[2])
          length_ate <- ci_ate[2] - ci_ate[1]
          
          # Save simulation results for this iteration.
          results_list[[idx]] <- data.frame(
            rep = rep, n = n,
            mu_type = mu_type, tau_type = tau_type,
            rmse_ate = rmse_ate,
            rmse_cate = rmse_cate,
            coverage_ate = coverage_ate,
            coverage_cate = coverage_cate,
            length_ate = length_ate,
            length_cate = length_cate
          )
          idx <- idx + 1
          
        }, error = function(e) {
          cat("?? Error in iteration. Skipping...\n")
          # Uncomment the next line to print error details for debugging:
          # cat("Error details:", e$message, "\n")
        })
      }
    }
  }
}

# Combine all results into one data frame.
results <- do.call(rbind, results_list)

# Create a summary table by grouping over the simulation settings.
summary_results <- results %>% 
  group_by(n, mu_type, tau_type) %>% 
  summarise(
    ATE_RMSE    = mean(rmse_ate),
    ATE_Coverage = mean(coverage_ate),
    ATE_Length  = mean(length_ate),
    CATE_RMSE   = mean(rmse_cate),
    CATE_Coverage = mean(coverage_cate),
    CATE_Length   = mean(length_cate),
    .groups = "drop"
  )

print(summary_results)
save(summary_results, file = "summary_results.RData")
