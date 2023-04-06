###############################################################################
# Runner script for a simulation study which accepts command line arguments
###############################################################################

# Load libraries
library(here)
library(foreach)
library(doParallel)

# Unpack command line arguments
args <- commandArgs(trailingOnly = T)
if (length(args) > 0){
  dgp_num <- as.integer(args[1])
  sample_size <- as.integer(args[2])
  num_variables <- as.integer(args[3])
  kappa <- as.numeric(args[4])
  number_simulations <- as.integer(args[5])
  ate_true <- as.numeric(args[6])
  n_prop_submodels <- as.numeric(args[7])
  on_remote <- as.integer(args[8])
  datestamp <- args[9]
  run_number <- as.integer(args[10])
  estimated_propensities <- as.integer(args[11])
  residualize_xbcf <- as.integer(args[12])
  project_pi_yhat <- as.integer(args[13])
  script_iteration <- as.integer(args[14])
  n_yhat_submodels <- as.numeric(args[15])
  use_yhat_covariate <- as.integer(args[16])
  grf_default <- as.integer(args[17])
} else{
  # Default arguments
  dgp_num <- 1
  sample_size <- 500
  num_variables <- 50
  kappa <- 2.
  number_simulations <- 20
  ate_true <- 0.5
  n_prop_submodels <- 19
  on_remote <- 0
  datestamp <- format(Sys.time(), "%Y%m%d%H%M")
  run_number <- 1
  estimated_propensities <- 1
  residualize_xbcf <- 0
  project_pi_yhat <- 1
  script_iteration <- 1
  n_yhat_submodels <- 9
  use_yhat_covariate <- 1
  grf_default <- 1
}

# Define a function that runs one iteration of the simulation
simulation_function = function(n, p, ate_true, dgp_num, snr, n_prop_submodels, estimated_propensities, residualize_xbcf, project_pi_yhat, print_during_sim=F){
    ate_vector <- rep(NA, 7)
    rmse_vector <- rep(NA, 7)
    runtime_vector <- rep(NA, 7)

    # DGP-specific random components of the covariates, outcome and treatment
    if (dgp_num == 1){
        # Covariates
        x1 <- rnorm(n)
        x2 <- rbinom(n, 1, 0.5)
        x3 <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.3, 0.3))
        x4 <- rnorm(n)
        x5 <- rbinom(n, 1, 0.5)
        
        if (p > 5){
            x_rest <- matrix(rnorm(n*(p-5)), ncol = p-5)
            colnames(x_rest) <- paste0("x", 6:p)
            x <- cbind(x1, x2, x3, x4, x5, x_rest)
        } else {
            x <- cbind(x1, x2, x3, x4, x5)
        }
        
        # Prognostic function
        mu <- function(x) {
            lev <- c(-4, 4, 0)
            result <- 1 + 2. * x[, 1] * (2 * x[, 2] - 2 * (1 - x[, 2])) + lev[x[,3]]
            return(result)
        }
        
        # Treatment effect
        tau <- ate_true + 0.5 * x[, 1] + 0.5 * (2 * x[, 2] - 1)
        
        # Propensity score
        pi <- pnorm(-mean(mu(x)) + mu(x) - 2. * (2*x[, 5] - 1) + 2. * x[, 4], 0, 3)
        
        if (estimated_propensities == 0){
            # Propensity submodels
            # 1. Oracle instrument-free propensity
            grid_size <- 1000
            pi_no_instr <- rep(0, nrow(x))
            lev <- c(-4, 4, 0)
            for (k in 1:grid_size){
                temp_x <- x
                temp_x[,4] <- rnorm(n)
                temp_x[,5] <- rbinom(n,1,0.5)
                pi_no_instr <- pi_no_instr + pnorm(-mean(mu(temp_x)) + mu(temp_x) - 2. * (2*temp_x[, 5] - 1) + 2. * temp_x[, 4], 0, 3)
            }
            pi_no_instr <- pi_no_instr/grid_size
            
            # 2. Randomly-chosen submodels
            pi.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
            for (j in 1:n_prop_submodels){
                # Select a subset of covariates
                # First, choose the number of covariates to draw
                num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
                # Then, choose the covariates
                covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
                # Integrate out all variables not selected
                pi_temp <- rep(0, nrow(x))
                lev <- c(-4, 4, 0)
                for (k in 1:grid_size){
                    temp_x <- x
                    if (!(1 %in% covariates_sampled)){
                        temp_x[,1] <- rnorm(n)
                    }
                    if (!(2 %in% covariates_sampled)){
                        temp_x[,2] <- rbinom(n, 1, 0.5)
                    }
                    if (!(3 %in% covariates_sampled)){
                        temp_x[,3] <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.3, 0.3))
                    }
                    if (!(4 %in% covariates_sampled)){
                        temp_x[,4] <- rnorm(n)
                    }
                    if (!(5 %in% covariates_sampled)){
                        temp_x[,5] <- rbinom(n, 1, 0.5)
                    }
                    pi_temp <- pi_temp + pnorm(-mean(mu(temp_x)) + mu(temp_x) - 2. * (2*temp_x[, 5] - 1) + 2. * temp_x[, 4], 0, 3)
                }
                pi_temp <- pi_temp/grid_size
                pi.subset[,j] <- pi_temp
            }
        }
        
        # Arrange covariates according to the way XBCF distinguishes continuous and categorical features
        p_categorical_x <- 5
        x_orig <- x
        x_transformed <- data.frame(x)
        x_transformed[, 3] <- as.factor(x_transformed[, 3])
        x_transformed <- makeModelMatrixFromDataFrame(data.frame(x_transformed))
        cat_inds <- c(2, 3, 4, 5, 7)
        x_transformed <- cbind(x_transformed[, -cat_inds], x_transformed[, cat_inds])
        categorical_index <- c(rep(0, ncol(x_transformed)-length(cat_inds)), rep(1, length(cat_inds)))
    } else if (dgp_num == 2) {
        # Covariates
        x1 <- runif(n)
        x2 <- runif(n)
        x3 <- runif(n)
        x4 <- runif(n)
        x5 <- runif(n)
        x6 <- runif(n)
        x7 <- runif(n)
        x8 <- runif(n)
        
        if (p > 8){
            x_rest <- matrix(runif(n*(p-8)), ncol = p-8)
            colnames(x_rest) <- paste0("x", 9:p)
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x_rest)
        } else {
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8)
        }
        
        # Prognostic function
        mu <- function(x) {
            3*sin(2*pi*x[,1]) + 6*sin(2*pi*x[,2]) + 6*sin(2*pi*x[,3])
        }
        
        # Treatment effect
        tau <- ate_true + 2 * (x[, 4] - 0.5) + 0.5 * (x[, 5] - 0.5)
        
        # Propensity score
        pi <- pnorm(1*(x[,1]-0.5) + 1*(x[,4]-0.5) + 5*(x[,6]-0.5) + 5*(x[,7]-0.5) + 5*(x[,8]-0.5))
        
        if (estimated_propensities == 0){
            # Propensity submodels
            # 1. Oracle instrument-free propensity
            grid_size <- 1000
            pi_no_instr <- rep(0, nrow(x))
            for (k in 1:grid_size){
                temp_x <- x
                temp_x[,6] <- runif(n)
                temp_x[,7] <- runif(n)
                temp_x[,8] <- runif(n)
                pi_no_instr <- pi_no_instr + pnorm(1*(temp_x[,1]-0.5) + 1*(temp_x[,4]-0.5) + 5*(temp_x[,6]-0.5) + 5*(temp_x[,7]-0.5) + 5*(temp_x[,8]-0.5))
            }
            pi_no_instr <- pi_no_instr/grid_size
            
            # 2. Randomly-chosen submodels
            pi.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
            for (j in 1:n_prop_submodels){
                # Select a subset of covariates
                # First, choose the number of covariates to draw
                num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
                # Then, choose the covariates
                covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
                # Integrate out all variables not selected
                pi_temp <- rep(0, nrow(x))
                for (k in 1:grid_size){
                    temp_x <- x
                    if (!(1 %in% covariates_sampled)){
                        temp_x[,1] <- runif(n)
                    }
                    if (!(2 %in% covariates_sampled)){
                        temp_x[,2] <- runif(n)
                    }
                    if (!(3 %in% covariates_sampled)){
                        temp_x[,3] <- runif(n)
                    }
                    if (!(4 %in% covariates_sampled)){
                        temp_x[,4] <- runif(n)
                    }
                    if (!(5 %in% covariates_sampled)){
                        temp_x[,5] <- runif(n)
                    }
                    if (!(6 %in% covariates_sampled)){
                        temp_x[,6] <- runif(n)
                    }
                    if (!(7 %in% covariates_sampled)){
                        temp_x[,7] <- runif(n)
                    }
                    if (!(8 %in% covariates_sampled)){
                        temp_x[,8] <- runif(n)
                    }
                    pi_temp <- pi_temp + pnorm(1*(temp_x[,1]-0.5) + 1*(temp_x[,4]-0.5) + 5*(temp_x[,6]-0.5) + 5*(temp_x[,7]-0.5) + 5*(temp_x[,8]-0.5))
                }
                pi_temp <- pi_temp/grid_size
                pi.subset[,j] <- pi_temp
            }
        }
        
        # Arrange covariates according to the way XBCF distinguishes continuous and categorical features
        p_categorical_x <- 0
        x_orig <- x
        x_transformed <- x
        categorical_index <- rep(0, ncol(x_transformed))
    } else if (dgp_num == 3) {
        # Covariates
        x1 <- runif(n)
        x2 <- runif(n)
        x3 <- runif(n)
        x4 <- runif(n)
        x5 <- runif(n)
        x6 <- runif(n)
        x7 <- runif(n)
        x8 <- runif(n)
        
        if (p > 8){
            x_rest <- matrix(runif(n*(p-8)), ncol = p-8)
            colnames(x_rest) <- paste0("x", 9:p)
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x_rest)
        } else {
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8)
        }
        
        # Prognostic function
        mu <- function(x) {
            3*x[,1] - 3*x[,2]
        }
        
        # Treatment effect
        tau <- ate_true
        
        # Propensity score
        pi <- (1 + 2*x[,1])/4
        
        if (estimated_propensities == 0){
            # Propensity submodels
            # 1. Oracle instrument-free propensity
            pi_no_instr = pi
            
            # 2. Randomly-chosen submodels
            pi.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
            for (j in 1:n_prop_submodels){
                # Set grid size
                grid_size <- 1000
                # Select a subset of covariates
                # First, choose the number of covariates to draw
                num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
                # Then, choose the covariates
                covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
                # Integrate out all variables not selected
                pi_temp <- rep(0, nrow(x))
                for (k in 1:grid_size){
                    temp_x <- x
                    if (!(1 %in% covariates_sampled)){
                        temp_x[,1] <- runif(n)
                    }
                    if (!(2 %in% covariates_sampled)){
                        temp_x[,2] <- runif(n)
                    }
                    if (!(3 %in% covariates_sampled)){
                        temp_x[,3] <- runif(n)
                    }
                    if (!(4 %in% covariates_sampled)){
                        temp_x[,4] <- runif(n)
                    }
                    if (!(5 %in% covariates_sampled)){
                        temp_x[,5] <- runif(n)
                    }
                    if (!(6 %in% covariates_sampled)){
                        temp_x[,6] <- runif(n)
                    }
                    if (!(7 %in% covariates_sampled)){
                        temp_x[,7] <- runif(n)
                    }
                    if (!(8 %in% covariates_sampled)){
                        temp_x[,8] <- runif(n)
                    }
                    pi_temp <- pi_temp + (1 + 2*temp_x[,1])/4
                }
                pi_temp <- pi_temp/grid_size
                pi.subset[,j] <- pi_temp
            }
        }
        
        # Arrange covariates according to the way XBCF distinguishes continuous and categorical features
        p_categorical_x <- 0
        x_orig <- x
        x_transformed <- x
        categorical_index <- rep(0, ncol(x_transformed))
    } else if (dgp_num == 4){
        # Covariates
        x1 <- rnorm(n)
        x2 <- rbinom(n, 1, 0.5)
        x3 <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.3, 0.3))
        x4 <- rnorm(n)
        x5 <- rbinom(n, 1, 0.5)
        
        if (p > 5){
            x_rest <- matrix(rnorm(n*(p-5)), ncol = p-5)
            colnames(x_rest) <- paste0("x", 6:p)
            x <- cbind(x1, x2, x3, x4, x5, x_rest)
        } else {
            x <- cbind(x1, x2, x3, x4, x5)
        }
        
        # Prognostic function
        mu <- function(x) {
            lev <- c(-4, 4, 0)
            result <- 1 + 2. * x[, 1] * (2 * x[, 2] - 2 * (1 - x[, 2])) + lev[x[,3]]
            return(result)
        }
        
        # Treatment effect
        tau <- rep(ate_true, n)
        
        # Propensity score
        pi <- pnorm(-mean(mu(x)) + mu(x) - 2. * (2*x[, 5] - 1) + 2. * x[, 4], 0, 3)
        
        if (estimated_propensities == 0){
            # Propensity submodels
            # 1. Oracle instrument-free propensity
            grid_size <- 1000
            pi_no_instr <- rep(0, nrow(x))
            lev <- c(-4, 4, 0)
            for (k in 1:grid_size){
                temp_x <- x
                temp_x[,4] <- rnorm(n)
                temp_x[,5] <- rbinom(n,1,0.5)
                pi_no_instr <- pi_no_instr + pnorm(-mean(mu(temp_x)) + mu(temp_x) - 2. * (2*temp_x[, 5] - 1) + 2. * temp_x[, 4], 0, 3)
            }
            pi_no_instr <- pi_no_instr/grid_size
            
            # 2. Randomly-chosen submodels
            pi.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
            for (j in 1:n_prop_submodels){
                # Select a subset of covariates
                # First, choose the number of covariates to draw
                num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
                # Then, choose the covariates
                covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
                # Integrate out all variables not selected
                pi_temp <- rep(0, nrow(x))
                lev <- c(-4, 4, 0)
                for (k in 1:grid_size){
                    temp_x <- x
                    if (!(1 %in% covariates_sampled)){
                        temp_x[,1] <- rnorm(n)
                    }
                    if (!(2 %in% covariates_sampled)){
                        temp_x[,2] <- rbinom(n, 1, 0.5)
                    }
                    if (!(3 %in% covariates_sampled)){
                        temp_x[,3] <- sample(1:3, n, replace = TRUE, prob = c(0.3, 0.3, 0.3))
                    }
                    if (!(4 %in% covariates_sampled)){
                        temp_x[,4] <- rnorm(n)
                    }
                    if (!(5 %in% covariates_sampled)){
                        temp_x[,5] <- rbinom(n, 1, 0.5)
                    }
                    pi_temp <- pi_temp + pnorm(-mean(mu(temp_x)) + mu(temp_x) - 2. * (2*temp_x[, 5] - 1) + 2. * temp_x[, 4], 0, 3)
                }
                pi_temp <- pi_temp/grid_size
                pi.subset[,j] <- pi_temp
            }
        }
        
        # Arrange covariates according to the way XBCF distinguishes continuous and categorical features
        p_categorical_x <- 5
        x_orig <- x
        x_transformed <- data.frame(x)
        x_transformed[, 3] <- as.factor(x_transformed[, 3])
        x_transformed <- makeModelMatrixFromDataFrame(data.frame(x_transformed))
        cat_inds <- c(2, 3, 4, 5, 7)
        x_transformed <- cbind(x_transformed[, -cat_inds], x_transformed[, cat_inds])
        categorical_index <- c(rep(0, ncol(x_transformed)-length(cat_inds)), rep(1, length(cat_inds)))
    } else if (dgp_num == 5) {
        # Covariates
        x1 <- runif(n)
        x2 <- runif(n)
        x3 <- runif(n)
        x4 <- runif(n)
        x5 <- runif(n)
        x6 <- runif(n)
        x7 <- runif(n)
        x8 <- runif(n)
        
        if (p > 8){
            x_rest <- matrix(runif(n*(p-8)), ncol = p-8)
            colnames(x_rest) <- paste0("x", 9:p)
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x_rest)
        } else {
            x <- cbind(x1, x2, x3, x4, x5, x6, x7, x8)
        }
        
        # Prognostic function
        mu <- function(x) {
            3*sin(2*pi*x[,1]) + 6*sin(2*pi*x[,2]) + 6*sin(2*pi*x[,3])
        }
        
        # Treatment effect
        tau <- rep(ate_true, n)
        
        # Propensity score
        pi <- pnorm(1*(x[,1]-0.5) + 1*(x[,4]-0.5) + 5*(x[,6]-0.5) + 5*(x[,7]-0.5) + 5*(x[,8]-0.5))
        
        if (estimated_propensities == 0){
            # Propensity submodels
            # 1. Oracle instrument-free propensity
            grid_size <- 1000
            pi_no_instr <- rep(0, nrow(x))
            for (k in 1:grid_size){
                temp_x <- x
                temp_x[,6] <- runif(n)
                temp_x[,7] <- runif(n)
                temp_x[,8] <- runif(n)
                pi_no_instr <- pi_no_instr + pnorm(1*(temp_x[,1]-0.5) + 1*(temp_x[,4]-0.5) + 5*(temp_x[,6]-0.5) + 5*(temp_x[,7]-0.5) + 5*(temp_x[,8]-0.5))
            }
            pi_no_instr <- pi_no_instr/grid_size
            
            # 2. Randomly-chosen submodels
            pi.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
            for (j in 1:n_prop_submodels){
                # Select a subset of covariates
                # First, choose the number of covariates to draw
                num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
                # Then, choose the covariates
                covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
                # Integrate out all variables not selected
                pi_temp <- rep(0, nrow(x))
                for (k in 1:grid_size){
                    temp_x <- x
                    if (!(1 %in% covariates_sampled)){
                        temp_x[,1] <- runif(n)
                    }
                    if (!(2 %in% covariates_sampled)){
                        temp_x[,2] <- runif(n)
                    }
                    if (!(3 %in% covariates_sampled)){
                        temp_x[,3] <- runif(n)
                    }
                    if (!(4 %in% covariates_sampled)){
                        temp_x[,4] <- runif(n)
                    }
                    if (!(5 %in% covariates_sampled)){
                        temp_x[,5] <- runif(n)
                    }
                    if (!(6 %in% covariates_sampled)){
                        temp_x[,6] <- runif(n)
                    }
                    if (!(7 %in% covariates_sampled)){
                        temp_x[,7] <- runif(n)
                    }
                    if (!(8 %in% covariates_sampled)){
                        temp_x[,8] <- runif(n)
                    }
                    pi_temp <- pi_temp + pnorm(1*(temp_x[,1]-0.5) + 1*(temp_x[,4]-0.5) + 5*(temp_x[,6]-0.5) + 5*(temp_x[,7]-0.5) + 5*(temp_x[,8]-0.5))
                }
                pi_temp <- pi_temp/grid_size
                pi.subset[,j] <- pi_temp
            }
        }
        
        # Arrange covariates according to the way XBCF distinguishes continuous and categorical features
        p_categorical_x <- 0
        x_orig <- x
        x_transformed <- x
        categorical_index <- rep(0, ncol(x_transformed))
    }
    
    # Treatment assignment
    # hist(pi,100)
    z <- rbinom(n, 1, pi)
    
    # Outcome variable
    mu_x <- mu(x)
    Ey <- mu_x + tau * z
    sig <- snr * sd(Ey)
    y <- Ey + sig * rnorm(n)
    
    if (estimated_propensities == 1){
        # Fit a "full" propensity model using xgboost
        pihat <- estimate_pi_x(x, z)
        # plot(pi, pihat)
        
        # Randomly-chosen propensity submodels
        pihat.subset <- matrix(NA, nrow = n, ncol = n_prop_submodels)
        for (j in 1:n_prop_submodels){
            # Select a subset of covariates
            # First, choose the number of covariates to draw
            num_covariates_subset <- sample(1:(ncol(x)-1), size = 1)
            # Then, choose the covariates
            covariates_sampled <- sort(sample(1:ncol(x), size = num_covariates_subset, replace = F))
            # Fit a "subset" propensity model
            pihat.subset.temp <- estimate_pi_x(x[,covariates_sampled,drop=F], z)
            # plot(pi, pihat.subset.temp)
            pihat.subset[,j] <- pihat.subset.temp
        }
        
        # Define pi(X) covariates for prognostic and treatment models (the same in this case)
        pi_x_con = cbind(pihat, pihat.subset)
        pi_x_mod <- cbind(pihat, pihat.subset)
        pi_in_use <- pihat
        pi_in_use_xbcf <- pihat
    } else {
        # Define pi(X) covariates for prognostic and treatment models (the same in this case)
        pi_x_con = cbind(pi, pi.subset)
        pi_x_mod <- cbind(pi, pi.subset)
        pi_in_use <- pi
        pi_in_use_xbcf <- pi
    }
    
    # Define covariates for prognostic and treatment models (the same in this case)
    x_con <- x_transformed
    x_mod <- x_transformed
    
    #### 2. Model Fitting and Estimation
    
    #### Initial estimate of yhat for XBCF residualization purposes
    yhat.xbart.marginal.full <- estimate_yhat(
        x_con, y, p_categorical_x, num_sweeps = 60, burnin = 30, num_trees = 40
    )
    y_tilde <- y - yhat.xbart.marginal.full
    
    if (use_yhat_covariate == 1){
        #### Iterate through yhat submodels
        # Randomly-chosen propensity submodels
        yhat.subset <- matrix(NA, nrow = n, ncol = n_yhat_submodels)
        for (j in 1:n_yhat_submodels){
            # Select a subset of covariates
            # First, choose the number of covariates to draw
            num_covariates_subset <- sample(1:(ncol(x_con)-1), size = 1)
            # Then, choose the covariates
            covariates_sampled <- sort(sample(1:ncol(x_con), size = num_covariates_subset, replace = F))
            # Fit a "subset" propensity model
            p_categorical_subset <- sum(covariates_sampled > (ncol(x_con) - p_categorical_x))
            yhat.subset.temp <- estimate_yhat(
                x_con[,covariates_sampled,drop=F], y, p_categorical_subset, 
                num_sweeps = 60, burnin = 30, num_trees = 40
            )
            # plot(pi, pihat.subset.temp)
            yhat.subset[,j] <- yhat.subset.temp
        }
        
        # Add yhats to the covariates
        x_con <- cbind(x_con, yhat.xbart.marginal.full, yhat.subset)
        x_mod <- cbind(x_mod, yhat.xbart.marginal.full, yhat.subset)
    }
    
    #### Create a new outcome with yhat residualized out
    if (residualize_xbcf == 1){
        y_xbcf <- y_tilde
    } else {
        y_xbcf <- y
    }
    
    #### Project pi onto yhat
    if (project_pi_yhat == 1){
        project_model <- rpart(pi_in_use ~ yhat.xbart.marginal.full, control = rpart.control(minbucket = 3, cp=1e-4))
        pi_in_use_xbcf <- predict(project_model)
        pi_x_con[,1] = pi_in_use_xbcf
        pi_x_mod[,1] <- pi_in_use_xbcf
    }
    
    #### 2.a. Default XBCF
    
    # Run `num_sweeps` of the algorithm
    t1 = proc.time()
    xbcf.fit <- XBART::XBCF.discrete(
        y = y_xbcf, Z = z, X_con = x_con, X_mod = x_mod, pihat = pi_in_use_xbcf,
        p_categorical_con = p_categorical_x, p_categorical_mod = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, num_cutpoints = n, 
        include_pi = "both", verbose = F, num_trees_con = 30, num_trees_mod = 10, 
        alpha_con = 0.95, beta_con = 1.25, alpha_mod = 0.25, beta_mod = 3, 
        a_scaling = T, b_scaling = T, Nmin = 1
    )
    
    # Compute tauhat(X)
    pred_xbcf <- predict(xbcf.fit, X_con = x_con, X_mod = x_mod, Z = z,
                         pihat = pi_in_use, burnin = 30,
                         include_pi = "both")
    tauhats_xbcf <- pred_xbcf$tau.adj.mean

    # Evaluate RMSE and runtime
    ate_xbcf <- mean(tauhats_xbcf)
    rmse_xbcf <- sqrt(mean((tauhats_xbcf - tau)^2))
    runtime_xbcf <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBCF RMSE: ", rmse_xbcf))
        print(paste0("XBCF Runtime: ", runtime_xbcf, " seconds"))
    }

    rmse_vector[1] <- rmse_xbcf
    ate_vector[1] <- ate_xbcf
    runtime_vector[1] <- runtime_xbcf
    
    # Free memory from the XBCF fit
    rm(xbcf.fit)
    
    #### 2.b. XBCF with Multiple Propensities
    
    # Run `num_sweeps` of the algorithm
    t1 = proc.time()
    xbcf.fit.mps <- XBART::XBCF.discrete(
        y = y_xbcf, Z = z, X_con = x_con, X_mod = x_mod, pihat = pi_x_con,
        p_categorical_con = p_categorical_x, p_categorical_mod = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, num_cutpoints = n, 
        include_pi = "both", verbose = F, num_trees_con = 30, num_trees_mod = 10, 
        alpha_con = 0.95, beta_con = 1.25, alpha_mod = 0.25, beta_mod = 3, 
        a_scaling = T, b_scaling = T, Nmin = 1
    )
    t1 = proc.time() - t1
    
    # Compute tauhat(X)
    pred_xbcf_mps <- predict(xbcf.fit.mps, X_con = x_con, X_mod = x_mod, Z = z,
                             pihat = pi_x_con, burnin = 30,
                             include_pi = "both")
    tauhats_xbcf_mps <- pred_xbcf_mps$tau.adj.mean

    # Evaluate RMSE and runtime
    ate_xbcf_mps <- mean(tauhats_xbcf_mps)
    rmse_xbcf_mps <- sqrt(mean((tauhats_xbcf_mps - tau)^2))
    runtime_xbcf_mps <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("Multiple Propensity XBCF RMSE: ", rmse_xbcf_mps))
        print(paste0("Multiple Propensity XBCF Runtime: ", runtime_xbcf_mps, " seconds"))
    }
    
    rmse_vector[2] <- rmse_xbcf_mps
    ate_vector[2] <- ate_xbcf_mps
    runtime_vector[2] <- runtime_xbcf_mps
    
    # Free memory from the XBCF fit
    rm(xbcf.fit.mps)
    
    #### 2.c. GRF with known propensity and margin estimated by XBART
    
    # Fit causal forest with yhat estimated using XBART and pihat = pi_in_use
    t1 = proc.time()
    if ((grf_default == 1) & (estimated_propensities == 1)){
        causal.forest.grf <- causal_forest(
            X = x, Y = y, W = z, 
            Y.hat = NULL, 
            W.hat = NULL
        )
    } else if ((grf_default == 1) & (estimated_propensities == 0)) {
        causal.forest.grf <- causal_forest(
            X = x, Y = y, W = z, 
            Y.hat = NULL, 
            W.hat = pi_in_use
        )
    } else if (grf_default == 0){
        causal.forest.grf <- causal_forest(
            X = x, Y = y, W = z, 
            Y.hat = yhat.xbart.marginal.full, 
            W.hat = pi_in_use
        )
    }
    t1 = proc.time() - t1
    
    # Compute tauhat(X)
    tauhats_xbart_grf <- predict(causal.forest.grf)$predictions

    # Evaluate RMSE and runtime
    ate_xbart_grf <- mean(tauhats_xbart_grf)
    rmse_xbart_grf <- sqrt(mean((tauhats_xbart_grf - tau)^2))
    runtime_xbart_grf <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBART GRF RMSE: ", rmse_xbart_grf))
        print(paste0("XBART GRF Runtime: ", runtime_xbart_grf, " seconds"))
    }
    
    rmse_vector[3] <- rmse_xbart_grf
    ate_vector[3] <- ate_xbart_grf
    runtime_vector[3] <- runtime_xbart_grf

    # Free memory from the XBCF fit
    rm(causal.forest.grf)
    
    #### 2.d. XBART DR-Learner (Kennedy 2022)
    
    # Sample splitting
    fold_assigned <- sample(1:2, size = n, replace = T)
    fold_1_inds <- (1:n)[fold_assigned==1]
    fold_2_inds <- (1:n)[fold_assigned==2]
    x_1 <- x_transformed[fold_1_inds,,drop=F]
    x_2 <- x_transformed[fold_2_inds,,drop=F]
    x_untransformed_1 <- x[fold_1_inds,,drop=F]
    x_untransformed_2 <- x[fold_2_inds,,drop=F]
    z_1 <- z[fold_1_inds]
    z_2 <- z[fold_2_inds]
    y_1 <- y[fold_1_inds]
    y_2 <- y[fold_2_inds]
    pi_1 <- pi_in_use[fold_1_inds]
    pi_2 <- pi_in_use[fold_2_inds]
    
    # Stage 1: Estimate a propensity model, 
    # treated and control model on fold 1
    t1 = proc.time()
    
    # 1(a): Propensity score (xgboost)
    # Extract data info
    n_1 <- nrow(x_untransformed_1)
    data_inds_1 <- 1:n_1
    
    # Construct simple test train split
    train_inds_1 <- sort(sample(data_inds_1, floor(n_1*0.8), replace = F))
    test_inds_1 <- (data_inds_1)[!((data_inds_1) %in% train_inds_1)]
    dtrain <- xgb.DMatrix(x_untransformed_1[train_inds_1,,drop=F], label = z_1[train_inds_1], nthread = 1)
    dtest <- xgb.DMatrix(x_untransformed_1[test_inds_1,,drop=F], label = z_1[test_inds_1], nthread = 1)
    
    # Train model, with early stopping on the "test" set
    watchlist <- list(eval = dtest, train = dtrain)
    param <- list(max_depth = 2, eta = 0.05, nthread = 1, gamma = 1e-10)
    pihat.model.1 <- xgb.train(param, dtrain, nrounds = 50, watchlist, objective = "binary:logistic", early_stopping_rounds = 10)
    
    # Predict on fold 2
    pihat_fold2 <- predict(pihat.model.1, x_untransformed_2)
    
    # 1(b): Control and treated outcome model (XBART)
    xbart.drlearner.mu0.fold1 <- XBART::XBART(
        y = y_1[z_1==0], X = x_1[z_1==0,,drop=F], p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, 
        num_trees = 40, verbose = F
    )
    xbart.drlearner.mu1.fold1 <- XBART::XBART(
        y = y_1[z_1==1], X = x_1[z_1==1,,drop=F], p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, 
        num_trees = 40, verbose = F
    )
    
    # Predict on fold 2
    mu0_hat_drlearner_xbart_fold2 <- rowMeans(predict(xbart.drlearner.mu0.fold1, X = x_2, burnin = 30))
    mu1_hat_drlearner_xbart_fold2 <- rowMeans(predict(xbart.drlearner.mu1.fold1, X = x_2, burnin = 30))
    
    # Stage 2: construct the "pseudo-outcome" and estimate the CATE on fold 2
    treatment_component <- ((z_2-pihat_fold2)/(pihat_fold2*(1-pihat_fold2)))
    outcome_resid <- y_2 - z_2*mu1_hat_drlearner_xbart_fold2 - (1-z_2)*mu0_hat_drlearner_xbart_fold2
    counterfactual_diff <- mu1_hat_drlearner_xbart_fold2 - mu0_hat_drlearner_xbart_fold2
    pseudo_outcome_fold2 <- (
        treatment_component*outcome_resid + counterfactual_diff
    )
    
    # CATE model
    xbart.drlearner.cate.fold2 <- XBART::XBART(
        y = pseudo_outcome_fold2, X = x_2, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, 
        num_trees = 40, verbose = F
    )
    t1 = proc.time() - t1
    
    # Predict the CATE on the entire sample and store it
    tauhat_xbart_drlearner <- rowMeans(predict(xbart.drlearner.cate.fold2, X = x_transformed, burnin = 30))
    
    # Estimate ATE from the CATE
    ate_xbart_drlearner <- mean(tauhat_xbart_drlearner)
    
    # Evaluate RMSE and runtime
    rmse_xbart_drlearner <- sqrt(mean((tauhat_xbart_drlearner - tau)^2))
    runtime_xbart_drlearner <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBART DR-Learner RMSE: ", rmse_xbart_drlearner))
        print(paste0("XBART DR-Learner Runtime: ", runtime_xbart_drlearner, " seconds"))
    }
    
    rmse_vector[4] <- rmse_xbart_drlearner
    ate_vector[4] <- ate_xbart_drlearner
    runtime_vector[4] <- runtime_xbart_drlearner
    
    #### 2.e. XBART S-Learner
    
    # Run `num_sweeps` of the algorithm
    t1 = proc.time()
    x_combined = cbind(x_transformed, z)
    x_combined_trt = cbind(x_transformed, 1)
    x_combined_ctrl = cbind(x_transformed, 0)
    xbart.s.learner <- XBART::XBART(
        y = y, X = x_combined, p_categorical = p_categorical_x+1,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    t1 = proc.time() - t1
    
    # Compute tauhat(X)
    pred_xbart_slearner_1 <- predict(xbart.s.learner, X = x_combined_trt, burnin = 30)
    pred_xbart_slearner_0 <- predict(xbart.s.learner, X = x_combined_ctrl, burnin = 30)
    tauhats_xbart_slearner <- rowMeans(pred_xbart_slearner_1) - rowMeans(pred_xbart_slearner_0)
    
    # Evaluate RMSE and runtime
    ate_xbart_slearner <- mean(tauhats_xbart_slearner)
    rmse_xbart_slearner <- sqrt(mean((tauhats_xbart_slearner - tau)^2))
    runtime_xbart_slearner <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBART S-Learner RMSE: ", rmse_xbart_slearner))
        print(paste0("XBART S-Learner Runtime: ", runtime_xbart_slearner, " seconds"))
    }
    
    rmse_vector[5] <- rmse_xbart_slearner
    ate_vector[5] <- ate_xbart_slearner
    runtime_vector[5] <- runtime_xbart_slearner

    # Free memory from the XBCF fit
    rm(xbart.s.learner)
    
    #### 2.f. XBART T-Learner
    
    # Run `num_sweeps` of the algorithm
    t1 = proc.time()
    y_trt = y[z==1]
    y_ctrl = y[z==0]
    x_trt = x_transformed[z==1,]
    x_ctrl = x_transformed[z==0,]
    xbart.tlearner.f1 <- XBART::XBART(
        y = y_trt, X = x_trt, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    xbart.tlearner.f0 <- XBART::XBART(
        y = y_ctrl, X = x_ctrl, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    t1 = proc.time() - t1
    
    # Compute tauhat(X)
    pred_xbart_tlearner_f1 <- predict(xbart.tlearner.f1, X = x_transformed, burnin = 30)
    pred_xbart_tlearner_f0 <- predict(xbart.tlearner.f0, X = x_transformed, burnin = 30)
    tauhats_xbart_tlearner <- rowMeans(pred_xbart_tlearner_f1) - rowMeans(pred_xbart_tlearner_f0)
    
    # Evaluate RMSE and runtime
    ate_xbart_tlearner <- mean(tauhats_xbart_tlearner)
    rmse_xbart_tlearner <- sqrt(mean((tauhats_xbart_tlearner - tau)^2))
    runtime_xbart_tlearner <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBART T-Learner RMSE: ", rmse_xbart_tlearner))
        print(paste0("XBART T-Learner Runtime: ", runtime_xbart_tlearner, " seconds"))
    }
    
    rmse_vector[6] <- rmse_xbart_tlearner
    ate_vector[6] <- ate_xbart_tlearner
    runtime_vector[6] <- runtime_xbart_tlearner

    # Free memory from the XBCF fit
    rm(xbart.tlearner.f1, xbart.tlearner.f0)
    
    #### 2.g. XBART X-Learner
    
    # Run `num_sweeps` of the algorithm
    t1 = proc.time()
    y_trt = y[z==1]
    y_ctrl = y[z==0]
    x_trt = x_transformed[z==1,]
    x_ctrl = x_transformed[z==0,]
    xbart.xlearner.f1 <- XBART::XBART(
        y = y_trt, X = x_trt, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    xbart.xlearner.f0 <- XBART::XBART(
        y = y_ctrl, X = x_ctrl, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    
    # Predict f1 on control data and f0 on treated data
    pred_xbart_xlearner_f1 <- rowMeans(predict(xbart.xlearner.f1, X = x_ctrl, burnin = 30))
    pred_xbart_xlearner_f0 <- rowMeans(predict(xbart.xlearner.f0, X = x_trt, burnin = 30))
    
    # Model the residual predictions to "boost" estimates of tau(X) on treated and control groups
    xbart.xlearner.h1 <- XBART::XBART(
        y = y_trt-pred_xbart_xlearner_f0, X = x_trt, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    xbart.xlearner.h0 <- XBART::XBART(
        y = pred_xbart_xlearner_f1-y_ctrl, X = x_ctrl, p_categorical = p_categorical_x,
        num_sweeps = 60, burnin = 30, parallel = F, random_seed = 1234, 
        num_trees = 40, verbose = F
    )
    t1 = proc.time() - t1
    
    # Predict h1(X) and h0(X) on full dataset
    pred_xbart_xlearner_h1 <- predict(xbart.xlearner.h1, X = x_transformed, burnin = 30)
    pred_xbart_xlearner_h0 <- predict(xbart.xlearner.h0, X = x_transformed, burnin = 30)
    tauhats_xbart_xlearner <- pi_in_use*rowMeans(pred_xbart_xlearner_h0) + (1-pi_in_use)*rowMeans(pred_xbart_xlearner_h1)
    
    # Evaluate RMSE and runtime
    ate_xbart_xlearner <- mean(tauhats_xbart_xlearner)
    rmse_xbart_xlearner <- sqrt(mean((tauhats_xbart_xlearner - tau)^2))
    runtime_xbart_xlearner <- round(as.list(t1)$elapsed, 2)
    if (print_during_sim){
        print(paste0("XBART X-Learner RMSE: ", rmse_xbart_xlearner))
        print(paste0("XBART X-Learner Runtime: ", runtime_xbart_xlearner, " seconds"))
    }
    
    rmse_vector[7] <- rmse_xbart_xlearner
    ate_vector[7] <- ate_xbart_xlearner
    runtime_vector[7] <- runtime_xbart_xlearner

    # Free memory from the XBCF fit
    rm(xbart.xlearner.f1, xbart.xlearner.f0, xbart.xlearner.h1, xbart.xlearner.h0)
    
    # Output results as a row vector per simulation
    sim_vector <- c(rmse_vector, ate_vector, runtime_vector)
    return(sim_vector)
}

# Function to estimate pihat with hyperparameter tuning
estimate_pi_x <- function(X, Z, num_random_eval = 50, num_folds = 5){
    # Extract data info
    n <- nrow(X)
    data_inds <- 1:n
    
    # Construct simple test train split
    train_inds <- sort(sample(data_inds, floor(n*0.8), replace = F))
    test_inds <- (data_inds)[!((data_inds) %in% train_inds)]
    dtrain <- xgb.DMatrix(X[train_inds,,drop=F], label = Z[train_inds], nthread = 1)
    dtest <- xgb.DMatrix(X[test_inds,,drop=F], label = Z[test_inds], nthread = 1)
    
    # Prepare to train model
    watchlist <- list(eval = dtest, train = dtrain)
    param <- list(max_depth = 2, eta = 0.05, nthread = 1, gamma = 1e-10)
    pihat.model <- xgb.train(param, dtrain, nrounds = 50, watchlist, objective = "binary:logistic", early_stopping_rounds = 10)
    return(predict(pihat.model, X))
}

# Function to estimate yhat using XBART
estimate_yhat <- function(X, y, p_categorical_x, num_sweeps = 60, 
                          burnin = 30, num_trees = 40){
    xbart.marginal.y <- XBART::XBART(
        y = y, X = X, p_categorical = p_categorical_x,
        num_sweeps = num_sweeps, burnin = burnin, parallel = F, 
        num_trees = num_trees, verbose = F
    )
    yhat.xbart.marginal <- rowMeans(predict(xbart.marginal.y, X = X, burnin = burnin))
    return(yhat.xbart.marginal)
}

# Setup a parallel cluster
if (on_remote == 0){
  mc <- detectCores()
  cl <- makeCluster(mc)
} else if (on_remote == 1){
  cl <- makeCluster(strtoi(system("nproc",intern=TRUE)))
}
registerDoParallel(cl)

# Capture the project directory using "here"
project_dir = here()

# Create "outputs / snapshots / datestamp" subdirectory, if doesn't exist
snapshots_subfolder = file.path(project_dir, "outputs", "snapshots", datestamp)
ifelse(!dir.exists(file.path(project_dir, "outputs")), dir.create(file.path(project_dir, "outputs")), FALSE)
ifelse(!dir.exists(file.path(project_dir, "outputs", "snapshots")), dir.create(file.path(project_dir, "outputs", "snapshots")), FALSE)
ifelse(!dir.exists(snapshots_subfolder), dir.create(snapshots_subfolder), FALSE)
output_csv <- file.path(snapshots_subfolder, paste0("simulation_results_",run_number,"_",script_iteration,".csv"))

# Run the simulation
sim_results = foreach(i = 1:number_simulations, .combine = rbind) %dopar% {
    if (on_remote == 0){
        library("XBART")
        library("dbarts")
        library("grf")
        library("xgboost")
        library("rpart")
    } else if (on_remote == 1){
        library("XBART", lib="~/local/R_libs")
        library("dbarts", lib="~/local/R_libs")
        library("grf", lib="~/local/R_libs")
        library("xgboost", lib="~/local/R_libs")
        library("rpart", lib="~/local/R_libs")
    }

    simulation_function(sample_size, num_variables, ate_true, dgp_num, kappa, n_prop_submodels, estimated_propensities, residualize_xbcf, project_pi_yhat, F)
}
stopCluster(cl)

# Save the results
augmented_results <- data.frame(cbind(
    dgp_num, sample_size, num_variables, kappa, ate_true, 
    n_prop_submodels, estimated_propensities, residualize_xbcf, project_pi_yhat, 
    n_yhat_submodels, use_yhat_covariate, grf_default, sim_results
))
colnames(augmented_results) <- c(
  "DGP Number", "Sample Size", "Number of Variables", "Noise to Signal Ratio",
  "True ATE", "Number of Propensity Submodels", "Estimated Propensities", "Residualized XBCF", 
  "Project Pi Yhat", "Number of Outcome Submodels", "Use Yhat Covariate", "GRF Default",
  "CATE RMSE (XBCF)", "CATE RMSE (Multiple Propensity XBCF)", 
  "CATE RMSE (GRF)", "CATE RMSE (DR-Learner)", 
  "CATE RMSE (S-Learner)", "CATE RMSE (T-Learner)", "CATE RMSE (X-Learner)", 
  "ATE (XBCF)", "ATE (Multiple Propensity XBCF)", 
  "ATE (GRF)", "ATE (DR-Learner)", 
  "ATE (S-Learner)", "ATE (T-Learner)", "ATE (X-Learner)", 
  "Runtime (XBCF)", "Runtime (Multiple Propensity XBCF)", 
  "Runtime (GRF)", "Runtime (DR-Learner)", 
  "Runtime (S-Learner)", "Runtime (T-Learner)", "Runtime (X-Learner)"
)
write.csv(augmented_results, output_csv, row.names = F)
