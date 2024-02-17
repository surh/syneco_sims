# setwd("/Users/sur/lab/exp/2024/today2")
# Simulate community reaction norms under mixed model

library(tidyverse)
library(lme4)


sim_syncom_f2 <- function(n_sims = 10,
                          beta_env = NULL,
                          b_strains = NULL,
                          b_envstrain = NULL){
  
  # Simulate params if not passed
  # Missing checks
  if( is.null(beta_env) ){
    beta_env <- rnorm(n = 1, mean = 0, sd = 1)
  }
  if( is.null(b_strains) ){
    b_strains <- rnorm(n = 5, mean = 0, sd = 1)
  }
  if( is.null(b_envstrain) ){
    b_envstrain <- rnorm(n = 5, mean = 0, sd = 1)
  }
  
  
  # Generate design
  Dat <- tibble(CCB1 = rep(c(0,1), each = 2),
                CCB2 = rep(c(0,1), times = 2)) %>%
    bind_rows(.,.) %>%
    mutate(CCB3 = rep(c(0,1), each = 4)) %>%
    bind_rows(.,.) %>%
    mutate(CCB4 = rep(c(0,1), each = 8)) %>%
    bind_rows(.,.) %>%
    mutate(CCB5 = rep(c(0,1), each = 16))
  # Dat %>% print(n = 100)
  
  Dat <- bind_rows(Dat, Dat, Dat) %>%
    mutate(env = rep(-1:1, each = 32),
           com = rep(paste0("c",1:32), times = 3))
  # Dat %>% print(n = 100)
  
  # Fixed and random coefficient matrices
  Beta <- matrix(c(0, beta_env), ncol = 1)
  B <- matrix(c(b_strains[1], b_envstrain[1],
                b_strains[2], b_envstrain[2],
                b_strains[3], b_envstrain[3],
                b_strains[4], b_envstrain[4],
                b_strains[5], b_envstrain[5]),
              ncol = 1)
  
  # Fixed effect design matrix
  X <- cbind(1, rep(Dat$env, times = 3))
  # X %*% Beta
  
  # Random effects design matrix
  Z <- cbind(Dat[,1], Dat[,1] * Dat[, 6],
             Dat[,2], Dat[,2] * Dat[, 6],
             Dat[,3], Dat[,3] * Dat[, 6],
             Dat[,4], Dat[,4] * Dat[, 6],
             Dat[,5], Dat[,5] * Dat[, 6]) 
  Z <- rbind(Z, Z, Z) %>% as.matrix()
  # Z %*% B
  
  
  # Make simulation
  Sims <- as.numeric(X %*% Beta + Z %*% B) + sapply(rep(nrow(X), times = n_sims), rnorm, mean = 0, sd = 1)
  colnames(Sims) <- paste0("sim_", 1:n_sims)
  # Sims
  
  # Return results
  list(Dat = rbind(Dat, Dat, Dat),
       Sims = Sims,
       X = X,
       Z = Z,
       beta_env = beta_env,
       b_strains = b_strains,
       b_envstrain = b_envstrain)
  
}

sim_syncom_f2()














list()