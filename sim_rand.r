# setwd("~/syneco/exp/2024/today/")
# Simulate community reaction norms under mixed model from random communities

library(tidyverse)
library(lme4)

get_b <- function(x){
  c(x[2,1] - x[1,1],
    x[2,2] - x[1,2])
}


sim_syncom_rand <- function(n_sims = 10,
                            n_coms = 32){
  
  beta_env <- rnorm(n = 1, mean = 0, sd = 1)
  b_strains <- rnorm(n = 10, mean = 0, sd = 1)
  b_envstrain <- rnorm(n = 10, mean = 0, sd = 1)
  
  X <- matrix(nrow = 0, ncol = 10)
  # Sim communties
  while(nrow(X) < 32){
    com <- rep(0, 10)
    # ii <- sample(1:10, size = 5, replace = FALSE)
    ii <- which(sample(0:1, size = 10, replace = TRUE) == 1)
    com[ii] <- 1
    X <- rbind(X, com)
    
    if(nrow(X) == n_coms){
      X <- X[ !duplicated(X), ]
    }
  }
  
  colnames(X) <- paste0("b_strain_", 1:10)
  Dat <- as_tibble(X) %>%
    mutate(com = paste0("com_", 1:n_coms))
  Dat <- bind_rows(Dat, Dat, Dat) %>%
    mutate(env = rep(-1:1, each = n_coms))
  
  
  
  # Fixed and random coefficient matrices
  Beta <- matrix(c(0, beta_env), ncol = 1)
  B <- matrix(c(b_strains[1], b_envstrain[1],
                b_strains[2], b_envstrain[2],
                b_strains[3], b_envstrain[3],
                b_strains[4], b_envstrain[4],
                b_strains[5], b_envstrain[5],
                b_strains[6], b_envstrain[6],
                b_strains[7], b_envstrain[7],
                b_strains[8], b_envstrain[8],
                b_strains[9], b_envstrain[9],
                b_strains[10], b_envstrain[10]),
              ncol = 1)
  
  # Fixed effect design matrix
  X <- cbind(1, Dat$env)
  # X %*% Beta
  
  
  # Random effects design matrix
  Z <- cbind(Dat[,1], Dat[,1] * Dat[, 12],
             Dat[,2], Dat[,2] * Dat[, 12],
             Dat[,3], Dat[,3] * Dat[, 12],
             Dat[,4], Dat[,4] * Dat[, 12],
             Dat[,5], Dat[,5] * Dat[, 12],
             Dat[,6], Dat[,6] * Dat[, 12],
             Dat[,7], Dat[,7] * Dat[, 12],
             Dat[,8], Dat[,8] * Dat[, 12],
             Dat[,9], Dat[,9] * Dat[, 12],
             Dat[,10], Dat[,10] * Dat[, 12]) %>%
    as.matrix()
  # Z %*% B
  
  Sims <- as.numeric(X %*% Beta + Z %*% B) + sapply(rep(nrow(X), times = n_sims), rnorm, mean = 0, sd = 1)
  colnames(Sims) <- paste0("sim_", 1:n_sims)
  
  # Return results
  list(Dat = rbind(Dat),
       Sims = Sims,
       X = X,
       Z = Z,
       beta_env = beta_env,
       b_strains = b_strains,
       b_envstrain = b_envstrain)
}

extract_estimates_rand <- function(Sims){
  Res <- NULL
  for(i in 1:ncol(Sims$Sims)){
    
    Dat <- Sims$Dat
    Dat$y <- Sims$Sims[,i]
    
    m1 <- lmer(y ~ env + (env | b_strain_1) +
                 (env | b_strain_2) +
                 (env | b_strain_3) +
                 (env | b_strain_4) +
                 (env | b_strain_5) +
                 (env | b_strain_6) +
                 (env | b_strain_7) +
                 (env | b_strain_8) +
                 (env | b_strain_9) +
                 (env | b_strain_10),
               data = Dat)
    # summary(m1)
    # ranef(m1)$b_strain_1
    
    B_est <- rbind( get_b(ranef(m1)$b_strain_1),
                    get_b(ranef(m1)$b_strain_2),
                    get_b(ranef(m1)$b_strain_3),
                    get_b(ranef(m1)$b_strain_4),
                    get_b(ranef(m1)$b_strain_5),
                    get_b(ranef(m1)$b_strain_6),
                    get_b(ranef(m1)$b_strain_7),
                    get_b(ranef(m1)$b_strain_8),
                    get_b(ranef(m1)$b_strain_9),
                    get_b(ranef(m1)$b_strain_10))
    
    Res <- Res %>%
      bind_rows(bind_cols(b_strains_est = B_est[,1],
                          b_strains_true = Sims$b_strains,
                          b_envstrain_est = B_est[,2],
                          b_envstrain_true = Sims$b_envstrain) %>%
                  mutate(sim_id = colnames(Sims$Sims)[i]))
    
  }
  Res
    
}


n_sets <- 100
set.seed(39)
# Sims <- sim_syncom_rand(n_sims = 10, n_coms = 32)
# Res <- extract_estimates_rand(Sims)

Res <- NULL
for(i in 1:n_sets){
  Sims <- sim_syncom_rand(n_sims = 10)
  
  Res <- Res %>%
    bind_rows(extract_estimates_rand(Sims) %>%
                mutate(set_id = paste0("set_", i))) 
}
Res
write_tsv(Res, file = "rand_sims_est_vs_true.tsv")


p1 <- Res %>%
  ggplot(aes(x = b_strains_true, y = b_strains_est)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "blue", size = 2) +
  AMOR::theme_blackbox()
p1 
ggsave("rand_b_strains.png", width = 5, height = 5)

p1 <- Res %>%
  ggplot(aes(x = b_envstrain_true, y = b_envstrain_est)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "blue", size = 2) +
  AMOR::theme_blackbox()
p1
ggsave("rand_b_envstrain.png", width = 5, height = 5)


p1 <- Res %>%
  ggplot(aes(x = abs(b_strains_true - b_strains_est), group = set_id)) +
  geom_density() +
  AMOR::theme_blackbox()
p1
ggsave("rand_b_strain_error.png", width = 5, height = 5)

p1 <- Res %>%
  ggplot(aes(x = abs(b_envstrain_true - b_envstrain_est), group = set_id)) +
  geom_density() +
  AMOR::theme_blackbox()
p1
ggsave("rand_b_envstrain_error.png", width = 5, height = 5)









