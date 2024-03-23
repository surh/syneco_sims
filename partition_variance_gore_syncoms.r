# setwd("/Users/sur/lab/exp/2024/today2")
library(tidyverse)
library(lme4)
# library(brms)

#' Code to partition the variance into VPlas and VGen. Code from
#' https://github.com/devillemereuil/CodePartReacNorm from
#' 
#'   Villemereuil, Pierre de, and Luis-Miguel Chevin.
#'   “Partitioning the Phenotypic Variance of Reaction Norms,” September 1,
#'   2023. https://ecoevorxiv.org/repository/view/5891/.

#' Set global parameters
args <- list()
args$input <- "data_from_gore_for_models/day1_data.tsv"


#' Read data and plot
Dat <- read_tsv(args$input,
                col_types = cols(Day = col_number(),
                                 Rep = col_character(),
                                 OD = col_number(), 
                                 S0 = col_number(),
                                 ID = col_character(),
                                 Nutrient = col_number(),
                                 Community = col_character()))
Dat 

p1 <- Dat %>%
  ggplot(aes(x = Nutrient, y = OD,
             group = interaction(Community,Rep, sep = "_")) ) +
  geom_line(aes(col = factor(S0))) +
  scale_color_brewer(type = "seq", palette = "YlOrRd") + 
  AMOR::theme_blackbox() +
  theme(panel.background = element_rect(fill = "grey"))
p1

#' We also create a nutrient squared column for the polynomial models and
#' a factor version of nutrient for the character state model
Dat <- Dat %>%
  mutate(Nutrient_sq = Nutrient * Nutrient,
         Nutrient_f = factor(Nutrient))

#' # Model
#' We are going to try 3 basic polynomial models with degrees 0-2

mp0 <- lmer(OD ~ 1 + (1 | Community), data = Dat)
summary(mp0)

mp1 <- lmer(OD ~ 1 + Nutrient + (1 + Nutrient | Community), data = Dat)
summary(mp1)

mp2 <- lmer(OD ~ 1 + Nutrient + Nutrient_sq + (1 + Nutrient + Nutrient_sq | Community), data = Dat)
summary(mp2)

#' We also fit a character state model
mc <- lmer(OD ~ 0 + Nutrient_f + (0 + Nutrient_f | Community),
           data = Dat)
summary(mc)

#' Compare all modelos with AIC, BIC
AIC(mp0, mp1, mp2, mc)
BIC(mp0, mp1, mp2, mc)

#' ## Extract params from polynomial models
mp <- mp2

#' Design matrix
design_mat <- model.matrix(mp)

#' Variance-covariance matrix
varcov <- VarCorr(mp)[["Community"]]
attr(varcov, "stddev") <- attr(varcov, "correlation") <- NULL

# Average (fixed) effects
out <- summary(mp)[["coefficients"]]
coefs <- out[ , "Estimate"]
se    <- out[ , "Std. Error"]
  
# Residual variance
V_Res <- attr(VarCorr(mp), "sc")^2

# Variance-covariance of environment
varcov_design <- cov(design_mat)
  
# # Returning the list of parameters
# tibble(a    = coefs[1],
#        b    = coefs[2],
#        c    = coefs[3],
#        a_se = se[1],
#        b_se = se[2],
#        c_se = se[3],
#        Va   = varcov[1, 1],
#        Vb   = varcov[2, 2],
#        Vc   = varcov[3, 3],
#        Cab  = varcov[1, 2],
#        Cac  = varcov[1, 3],
#        Cbc  = varcov[2, 3],
#        V_R   = V_R) 


#'  ## Variance partitioning
#'  For a polynomial model we have an unbiased estimator for V_plas
#'  \hat{V}_{Plas} = \bar\theta^TX\theta - Tr(S_\theta X)

V_Plas <- coefs %*% varcov_design %*% coefs - se %*% varcov_design %*% se
V_Plas <- as.numeric(V_Plas)

#' For V_gen we have:
#' V_{Gen} = E_\epsilon(x^T\Theta x) = \bar{x}^T\Theta\bar{x} + Tr(\Theta X)
#' Is this equivalent?
V_Gen <- sum( diag( (1/nrow(design_mat)) * (t(design_mat) %*% design_mat) %*% varcov ) )

V_Tot <- V_Plas + V_Gen + V_Res
V_Phen <- var(Dat$OD)
V_Plas
V_Gen
V_Res
V_Tot
V_Phen


Pi <- (coefs^2 * diag(varcov_design) - se^2) / V_Plas
Pi
Gamma <- 2 * varcov[1, 3] * mean(design_mat[,3]) / V_Gen
Gamma











