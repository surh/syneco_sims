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











