library(tidyverse)
library(lme4)




Dat <- tibble(CCB1 = rep(c(0,1), each = 2),
              CCB2 = rep(c(0,1), times = 2)) %>%
  bind_rows(.,.) %>%
  mutate(CCB3 = rep(c(0,1), each = 4)) %>%
  bind_rows(.,.) %>%
  mutate(CCB4 = rep(c(0,1), each = 8)) %>%
  bind_rows(.,.) %>%
  mutate(CCB5 = rep(c(0,1), each = 16))
Dat %>% print(n = 100)

Dat <- bind_rows(Dat, Dat, Dat) %>%
  mutate(env = rep(-1:1, each = 32),
         com = rep(paste0("c",1:32), times = 3))
Dat %>% print(n = 100)

set.seed(942)
Dat$y <- rnorm(96)
Dat


m1 <- lmer(y ~ env + (1 + env|com), data = Dat)
summary(m1)
ranef(m1, which = "env")
plot(ranef(m1))
# m2 <- lmer(y ~ env + (1 + CCB1 + CCB2 | CCB1:CCB2), data = Dat)
# summary(m2)

# m2 <- lmer(y ~ env + (1 + CCB1 | CCB1) +
#              (1 + CCB2|CCB2) +
#              (1 + CCB3|CCB3) +
#              (1 + CCB4|CCB4) +
#              (1 + CCB5|CCB5), data = Dat)
# summary(m2)





Beta <- matrix(c(2, 2), ncol = 1)
B <- matrix(c(1, 0.8,
              -1.4, -0.9,
              1.2, 1.1,
              -2, -1.3,
              1.5, 0), ncol = 1)

# B <- matrix(c(1, 1.8,
#               -1.4, -1.9,
#               1.2, 3.1,
#               -2, -1.3,
#               1.5, 0), ncol = 1)
# B <- matrix(c(1, 0,
#               -1.4, 0,
#               1.2, 0,
#               -2, 0,
#               1.5, 0), ncol = 1)


X <- cbind(1, rep(Dat$env, times = 3))
# X %*% Beta

Z <- cbind(Dat[,1], Dat[,1] * Dat[, 6],
           Dat[,2], Dat[,2] * Dat[, 6],
           Dat[,3], Dat[,3] * Dat[, 6],
           Dat[,4], Dat[,4] * Dat[, 6],
           Dat[,5], Dat[,5] * Dat[, 6]) 
Z <- rbind(Z, Z, Z) %>% as.matrix()
# Z %*% B


set.seed(1329)
Dat2 <- rbind(Dat, Dat, Dat)
Dat2$y2 <- as.numeric(X %*% Beta + Z %*% B + rnorm(n = nrow(Dat2)))
Dat2

# ggplot(Dat2, aes(x = env, y = y2)) +
#   facet_wrap(~ com, ncol = 4) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ x) +
#   theme_classic()


p1 <- Dat2 %>%
  group_by(com, env) %>%
  summarise(y2 = mean(y2),
            .groups = 'drop') %>%
  ggplot(aes(x = env, y = y2, group = com)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              col = "blue", alpha = 0.2) +
  ylab(label = "phenotype") +
  theme_classic()
p1
ggsave("reaction_norm.png", p1, width = 6, height = 4)
ggsave("reaction_norm.svg", p1, width = 6, height = 4)

m1 <- lmer(y2 ~ env + ( 1 |com),
           data = Dat2)
summary(m1)
ranef(m1)

m2 <- lmer(y2 ~ env + ( env |com),
           data = Dat2)
summary(m2)
ranef(m2)

m3 <- lmer(y2 ~ env + (1 | CCB1) +
             (1 | CCB2) +
             (1 | CCB3) +
             (1 | CCB4) +
             (1 | CCB5), data = Dat2)
summary(m3)
ranef(m3)

m4 <- lmer(y2 ~ env + (env | CCB1) +
             (env | CCB2) +
             (env | CCB3) +
             (env | CCB4) +
             (env | CCB5), data = Dat2)
summary(m4)
ranef(m4)



m5 <- lmer(y2 ~ env + (env | CCB1) +
             (env | CCB2) +
             (env | CCB3) +
             (env | CCB4) +
             (1 | CCB5), data = Dat2)
summary(m5)
ranef(m5)


Zt <- attributes(m2)$pp$Zt
t(Zt)
head(Z)


AIC(m1,m2,m3,m4, m5)

# m2 <- nlme::lme(y2 ~ env, random = ~CCB1|CCB1, data = Dat2)
# summary(m2)
# ranef(m2)



# library(brms)
# 
# m3 <- brm(y2 ~ env + (CCB1 | CCB1) +
#             (CCB2 | CCB2) +
#             (CCB3 | CCB3) + 
#             (CCB4 | CCB4) +
#             (CCB5 | CCB5), data = Dat2, cores = 4)
# summary(m3)
# ranef(m3)


