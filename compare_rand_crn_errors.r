library(tidyverse)

Rand <- read_tsv("rand_sims_est_vs_true.tsv")
Rand

Crn <- read_table(file = "today2/crn_sims_est_vs_true (copia).tsv")
Crn



p1 <- Rand %>%
  mutate(type = "Random") %>%
  bind_rows(Crn %>% mutate(type = "F2-like")) %>%
  ggplot(aes(x = abs(b_strains_true - b_strains_est),
             group = interaction(type, set_id), col = type)) +
  geom_density() +
  xlab(label = "Absolute error in specific\nstrain estimates") +
  AMOR::theme_blackbox() +
  theme(legend.position = "bottom")
p1
ggsave("b_strain_error_comparison.png", width = 5, height = 5)



p1 <- Rand %>%
  mutate(type = "Random") %>%
  bind_rows(Crn %>% mutate(type = "F2-like")) %>%
  ggplot(aes(x = abs(b_envstrain_true - b_envstrain_est),
             group = interaction(type, set_id), col = type)) +
  geom_density() +
  xlab(label = "Absolute error in specific\nstrain by environment estimates") +
  AMOR::theme_blackbox() +
  theme(legend.position = "bottom")
p1
ggsave("b_env_error_comparison.png", width = 5, height = 5)

