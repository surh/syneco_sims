library(tidyverse)
apply(combn(x = 10, m = 5),2,function(ii){mat <- matrix(0,5,2); mat[ii] <- 1; mat})


dat <- combn(x = 10, m = 5) %>%
  as_tibble() %>%
  as.list %>%
  map_df(function(ii){
    mat <- matrix(0,5,2)
    mat[ii] <- 1
    mat <- as_tibble(mat)
    
    if(all(rowSums(mat) == 1)){
      mat$taxon <- 1:5
      mat$n_ccb <- sum(mat$V1)
      return(mat)
    }else{
      return(NULL)
    }
  }, .id = "id") 
dat

dat <- dat %>%
  pivot_longer(c(-id, -taxon, -n_ccb), names_to = "origin", values_to = "presence") %>%
  mutate(origin = replace(origin, origin == "V1", "CCB")) %>%
  mutate(origin = replace(origin, origin == "V2", "REF")) %>%
  mutate(id = factor(id, levels = unique(id[ order(n_ccb) ])))
dat$presence[ dat$presence == 1 & dat$origin == "REF" ] <- 2

p1 <- dat %>%
  ggplot(aes(x = origin, y = taxon)) +
  facet_wrap(~id, ncol = 8) + 
  geom_tile(aes(fill = as.factor(presence)), col = "black") + 
  scale_fill_manual(values = c("white"," red", "blue"),
                    name = "presence",
                    labels = c("none", "CCB", "REF")) +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 90))
p1
ggsave("all_swap_coms.png", p1, width = 5, height = 4)
ggsave("all_swap_coms.svg", p1, width = 5, height = 4)



p1 <- dat %>%
  split(.$id) %>%
  map_df(function(d){
    if(d$presence[ d$taxon == 5 & d$origin == "CCB" ] == 0){
      return(NULL)
    }else(
      return(d)
    )
  }) %>%
  ggplot(aes(x = origin, y = taxon)) +
  facet_wrap(~id, ncol = 8) + 
  geom_tile(aes(fill = as.factor(presence)), col = "black") + 
  scale_fill_manual(values = c("white"," red", "blue"),
                    name = "presence",
                    labels = c("none", "CCB", "REF")) +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.text.x = element_text(angle = 90))
p1
ggsave("CCB5_swap_coms.png", p1, width = 6, height = 3)
ggsave("CCB5_swap_coms.svg", p1, width = 6, height = 3)


    