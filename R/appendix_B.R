# required packages
library(tidyverse)
library(tidytext)
library(data.table)
library(ggpubr)
library(facetscales)
library(ggh4x)
library(grid)
library(sf)
library(doParallel)
library(graphics)
library(scales)

select <- dplyr::select

################
# Plot results #
################
results <- readRDS("data/simresults.rds")

## diagnostics -----------------------------------------------------------------
# divergent chains
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  ggplot() +
  geom_boxplot(aes(x = hrc_hra_lbl,
                   y = divergent/4000, 
                   col = N_lbl)) + 
  geom_hline(aes(yintercept = 0), linetype = "dotted", size = 1) +
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = "Proportion of divergent chains", col = "", x = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/pct_divergent_chains.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

# rhat
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  ggplot() +
  geom_boxplot(aes(x = hrc_hra_lbl,
                   y = rhat, 
                   col = N_lbl)) + 
  geom_hline(aes(yintercept = 1), linetype = "dotted", size = 1) +
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = "Rhat values", col = "", x = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/diagnostics/rhat.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

# ESS-tail and ESS-bulk
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  ggplot() +
  geom_boxplot(aes(x = hrc_hra_lbl,
                   y = ess_bulk, 
                   col = N_lbl)) + 
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = "ESS bulk", col = "", x = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/ess_bulk.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  ggplot() +
  geom_boxplot(aes(x = hrc_hra_lbl,
                   y = ess_tail, 
                   col = N_lbl)) + 
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = "ESS tail", col = "", x = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/ess_tail.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

# LOO-statistic
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  group_by(hrc_hra_lbl, N_lbl, speed, model) %>%
  summarize(loo_mu = mean(loo, na.rm = T),
            loo_se = mean(loo_se, na.rm = T), .groups = "keep") %>%
  group_by(model) %>%
  mutate(loo_avg = mean(loo_mu, na.rm = T)) %>%
  
  ggplot(aes(x = loo_mu, xmin = loo_mu - loo_se, xmax = loo_mu + loo_se, 
             y = N_lbl, group = model, col = model)) +
  geom_point(position = position_dodge(0.5), size = 3) + 
  geom_linerange(position = position_dodge(0.5), size = 1) +
  geom_vline(aes(xintercept = loo_avg, col = model), linetype = "dashed", size = 1) +
  facet_grid(hrc_hra_lbl~speed, labeller = label_parsed, scales = "free_y") +
  labs(x = "LOO", col = "Model", y = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  ) 

ggsave("fig/appendix_B/loo.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

# Bayesian p-value
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  ggplot(aes(x = hrc_hra_lbl, y = pval, col = N_lbl)) +
  geom_boxplot() + 
  geom_point(alpha = 0.1, position = position_jitterdodge(jitter.width = 2)) + 
  geom_hline(yintercept = 0.05, linetype = "dotted", size = 1) +
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = expression(bold(P[B])), col = "", x = "") +  
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/bayesian_pval.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  group_by(hrc_hra_lbl, N_lbl, speed, model) %>%
  summarize(prop_sig = sum(pval < 0.05)/ n()) %>%
  group_by(model) %>%
  mutate(prop_sig_avg = mean(prop_sig, na.rm = T)) %>%
  ggplot(aes(x = prop_sig, y = N_lbl, col = model)) +
  geom_point(position = position_dodge(0.5)) +
  geom_linerange(aes(xmax = prop_sig, col = model), xmin = 0, position = position_dodge(0.5)) +
  geom_vline(aes(xintercept = prop_sig_avg, col = model), linetype = "dashed", size = 1) +
  facet_grid(hrc_hra_lbl~speed, labeller = label_parsed, scales = "free_y") +
  labs(x = expression(bold(Proportion~of~simulations~with~P[B]<0.05)), col = "", y = "") +  
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  ) 

ggsave("fig/appendix_B/prop_pb_sig.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

# Chi-square discrepancy
results %>% 
  filter(method == "Stan", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  mutate(chat = ifelse(chat < quantile(results$chat, 0.85, na.rm = T), 
                       chat, quantile(results$chat, 0.85, na.rm = T))) %>%
  ggplot(aes(x = hrc_hra_lbl, y = chat, col = N_lbl)) +
  geom_boxplot() + 
  geom_point(alpha = 0.1, position = position_jitterdodge(jitter.width = 2)) + 
  geom_hline(yintercept = 1, linetype = "dotted", size = 1) +
  scale_x_discrete(labels = function(l) parse(text=l)) +
  facet_grid(model~speed) +
  labs(y = expression(bold(T^(s))), 
       col = "", x = "") +  
  scale_color_discrete(na.value = "grey") +
  #scale_color_grey() +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold", angle = 45, vjust = 0.5),
    strip.text = element_text(size = 14),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )

ggsave("fig/appendix_B/chat.jpg", 
       width = 24, height = 24, unit = "cm", dpi = 300, device = "jpg")

table_gof <- 
  results %>% 
  filter(N %in% c(10, 30, 50, 150, 300, 600), method == "Stan") %>%
  group_by(factor(hrc_hra), N = desc(N), speed, model) %>%
  summarize(pval = sum(pval < 0.05)/ n(),
            loo = mean(loo)) %>%
  pivot_wider(names_from = model, values_from = c(pval, loo)) %>%
  ungroup() %>% 
  bind_rows(summarise_all(., ~ if(is.numeric(.)) mean(.) else "Average"))

write_excel_csv(table_gof, "tab/gof.csv")
