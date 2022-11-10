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

## Simulated ground truths -----------------------------------------------------
# tabular format - detection parameters
truth_table <- results %>% filter(model %in% c("BP", "PP")) %>%
  filter(par == "p", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  group_by(hrc_hra = factor(hrc_hra_lbl), desc(N), speed, model) %>%
  summarise(p = mean(truth), min_p = min(truth), max_p = max(truth)) %>%
  pivot_wider(values_from = c(p, min_p, max_p), names_from = model) %>%
  mutate(p = paste0(round(min_p_BP, 4), " - ", round(max_p_BP, 4)),
         mu = paste0(round(min_p_PP, 4), " - ", round(max_p_PP, 4))) %>%
  select(hrc_hra:speed, p, mu)

write_excel_csv(truth_table, "tab/ground_truth_p.csv")

# tabular format - abundance and sit_use frequency
truth_table <- results %>% 
  filter(par == "lambda", model == "BernP", N %in% c(10, 30, 50, 150, 300, 600)) %>%
  group_by(factor(hrc_hra_lbl), N_lbl) %>%
  summarise(abundance = mean(truth), 
            min_su_freq = min(truth2), 
            max_su_freq = max(truth2))

write_excel_csv(truth_table, "tab/ground_truth_lambda.csv")

### Bias (figures) -------------------------------------------------------------

# GENERAL FUNCTION
# @params: dat       data.frame containing variables to plot
# @params: param     a string, containing the variable for which to plot the bias
# @params: bias_col  a string, containing the name of the column containing the bias
# @params: drop_hra  drop the largest home range area (TRUE)
#
# @returns: p        a ggplot, with facets displaying the bias for each simulation scenario

bias_facet_plot <- function(dat, param, bias_col, drop_hra = T) {
  
  fscale <- ifelse(param == "p", "free_y", "free")
  if(drop_hra) { dat <- dat %>% filter(hra != 91.61) }
  
  p <- dat %>% 
    
    # Filter predicates
    filter(par == param, ncam == 36,
           N %in% c(10, 30, 50, 150, 300, 600),
           truth > 0, divergent < 200, rhat < 1.02) %>%
    
    # Convert N to a factor
    mutate(N = factor(N, levels = rev(levels(factor(N))))) %>%
    
    # Calculate the average bias per model over all scenario's
    group_by(model) %>%
    mutate(avg_bias = mean(get(bias_col), na.rm = T)) %>%
    ungroup() %>%
    
    # Calculate mean and sd of the bias for all scenarios
    group_by(model, hrc, hra, hrc_hra_lbl, N, speed, avg_bias) %>%
    summarize(bias_mu = mean(get(bias_col), na.rm = T),
              bias_ll = quantile(get(bias_col), 0.025),
              bias_ul = quantile(get(bias_col), 0.975),
              .groups = "keep") %>%
    ungroup() %>%
    
    # Plot!
    ggplot(aes(x = bias_mu, y = N, col = bias_mu)) +
    geom_vline(xintercept = 0, linetype = "dotted", size = 0.75) +
    geom_vline(aes(xintercept = avg_bias, col = avg_bias), 
               linetype = "dashed", size = 0.75) +
    geom_point(aes(shape = speed), position = position_dodge(0.5), size = 2) +
    geom_linerange(aes(xmin = bias_ll, xmax = bias_ul, group = speed), 
                   size = 0.75, position = position_dodge(0.5)) +
    # scales
    scale_y_discrete(expand = expansion(add = c(-0.5, -0.5))) +
    scale_shape_manual(values = c(16, 17)) +
    # facets
    facet_nested(hrc + hra + N ~ model, drop = T, 
                 scales = fscale, labeller = label_parsed, switch = "y") +
    # theme
    guides(col = "none", shape = "none") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 11, face = "bold"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5, 
                                 margin = margin(5,40,5,40)),
      strip.text = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      panel.grid = element_blank()
    )
  
  if(param == "p") {
    p <- p + labs(y = "", x = expression(bold(Bias(theta[det])))) +
      scale_color_gradient2(low = "steelblue", mid = "grey", high = "firebrick")
  }
  
  if(param == "lambda") {
    
    scales_x <- list(
      `BernP` = scale_x_continuous(limits = c(-20, 10), breaks = seq(-20, 20, 10)),
      `BP` = scale_x_continuous(limits = c(-20, 65), breaks = seq(-20, 60, 20)),
      `PP` = scale_x_continuous(limits = c(-20, 10), breaks = seq(-20, 20, 10)),
      `BZIP` = scale_x_continuous(limits = c(-25, 105), breaks = seq(0, 100, 50)),
      `PZIP` = scale_x_continuous(limits = c(-25, 105), breaks = seq(0, 100, 50))
    )
    
    if(bias_col == "bias"){
      
      p <- p + labs(y = "", x = expression(bold(Bias(lambda)))) +
        scale_color_gradient2(low = "steelblue", mid = "grey", high = "firebrick",
                              limits = c(-20, 20), na.value = "firebrick") +
        facetted_pos_scales(x = scales_x)
      
    } else {
      
      p <- p + labs(y = "", x = expression(bold(Bias(lambda[use])))) +
        scale_color_gradient2(low = "steelblue", mid = "grey", high = "firebrick",
                              limits = c(-20, 20), na.value = "firebrick") +
        facetted_pos_scales(x = scales_x)
    }
  }
  return(p)
}

## Detection probability 
det <- bias_facet_plot(dat = results %>% filter(!grepl("ZIP", model)), 
                        param = "p", bias_col = "bias")

## Abundance
abun <- bias_facet_plot(dat = results %>% filter(!grepl("ZIP", model)),
                         param = "lambda", bias_col = "bias")

use_freq <- bias_facet_plot(dat = results %>% filter(!grepl("ZIP", model)),
                             param = "lambda", bias_col = "bias2")

## Relative abundance (trend analysis)
N_all <- unique(results$N)
trend_results <-
  foreach(i = seq(1,length(N_all), by = 3), .combine = "rbind") %do% {
    
    results_N <- results %>% 
      # filter per density
      filter(par == "lambda" & model != "BB",
             truth > 0, divergent < 200, rhat < 1.02,
             N %in% N_all[i:(i+2)]) %>%
      # drop unneeded columns
      select(-c(N_lbl, n1:p1, naive:CI_coverage2)) %>% 
      # pivot
      pivot_wider(-id, values_from = c(truth, estimate), 
                  names_from = N) %>%
      ungroup()
    
    meta <- results_N %>% select(method:sim_id)
    table_N <- results_N %>% select(-c(method:sim_id))
    
    table_N <-
      cbind(meta,
            foreach(j = c(1,4), .combine = "cbind") %do% {
              as_tibble(sapply(1:3, function(k){
                round(1 - table_N[,j:(j+2)][,k]/table_N[,j:(j+2)][,1], 2)})
              )
            }
      )
    
    table_N <- table_N[,-c(19,22)]
    colnames(table_N)[19:22] <- rep(c("trend10", "trend20"), 2)
    
    trends <- 
      cbind(
        pivot_longer(table_N[,19:20], cols = trend10:trend20, 
                     values_to = "truth", names_to = "trend"),
        pivot_longer(table_N[,21:22], cols = trend10:trend20, 
                     values_to = "estimate", names_to = "trend")[,2]
      )
    
    rep_row <- function(x) x[rep(seq_len(nrow(x)), each = 2), ]
    out <- cbind(rep_row(table_N[,-c(19:22)]),trends) 
    out <- out %>% 
      mutate(bias = estimate - truth,
             rel_bias = (estimate - truth)/ truth,
             Nref = N_all[i],
             Nref_lbl = paste0("N = ", Nref)) %>%
      select(method:sim_id, Nref, Nref_lbl, everything())
    
    
    return(out)
  }

# reorder factor levels
trend_results <- trend_results %>%
  mutate(Nref_lbl = factor(Nref_lbl, levels = paste0("N = ", sort(N_all))))

# facet labels
trend_labs <- c("10 % decrease", "20 % decrease")
names(trend_labs) <- c("trend10", "trend20")

trend_plot <- trend_results %>% 
  
  # Filter predicates
  filter(ncam == 36, Nref %in% c(10, 30, 50, 150, 300, 600),
         !grepl("ZIP", model), hra != 91.61) %>%
  
  # Convert N to a factor
  mutate(Nref = factor(Nref, levels = rev(levels(factor(Nref)))),
         trend = factor(trend, levels = c("trend10", "trend20"), 
                        labels = c("10%", "20%"))) %>%
  
  # Calculate the average bias per model over all scenario's
  group_by(model) %>%
  mutate(avg_bias = mean(bias, na.rm = T)) %>%
  ungroup() %>%
  
  # Calculate mean and sd of the bias for all scenarios
  group_by(model, hrc, hra, hrc_hra_lbl, Nref, speed, avg_bias) %>%
  summarize(bias_mu = mean(bias, na.rm = T),
            bias_ll = quantile(bias, 0.025, na.rm = T),
            bias_ul = quantile(bias, 0.975, na.rm = T),
            .groups = "keep") %>%
  ungroup() %>%
  
  # Plot!
  ggplot(aes(x = bias_mu, y = Nref, col = bias_mu)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.75) +
  geom_vline(aes(xintercept = avg_bias, col = avg_bias), 
             linetype = "dashed", size = 0.75) +
  geom_point(aes(shape = speed), position = position_dodge(0.5), size = 2) +
  geom_linerange(aes(xmin = bias_ll, xmax = bias_ul, group = speed), 
                 size = 0.75, position = position_dodge(0.5)) +
  labs(x = expression(Bias~(lambda-trend))) +
  # scales
  scale_x_continuous(limits = c(-0.2, 0.4), breaks = seq(-0.2, 0.4, 0.2)) +
  scale_y_discrete(expand = expansion(add = c(-0.5, -0.5))) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_gradient2(low = "steelblue", mid = "grey", high = "firebrick") +
  # facets
  facet_nested(hrc + hra + Nref ~ model, drop = T, 
               scales = "free", labeller = label_parsed, switch = "y") +
  # theme
  guides(col = "none", shape = "none") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 11, angle = 45, vjust = 0.5, 
                               margin = margin(5,40,5,40)),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.grid = element_blank()
  )

ggarrange(
  ggarrange(det + theme(axis.title.x = element_blank(),
                         plot.margin = unit(c(1,rep(0.2,3)), "cm")), 
            use_freq + theme(axis.title.x = element_blank(),
                              plot.margin = unit(c(1,rep(0.2,3)), "cm")),
            align = "v", ncol = 1, labels = c("(a)", "(c)")),
  ggarrange(abun + theme(axis.title.x = element_blank(),
                          strip.background.y = element_blank(),
                          strip.text.y = element_blank(),
                          plot.margin = unit(c(1,rep(0.2,3)), "cm")),
            trend_plot + theme(axis.title.x = element_blank(),
                               strip.background.y = element_blank(),
                               strip.text.y = element_blank(),
                               plot.margin = unit(c(1,rep(0.2,3)), "cm")),
            align = "v", ncol = 1, labels = c("(b)", "(d)")),
  align = "hv", widths = c(1.33, 1))


ggsave("fig/main/bias_facet_plot.jpg", 
       width = 16, height = 24, unit = "cm", dpi = 300, device = "jpg")

### Bias (tables) --------------------------------------------------------------

get_summary_table <- function(dat, param, use = F) {
  
  vars <- c("CI_coverage", "rmse", "rel_bias")
  if (use) { vars <- paste0(vars, "2") }
  
  dat %>% 
    filter(par == param, ncam == 36, N %in% c(10, 30, 50, 150, 300, 600),
           truth > 0, divergent < 200, rhat < 1.02) %>%
    group_by(hrc = factor(hrc), hra = factor(hra), N = desc(N), speed, model) %>%
    summarize(CI_coverage = ifelse(is.na(mean(get(vars[1]))), 
                                   0, mean(get(vars[1])*100)),
              rmse = mean(get(vars[2])),
              pctOK = sum(abs(get(vars[3])) <= 0.5)/ n(),
              .groups = "keep") %>%
    mutate(N = factor(abs(as.numeric(N)))) %>%
    pivot_wider(names_from = model, values_from = c(CI_coverage, rmse, pctOK)) %>%
    ungroup() %>% 
    bind_rows(summarise_all(., ~ if(is.numeric(.)) median(.) else "Median"))
}

## Detection probability 
table_p <- get_summary_table(results, "p")
write_excel_csv(table_p, "tab/summary_table_p.csv")

## Abundance
# traditional abundance
table_abun <- get_summary_table(results, "lambda")
write_excel_csv(table_abun, "tab/summary_table_abundance.csv")

# site use frequency 
table_use_freq <- get_summary_table(results, "lambda", use = T)
write_excel_csv(table_use_freq,"tab/summary_table_site-use_frequency.csv")

bias_table <- bind_cols(
  table_p %>% select(hrc:speed, starts_with("pctOK"), -ends_with("ZIP")),
  table_abun %>% select(starts_with("pctOK"), -ends_with("ZIP")),
  table_use_freq %>% select(starts_with("pctOK"), -ends_with("ZIP"))
  )
write_excel_csv(bias_table, "tab/summary_table_relbias.csv")

## Relative abundance (trend analysis)
table_trend <- 
  trend_results %>% 
  filter(par == "lambda", Nref %in% c(10, 30, 50, 150, 300, 600), 
         method == "Stan", !model %in% c("PZIP", "BZIP")) %>%
  group_by(hrc_hra = factor(hrc_hra), Nref = desc(Nref), speed, model, trend) %>%
  summarize(pctOK25 = sum(abs(rel_bias) <= 0.25)/ n(),
            pctOK50 = sum(abs(rel_bias) <= 0.50)/ n()) %>%
  mutate(Nref = factor(abs(as.numeric(Nref)))) %>%
  pivot_wider(names_from = c(model, trend), values_from = c(pctOK25, pctOK50)) %>%
  ungroup() %>% 
  select(hrc_hra:speed, ends_with("trend10"), ends_with("trend20")) %>%
  bind_rows(summarise_all(., ~ if(is.numeric(.)) median(.) else "Median"))

write_excel_csv(table_trend, "tab/summary_table_trend.csv")

### END Bias (tables) ----------------------------------------------------------