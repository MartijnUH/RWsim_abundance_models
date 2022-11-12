# required packages
library(tidyverse)
library(data.table)
library(sf)
source("R/functions.R")
select <- dplyr::select

## set key
key <- 
  expand.grid(
    sim = 40,
    days = 25,
    hsteps = 6,
    r_hrc = c(F, T),
    hr_r = c(450, 900, 1800, 5400),
    speed = c("slow", "fast"),
    samp_prc = c(100, 50, 25),
    model = c("RN", "NMM", "PPM", "ZIPNM", "ZIPPM")
  ) %>% filter(r_hrc == F & hr_r == 450 | r_hrc == T) %>%
  mutate(
    Nmax = case_when(
      hr_r == 450 ~ 600,
      hr_r == 900 ~ 300,
      hr_r == 1800 ~ 150,
      hr_r == 5400 ~ 30),
    
    div = case_when(
      Nmax %in% c(600, 300) ~ 2, 
      Nmax %in% c(150, 30) ~ 3),
    
    Nmax0.9 = Nmax * 0.9, 
    Nmax0.8 = Nmax * 0.8,
    Nmax_div = Nmax / div,
    Nmax_div0.9 = (Nmax / div) * 0.9,
    Nmax_div0.8 = (Nmax / div) * 0.8) %>%
  
  select(-div) %>%
  pivot_longer(Nmax:Nmax_div0.8, names_to = "Nstep", values_to = "N") %>%
  mutate(id = 1:nrow(.)) %>%
  select(-Nstep) %>% select(id, sim:hr_r, N, everything())

# read survey data
ll <-  readRDS("data/rw/survey_infolist.rds")

# initialize output list
out <- list()

for(i in 1:nrow(key)){
  # print progress in key
  print(paste0("situation: ", i, " out of ", nrow(key), " situations"))
  
  if(file.exists(paste0("data/model_output/est_list", i, ".rds"))){
    ## read data
    
    # count histories + survey info
    n1 <- paste0("rHRC = ", key$r_hrc[i])
    n2 <- paste0("HR_r = ", key$hr_r[i])
    n3 <- paste0("speed = ", key$speed[i])
    n4 <- paste0("N = ", key$N[i])
    n5 <- paste0("sampling_prc = ", key$samp_prc[i])
    
    datl <- ll[[paste0(n1, "_", n2, "_", n3)]][[n4]][[n5]]
    
    # estimates retrieved from stanfits
    estl <- readRDS(paste0("data/model_output/est_list", i, ".rds"))
    
    ## store results for every simulation
    pars <- c("p", "psi", "lambda")
    npars <- ifelse(key$model[i] == "OCC", 2, 3)
    pars_array <- array(NA, dim = c(key$sim[i], 23, 3))
    
    for(s in 1:key$sim[i]){
      
      ctdat <- datl[[s]]$counts$unmarked
      # proportion of 1's
      prop1s <- mean(as.logical(ctdat))
      # number of 1's 
      n1s <- length(which(ctdat>0))
      
      for(p in 1:npars){
        
        ## add results to pars_array
        # simulation id
        pars_array[s, 1, p] <- s 
        
        # number of 1's
        pars_array[s, 2, p] <- length(which(ctdat>0))
        
        # proportion of 1's
        pars_array[s, 3, p] <- mean(as.logical(ctdat))
        
        # true parameter values
        pars_array[s, 4, p] <- 
          switch(p, `1` = ifelse(npars == 2, 
                                 datl[[s]]$p, 
                                 switch(as.character(key$model[i]),
                                        "RN" = datl[[s]]$r2,
                                        "NMM" = datl[[s]]$r2,
                                        "ZIPNM" = datl[[s]]$r2,
                                        "PPM" = datl[[s]]$mu,
                                        "ZIPPM" = datl[[s]]$mu)), 
                                      `2` = datl[[s]]$psi,
                                      `3` = datl[[s]]$lambda$marked)
        # estimated parameter values
        pars_array[s, 5, p] <- 
          #estl[[s]] %>% filter(parameter == pars[p]) %>% pull(mean) 
          estl[[s]]$est %>% filter(variable == pars[p]) %>% pull(mean) 
        
        # naive estimators
        pars_array[s, 6, p] <- get_naive(ctdat, pars[p])
        
        # 90% ci
        pars_array[s, 7:8, p] <- 
          c(estl[[s]]$est %>% filter(variable == pars[p]) %>% pull(q5),
            estl[[s]]$est %>% filter(variable == pars[p]) %>% pull(q95))    
        # bias
        pars_array[s, 10, p] <- pars_array[s, 5, p] - pars_array[s, 4, p] 
        # relative bias
        pars_array[s, 11, p] <- (pars_array[s, 10, p] / pars_array[s, 4, p])
        # squared errors
        pars_array[s, 12, p] <- (pars_array[s, 5, p] - pars_array[s, 4, p])^2 
        # HMC convergence diagnostics
        HMC_dat <- estl[[s]]$est %>% 
          filter(variable == pars[p]) %>% 
          select(rhat:walltime)
        
        pars_array[s, 14:19, p] <- 
          c(HMC_dat$ess_bulk, HMC_dat$ess_tail, HMC_dat$rhat, 
            HMC_dat$divergent, HMC_dat$max_treedepth, HMC_dat$walltime)
        
        # GOF
        pars_array[s, 20:21, p] <- 
          estl[[s]]$est %>% filter(variable %in% c("chat", "pval")) %>% pull(mean)
        
        # LOO
        pars_array[s, 22:23, p] <-
          estl[[s]]$loo %>% as.table() %>% as.data.frame() %>% 
          filter(Var1 == "elpd_loo") %>%
          pivot_wider(names_from = c("Var1", "Var2"), values_from = "Freq") %>%
          as_vector()
        
      }
    }
    for(p in 1:npars){
      # 95% CI coverage
      pars_array[, 9, p] <- mean((pars_array[, 7, p] <= pars_array[, 4, p]) & 
                                   (pars_array[, 8, p] >= pars_array[, 4, p]), 
                                 na.rm = T)               
      # root mean square error (RMSE)
      pars_array[, 13, p] <- sqrt(mean(pars_array[, 12, p]))
    }
    # round numbers
    pars_array <- apply(pars_array, c(1,2,3), round, 4) 
    
    results.df <- 
      rbind.data.frame(data.frame(par = rep("p", length(estl)), pars_array[,,1]),
                       data.frame(par = rep("psi", length(estl)), pars_array[,,2]),
                       data.frame(par = rep("lambda", length(estl)), pars_array[,,3])
      ) %>%
      mutate(id = i,
             method = "Stan", 
             N = key[i,]$N,
             speed = key[i,]$speed,
             days = key[i,]$days,
             hourly_steps = key[i,]$hsteps,
             random_hr_center = key[i,]$r_hrc,
             hr_r = key[i,]$hr_r, 
             cam_grid = "g12",
             cam_sites = 144,
             ncam = (key[i,]$samp_prc/100) * 144, 
             model = key[i,]$model,
             estimate = "mean") %>% 
      select(id, method, N, speed, days, hourly_steps, random_hr_center, hr_r, 
             cam_grid, cam_sites, ncam, model, estimate, everything())
    
    colnames(results.df) <- 
      c("id", "method", "N", "speed", "days", "hourly_steps", "random_hr_center", 
        "hr_radius", "cam_grid", "cam_sites", "ncam", "model", "point_estimate", 
        "par", "sim_id", "n1", "p1", "truth", "estimate", "naive", "5%", "95%", 
        "CI_coverage", "bias", "rel_bias", "se", "rmse", "ess_bulk", "ess_tail", 
        "rhat", "divergent", "max_treedepth", "walltime", "chat", "pval", "loo", 
        "loo_se")
    
    out[[i]] <- results.df
  }
}

results <- rbindlist(out)

# add variables of interest
results <- results %>%
  mutate(
    
    par = factor(par, levels = c("p", "psi", "lambda")),
    model = factor(model, levels = c("RN", "NMM", "PPM", "ZIPNM", "ZIPPM"),
                   labels = c("BernP", "BP", "PP", "BZIP", "PZIP")),
    N_lbl = factor(N, levels = as.numeric(levels(factor(N))), 
                   labels = paste0("N = ", levels(factor(N)))),
    hrc = factor(random_hr_center, levels = c(F, T), 
                 labels = c("closure", "non-closure")),
    hra = factor(round(hr_radius^2*pi/1000^2,2)),
    hra_lbl = factor(hra,
                     labels = sapply(1:length(levels(hra)), function(x){
                       
                       tmp <- as.character(levels(hra))[x]
                       deparse(substitute(HRA == hra ~ km^2, list(hra = tmp)))
                       
                     })),
    hrc_hra = interaction(hrc, hra, sep = "_", drop = T),
    hrc_hra_lbl = factor(hrc_hra, 
                         labels = sapply(1:length(levels(hrc_hra)), function(x){
                           
                           # split the hrc_hra string
                           tmp <- strsplit(as.character(levels(hrc_hra)), "_")
                           
                           # store the hra and hrc characters
                           hrc <- ifelse(tmp[[x]][1] == "centroid", 
                                         "(closure)", "(non-closure)")
                           hra <- tmp[[x]][2]
                           
                           # return plotmath expression with two lines
                           return(deparse(substitute(atop(HRA == hra ~ km^2, hrc), 
                                                     list(hra = hra, hrc = hrc))))
                         })
    )
  ) %>%
  select(id:N, N_lbl, speed:hr_radius, hrc:hra_lbl, hrc_hra:hrc_hra_lbl, everything())


## add asymptotic truth and calculate bias, rel. bias, se, rmse, CI coverage
asymptotic_truth <- lapply(1:nrow(key), function(i){
  
  # count histories + survey info
  n1 <- paste0("rHRC = ", key$r_hrc[i])
  n2 <- paste0("HR_r = ", key$hr_r[i])
  n3 <- paste0("speed = ", key$speed[i])
  n4 <- paste0("N = ", key$N[i])
  n5 <- paste0("sampling_prc = ", key$samp_prc[i])
  
  datl <- ll[[paste0(n1, "_", n2, "_", n3)]][[n4]][[n5]]
  lam <- sapply(1:key$sim[i], function(j) datl[[j]]$lambda$unmarked)
  data.frame(id = rep(i, key$sim[i]), sim_id = 1:key$sim[i], 
             par = "lambda", asymptotic_truth = lam)
}) %>% rbindlist()

results <- 
  left_join(results, asymptotic_truth, c("id", "sim_id", "par")) %>%
  rename(truth2 = asymptotic_truth) %>% 
  mutate(bias2 = estimate - truth2,
         rel_bias2 = (estimate - truth2)/ truth2,
         se2 = (estimate - truth2)^2,
         temp = ifelse(`5%` <= truth2 & truth2 <= `95%`, 1, 0)) %>%
  group_by(id, par) %>%
  mutate(rmse2 = sqrt(mean(se2)),
         CI_coverage2 = mean(temp)) %>%
  select(-temp)

saveRDS(results, "data/simresults.rds")
