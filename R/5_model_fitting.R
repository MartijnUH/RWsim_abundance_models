## Load R packages
library(readr) 
library(dplyr) 
library(tidyr)
library(posterior)
library(tidybayes)
library(cmdstanr)
library(parallel)
select <- dplyr::select
options(mc.cores=parallel::detectCores()-1)

## input arguments
#args <- commandArgs(TRUE) # Comment out if not working on HPC
#i <- as.integer(args[1])  # Comment out if not working on HPC

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

for (i in 1:nrow(key)) { # Use loop if not working on HPC

## input files
# list of count data
n1 <- paste0("rHRC = ", key$r_hrc[i])
n2 <- paste0("HR_r = ", key$hr_r[i])
n3 <- paste0("speed = ", key$speed[i])
n4 <- paste0("N = ", key$N[i])
n5 <- paste0("sampling_prc = ", key$samp_prc[i])

datalist <- readRDS(paste0("data/rw/hpc/chlist_", 
                           n1, "_", n2, "_", n3, "_", n4, "_", n5, ".rds"))

if(key$model[i] == "RN"){
  for(l in 1:length(datalist)){
    datalist[[l]][datalist[[l]]>0] <- 1
  }
}

# stan model
mod_file <- paste0("stan/", as.character(key$model[i]), "_multithreaded.stan")
mod <- cmdstan_model(stan_file = mod_file, 
                     cpp_options = list(stan_threads = TRUE))

## fit stan model to simulation data
time <- system.time({
  est.list <-
    lapply(1:40, function(s){

      # print progress in fitting to model for each simulation
      print(paste0("fit stan model for simulation round: ", 
                   s, " out of 40 simulations"))

      # Select appropriate detection history
      dat <- datalist[[s]]

      # drop rows with all NA's from the detection history
      dat <- dat[rowSums(is.na(dat)) != ncol(dat), ]

      # prepare input list to stanmodel
      standata <- list(
              R = nrow(dat),
              T = ncol(dat),
              y = dat)
      
      if(as.character(key$model[i]) == "RN"){
        standata$K <- 100
      }
            
      if(as.character(key$model[i]) %in% c("NMM", "ZIPNM", "PPM", "ZIPPM")){
        standata$K <- max(dat) + 100
      }
      
      fit <- mod$sample(data = standata, seed = 2022, 
                        chains = 1, parallel_chains = 1,
                        threads_per_chain = floor(detectCores()/4)-1,
                        iter_warmup = 100, iter_sampling = 100,
                        refresh = 5, adapt_delta = 0.90) 

      diagnostics <- as_draws_df(fit$sampler_diagnostics())
      loo <- fit$loo()$estimates           
      est <- fit$summary(c("p", "psi", "lambda", "chat", "pval")) %>%
        mutate(divergent = sum(diagnostics$divergent__),
               max_treedepth = sum(diagnostics$treedepth__ >= 10),
               walltime = fit$time()$total)
      # return final estimates table
      return(list(est = est, loo = loo))
    })
  })

# save estimates and time
saveRDS(est.list, paste0("data/model_output/est_list", i, ".rds"))
saveRDS(time, paste0("data/model_output/time", i, ".rds"))

} # Use loop if not working on HPC
