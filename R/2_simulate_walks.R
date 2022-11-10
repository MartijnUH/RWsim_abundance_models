# Required packages, data and data paths
library(MASS)
library(mvtnorm)
library(data.table)
library(raster)
library(sf)
library(lubridate)
library(tidyverse)
library(ggpubr)
library(activity)
source("R/functions.R")
memory.size(max = 5*10^6)
select <- dplyr::select

grid.list <- readRDS("data/grids.rds")
ct_data <- read.csv("data/NPHK_data.csv", sep = ";", stringsAsFactors = FALSE)

# output directory
odir <- "data/rw/"

# set situation
key <- 
  expand.grid(
    sim = 40,
    days = 25,
    hsteps = 6,
    r_hrc = c(F, T),
    hr_r = c(450, 900, 1800, 5400),
    speed = c("slow", "fast")
  )

key <- 
  key %>% mutate(
    hr_grid = case_when(
      hr_r == 450 ~ "g12",
      hr_r == 900 ~ "g6",
      hr_r == 1800 ~ "g3",
      hr_r == 5400 ~ "g2",
    ),
    Nmax = case_when(
      hr_r == 450 ~ 600,
      hr_r == 900 ~ 300,
      hr_r == 1800 ~ 150,
      hr_r == 5400 ~ 30,
    )
  ) %>% filter(r_hrc == F & hr_grid == "g12" | r_hrc == T)

###############
# Simulation
###############
## A. Generate temporal covariate for state probabilities
# trim ct data
ct_data <- ct_data %>%
  filter(!grepl("x2", deployment_sampling_point),  # exclude redundant deployments
         !deployment_sampling_point %in% c("JW_0361", "JW_0362"),  # exclude deployments manually
         !is.na(deployment_sampling_point)) %>%  # exclude deployments with NA value
  filter(animal_scientific_name == "Sus scrofa")  # subset data for species ...

groupsize <- ct_data %>% pull(animal_count)

# extract time columns for ct_data
ct_time_data <- ct_data %>% 
  separate(image_sequence_datetime, 
           into = c("image_sequence_date",
                    "image_sequence_timez"), 
           sep = "T", remove = FALSE) %>%
  separate(image_sequence_timez, 
           into = c("image_sequence_time","x"), 
           sep = "Z", remove = TRUE) %>%
  dplyr::select(animal_scientific_name, 
                image_sequence_datetime,
                image_sequence_date, 
                image_sequence_time) %>%
  mutate(image_sequence_time_rad = hms_to_rad(image_sequence_time))

# calculate activity density over radians in interval [0; 2*pi]
dens <- density2(ct_time_data$image_sequence_time_rad, from = 0, to = 2*pi)
plot(dens$x, dens$y, type="l")

# store densities in dataframe activity.df
activity.df <- tibble(time_rad = dens$x, density = dens$y) %>% 
  mutate(time_rad_bin = cut(time_rad, 24),
         time = rad_to_hms(time_rad),
         hour = lubridate::hour(time))
plot(activity.df$time_rad, activity.df$density)

# extract average densities for each time bin
z_cov <- activity.df %>%
  group_by(time_rad_bin) %>%
  dplyr::summarise(density = mean(density)) %>%
  mutate(step_id = 1:nrow(.)) %>%
  dplyr::select(step_id, everything()) %>%
  pull(density)

# scale activity covariate to [0,1]
z_cov_scaled <- z_cov/max(z_cov)

for (j in 1:nrow(key)) {
  
  # print progress in key
  print(paste0("situation: ", j, " out of ", nrow(key), " situations"))
  
  ## Simulation
  # simulate nsim random walks
  # set seed to generate nsim seed numbers
  set.seed(2021)
  
  # number of simulations
  nsim <- key$sim[j]
  
  # seed number for each simulation round
  seed_numbers <- sample(1:nsim, nsim)
  
  # number of individuals N to consider
  div <- switch(as.character(key$Nmax[j]), 
                "600" = 2, "300" = 2, "150" = 3, "30" = 3)
  
  Ntot <- key$Nmax[j]
  Nsteps <- c(key$Nmax[j]*c(1,0.9,0.8), (key$Nmax[j]/div) *c(1,0.9,0.8))
  
  # number of days
  ndays <- key$days[j]
  
  # bounding box of the study area
  bbox <- matrix(c(0,0,10800,10800), ncol = 1)
  rownames(bbox) <- c("xmin", "ymin", "xmax", "ymax")
  
  # camera selection grid
  cam_grid <- grid.list[["g12"]]
  
  # save movement data?
  save <- T
  # also store sf data?
  sf <- F
  # list to store simulation data
  rw.list <- list()
  true_states.list <- list()
  Ndata.list <- list()
  
  
  for (s in 1:nsim) {
    # print progress in random walks for each simulation
    print(paste0("simulate random walk: ", 
                 s, " out of ", key$sim[j], " simulations"))
    
    # set.seed for each simulation round
    set.seed(seed_numbers[s])
    
    # assign 
    N <- 0; Ntot <- key$Nmax[j]
    while(sum(N) != Ntot){
      if(sum(N) < Ntot){
        N <- c(N, sample(groupsize, 1))
      } else {
        N <- N[-length(N)]
        N <- c(N, sample(groupsize, 1))
      }
    }
    N <- N[-1]
    ngroup <- length(N)
    
    Nnames <- sapply(1:length(Nsteps), function(x) paste0("N", Nsteps[x]))
    Ndata <- data.frame(matrix(rep(N, length(Nsteps)), ncol = length(Nsteps)))
    colnames(Ndata) <- Nnames
    
   
    ##############################
    # activity-adjusted movement #
    ##############################
    # parameters
    #-----------#
    step.mu <- switch(as.character(key$speed[j]), 
                      "fast" = 380, "slow" = 100) * z_cov_scaled
    step.sd <- switch(as.character(key$speed[j]), 
                      "fast" = 140, "slow" = 35) * z_cov_scaled

    #############################
    # simulate random walk data #
    #############################
    if(key$r_hrc[j]){
        
      # sample random home range centers
      hrc <- round( matrix(c(
        runif(ngroup, bbox["xmin",], bbox["xmax",]),
        runif(ngroup, bbox["ymin",], bbox["ymax",])
      ), ncol = 2), 2)
        
    } else {
        
      # sample home range centers at sampling grid centroids
      idx <- sample(cam_grid$site_id, ngroup, replace = T)
        
      hrc <- 
        t(sapply(1:ngroup, function(i) {
          cam_grid %>% 
            filter(site_id %in% idx[i]) %>%
            st_centroid() %>%
            st_coordinates()}))
        
    }
      
    # compute home range areas
    hra <- st_buffer(
      st_cast(st_sfc(st_multipoint(hrc)), "POINT"), key$hr_r[j])
      
    # randomly sample initial positions within home range areas
    step0 <- 
      round(matrix(
        unlist(lapply(1:length(hra), function(i) {
          st_coordinates(st_sample(hra[[i]], 1, type = "random"))})), 
        byrow = T, ncol = 2), 2)
      
    if(0){
      grid_shift <- grid.list[["g12"]]
      grid_shift$x  <- grid_shift$x + c(50, 50)
      gp <- ggplot() +
        geom_sf(aes(col = "black"), fill = NA, linetype = "dotted", data = grid_shift) +
        geom_sf(aes(col = "black"), fill = NA, data = cam_grid) +
        geom_sf(aes(col = NA), data = st_buffer(st_multipoint(hrc), key$hr_r[j]), 
                alpha = 0.1, fill = "firebrick") +
        geom_sf(aes(col = "firebrick"), data = st_multipoint(hrc)) +
        geom_point(aes(x = step0[,1], y = step0[,2]), shape = 8) +
        #geom_sf(aes(col = "black"), data = init_pnts) +
        #theme_void() +
        scale_color_manual(name = NULL, values = c("black", "firebrick"),
                           labels = c("Intial location", 
                                      "Home range centre + Home range area")) +
        theme(legend.position = "none") 
        
      # plot home range
      gp + 
        labs(y = "Northing (m)", x = "Easting (m)") + 
        theme_bw() +
        theme(
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
        base_breaks_x(seq(0,10800,by=1800)) +
        base_breaks_y(seq(0,10800,by=1800))
        
      ggsave("fig/main/design_example.jpg", 
             width = 16, height = 16, unit = "cm", device = "jpg", dpi = 300)
    }
      
      # simulate a random walk, constrained by individual's home ranges
      rw <- random_walk_hrc(ndays = key$days[j], ngroup = ngroup, groupsize = N,
                             hrc = hrc, r = key$hr_r[j], step0 = step0,
                             step.mu = step.mu, step.sd = step.sd,
                             plot = F, seed = seed_numbers[s])
      
      if(key$hsteps[j] > 1){
        # add brownian motion to stored random walk if hsteps > 1
        rw <- add_brownian(rw, hsteps = key$hsteps[j], 
                           z = z_cov_scaled, plot = F)
        
      }
      
    tmp_states.list <- list()
    
    for(step in 1:length(Nsteps)){
      
      while(sum(N) != Nsteps[step]){
        id <- sample(which(N>0), 1)
        N[id] <- N[id] - 1
      }
      Ndata[,step] <- N
      
      # assign true latent states
      locs <- 
        cbind(N = Ndata[,step], 
              st_as_sf(
                st_buffer(st_sfc(st_multipoint(hrc)) %>% 
                            st_cast(.,"POINT"), key$hr_r[j]))
        ) 
      
      truth <- 
        st_intersection(cam_grid, st_as_sf(locs[locs$N>0,])) %>% 
        filter(st_is(.,"POLYGON")) %>%
        group_by(site_id) %>%
        summarize(ngroup = n(), N = sum(N)) %>%
        st_drop_geometry()
      
      tmp_states.list[[step]] <- list(psi = nrow(truth)/nrow(cam_grid), 
                                      lambda_M = Nsteps[step]/nrow(cam_grid), 
                                      lambda_unM = sum(truth$N)/nrow(cam_grid), 
                                      data = truth)
      
    }
    
    # assign rw_data in each simulation to a list
    rw.list[[s]] <- rw
    
    # assign true states in each simulation to a list
    names(tmp_states.list) <- Nsteps
    true_states.list[[s]] <- tmp_states.list
    
    # assign Ndata in each simulation to a list
    Ndata.list[[s]] <- Ndata
    
  }
  
  if(save){
    
    # save list (of lists) of true states
    saveRDS(true_states.list, paste0(odir, "SList_",
                                     "s", key$sim[j], 
                                     "_N", key$Nmax[j],
                                     "_d", key$days[j], 
                                     "_hs", key$hsteps[j],  
                                     "_hrr", key$hr_r[j], 
                                     "_rhrc", key$r_hrc[j],
                                     "_", key$speed[j],
                                     ".rds"))
    
    # save list (of lists) of Ndata
    saveRDS(Ndata.list, paste0(odir, "NList_",
                               "s", key$sim[j], 
                               "_N", key$Nmax[j],
                               "_d", key$days[j], 
                               "_hs", key$hsteps[j],  
                               "_hrr", key$hr_r[j], 
                               "_rhrc", key$r_hrc[j],
                               "_", key$speed[j],
                               ".rds"))
    
    # save random walk trajectory datalist
    saveRDS(rw.list, paste0(odir, "Wlist_",
                            "s", key$sim[j], 
                            "_N", key$Nmax[j],
                            "_d", key$days[j], 
                            "_hs", key$hsteps[j],  
                            "_hrr", key$hr_r[j], 
                            "_rhrc", key$r_hrc[j],
                            "_", key$speed[j],
                            ".rds"))
    
  }
}
