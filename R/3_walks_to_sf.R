#################################
# RW trajectory data to sf data #
#################################

# This script takes a list of data frames generated 
# by Random Walks in the RW_simulation.R script
# and converts this to a list of Sf data frames.
# 
# Note: this script makes use of the build-in function foreach,
# which runs more efficiently on Linux! 


## Initialize
# required R-packages
library(dplyr)
library(data.table)
library(sf)
library(parallel)
library(doParallel)
select <- dplyr::select

# increase memory allocation
memory.size(max = 5*10^6)

# file directory
fdir <- "data/rw/"

# set key
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

# save as .RDS (Robject), else as .SHP (shapefile)
as_rds <- T

for(i in 1:nrow(key)){
  
  # 0. Print progress
  print(paste0("situation: ", i, " out of ", nrow(key), " situations"))
  
  # 1. Read list of data frames containing random walk trajectories
  rw.list <- readRDS(paste0(fdir, "Wlist_",
                            "s", key$sim[i], 
                            "_N", key$Nmax[i],
                            "_d", key$days[i], 
                            "_hs", key$hsteps[i],  
                            "_hrr", key$hr_r[i], 
                            "_rhrc", as.logical(key$r_hrc[i]),
                            "_", key$speed[i],
                            ".rds"))
  
  # 2. Convert to list of sf data frames
  rw.sf <- 
      lapply(1:length(rw.list), function(i){
        # print progress
        print(paste0(i, " / ", length(rw.list)))
        # convert to data frame and drop step 0, 
        # except for the initial point location (i.e. at day 1)
        dat <- as_tibble(rw.list[[1]]) %>% 
          filter(step_id > 0 | day_id == 1)
        
        # split data 
        coords <- dat %>% select(group_id, step_x:step_y)
        cdat <- dat %>% select(-group_id, -c(step_x:step_y))
        
        # set up PSOCK clusters
        cl <- makeCluster(detectCores()-2, type = "PSOCK") 
        registerDoParallel(cl)
        
        # convert point coordinates into sf linestrings
        system.time({
          coords <- 
          foreach(g = 1:max(coords$group_id), .combine = "rbind") %dopar% {
            # run garbage collector
            gc()
            # set up the environment for parallel R-sessions
            library(sf); library(dplyr); library(purrr)
            
            st_segment = function(step_x, step_y, step_x2, step_y2) {
              if(is.na(step_x2)|is.na(step_y2)){
                st_point(c(step_x,step_y))
              } else {
                st_linestring(matrix(c(step_x, step_x2, step_y, step_y2), 2, 2))
              }
            }
            
            out <- as_tibble(coords) %>%
              filter(group_id == g) %>%
              mutate(step_x2 = step_x[row_number() + 1],
                     step_y2 = step_y[row_number() + 1]) 
            
            geom <- out %>%
              select(starts_with("step")) %>%
              pmap(st_segment)
            
            out$geometry <- geom
            return(out %>% select(-step_x2, -step_y2))
            
          } 
        })
        
        # stop the PSOCK clusters
        stopCluster(cl)
        
        return(cbind.data.frame(cdat, coords) %>% 
                 select(day_id, group_id, everything()) %>% 
                 filter(step_id > 0) 
        )
      }) %>% rbindlist(idcol = "sim_id") %>% st_as_sf()
  
  # 3. Save
  if(as_rds){
    # 3.1. Save binary geometries + standard dataframe
    saveRDS(st_as_binary(st_geometry(rw.sf)), 
            paste0(fdir, "shapefiles/Wsf_hr", key$hr_grid[i], 
                   "_rhrc", key$r_hrc[i], "_", key$speed[i], "_geom.rds"))
    
    saveRDS(st_drop_geometry(rw.sf), 
            paste0(fdir, "shapefiles/Wsf_hr", key$hr_grid[i], 
                   "_rhrc", key$r_hrc[i], "_", key$speed[i], ".rds"))
    
  }
  else {
    # 3.2. Save shapefile
    st_write(rw.sf, paste0(fdir, "Wsf", "_hr", key$hr_grid[i], ".shp"))
  }
}
  
    