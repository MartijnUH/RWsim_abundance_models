# This script contains functions that are needed to simulate animal trajectories

## Functions to obtain r baseplot-style x- and y-axis in ggplot
base_breaks_x <- function(b){
  d <- data.frame(y=-Inf, yend=-Inf, x=min(b), xend=max(b))
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), 
                    inherit.aes=FALSE),
       scale_x_continuous(breaks=b))
}
base_breaks_y <- function(b){
  d <- data.frame(x=-Inf, xend=-Inf, y=min(b), yend=max(b))
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), 
                    inherit.aes=FALSE),
       scale_y_continuous(breaks=b))
}

## Function to construct a spatial sf dataframe grid layer
sf_make_grid <- function(n_sites,
                         scale,
                         square = T,
                         cov_raster){
    if(is.null(class(cov_raster))){
      grid <- 
        st_as_sf(
          st_make_grid(
            as(
              as(raster::extent(0, 1 * 10^scale, 0, 1 * 10^scale), "SpatialPolygons")
              , "sf"), 
            square = square, n = sqrt(n_sites)),
          site_id = as.integer(factor(1:n_sites))
        )
    } else {
      grid <- 
        st_as_sf(
          st_make_grid(
            as(
              as(raster::extent(0, 1 * 10^scale, 0, 1 * 10^scale), "SpatialPolygons")
              , "sf"), 
            square = square, n = sqrt(n_sites)),
          site_id = as.integer(factor(1:n_sites))
        )
        
        grid <- 
          # add covariate info
          grid %>%
          mutate(cov = unlist(lapply(raster::extract(cov_raster, grid), mean)),
                 cov_cat = ifelse(cov >= quantile(cov, 0.50), "high", "low"))
    }
  return(grid)
}

## Function to sample a sf grid layer, according to various sampling schemes
sample_grid <- function(grids.sf,   # sf data.frame with grid layers to sample from
                        n,          # sample size
                        type,       # type of sampling: 'random', 'grts', 'regular' or 'clusterd'
                        cov,        # covariate values [0 - 1] for 'preferential' probability
                        seed        # seed number 
) {
  set.seed(seed)
  if(type == "grts") {
    # grid random tesselation sampling scheme
    sample <- grts.polygon(as_Spatial(grids.sf), n = n)
    out.sf <- grids.sf[st_as_sf(sample), ] %>% st_centroid()
  }
  
  else if(type == "preferential") {
    
    sample <- sample(grids.sf$site_id, n, 
                     prob = plogis(2 * unlist(st_drop_geometry(grids.sf[, cov]))))
    out.sf <- grids.sf[sample,] %>% st_centroid()
    
  } else {
    # random, regular or clustered sampling scheme 
    sample <- spsample(as_Spatial(grids.sf), n = n, type = type)
    out.sf <- grids.sf[st_as_sf(sample), ] %>% st_centroid()
  }
  return(out.sf)
}

## Function to return a point on the arc of a circle, from:
get_arcpoint <- function(x1, y1,  # x and y- coords of the startpoint
                         theta,   # angle in rad
                         l        # length of the line connecting the points
) {
  return( c(x1 + l * cos(theta), y1 + l * sin(theta)) )
}

## Function that returns the camera field of view from: 
get_circ_segement <- function(sf,            # sf dataframe with deployments
                              theta,         # compas angle in degrees
                              r              # radius 
){
  # convert half of the fov angle (in degrees) to radians
  theta = theta * (pi/180)
  theta_min <- (pi/2) - (theta/2)
  theta_plus <- (pi/2) + (theta/2)
  
  # point 1
  p1 <-  st_sf(st_sfc(
    lapply(1:nrow(sf), function(i) {
      st_point(
        get_arcpoint(x1 = st_coordinates(sf)[i,1], 
                     y1 = st_coordinates(sf)[i,2], 
                     theta = (theta_min), l = r * 1.25)
      )
    })
  ))
  
  # point 2
  p2 <-  st_sf(st_sfc(
    lapply(1:nrow(sf), function(i) {
      st_point(
        get_arcpoint(x1 = st_coordinates(sf)[i,1], 
                     y1 = st_coordinates(sf)[i,2], 
                     theta = (theta_plus), l = r * 1.25)
      )
    })
  ))
  
  # polygon to clip, from deployment point, point 1 and point 2
  clip <- st_sf(st_polygonize(
    st_cast(st_as_sfc(
      lapply(1:nrow(sf), function(i) {
        st_linestring(rbind(st_coordinates(p1)[i,],
                            st_coordinates(sf)[i,],
                            st_coordinates(p2)[i,],
                            st_coordinates(p1)[i,]))
      })), 
      "MULTILINESTRING")
  ))
  
  # buffer layer to be clipped
  circles <- sf %>% st_buffer(dist = r)
  # clipped circlular segment (camera's field of view)
  camera_fov <- st_as_sf(site_id = sf$site_id, st_as_sfc(
    lapply(1:nrow(circles), function(i) {
      st_geometry(st_intersection(clip[i,], circles[i,]))[[1]]
    } )
  ))
  
  return(camera_fov)
}

## Function to calculate the compass bearing of the line between two points
calc_angle <- function(
  xx,            # differences in x coordinates between two points
  yy,            # differences in y coordinates between two points
  bearing = TRUE,  # bearing = FALSE returns +/- pi instead of 0:2*pi
  as.deg = FALSE   # as.deg = TRUE returns degrees instead of radians
) {
  c = 1
  if (as.deg){
    c = 180/pi
  }
  if(xx + yy != 0){
    b <- sign(xx)
    b[b == 0] <- 1  #corrects for the fact that sign(0) == 0
    tempangle = b*(yy < 0) * pi + atan(xx/yy)
    if(bearing){
      #return a compass bearing 0 to 2pi
      #if bearing==FALSE then a heading (+/- pi) is returned
      tempangle[tempangle < 0] <- tempangle[tempangle < 0] + 2*pi
    }
  } else {
    tempangle = 0
  }
  return(tempangle * c)
}

## Function to alculates the bearing and length of the two lines formed by three points
bearing.ta <- function(
  loc1,loc2,loc3,  # locations are assumed to be in (X,Y) format
  as.deg=FALSE     # as.deg = TRUE returns degrees instead of radians
  ){
  
  if (length(loc1) != 2 | length(loc2) != 2 | length(loc3) !=2){
    print("Locations must consist of either three vectors, length == 2,
or three two-column dataframes")
    return(NaN)
  }
  c = 1
  if (as.deg){
    c = 180/pi
  }
  
  locdiff1 <- loc2 - loc1
  locdiff2 <- loc3 - loc2
  bearing1 <- anglefun(locdiff1[1], locdiff1[2], bearing=F)
  bearing2 <- anglefun(locdiff2[1], locdiff2[2], bearing=F)
  if(is.data.frame(locdiff1)){
    dist1 <- sqrt(rowSums(locdiff1^2))
    dist2 <- sqrt(rowSums(locdiff2^2))
  }else{
    dist1 <- sqrt(sum(locdiff1^2))
    dist2 <- sqrt(sum(locdiff2^2))
  }
  
  ta = (bearing2 - bearing1)
  
  ta[ta < -pi] = ta[ta < -pi] + 2*pi
  ta[ta > pi] = ta[ta > pi] - 2*pi
  return(list(bearing1 = unlist(bearing1*c), bearing2 = unlist(bearing2*c),
              ta = unlist(ta*c), dist1 = unlist(dist1), dist2 = unlist(dist2)))
}

## Function to convert hours, minutes and seconds to radians
hms_to_rad <- function(
  time_col,  # column with image_sequence_time as a character
  tz = "UTC" # time zone, default is UTC
) {
  
  return(
    (as.numeric(
      as.POSIXct(strptime(time_col, format = "%H:%M:%S", tz = tz))) -
        as.numeric(as.POSIXct(strptime("0", format = "%S", tz = tz)))) / 3600 * (pi/12)
  )
}

## Function to convert radians into hours, minutes and seconds
rad_to_hms <- function(
  rad_col,   # column with radians
  tz = "UTC" # time zone, default is UTC
) {
  
  degrees <- rad_col * (180/pi)
  
  deg <- floor(degrees)
  arcmin <- floor((degrees - floor(degrees)) * 60)
  arcsec <- round(((degrees - floor(degrees)) * 60 - floor((degrees - floor(degrees)) * 60))*60, 0)
  
  hours <- as.character(floor(deg*(24/360)))
  minutes <- as.character(floor(arcmin * (15/60)))
  seconds <- as.character(floor(arcsec * (15/60)))
  
  return(
    hms( paste0(ifelse(nchar(hours) != 2, paste0(0, hours), hours), ":", 
                ifelse(nchar(minutes) != 2, paste0(0, minutes), minutes), ":", 
                ifelse(nchar(seconds) != 2, paste0(0, seconds), seconds)) )
  )
}

## Function to convert minutes to hms
min_to_hms <- function(x) { 
  require(lubridate)
  h = floor(x/60)
  m = floor(x %% 60)
  if(h < 10){
    return(hms(paste0("0", h, ":", m, ":00")))
  } else {
    return(hms(paste0(h, ":", m, ":00")))
  }
}

## Function to convert hours to hms
hour_to_hms <- function(x){
  require(lubridate)
  h <- floor(x)
  m <- round((x %% 1) * 60 ,2)
    
  if(h < 10  & m < 10){
    return(hms(paste0("0", h, ":", "0", m, ":00")))
  } else  if( h < 10 & m > 10){
    return(hms(paste0("0", h, ":", m, ":00")))
  } else if( h > 10 & m < 10){
    return(hms(paste0(h, ":", "0", m, ":00")))
  } else(
    return(hms(paste0(h, ":", m, ":00")))
  )
}

# Function to sample a point within a home range, given the radius and latent center
# and within the study field
sample_pnt_from_circ <- function(latent_centers, radius, bbox, seed){
  set.seed(seed)
  x <- y <- - 100
  while(
    any(c(x,y) < bbox[c("xmin","ymin"),] | c(x,y) > bbox[c("xmax","ymax"),])
  ){
    r <- radius 
    theta <- runif(1, 0, 1) * 2 * pi
    x <- r * cos(theta) + latent_centers[1]
    y <- r * sin(theta) + latent_centers[2]
  }
  return(c(x, y))
}

# Function to simulate a random walk within a home range, based on a binary activity state z
random_walk_hrc <- 
  function(
    ndays,          # number of days (numeric)
    ngroup,         # number of individuals (numeric)
    groupsize,      # number of animals in each individual group
    hrc,            # a (ngroup x 2) matrix with coordinates for latent home range centres of the group
    r,              # numeric, size of the home range
    step0,          # a (ngroup x 2) matrix with coordinates at step0 for each group
    step.mu,        # max. step length in x and y direction
    step.sd,        # max. standard deviation of step length in x and y direction
    seed,           # seed number
    plot = TRUE
  ){ 
    
    set.seed(seed)
    require(ggplot2)
    
    # Total number of steps per group
    nstep <- ndays*24
    step.mu <- rep(step.mu, ndays)
    step.sd <- rep(step.sd, ndays)

    # Initialize output 
    out <- matrix(0, ncol = 8)
    colnames(out) <- c("group_id", "nind", "day_id", "step_id", 
                       "step_x", "step_y", "step_l", "step_a")
    
    for(n in 1:ngroup){
      # Initialize empty vectors
      Xt <- Yt <- dist <- ang <- rep(-1e6, nstep)
      Xt[1] <- step0[n,1]; Yt[1] <- step0[n,2]

      for (i in 2:nstep) {
        # Constrain walk to the home range
        while((Xt[i] - hrc[n, 1])^2 + (Yt[i] - hrc[n, 2])^2 > r^2) {
          
          # Consecutive steps in the random walk
          step <- replicate(2, rnorm(1, step.mu[i], step.sd[i]))   
          Xt[i] <- Xt[i-1] + ifelse(runif(1)>0.5, step[1], -step[1])
          Yt[i] <- Yt[i-1] + ifelse(runif(1)>0.5, step[2], -step[2])
          
          # Calculate distance moved in 2d plane
          dist[i] <- pointDistance(c(Xt[i-1], Yt[i-1]), c(Xt[i], Yt[i]), longlat = F)
          
          # Calculate angle of movement in 2d plane
          ang[i] <- calc_angle(Xt[i-1] - Xt[i], Yt[i-1] - Yt[i], bearing = F)
          
        }
      }
      # Bind together walks of all groups
      out <- rbind(out, 
                   as.matrix(data.frame(group_id = rep(n, nstep), 
                                        nind = rep(groupsize[n], nstep), 
                                        day_id = rep(1:ndays, each=24), 
                                        hour_id = rep(1:24, ndays),
                                        step_x = Xt, step_y = Yt,
                                        step_l = dist, step_a = ang)))
    }
    
    return(out[-1,])
  }

# Function to simulate additional Brownian motion between consecutive locations
# generated from a random walk with the FUNCTION 'random_walk_hrc'
add_brownian <- 
  function(
    data,    # discrete trajectory data - output from FUNCTION 'random_walk_hrc'
    hsteps,  # number of steps per hour,
    z,       # activity weights per hour
    plot = TRUE
  ) {
    # initialize
    data <- data.frame(data)
    cols <- ncol(data)+1
    out <- matrix(ncol = cols)
    ngroup <- max(data$group_id)
    ndays <- max(data$day_id)

    for(i in 1:ngroup) {

      # arrange by day id and step id
      data_i <- data %>% 
        filter(group_id == i) %>%
        arrange(day_id, step_id)
      
      # initialize container for Brownian bridges
      brownian <- as.numeric(data_i[1, c("step_x","step_y")])
      l <- list()

      for(j in 2:nrow(data_i)) {
        
        # total time in seconds
        T <- 60*60
        # time increment
        Dt <- T/hsteps
        # time points
        t <- seq(0,T,length=hsteps+1)
        
        # initialization of the matrix W
        W <- matrix(rep(numeric(hsteps+1),2),ncol=2) 
        
        # weighing factor for amount of drift
        wd <- z[data_i$step_id[j-1]]
        
        # simulate Brownian motion
        for(hs in 2:(hsteps+1))
          W[hs,]<-W[hs-1,]+rnorm(2)*sqrt(Dt)*wd

        # start and end point
        x <- as.numeric(data_i[j-1, c("step_x","step_y")])
        y <- as.numeric(data_i[j, c("step_x","step_y")])
        
        # simulate Brownian bridge
        BB <- sapply(1:2, function(k){
          x[k] + W[,k] - t/T * (W[hsteps+1,k] - y[k] + x[k])})
      
        # add to container
        brownian <- rbind(brownian, BB[-1,]) 
      }
      # create new trajectory data for group i
      new_step_id <- rep(data_i$step_id, each=6)
      new_step_id <- new_step_id[-(length(new_step_id)-4):-length(new_step_id)]
      
      data_i <- data_i %>%
        slice(rep(1:(n()), each = hsteps)) %>%
        slice(1:(n()-5)) %>%
        rename(hour_id = step_id) %>%
        mutate(step_x = brownian[,1], step_y = brownian[,2]) %>%
        mutate(step_id = new_step_id) %>%
        select(day_id:hour_id, step_id, everything())
      
      # combine all group-specific trajectories
      out <- rbind(out, as.matrix(data_i[,1:cols]))
    }
    
    if(plot){
      require(ggplot2); require(tibble)
      p <- ggplot(as_tibble(out[-1,]), 
                  aes(x = step_x, y = step_y, colour = as.factor(group_id))
      )
      p <- p + geom_path()
      print(p)
    }
    return(out[-1,])
  }

# Function to calculate the step length and
add_slength_sangle <- function(data){
  
  data <- data.frame(data)
  out <- matrix(ncol = ncol(data))
  
  for(i in 1:max(data$group_id)){
    data_i <- data %>% filter(group_id == i)
    data_i$step_length <- 0
    data_i$step_angle <- 0
    
    for(r in 2:nrow(data_i)){
      # Calculate distance moved in 2d plane
      data_i$step_length[r] <- pointDistance(
        data_i[r-1, c("step_x","step_y")],
        data_i[r, c("step_x","step_y")],
        longlat = F)
      
      # Calculate angle of movement in 2d plane
      data_i$step_angle[r] <- calc_angle(
        data_i[r-1, "step_x"] - data_i[r, "step_x"], 
        data_i[r-1, "step_y"] - data_i[r, "step_y"],
        bearing = F)
    }
    out <- rbind(out, as.matrix(data_i))
  }
  return(out[-1,])
}


## Function to convert point geometries of consecutive steps to a collection of linestrings    
points_to_linestrings <- 
  function(
    dat  # a sf data frame, with point geometries
  ){
    
    # convert to sf dataframe
    if(class(dat)[1] != "sf") {
      dat <- st_as_sf(dat, coords = c("step_x", "step_y"))
    }
    
    # extract coordinates 
    coords <- st_coordinates(dat)
               
    # initialize list to store new linestring geometries
    line <- list()
    line[[1]] <- st_geometry(dat)[[1]]
                       
    for(i in 2:nrow(coords)){
      # generate linestrings from consecutive point locations
      line[[i]] <- st_linestring(rbind(coords[i-1,], coords[i,]))
    }
    
    # replace point geometries by linestrings
    st_geometry(dat) <- st_as_sfc(line)
               
    # remove initial step from the data (remaining point geometry)
    return(dat)
  }


## Function to convert ordered point geometries to a collection of linestrings 
st_segment = function(step_x, step_y, step_x2, step_y2) {
  if(is.na(step_x2)|is.na(step_y2)){
    st_point(c(step_x, step_y))
  } else {
    st_linestring(matrix(c(step_x, step_x2, step_y, step_y2), 2, 2))
  }
}

## Function to convert trajectories' data to a data frame,
## that contains true N's and counts.
trajectory_to_surveydata <- 
  function(
    data.sf,
    grid.sf,
    cam.sf,
    
    p = 1, 
    r = 1,
    arrays = F
  ){
      # data, drop intial point location
      data.sf <- data.sf[-1,]
      # metadata
      meta.df <- 
        expand.grid(
          site_id = 1:max(grid.sf$site_id), 
          day_id = 1:max(data.sf$day_id)
        )
      
      # Derive No. of animal present per grid per day 
      # from trajectory sf dataframe 
      in_grid <- st_intersection(data.sf, grid.sf) 
      
      n.df <- in_grid %>%
        st_drop_geometry() %>%
        group_by(group_id, site_id, .add = T) %>%
        summarize(n = mean(n)) %>%
        group_by(site_id) %>%
        summarize(n = sum(n))
      
      # Derive No. of animal passing the camera's fov per day 
      # from trajectory sf dataframe 
      in_fov <- st_intersection(in_grid, cam.sf) 
      
      n_fov.df <- in_fov <- 
        st_intersection(data.sf, cam.sf) %>%
        st_drop_geometry() %>%
        
        # 1.1. individuals detected per group, per site, per day
        # aggregate n over 24 hours: unique() --> marked; sum() --> unmarked
        group_by(group_id, site_id, day_id) %>%
        summarize(n_M = unique(n), n_unM = sum(n), .groups = "keep") %>%
        
        # 1.2. individuals detected per site, per day
        # aggregate n over groups: sum()
        group_by(site_id, day_id) %>%
        summarize(n_fov_M = sum(n_M), n_fov_unM = sum(n_unM), .groups = "keep")
        
      
      ## joined data frame: per site_id and day_id:
      # n_unM: No. of animals present, if individuals were unmarked (redundant counts)
      # n_fov_unM: No. of animals in camera field of view, if individuals were unmarked (redundant counts)
      # n_M: No. of animal present per grid per day, if individuals were marked (no redundant counts)
      # n_fov_M: No. of animals in camera field of view, if individuals were marked (no redundant counts)
      # occu: occupancy status (0: absent, 1: present)
      # sampled: site sampled by a CT (0: no, 1: yes)
      # in_fov: at least one individual in camera fov (0: no, 1: yes)
      # count_unM: number of animals detected from n_unM animals present
      # count_M: number of animals detected from n_M animals present
      # det: at least one detected (0: undetected, 1: detected)
      
      out.df <- meta.df %>%
        left_join(n.df, "site_id") %>%
        left_join(n_fov.df, c("site_id", "day_id")) %>%
        group_by(site_id, day_id) %>%
        
        mutate(occu = ifelse(n > 0, 1, 0),
               sampled = ifelse(site_id %in% cam.sf$site_id, 1, 0),
               in_fov = ifelse(is.na(n_fov_M), 0, 1),
          
               count_M = ifelse(sampled & in_fov, 
                                rbinom(1, n_fov_M, prob = r), 
                                NA),
               
               count_unM = ifelse(sampled & in_fov, 
                                  rbinom(1, n_fov_unM, prob = r), 
                                  NA),
               
               det = ifelse(sampled & in_fov, 
                            rbinom(1, 1, prob = p), 
                            NA)) %>%
        mutate_all(~replace(., is.na(.), 0))
      
     return(st_as_sf(left_join(out.df, grid.sf, by = "site_id")))
  }

## Function to calculate the true occupancy (psi), 
## detection probability (p), 
## daily detection probability (p_daily) and
## detection history (det_hist) for a list of data
list_survey_info <- 
  function(data.sf,
           grid.sf,
           cam.sf,
           deploy,
           true_states,
           p = 1,
           r = 1){
    
    ## 1. no. of individuals passing the camera's fov per day 
    in_fov <- 
      st_intersection(data.sf, cam.sf) %>%
      st_drop_geometry() %>%
      
      # 1.1. individuals detected per group, per site, per day
      # aggregate n over 24 hours: unique() --> marked; sum() --> unmarked
      group_by(group_id, site_id, day_id) %>%
      summarize(n_M = unique(n), n_unM = sum(n), .groups = "keep") %>%
      
      # 1.2. individuals detected per site, per day
      # aggregate n over groups: sum()
      group_by(site_id, day_id) %>%
      summarize(n_fov_M = sum(n_M), n_fov_unM = sum(n_unM), .groups = "keep")

    
    ## 2. metadata to perform left join
    meta <- 
      expand.grid(
        site_id = 1:max(grid.sf$site_id), 
        day_id = 1:max(data.sf$day_id)
      )
    
    ## 3. marked/ unmarked counts and detections binned per day
    counts <- meta %>%
      left_join(in_fov, c("site_id", "day_id")) %>%
      group_by(site_id, day_id) %>%
      
      mutate(sampled = ifelse(site_id %in% cam.sf$site_id, 1, 0),
             in_fov = ifelse(is.na(n_fov_M), 0, 1),
             
             count_M = ifelse(sampled & in_fov, 
                              rbinom(1, n_fov_M, prob = r), 
                              NA),
             
             count_unM = ifelse(sampled & in_fov, 
                                rbinom(1, n_fov_unM, prob = r), 
                                NA),
             
             det = ifelse(sampled & in_fov, 
                          rbinom(1, 1, prob = p), 
                          NA)) %>%
      ungroup() %>%
      mutate_all(~replace(., is.na(.), 0))
    
    # arrange data by day
    counts <- counts %>% arrange(day_id)
      
    # sampling effort
    nsites <- max(counts$site_id)
    ndays <- max(counts$day_id)

    ## 4. Count history
    # generate count history as a list of two matrices:
    # for respectively unmarked counts and marked counts
    count_hist.list <-
      lapply(1:2, function(m){
        
        col <- ifelse(m == 1, "count_unM", "count_M")
        count_long <- counts[[col]]
        count_hist <- matrix(count_long, nrow = nsites, ncol = ndays)
        count_hist[-cam.sf$site_id,] <- NA
        
        return(count_hist)
        
        })
    
    names(count_hist.list) <- c("unmarked", "marked")
    
    ## 5. Sampled sites
    sampled <- deploy
      
    ## 5. true states
    # 5.1. mean site-occupancy 
    psi <- true_states$psi
      
    # 5.2. mean site-abundance
    lambda <- 
      lapply(1:2, function(m){
        
        col <- ifelse(m == 1, "lambda_unM", "lambda_M")
        return(true_states[[col]])
        
      })
      
    names(lambda) <- c("unmarked", "marked")

    ## 6. true detectability
    # 6.1. per species detectability
    z <- rep(F,nsites); z[true_states$data$site_id] <- T
    p <- 
      lapply(1:length(deploy), function(x){
        d <- rep(F,144); d[deploy[[x]]] <- T
        mean(as.numeric(as.logical(count_hist.list$unmarked[z & d,])))
      })
    names(p) <- names(deploy)
    # 6.2. per capita detectability
    N <- rep(0,144); N[true_states$data$site_id] <- true_states$data$N
    r <-
      lapply(1:length(deploy), function(x){
        d <- rep(F,144); d[deploy[[x]]] <- T
        mean(1 - (1-p[[x]])^(1/N[d & z]))
      })
    names(r) <- names(deploy)
    
    r2 <- 
      lapply(1:length(deploy), function(x){
        d <- rep(F,144); d[deploy[[x]]] <- T
        mean(count_hist.list$marked[z & d,]/ N[d & z])
      })
    names(r2) <- names(deploy)
   # 6.3. per capita encounter rate
   mu <- 
      lapply(1:length(deploy), function(x){
        d <- rep(F,144); d[deploy[[x]]] <- T
        mean(count_hist.list$unmarked[z & d,]/ N[d & z])
      })
   names(mu) <- names(deploy)
   
   

    return(list(count_hist = count_hist.list, psi = psi, lambda = lambda, 
                p = p, r = r, r2 = r2, mu = mu, sampled = sampled))
  }


## Function to obtain naive estimates from counts
get_naive <- function(x, par){
  if(par == "psi"){
    mean(as.numeric(as.logical(rowSums(x))))
  }
  else if(par == "lambda") {
    mean(apply(x, 1, max, na.rm=T))
  }
  else if (par == "p"){
    mean(as.numeric(as.logical(x[rowSums(x)>1, na.rm=T])))
  }
  else {
    print("Choose either 'p', 'psi' or 'lambda' as par!")
  }
}

## Function to make a panel plot summarizing the ecological survey
## Panel A: group-specific trajectories
## Panel B: occupancy and detection status over all sampling days
## Panel C: True number of animals per site
## Panel D: Unmarked animal counts per site
summary_plot <- 
  function(data,      # detection/count data per sampling grid (sf)
           hr_grid,   # home range grid (sf),
           cam_grid,  # camera grid (sf),
           deploy,    # camera locations (sf)
           traj       # trajectory data (sf)
  ){
    gp1 <- ggplot() +
      #ggtitle("A. Trajectories") +
      geom_sf(data = hr_grid, aes(), fill = NA, linetype = "dotted") +
      geom_sf(data = cam_grid, aes(), fill = NA) +
      geom_sf(data = traj, col = "firebrick") +
      geom_sf(data = deploy, col = "black") + 
      labs(y = "Northing (m)", x = "Easting (m)") +
      base_breaks_x(seq(0,10800,by=1800)) +
      base_breaks_y(seq(0,10800,by=1800)) +
      theme_bw() +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
    
    f <- factor(factor(data$occu):factor(data$det))
    gp2 <- ggplot() +
      geom_sf(data = data, aes(fill = f)) +
      geom_sf(data = deploy, aes()) +
      theme_void() +
      theme(legend.title = element_blank(),
            plot.margin = margin(0.5,0.2,0.02,1, "cm")) +
      scale_fill_manual(values = c("grey90", "grey50", "green"),
                        labels = c("Unoccupied", "Occupied", "Detected"),
                        drop = F)
    
    gp3 <- ggplot() +
      geom_sf(data = data, aes(fill = factor(n))) +
      geom_sf(data = deploy, aes()) +
      theme_void() +
      theme(legend.title = element_blank(),
            plot.margin = margin(0.5,0.2,0.02,1, "cm")) +
      scale_fill_manual("LA", values = c("white", brewer.pal(length(unique(data$n))-1, "Reds")))
   

    gp4 <- ggplot() +
      geom_sf(data = data, aes(fill = factor(count_unM))) +
      geom_sf(data = deploy, aes()) +
      theme_void() +
      theme(legend.title = element_blank(),
            plot.margin = margin(0.5,0.2,0.02,1, "cm")) +
      scale_fill_manual("Counts", values = c("white", brewer.pal(length(unique(data$n))-1, "Reds")))
    
    return(
      ggarrange(
        ggarrange(gp1, gp2, ncol = 2,legend = "bottom",
                  labels = c("(a)", "(b)"), hjust = - 0.75),
        ggarrange(gp3, gp4, ncol = 2, legend = "bottom",
                  labels = c("(c)", "(d)"), hjust = - 0.75), 
        nrow = 2, align = "hv")
    )
  }
