# This code has been adapted from Matt Craddock's code.
# (see https://craddm.github.io/blog/2017/02/25/EEG-topography)

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(magrittr)

# ==== Load data ====
raw.data <- read.csv("filename") %>%
  as.tbl() %>% 
  mutate(Condition   = as.factor(Condition),
         # convert Theta and Radius to (x,y)
         radianTheta = pi/180*Theta,
         x           = Radius*sin(radianTheta),
         y           = Radius*cos(radianTheta),
         r           = sqrt(x^2 + y^2))

# Support functions -------------------------------------------------------

expand.topo <- function(data, gridRes = NULL) {
  library(mgcv)
  if (is.null(gridRes))
    gridRes <- data[[1]] %$% unique(Channel) %>% {length(.)/2}
  
  res <- list()
  pb <- txtProgressBar(1,length(data)+1)
  for (c in 1:length(data)) {
    LL <- data[[c]]
    
    minx <- min(LL$x)*1.1
    miny <- min(LL$y)*1.1
    maxx <- max(LL$x)*1.1
    maxy <- max(LL$y)*1.1
    
    tmp.GAM <- expand.grid(x = seq(minx*1.05,maxx*1.05,length = gridRes),
                           y = seq(miny*1.05,maxy*1.05,length = gridRes)) %>%
      data.frame()
    
    tmp.GAM$Amp <- LL %>%
      gam(Amp ~ s(x, y, bs = 'ts'), data = .) %>% 
      predict(tmp.GAM, type = "response")
    res[[c]] <- tmp.GAM %>% 
      filter(sqrt(x^2 + y^2) <= max(diff(range(x)),diff(range(y)))/2)
    rm(tmp.GAM)
    setTxtProgressBar(pb, c)
  }
  return(res)
}

circleFun <- function(diameter = 1, center = c(0,0), npoints = 100) {
  r  <- diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Topo theme
theme_topo <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(rect       = element_blank(),
          line       = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank()
    )
}

# Select Channel Groups
group_channels <- function(chanlocs,...){
  # '...' is a series of names vectors listing the channels in each group 
  channel.groups <- list(...)
  lapply(1:length(channel.groups), function(x) {
    y <- chanlocs[chanlocs$Channel %in% channel.groups[[x]],]
    y$Chan.Group <- names(channel.groups)[x]
    return(y)
  }) %>% 
    tibble(data = .) %>% 
    unnest() %>% 
    mutate(Chan.Group = factor(Chan.Group)) %>% 
    return()
}

# Prepare Data for Plotting -----------------------------------------------

# Channel Data
 <- raw.data %>% 
  select(Channel, x, y) %>%
  distinct()

## Channel Groups for plotting later
# channel.groups <- chanlocs %>% 
#   group_channels()
## e.g.: group_channels(Frontal = c("E15","E23","E18","E16"))

# Head Shapes
diam      <- chanlocs %$%
  max(diff(range(x)),diff(range(y))) # can be adjusted
head <- list()
head$shape <- circleFun(diam*0.75)
head$mask  <- circleFun(diam*1.15)
head$nose  <- data.frame(x = c(-.075,0,.075),
                         y = c(-.005,.075,-.005) + max(head$shape$y))

# Topo Data
topo.data <- raw.data %>% 
  group_by(Condition,TimePnt) %>% 
  select(Condition,TimePnt,Channel,x,y,Amp) %>% 
  nest() %>%
  mutate(data = expand.topo(data)) %>% 
  unnest() # %>% 
  # If you only want whats inside the head
  # filter(!sqrt(x^2 + y^2) > diff(range(head$shape$x))/2)


# Plot --------------------------------------------------------------------

# Plot Channels
ggplot(mapping = aes(x,y)) +
  geom_path(data = head$shape) +
  geom_line(data = head$nose) +
  geom_text(data = chanlocs, aes(label = Channel)) +
  theme_topo() +
  coord_equal()

# Plot topos!
topo.plot <- ggplot(topo.data,aes(x, y)) +
  ## The Topos
  geom_raster(aes(fill = Amp))+
  stat_contour(aes(z = Amp), binwidth = 0.5, color = "gray")+
  facet_grid(Condition~TimePnt) +
  ## The Head
  geom_path(data = head$mask, color = "white",size = 4) +
  geom_path(data = head$nose, size = 1.5) +
  geom_path(data = head$shape, size = 1.5) +
  ## The Channels
  geom_point(data = chanlocs) +
  geom_text(data = chanlocs, aes(label = Channel)) +
  # stat_ellipse(data = channel.groups, aes(group = Chan.Group)) +
  # geom_point(data = channel.groups, aes(shape = Chan.Group)) +
  ## Design
  scale_fill_distiller(palette = "RdYlBu", direction = -1,
                       limits = range(topo.data$Amp), oob = squish) +
  labs(fill = "Amplitude") +
  theme_topo() +
  coord_equal()

topo.plot


