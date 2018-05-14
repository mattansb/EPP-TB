library(tidyverse)
library(scales)

# ==== Load data ====

df <- read_csv("@filename@")

raster2rect <- function(data) {
  freq.vals <- unique(data$Frequency)
  
  d.time <- diff(unique(data$Time))[1]/2
  d.frq <- diff(freq.vals)/2
  
  upper <- freq.vals + c(d.frq, d.frq[length(d.frq)])
  lower <- freq.vals - c(d.frq[1], d.frq) 
  ntime <- data$Time %>% unique() %>% length
  
  data$xmin <- data$Time - d.time
  data$xmax <- data$Time + d.time
  data$ymin <- map_dbl(lower, ~rep(.x, ntime))
  data$ymax <- map_dbl(upper, ~rep(.x, ntime))
  
  
  return(data)
}


# Plot --------------------------------------------------------------------


ERSP.plot <- ggplot(df,aes(Time, Frequency, fill = ersp)) +
  ## Plot data
  geom_raster(interpolate = T) + # plot y-axis scale as is
  # geom_rect(aes(xmin=xmin, xmax=xmax, # change scaling of y-axis
  #               ymin=ymin, ymax=ymax),
  #           data = raster2rect(df),
  #           color = NA) +
  facet_grid(.~Condition) +
  ## Theme, Scales and Labels
  # scale_y_log10() + # if Frequencies are log-spaced
  scale_fill_distiller(palette = "OrRd", limits = range(df$ersp), oob=squish)

ITC.plot <- ggplot(df,aes(Time, Frequency, fill = itc)) +
  ## Plot data
  geom_raster(interpolate = T) + # plot y-axis scale as is
  # geom_rect(aes(xmin=xmin, xmax=xmax, # change scaling of y-axis
  #               ymin=ymin, ymax=ymax),
  #           data = raster2rect(df),
  #           color = NA) +
  facet_grid(.~Condition) +
  ## Theme, Scales and Labels
  # scale_y_log10() + # if Frequencies are log-spaced
  scale_fill_distiller(palette = "OrRd", limits = range(df$itc), oob=squish)

# show the plots!
ERSP.plot
ITC.plot
