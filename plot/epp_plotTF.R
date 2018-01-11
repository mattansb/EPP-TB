library(tidyverse)
library(scales)


# ==== Load data ====

df <- read_csv("@filename@")

raster2rect <- function(data) {
  freq.vals <- unique(data$Frequencies)
  
  d.time <- diff(unique(data$Time))[1]/2
  d.frq <- diff(freq.vals)/2
  
  upper <- freq.vals + c(d.frq, d.frq[length(d.frq)])
  lower <- freq.vals - c(d.frq[1], d.frq) 
  ntime <- length(unique(data$Time))
  
  data$xmin <- data$Time - d.time
  data$xmax <- data$Time + d.time
  data$ymin <- unlist(lapply(lower, function(i) rep(i, ntime)))
  data$ymax <- unlist(lapply(upper, function(i) rep(i, ntime)))
  
  return(data)
}


# Plot --------------------------------------------------------------------

# Use geom_raster to plot y-axis scale as is
TF.plot <- ggplot(df,aes(Time, Frequencies, fill = ersp)) + # or fill = itc
  geom_raster(interpolate = T) +
  facet_grid(.~Condition) +
  # scale_y_log10() + # if Frequencies are log-spaced
  scale_fill_distiller(palette = "OrRd", limits = c(-1.5,1.5),oob=squish)

# Use geom_rect to change scaling of y-axis
TF.plot <- ggplot(raster2rect(df),aes(Time, Frequencies, fill = ersp)) +  # or fill = itc
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color = NA) +
  facet_grid(.~Condition) +
  # scale_y_log10() + # if Frequencies are log-spaced
  scale_fill_distiller(palette = "OrRd", limits = c(-1.5,1.5),oob=squish)

TF.plot # show the plot!
