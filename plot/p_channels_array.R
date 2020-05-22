library(tidyverse)

df <- read_csv("@filename@")


# Plot Grid ---------------------------------------------------------------

chan_plot <- ggplot(df,aes(Time,amp, color = Condition, group = Condition)) + 
  ## Axes
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ## Plot ERPs
  geom_line() + 
  ## Theme, Scales and Labels
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
  labs(x = "Time", y = "Value") +
  facet_wrap(~Channel)

chan_plot


# Plot With Locations -----------------------------------------------------

by_chanel <- split(df,df$Channel)

annotate_channel <- function(df, w = 0.05, h = 0.08){
  p <- ggplot(df,aes(Time,amp, color = Condition, group = Condition)) + 
    ## Axes
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ## Plot ERPs
    geom_line() + 
    ## Theme, Scales and Labels
    # scale_y_reverse() + # minus up?
    scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
    theme_minimal() +
    guides(color=FALSE) + 
    theme_void() + 
    labs(title = df$Channel[1])
  
  x_min <- df$x[1] - w
  y_min <- df$y[1] - h
  x_max <- df$x[1] + w
  y_max <- df$y[1] + h
  return(annotation_custom(grob = ggplotGrob(p), xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max))
}

chan_plot <- ggplot(mapping = aes(x = 0,y = 0)) + 
  map(by_chanel,annotate_channel) + 
  theme_void()

chan_plot
