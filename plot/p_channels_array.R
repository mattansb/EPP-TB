library(tidyverse)

chanarray_data <- read_csv("@filename@")

# Plot With Locations -----------------------------------------------------

annotate_channel <- function(data, w = 0.05, h = 0.08){
  .p <- ggplot(data,aes(Time,amp, color = Condition, group = Condition)) + 
    ## Axes
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ## Plot ERPs
    geom_line() + 
    ## Theme, Scales and Labels
    # scale_y_reverse() + # minus up?
    scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
    theme_minimal() +
    guides(color = FALSE) + 
    theme_void() + 
    labs(title = data$Channel[1])
  
  .g <- annotation_custom(
    grob = ggplotGrob(.p),
    xmin = data$x[1] - w,
    xmax = data$x[1] + w,
    ymin = data$y[1] - h,
    ymax = data$y[1] + h
  )
  return(.g)
}



p_chanarray <- ggplot(mapping = aes(x = 0, y = 0)) + {
  chanarray_data %>% 
    split(.$Channel) %>% 
    map(annotate_channel)
} +
  theme_void()

p_chanarray
