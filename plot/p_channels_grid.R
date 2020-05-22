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

