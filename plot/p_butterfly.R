library(tidyverse)

butterfly_data <- read_csv("@filename@") %>% # load the data file
  group_by(Condition, Time) %>%
  mutate(meanAmp = mean(amp)) %>%
  ungroup()

p_butterfly <- butterfly_data %>% 
  ggplot(aes(x = Time, y = amp,
             color = Condition,
             group = ID)) + 
  facet_wrap( ~ Condition) +
  ## Axes
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ## Plot ERPs
  geom_line(alpha = 0.8) + # single subject ERPs
  geom_line(aes(y = meanAmp), color = "black", size = 1) + # grand
  ## Theme, Scales and Labels
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
  labs(x = "Time", y = "Value")

p_butterfly
