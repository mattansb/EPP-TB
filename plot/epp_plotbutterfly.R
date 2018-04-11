library(tidyverse)

butterfly.data <- read_csv("@filename@") %>% # load the data file
  group_by(Condition,Time) %>% 
  mutate(meanAmp = mean(amp)) %>% 
  ungroup()

butterfly_plot <- butterfly.data %>% 
  ggplot(aes(x = Time, y = amp,
             color = Condition,
             group = ID)) + 
  facet_wrap(~Condition)

butterfly_plot +
  geom_vline(xintercept = 0) + # add x-axis at amp = 0
  geom_hline(yintercept = 0) + # add y-axis lines at time = 0
  geom_line(alpha = 0.8) + # plot single subject ERPs
  geom_line(aes(y = meanAmp), color = "black", size = 1) + # plot grand ERPs
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) +   # remove padding around time line axis
  labs(x = "Time", y = expression(paste(mu,"V"))) # add axis labels
