library(tidyverse) # required to work


# load the data file
grand_data <- read_csv("@filename@") %>%
  mutate(
    across(c(Time, mean, @error@), as.numeric),
    Condition = factor(Condition, levels = unique(Condition))
  )

p_grand <- grand_data %>% 
  ggplot(aes(x = Time, y = mean,
             color = Condition, fill = Condition,
             group = Condition)) +
  ## Axes
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ## Plot ERPs
  geom_ribbon(aes(ymin = mean - @error@, ymax = mean + @error@), # error 
              alpha = 0.25, color = NA) + 
  geom_line() + # plot grand ERPs
  ## Theme, Scales and Labels
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
  labs(x = "Time", y = "Value")

p_grand
