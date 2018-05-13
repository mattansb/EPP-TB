library(tidyverse) # required to work

alpha <- 0.95 # confidence level

erp.plot.data <- read_csv("@filename@") %>%    # load the data file
  as.tbl() %>% 
  mutate_at(.vars = c("Time","mean","sd"), .funs = as.numeric) %>% 
  mutate(Condition = as.factor(Condition),     # set as factor
         N  = as.integer(N),                   # set as integer
         se = sd/sqrt(N),                      # calculate SE
         ci = se * qt(0.5 + alpha/2 ,N-1))     # calculate CI width

ERPplot <- erp.plot.data %>% 
  ggplot(aes(x = Time, y = mean,
             color = Condition, fill = Condition,
             group = Condition)) +
  ## Axes
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ## Plot ERPs
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se), # error (can be changed to se, sd or ci)
              alpha = 0.25, color = NA) + 
  geom_line() + # plot grand ERPs
  ## Theme, Scales and Labels
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) + # remove padding around time line axis
  labs(x = "Time", y = expression(paste(mu,"V")))

ERPplot