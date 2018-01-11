library(tidyverse) # required to work

alpha <- 0.95 # For confidence interval

erp.plot.data <- read_csv("@filename@") %>%    # load the data file
  mutate_at(.vars = c("Time","mean","sd"), .funs = as.numeric) %>% 
  mutate(Condition = as.factor(Condition),     # set as factor
         N = as.integer(N),                    # set as integer
         se = sd/sqrt(N),                      # calculate SE
         ci = se * qt(0.5 + alpha/2 ,N-1))     # calculate CI width

ERPplot <- erp.plot.data %>% 
  ggplot(aes(x = Time, y = mean,
             color = Condition, fill = Condition,
             group = Condition))

ERPplot +
  geom_vline(xintercept = 0) + # add x-axis at amp = 0
  geom_hline(yintercept = 0) + # add y-axis lines at time = 0
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se),
              alpha = 0.25, color = NA) + # Plot SE (can be changed to SD or CI)
  geom_line() + # plot grand ERPs
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0, 0)) +   # remove padding around time line axis
  labs(x = "Time", y = expression(paste(mu,"V"))) # add axis labels
