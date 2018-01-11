library(tidyverse) # required to work

alpha <- 0.95 # For confidence interval

erp.plot.data <- read.csv("@filename@") %>%      # load the data file
  as.tbl() %>% 
  mutate(Condition = as.factor(Condition),       # set as factor
         Time      = as.double(Time),            # set as numeric
         N         = as.integer(N),              # set as integer
         mean      = as.double(mean),            # set as numeric
         sd        = as.double(sd),              # set as numeric
         se        = sd/sqrt(N),                 # calculate SE
         ci        = se * qt(0.5 + alpha/2 ,N-1) # calculate CI width
  )

ERPplot <- ggplot(erp.plot.data, aes(x = Time, y = mean,
                                     group = Condition,
                                     color = Condition)) +
  ## add axis lines at time = 0, and amp = 0
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ## plot the grands and errors
  geom_ribbon(aes(ymin = mean-se, ymax = mean+se, fill = Condition),
              alpha = 0.25, color = NA) +  
  geom_line() +
  ## Design
  # scale_y_reverse() + # minus up?
  scale_x_continuous(expand = c(0,0)) +
  labs(x = "Time",
       y = expression(paste(mu,"V")))

ERPplot # show the plot!
