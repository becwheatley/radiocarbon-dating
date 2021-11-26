library(rcarbon)

theme_set(theme_classic())

mt <- readRDS("uniform_50p_vs_uniform.Rds")

# extract and combine data
spd.sims <- mt$sim
colnames(spd.sims) <- paste0("s",1:mt$nsim)
spd.data <- spd.sims %>% as_tibble %>% select(-s100) %>% # use 99 sims for convenience
  mutate(calBP = mt$result$calBP, obs = mt$result$PrDens)
nsim <- mt$nsim -1 # to account for deleted sim
spd.data.long <- spd.data %>% pivot_longer(-calBP, names_to = "rep", values_to = "pd")

# calculate p-values
spd.data$calBP %>% head
p_sig <- 0.05
dat_sig_p <- 
  spd.data%>% 
  select(-calBP) %>% 
  apply(1, function(row) (row - mean(row))^2) %>% # squared area
  apply(1, cumsum) %>% # cumulative area, i.e., integral
  apply(1, function(row) ((nsim + 2) - rank(row))/(nsim + 1)) %>% # progressive p-value
  {tibble(calBP = spd.data$calBP, p =.["obs",], pd = spd.data$obs)} %>% 
  mutate(sig = p<=p_sig, grp = c(0,diff(sig))) %>% 
  filter(sig) %>% 
  mutate(sig = NULL, grp = cumsum(grp))

p_total <- dat_sig_p %>% slice(n()) %>% pull(p)

# compute total discrepancy
spd.data.long %>% group_by(calBP) %>% summarise(pd = mean(pd)) %>% 
  mutate(obs = spd.data$obs, 
         dis = discrepancy_fun(pd,obs), 
         dis =cumsum(dis)) %>% 
  slice(n()) %>% # take the last value to get the total
  pull(dis)

discrepancy_fun <- function(x, y, q = 1, p =2) (abs(x^q - y^q))^p

dat_sig_p %>% ggplot(aes(-calBP,p)) + geom_line() + ylim(0,0.06) + 
  geom_hline(aes(yintercept = p_sig), lty = "dashed") +
  labs(subtitle = "p-values computed cumulatively from left to right")

spd.data.long %>% 
  ggplot(aes(-calBP,pd)) +
  geom_line(aes(group = rep), alpha = 0.2, color = blues9[5], size = 0.5) +
  geom_line(col = "grey10", size = 1, lty = "dashed",
            data = ~.x %>% group_by(calBP) %>% summarise(pd = mean(pd))) +
  geom_line(col = blues9[9], size = 0.8,
            data = ~ .x %>% filter(rep == "obs")) +
  geom_line(aes(group = grp), col = "red", data = dat_sig_p, size = 0.9) +
  labs(title = "Summed probabilty distribution",
       subtitle= paste0("Solid line is the emprical curve (red if p < ", p_sig,", blue otherwise)",
                        "\n","p-value for the total data set is ", p_total))
#ggsave("spd_1.png")

