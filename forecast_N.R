
library(MCMCvis)
library(foreach)
library(tidyverse)
library(truncnorm)

data_pop <- readRDS("data/data_pop_empe.rds")
env <- readRDS("data/data_coupled_transformed.rds")

# Site id info
data_site1 <- readRDS("data/data_env_empe.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  arrange(site_id) %>%
  select(site_id) %>%
  unique()

data_site2 <- readRDS("data/data_env_empe.rds") %>%
  filter(year >= 2009 & year <= 2018) %>%
  arrange(site_id) %>% 
  select(site_id) %>%
  unique()

# Parameter chains 
res_sac <- read_rds("results/results_sac.rds")

param_chains <- MCMCchains(res_sac, params = c("alpha", "beta", "sigma"))

set.seed(19)
idx <- sample(1:nrow(param_chains), 100)

b0 <- param_chains[idx,"alpha"]
b1 <- param_chains[idx,"beta"]

# Site effects
eps <- MCMCchains(res_sac, params = c("eps"), exact = F)
eps <- eps[idx,]
colnames(eps) <- data_site1$site_id
 
## Generate 0 site effects for colonies not in the original data 
#s_s <- rtruncnorm(100, a = 0, mean = 0, sd = param_chains[idx,"sigma"])
#e_s <- rnorm(100*16, 0, s_s) %>% 
  #matrix(ncol = 16, nrow = 100)
e_s <- matrix(0, ncol = 16, nrow = 100)
#idx_site <- which(!unique(data_site2$site_id) %in% data_site1$site_id)
#colnames(e_s) <- unique(data_site2$site_id)[idx_site]

## Combine site effects
eps2 <- cbind(eps, e_s)
idx_site2 <- order(colnames(eps2)) 
eps2 <- eps2[,idx_site2]

# Forecast and hindcast
N <- foreach(i = 1:length(env)) %do% {
  
  x <- env[[i]]
  x[x == 0] <- 0.0001
  x <- log(x)
  
  foreach(h = 1:100) %:% 
    foreach(k = 1:nrow(x), .combine = "rbind") %do% {
      b0[h] + b1[h]*x[k,] + eps2[h,k]
    }
}

N_global <- foreach(i = 1:50) %:%
  foreach(k = 1:100) %do% {
    z <- exp(N[[i]][[k]])
    z[is.infinite(z)] <- 0
    apply(z, 2, sum)
  }

N_global <- map(N_global, function(x) do.call(rbind, x))

# Plots
plot_empe <- function(y_pred_chains) {
  
  chains_dat <- foreach(i = 1:50, .combine = "rbind") %do%
    y_pred_chains[[i]]
  
  y_max <- apply(chains_dat, 2, quantile, 0.995)
  y_min <- apply(chains_dat, 2, quantile, 0.005)
  
  dat_ci <- data.frame(year = 1909:2100,
                       max = y_max,
                       min = y_min)
  
  y_mean <- foreach(i = 1:50, .combine = "c") %do%
    apply(y_pred_chains[[i]], 2, mean)
  
  dat_mean <- data.frame(y_mean = y_mean,
                         year = rep(1909:2100, 50),
                         run = rep(as.factor(1:50), each = 192))
  
  theme_set(theme_bw())
  ggplot() +
    geom_ribbon(data = dat_ci, mapping = aes(x = year, ymin = min, ymax = max), 
                fill = "blue4", alpha = 0.3) +
    geom_line(data = dat_mean, mapping = aes(x = year, y = y_mean, col = run), 
              alpha = 0.8, linewidth = 0.2) +
    #geom_hline(yintercept = 167000, linetype = 2) +
    #geom_vline(xintercept = 2087, linetype = 2) +
    labs(y = "Abundance") +
    theme(legend.position = "none",
          panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 8)) +
    scale_color_hue(l = 50) +
    #scale_y_continuous(limits = c(50000, 350000)) +
    scale_x_continuous(breaks = c(1909, 1950, 2000, 2050, 2100))
}

plot_empe(N_global)

ggsave("fig_forecast.pdf", width = 6, height = 5, units = "in") 
