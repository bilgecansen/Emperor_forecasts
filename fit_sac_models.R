
library(foreach)
library(tidyverse)
library(MCMCvis)
library(rstan)

empe_sites <- read.csv("data/empe_sites_new.csv")
sites_update <- read.csv("data/colony_attributes_update.csv")

data_pop <- readRDS("data/data_pop_empe.rds")

theme_set(theme_bw())


# Load environmental data -------------------------------------------------

sites_all <- unique(readRDS("data/data_env_empe.rds")$site_id)

data_esm <- readRDS("data/data_env_empe.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("aice"))) %>%
  arrange(site_id)

data_env_avg <- data_esm %>%
  group_by(site_id) %>%
  summarise(across(-year, mean)) %>%
  ungroup()

env_mat <- data_env_avg[,-1] %>%
  as.matrix()

env_mat_std <- apply(env_mat, 2, function(x) (x - mean(x))/sd(x))


# Fit linear regression ---------------------------------------------------

#N_chains <- readRDS("data/N_chains.rds")
#N_mean_chains <- foreach(i = 1:50, .combine = "cbind") %do% {
  #z <- t(apply(N_chains[,seq(i, 500, 50)], 1, function(x) {
    #as.vector(mean(x[which(x != 0)]))
  #}))
  
  #as.vector(z)
#}

#N_mean <- apply(N_mean_chains, 2, mean)
#N_sd <- apply(N_mean_chains, 2, sd)
#N_min <- apply(N_mean_chains, 2, quantile, 0.05)
#N_max <- apply(N_mean_chains, 2, quantile, 0.95)

#N_data <- data.frame(mean = N_mean,
                     #sd = N_sd,
                     #min = N_min,
                     #max = N_max,
                     #site_id = sites_update$site_id)
#N_data <- N_data[order(N_data$site_id),]

N_data <- readRDS("data/N_data.rds")

x1 <- log(env_mat[,1])
x2 <- env_mat_std[,1]

dat_lm <- list(y = N_data$mean,
               y_sd = N_data$sd,
               N = 50,
               X = x1,
               X_pred = seq(min(x1), max(x1), length.out = 100),
               S = nrow(env_mat),
               L = 100)

res_lm <- stan(file = 'lm2.stan', 
                 data = dat_lm,
                 iter = 5000,
                 cores = 4,
                 control = list(adapt_delta = 0.9999, 
                                max_treedepth = 20))

if (!"results"  %in% list.files()) dir.create("results")
saveRDS(res_lm, "results/results_sac.rds")

# Create posterior csv files
param_chains <- MCMCchains(res_lm, params = c("alpha", "beta"))
write.csv(param_chains, "results/param_chains.csv")

## Site effects
eps <- MCMCchains(res_sac, params = c("eps"), exact = F)
colnames(eps) <- unique(data_esm$site_id)

e_s <- matrix(0, ncol = 16, nrow = nrow(eps))
idx_site <- which(!sites_all %in% unique(data_esm$site_id))
colnames(e_s) <- sites_all[idx_site]

## Combine site effects
eps2 <- cbind(eps, e_s)
idx_site2 <- order(colnames(eps2)) 
eps2 <- eps2[,idx_site2]
write.csv(eps2, "results/eps_chains.csv")

# Model with standardization
dat_lm2 <- list(y = N_data$mean,
               y_sd = N_data$sd,
               N = 50,
               X = x2,
               X_pred = seq(min(x2), max(x2), length.out = 100),
               S = nrow(env_mat),
               L = 100)

res_lm2 <- stan(file = 'lm2.stan', 
               data = dat_lm2,
               iter = 4000,
               cores = 4)

# Model fit
MCMCsummary(res_lm, params = "Rsq")

# Univariate plots
y1 <- MCMCsummary(res_lm, params = "mu_pred")
y1_mean <- y1[,1]
y1_min <- y1[,3]
y1_max <- y1[,5]

x_pred1 <- seq(min(x1), max(x1), length.out = 100)

ggplot() +
  geom_ribbon(aes(x = x_pred1, ymin = y1_min, ymax = y1_max), alpha = 0.8, 
              fill = "grey") +
  geom_segment(aes(x = env_mat[,1], xend = env_mat[,1], 
                   y = N_data$min, yend = N_data$max)) +
  geom_line(aes(x = x_pred1, y = y1_mean), col = "darkorange", size = 1.5, 
            linetype = 2) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), alpha= 0.8, 
             color = "darkorange", size = 3) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), shape = 1, size = 3, 
             stroke = 1.1) +
  labs(y = "Colony Abundance (log)", x = "SIC (Laying)") +
  theme(#axis.title.y = element_blank(),
    panel.border = element_blank(), 
    #panel.grid.major = element_blank()) +
    panel.grid.minor = element_blank())

y2 <- MCMCsummary(res_lm2, params = "mu_pred")
y2_mean <- y1[,1]
y2_min <- y1[,3]
y2_max <- y1[,5]

x_pred2 <- seq(min(x2), max(x2), length.out = 100)*sd(env_mat[,1]) + 
  mean(env_mat[,1])

ggplot() +
  geom_ribbon(aes(x = x_pred2, ymin = y2_min, ymax = y2_max), alpha = 0.8, 
              fill = "grey") +
  geom_segment(aes(x = env_mat[,1], xend = env_mat[,1], 
                   y = N_data$min, yend = N_data$max)) +
  geom_line(aes(x = x_pred2, y = y2_mean), col = "darkorange", size = 1.5, 
            linetype = 2) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), alpha= 0.8, 
             color = "darkorange", size = 3) +
  geom_point(aes(x = env_mat[,1], y = N_data$mean), shape = 1, size = 3, 
             stroke = 1.1) +
  labs(y = "Colony Abundance (log)", x = "SIC (Laying)") +
  theme(#axis.title.y = element_blank(),
    panel.border = element_blank(), 
    #panel.grid.major = element_blank()) +
    panel.grid.minor = element_blank())
