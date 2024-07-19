
library(foreach)
library(tidyverse)
library(abind)

empe_sites <- read.csv("data/empe_sites_new.csv")
data_pop <- readRDS("data/data_pop_empe.rds")


# Load environmental data -------------------------------------------------

data_esm <- readRDS("data/data_env_empe.rds") %>%
  filter(site_id %in% unique(data_pop$sat$site_id)) %>%
  filter(year >= 2009 & year <= 2018) %>%
  select(site_id, year, contains(c("aice"))) %>%
  arrange(site_id)

data_aice_coupled <- readRDS("data/data_aice_coupled.rds") %>%
  map(., function(x) arrange(x, site_id)) %>%
  map(., function(x) rename(x, year = season))

get_env_mat <- function(data_env, x) {
  select(data_env, site_id, year, x) %>%
    filter(site_id %in% unique(data_pop$sat$site_id)) %>%
    pivot_wider(names_from = year, values_from = x) %>%
    select(-site_id) %>%
    as.matrix()
}

env_mat <- get_env_mat(data_esm, "aice_laying")
env_mat_avg <- apply(env_mat, 1, mean)
x_mean <- mean(env_mat_avg)
x_sd <- sd(env_mat_avg)

env_mat_fore <- foreach(h = 1:length(data_aice_coupled)) %do% 
  get_env_mat(data_aice_coupled[[h]], "aice_laying")

# correct env forecast data
trans_dat <- function(x, y) {
  
  mu_obs <- mean(x)
  sigma_obs <- sd(x)
  sigma_model <- sd(y[110:119])
  mu_model <- mean(y[110:119])
  
  (sigma_obs / sigma_model) * (y - mu_model) + mu_obs
}

env_mat_fore <- foreach(i = 1:length(env_mat_fore)) %:%
  foreach(k = 1:nrow(env_mat_fore[[i]]), .combine = "rbind") %do% {
    trans_dat(env_mat[k,], env_mat_fore[[i]][k,])
  }

for (i in 1:50) {
  env_mat_fore[[i]][env_mat_fore[[i]] < 0] <- 0
  env_mat_fore[[i]][env_mat_fore[[i]] > 1] <- 1
}   

# moving window average
env_mat_fore_avg <- foreach (i = 1:50) %:% 
  foreach(k = 1:50, .combine = "rbind") %:%
  foreach (h = 1:192, .combine = "c") %do% {
    mean(env_mat_fore[[i]][k,h:(h+9)])
  }

# make sure forecasts are within data bounds
#for (i in 1:50) {
  #env_mat_fore_avg[[i]][env_mat_fore_avg[[i]] < min(env_mat_avg)] <- 
    #min(env_mat_avg)
  #env_mat_fore_avg[[i]][env_mat_fore_avg[[i]] > max(env_mat_avg)] <- 
    #max(env_mat_avg)
#} 

for (i in 1:50) {
  colnames(env_mat_fore_avg[[i]]) <- 1909:2100
  rownames(env_mat_fore_avg[[i]]) <- unique(data_esm$site_id)
}

saveRDS(env_mat_fore_avg, "data/data_coupled_transformed.rds")

# standardize the data
env_mat_fore_avg_std <- foreach (i = 1:50) %do% {
  (env_mat_fore_avg[[i]] - x_mean)/x_sd
}

saveRDS(env_mat_fore_avg_std, "data/data_coupled_normalized.rds")
