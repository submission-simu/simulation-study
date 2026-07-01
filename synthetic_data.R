rm(list = ls())
set.seed(2026)

K <- 5
cov_dim <- 6

beta_star <- c(0.35, -0.65, 0.35, 0.16, -0.14, -0.6)
# see Table 2 in the manuscript
site_size <- c(2852, 1033, 842, 461, 240)
mean_mat <- matrix(c(2637, 968, 798, 370, 204,
                     680, 234, 176, 181, 68,
                     377, 144, 130, 51, 25,
                     902, 383, 345, 145, 72,
                     1053, 588, 399, 101, 48,
                     185, 73, 58, 102, 17), nrow = cov_dim, ncol = K, byrow = TRUE)
mean_mat[2, ] <- site_size - mean_mat[2, ]


mean_vec <- c(0.9, 0.25, 0.15, 0.25, 0.3, 0.1)

n_cum <- c(0, cumsum(site_size))
shape_seq <- seq(2.5, 1.5, length.out = K)
scale_seq <- seq(8, 10, length.out = K)
censor_low <- rep(2, K)
censor_high <- rep(4.5, K)
data_list <- as.list(rep(NA, K))

for(k in 1:K){
  covariate_mat <- matrix(0, nrow = site_size[k], ncol = cov_dim)
  covariate_mat[sample(1:site_size[k], size = mean_mat[1, k]), 1] <- 1
  covariate_mat[sample(1:site_size[k], size = mean_mat[2, k]), 2] <- 1
  covariate_mat[sample(1:site_size[k], size = mean_mat[3, k]), 3] <- 1
  covariate_mat[sample(1:site_size[k], size = mean_mat[4, k]), 4] <- 1
  covariate_mat[sample(1:site_size[k], size = mean_mat[5, k]), 5] <- 1
  covariate_mat[sample(1:site_size[k], size = mean_mat[6, k]), 6] <- 1
  
  hazard_scaling <- c(covariate_mat %*% beta_star)
  baseline_t2e <- rweibull(site_size[k], shape = shape_seq[k], scale = scale_seq[k])
  Cox_t2e <- baseline_t2e * exp(-hazard_scaling / shape_seq[k])
  censor <- runif(site_size[k], min = censor_low[k], max = censor_high[k])
  
  
  
  local_data <- data.frame(ID = n_cum[k] + c(1:site_size[k]),
                           site = rep(k, site_size[k]),
                           t2e = pmin(Cox_t2e, censor),
                           censor_ind = as.numeric(Cox_t2e <= censor),
                           V1 = covariate_mat[, 1],
                           V2 = covariate_mat[, 2],
                           V3 = covariate_mat[, 3],
                           V4 = covariate_mat[, 4],
                           V5 = covariate_mat[, 5],
                           V6 = covariate_mat[, 6])
  
  data_list[[k]] <- local_data
}

save(data_list, file = "real_data.rdata")
