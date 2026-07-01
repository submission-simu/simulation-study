rm(list = ls())
library("survival")
library("ggplot2")
library("xtable")
load("real_data.rdata")

K <- 5
cov_dim <- 6

lower_a <- rep(-5, cov_dim)
upper_b <- rep(5, cov_dim)

data_process <- function(data_list, covariate_names){
  # for each site, pre-calculate the failure time points, the risk sets, and the respective weights
  # also, remove those sites with zero events
  site_to_remove <- c()
  K <- length(data_list)
  failure_num <- rep(NA, K)
  failure_times <- as.list(rep(NA, K))
  risk_sets <- as.list(rep(NA, K))
  covariate_list <- as.list(rep(NA, K))
  failure_position <- as.list(rep(NA, K))
  for(k in 1:K){
    # prepare a list for covariates in matrix format so as to speed up computation of log partial likelihood, gradient, and hessian
    covariate_list[[k]] <- as.matrix(data_list[[k]][, covariate_names, drop = FALSE]) 
    # find over which position lies the failure times
    failure_position[[k]] <- which(data_list[[k]]$censor_ind == 1)
    # find failure times
    failure_times[[k]] <- data_list[[k]]$t2e[failure_position[[k]]]
    # the number of failures
    failure_num[k] <- length(failure_times[[k]])
    
    if(failure_num[k] == 0){
      site_to_remove <- c(site_to_remove, k)
    }else{
      temp_risk <- as.list(rep(NA, failure_num[k]))
      for(j in 1:failure_num[k]){
        temp_risk[[j]] <- which(data_list[[k]]$t2e >= failure_times[[k]][j])
      }
      risk_sets[[k]] <- temp_risk
    }
  }
  
  if(length(site_to_remove) > 0){
    data_list <- data_list[-site_to_remove]
    covariate_list <- covariate_list[-site_to_remove]
    failure_num <- failure_num[-site_to_remove]
    failure_times <- failure_times[-site_to_remove]
    failure_position <- failure_position[-site_to_remove]
    risk_sets <- risk_sets[-site_to_remove]
    K <- K - length(site_to_remove)
  }
  return(list(data_list = data_list,
              covariate_list = covariate_list,
              failure_position = failure_position,
              failure_num = failure_num,
              risk_sets = risk_sets,
              K = K))
}


log_plk <- function(beta, covariate = NA, 
                    failure_position = NA,
                    failure_num = NA,
                    risk_sets = NA){
  
  eta <- drop(covariate %*% beta)
  res <- sum(eta[failure_position])
  
  exp_eta <- exp(eta)
  
  res - sum(vapply(risk_sets[seq_len(failure_num)],
                   function(idx) log(sum(exp_eta[idx])),
                   numeric(1)))
}


grad_plk <- function(beta, 
                     covariate = NA, 
                     failure_position = NA,
                     failure_num = NA,
                     risk_sets = NA){
  
  eta <- drop(covariate %*% beta)
  exp_eta <- exp(eta)
  
  # numerator part: sum of covariates at failures
  # always keep matrix to make colSums fast and safe
  res <- colSums(covariate[failure_position, , drop = FALSE])
  
  # subtract weighted mean covariates for each risk set
  for (j in seq_len(failure_num)) {
    idx <- risk_sets[[j]]
    if (length(idx) == 1L) {
      res <- res - covariate[idx, ]
    } else {
      w  <- exp_eta[idx]
      sw <- sum(w)
      # t(X_idx) %*% w  (p x 1), returned as numeric vector
      wx <- drop(crossprod(covariate[idx, , drop = FALSE], w))
      res <- res - wx / sw
    }
  }
  
  res
}


hess_plk <- function(beta,
                     p = NA,
                     covariate = NA, 
                     failure_num = NA,
                     risk_sets = NA){
  
  eta <- drop(covariate %*% beta)
  exp_eta <- exp(eta)
  
  res <- matrix(0.0, p, p)
  
  for (j in seq_len(failure_num)) {
    idx <- risk_sets[[j]]
    if (length(idx) > 1L) {
      w  <- exp_eta[idx]
      sw <- sum(w)
      
      Xidx <- covariate[idx, , drop = FALSE]
      
      # wx = sum_i w_i x_i  (p-vector)
      wx <- drop(crossprod(Xidx, w))
      
      # sum_i w_i x_i x_i^T  = crossprod( Xidx * sqrt(w) )
      Xw <- Xidx * sqrt(w)
      EXX <- crossprod(Xw) / sw
      
      # (sum w x)(sum w x)^T / (sum w)^2
      outer_term <- tcrossprod(wx) / (sw * sw)
      
      # your original form: outer_term - EXX  (negative covariance)
      res <- res + outer_term - EXX
    }
  }
  
  res
}

third_plk <- function(beta,
                      p,
                      covariate,
                      failure_num,
                      risk_sets) {
  eta <- drop(covariate %*% beta)
  exp_eta <- exp(eta)
  
  # third derivative tensor: p x p x p
  T3 <- array(0.0, dim = c(p, p, p))
  
  for (j in seq_len(failure_num)) {
    idx <- risk_sets[[j]]
    m <- length(idx)
    if (m <= 1L) next
    
    w <- exp_eta[idx]
    sw <- sum(w)
    prob <- w / sw  # softmax weights in this risk set
    
    X <- covariate[idx, , drop = FALSE]  # m x p
    
    # weighted mean mu (p-vector)
    mu <- drop(crossprod(prob, X))       # 1xp -> p
    
    # centered X
    Xc <- sweep(X, 2, mu, "-")           # m x p
    
    # Contribution of this risk set to ∇^3 log(sum exp) is E[(X-μ)^{⊗3}]
    # But for partial log-likelihood, we subtract it: T3 <- T3 - that
    for (k in seq_len(p)) {
      wk <- prob * Xc[, k]               # length m
      # crossprod: p x m  %*%  m x p  -> p x p
      # element (a,b) = sum_i Xc[i,a] * (Xc[i,b] * wk[i])
      Mk <- crossprod(Xc, Xc * wk)
      T3[, , k] <- T3[, , k] - Mk
    }
  }
  
  T3
}

mini_SL_fun <- function(beta, covariate = NA, 
                        failure_position = NA,
                        failure_num = NA,
                        risk_sets = NA,
                        hess_external = NA,
                        initial_beta = NA){
  
  res <- 0.5 * t(beta - initial_beta) %*% hess_external %*% (beta - initial_beta)
  res <- res + log_plk(beta, covariate = covariate,
                       failure_position = failure_position,
                       failure_num = failure_num,
                       risk_sets = risk_sets)
  return(res)
}

ODACH <- function(beta, covariate = NA, 
                  failure_position = NA,
                  failure_num = NA,
                  risk_sets = NA,
                  grad_external = NA,
                  hess_external = NA,
                  initial_beta = NA){
  res <- sum(grad_external*beta) + 0.5 * t(beta - initial_beta) %*% hess_external %*% (beta - initial_beta)
  res <- res + log_plk(beta, covariate = covariate,
                       failure_position = failure_position,
                       failure_num = failure_num,
                       risk_sets = risk_sets)
  return(res)
}

aug_SL <- function(beta, covariate = NA,
                   failure_position = NA,
                   failure_num = NA,
                   risk_sets = NA,
                   est_list = NA,
                   grad_list = NA,
                   hess_list = NA,
                   lead_site = NA){
  
  K <- length(est_list)
  idx <- setdiff(seq_len(K), lead_site)
  res <- 0
  for(int_i in idx){
    res <- res + 0.5 * t(beta - est_list[[int_i]]) %*% hess_list[[int_i]] %*% (beta - est_list[[int_i]]) + sum(grad_list[[int_i]] * (beta - est_list[[int_i]])) 
  }
  res <- c(res) + log_plk(beta, covariate = covariate,
                          failure_position = failure_position,
                          failure_num = failure_num,
                          risk_sets = risk_sets)
  return(res)
}

# T: array with dim (p, p, p)
# v: length-p vector
# returns: p x p matrix with M[j,k] = sum_l T[j,k,l] * v[l]
contract_3rd <- function(T, v) {
  p <- dim(T)[1L]
  # unfold (j,k,l) into ((j,k), l)
  T_unf <- matrix(aperm(T, c(1,2,3)), nrow = p*p, ncol = p)
  M_vec <- drop(T_unf %*% v)           # length p*p
  matrix(M_vec, nrow = p, ncol = p)   # back to p x p
}


one_step_polish <- function(current_est, grad_list = NA, hess_list = NA, tensor_list = NA, est_list = NA,
                            coordination = "L", covariate = NA, failure_position = NA, failure_num = NA,
                            risk_sets = NA){
  effective_K <- length(est_list)
  cov_dim <- length(est_list[[2]])
  
  if(coordination == "L"){
    term1 <- hess_plk(beta = current_est, 
                      p = cov_dim,
                      covariate = covariate,
                      failure_num = failure_num,
                      risk_sets = risk_sets)
    
    term2 <- grad_plk(beta = current_est, 
                      covariate = covariate,
                      failure_position = failure_position,
                      failure_num = failure_num,
                      risk_sets = risk_sets)
    
    for (my_i in 2:effective_K) {
      delta <- current_est - est_list[[my_i]]
      
      Tm <- contract_3rd(tensor_list[[my_i]], delta)  # p x p
      
      Hi <- hess_list[[my_i]]
      gi <- grad_list[[my_i]]
      
      term1 <- term1 + Hi + Tm
      term2 <- term2 + gi + Hi %*% delta + 0.5 * (Tm %*% delta)
    }
    
  }else{
    term1 <- matrix(0, cov_dim, cov_dim)
    term2 <- rep(0, cov_dim)
    
    for (my_i in 1:effective_K) {
      delta <- current_est - est_list[[my_i]]
      
      Tm <- contract_3rd(tensor_list[[my_i]], delta)  # p x p
      
      Hi <- hess_list[[my_i]]
      
      term1 <- term1 + Hi + Tm
      term2 <- term2 + Hi %*% delta + 0.5 * (Tm %*% delta)
    }
  }
  
  OS_est <- c(current_est - solve(term1) %*% term2)
  return(OS_est)
}



pre_processing <- data_process(data_list, covariate_names = c("V1", "V2", "V3", "V4", "V5", "V6"))

data_list <- pre_processing$data_list
effective_K <- pre_processing$K

pooled_df <- data_list[[1]]
for(k in 2:effective_K){
  pooled_df <- rbind(pooled_df, data_list[[k]])
}

fit_pool <- coxph(Surv(t2e, censor_ind) ~ V1 + V2 + V3 + V4 + V5 + V6 + strata(site), data = pooled_df, method = "breslow")
pool_est <- as.numeric(fit_pool$coefficients)
pool_var <- fit_pool$var


### proposed method coordinated by lead site
lead_est <- optim(par = rep(0, cov_dim), fn = log_plk, method = "L-BFGS-B", control = list(fnscale = -1),
                  covariate = pre_processing$covariate_list[[1]],
                  failure_position = pre_processing$failure_position[[1]],
                  failure_num = pre_processing$failure_num[[1]],
                  risk_sets = pre_processing$risk_sets[[1]], lower = lower_a, upper = upper_b)$par


lead_hess <- hess_plk(beta = lead_est, 
                      p = cov_dim,
                      covariate = pre_processing$covariate_list[[1]],
                      failure_num = pre_processing$failure_num[[1]],
                      risk_sets = pre_processing$risk_sets[[1]])

lead_var <- -solve(lead_hess)


mini_SL_list <- as.list(rep(NA, effective_K))
grad_mini_SL <- as.list(rep(NA, effective_K))
hess_mini_SL <- as.list(rep(NA, effective_K))
tensor_mini_SL <- as.list(rep(NA, effective_K))

for(k in 2:effective_K){
  mini_SL_list[[k]] <- optim(par = lead_est, fn = mini_SL_fun, 
                             control = list(fnscale = -1), method = "L-BFGS-B",
                             covariate = pre_processing$covariate_list[[k]],
                             failure_position = pre_processing$failure_position[[k]],
                             failure_num = pre_processing$failure_num[[k]],
                             risk_sets = pre_processing$risk_sets[[k]],
                             hess_external = lead_hess,
                             initial_beta = lead_est, lower = lower_a, upper = upper_b)$par
  
  grad_mini_SL[[k]] <- grad_plk(beta = mini_SL_list[[k]],
                                covariate = pre_processing$covariate_list[[k]],
                                failure_position = pre_processing$failure_position[[k]],
                                failure_num = pre_processing$failure_num[[k]],
                                risk_sets = pre_processing$risk_sets[[k]])
  
  
  hess_mini_SL[[k]] <- hess_plk(beta = mini_SL_list[[k]],
                                p = cov_dim,
                                covariate = pre_processing$covariate_list[[k]],
                                failure_num = pre_processing$failure_num[[k]],
                                risk_sets = pre_processing$risk_sets[[k]])

  tensor_mini_SL[[k]] <- third_plk(beta = mini_SL_list[[k]],
                                   p = cov_dim,
                                   covariate = pre_processing$covariate_list[[k]],
                                   failure_num = pre_processing$failure_num[[k]],
                                   risk_sets = pre_processing$risk_sets[[k]])
}

our_L <- optim(par = lead_est, fn = aug_SL, 
               control = list(fnscale = -1), method = "L-BFGS-B",
               covariate = pre_processing$covariate_list[[1]],
               failure_position = pre_processing$failure_position[[1]],
               failure_num = pre_processing$failure_num[[1]],
               risk_sets = pre_processing$risk_sets[[1]],
               grad_list = grad_mini_SL,
               hess_list = hess_mini_SL,
               est_list = mini_SL_list,
               lead_site = 1, lower = lower_a, upper = upper_b)$par

temp_hess <- hess_plk(beta = our_L,
                      p = cov_dim,
                      covariate = pre_processing$covariate_list[[1]],
                      failure_num = pre_processing$failure_num[[1]],
                      risk_sets = pre_processing$risk_sets[[1]])

our_L_var <- -solve((Reduce("+", hess_mini_SL[-1]) + temp_hess))


our_L_OS <- one_step_polish(current_est = our_L, grad_list = grad_mini_SL, hess_list = hess_mini_SL,
                            tensor_list = tensor_mini_SL, est_list = mini_SL_list, coordination = "L",
                            covariate = pre_processing$covariate_list[[1]],
                            failure_position = pre_processing$failure_position[[1]],
                            failure_num = pre_processing$failure_num[[1]],
                            risk_sets = pre_processing$risk_sets[[1]])

temp_hess <- hess_plk(beta = our_L_OS,
                      p = cov_dim,
                      covariate = pre_processing$covariate_list[[1]],
                      failure_num = pre_processing$failure_num[[1]],
                      risk_sets = pre_processing$risk_sets[[1]])

our_L_OS_var <- -solve((Reduce("+", hess_mini_SL[-1]) + temp_hess))


### proposed method coordinated by central server
local_est <- as.list(rep(NA, effective_K))
hess_local <- as.list(rep(NA, effective_K))
tensor_local <- as.list(rep(NA, effective_K))

for(k in 1:effective_K){
  local_est[[k]] <- optim(par = rep(0, cov_dim), fn = log_plk, method = "L-BFGS-B", control = list(fnscale = -1),
                          covariate = pre_processing$covariate_list[[k]],
                          failure_position = pre_processing$failure_position[[k]],
                          failure_num = pre_processing$failure_num[[k]],
                          risk_sets = pre_processing$risk_sets[[k]], lower = lower_a, upper = upper_b)$par
  
  hess_local[[k]] <- hess_plk(beta = local_est[[k]], 
                              p = cov_dim,
                              covariate = pre_processing$covariate_list[[k]],
                              failure_num = pre_processing$failure_num[[k]],
                              risk_sets = pre_processing$risk_sets[[k]])

  tensor_local[[k]] <- third_plk(beta = local_est[[k]], 
                                 p = cov_dim,
                                 covariate = pre_processing$covariate_list[[k]],
                                 failure_num = pre_processing$failure_num[[k]],
                                 risk_sets = pre_processing$risk_sets[[k]])
}

term1 <- matrix(0, cov_dim, cov_dim)
term2 <- rep(0, cov_dim)
for(my_k in 1:effective_K){
  term1 <- term1 + hess_local[[my_k]]
  term2 <- term2 + hess_local[[my_k]] %*% local_est[[my_k]]
}
our_C <- c(solve(term1) %*% term2)

our_C_var <- -solve(Reduce("+", hess_local))

our_C_OS <- one_step_polish(current_est = our_C, hess_list = hess_local, tensor_list = tensor_local, est_list =  local_est,
                            coordination = "C")
our_C_OS_var <- our_C_var




### ODACH
local_est <- as.list(rep(NA, effective_K))
local_var <- as.list(rep(NA, effective_K))
for(k in 1:effective_K){
  fit_local <- coxph(Surv(t2e, censor_ind) ~ V1 + V2 + V3 + V4 + V5 + V6, data = data_list[[k]], method = "breslow")
  local_est[[k]] <- fit_local$coefficients
  local_var[[k]] <- fit_local$var
}


term1 <- rep(0, cov_dim)
term2 <- matrix(0, cov_dim, cov_dim)
for(k in 1:effective_K){
  term1 <- term1 + solve(local_var[[k]]) %*% local_est[[k]]
  term2 <- term2 + solve(local_var[[k]])
}
meta_est <- c(solve(term2) %*% term1) 



grad_meta <- as.list(rep(NA, effective_K))
hess_meta <- as.list(rep(NA, effective_K))
for(k in 2:effective_K){
  grad_meta[[k]] <- grad_plk(beta = meta_est, 
                             covariate = pre_processing$covariate_list[[k]],
                             failure_position = pre_processing$failure_position[[k]],
                             failure_num = pre_processing$failure_num[[k]],
                             risk_sets = pre_processing$risk_sets[[k]])
  
  
  hess_meta[[k]] <- hess_plk(beta = meta_est, 
                             p = cov_dim,
                             covariate = pre_processing$covariate_list[[k]],
                             failure_num = pre_processing$failure_num[[k]],
                             risk_sets = pre_processing$risk_sets[[k]])
  
}

odach_est <- optim(par = meta_est, fn = ODACH, control = list(fnscale = -1), method = "L-BFGS-B",
                   covariate = pre_processing$covariate_list[[1]],
                   failure_position = pre_processing$failure_position[[1]],
                   failure_num = pre_processing$failure_num[[1]],
                   risk_sets = pre_processing$risk_sets[[1]],
                   grad_external = Reduce(`+`, grad_meta[-1]),
                   hess_external = Reduce(`+`, hess_meta[-1]),
                   initial_beta = meta_est, lower = lower_a, upper = upper_b)$par

temp_hess <- hess_plk(beta = odach_est,
                      p = cov_dim,
                      covariate = pre_processing$covariate_list[[1]],
                      failure_num = pre_processing$failure_num[[1]],
                      risk_sets = pre_processing$risk_sets[[1]])

odach_var <- -solve((Reduce("+", hess_meta[-1]) + temp_hess))


pooled_lower <- pool_est - qnorm(0.975) * sqrt(diag(pool_var))
pooled_upper <- pool_est + qnorm(0.975) * sqrt(diag(pool_var))

odach_lower <- odach_est - qnorm(0.975) * sqrt(diag(odach_var))
odach_upper <- odach_est + qnorm(0.975) * sqrt(diag(odach_var))

our_L_lower <- our_L - qnorm(0.975) * sqrt(diag(our_L_var))
our_L_upper <- our_L + qnorm(0.975) * sqrt(diag(our_L_var))

our_C_lower <- our_C - qnorm(0.975) * sqrt(diag(our_C_var))
our_C_upper <- our_C + qnorm(0.975) * sqrt(diag(our_C_var))

our_L_OS_lower <- our_L_OS - qnorm(0.975) * sqrt(diag(our_L_OS_var))
our_L_OS_upper <- our_L_OS + qnorm(0.975) * sqrt(diag(our_L_OS_var))

our_C_OS_lower <- our_C_OS - qnorm(0.975) * sqrt(diag(our_C_OS_var))
our_C_OS_upper <- our_C_OS + qnorm(0.975) * sqrt(diag(our_C_OS_var))

lead_lower <- lead_est - qnorm(0.975) * sqrt(diag(lead_var))
lead_upper <- lead_est + qnorm(0.975) * sqrt(diag(lead_var))

my_est <- matrix(NA, 7, cov_dim * 3)

my_est[1, ] <- c(pool_est, pooled_lower, pooled_upper)
my_est[2, ] <- c(odach_est, odach_lower, odach_upper)
my_est[3, ] <- c(our_L_OS, our_L_OS_lower, our_L_OS_upper)
my_est[4, ] <- c(our_L, our_L_lower, our_L_upper)
my_est[5, ] <- c(our_C_OS, our_C_OS_lower, our_C_OS_upper)
my_est[6, ] <- c(our_C, our_C_lower, our_C_upper)
my_est[7, ] <- c(lead_est, lead_lower, lead_upper)


method_labels    <- c("Pool", "ODACH", "hat_theta^L", 
                      "check_theta^L",
                      "hat_theta^C",
                      "check_theta^C",
                      "Single-site")

dimension_labels <- c("Treatment", "Male", 
                      "Age between 18 and 30",
                      "Age between 30 and 40",
                      "Black/African American",
                      "Other race")


# Construct the data frame row by row
n_methods <- length(method_labels)
n_dims <- length(dimension_labels)

est_df <- data.frame()

for (i in 1:n_methods) {
  for (j in 1:n_dims) {
    est_df <- rbind(est_df, data.frame(
      method = method_labels[i],
      dimension = dimension_labels[j],
      estimate = my_est[i, j],
      lower = my_est[i, j + 6],
      upper = my_est[i, j + 12]
    ))
  }
}

# Ensure correct factor levels
est_df$method    <- factor(est_df$method, levels = rev(method_labels))
est_df$dimension <- factor(est_df$dimension, levels = dimension_labels)

method_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                   "#66a61e", "#e6ab02", "#a6761d", "#666666")

pooled_vals <- my_est[1, 1:6]

vline_df <- data.frame(
  dimension = factor(dimension_labels, levels = dimension_labels),
  xintercept = pooled_vals
)


# Flag if CI covers zero
est_df$ci_cover_zero <- est_df$lower <= 0 & est_df$upper >= 0
est_df$ci_color <- ifelse(est_df$ci_cover_zero, "black", "red")

# Plot
p <- ggplot(est_df, aes(x = estimate, y = method)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = ci_color), height = 0.25) +
  geom_point(aes(color = ci_color), size = 2) +
  geom_vline(data = vline_df, aes(xintercept = xintercept),
             linetype = "dashed", color = "blue", inherit.aes = FALSE) +
  facet_wrap(~ dimension, ncol = 2, scales = "free_x") +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  labs(x = "", y = "Method") +
  theme_minimal() +
  theme(
    strip.text         = element_text(size = 15, face = "bold"),
    axis.text.y        = element_text(size = 15),
    axis.text.x        = element_text(size = 15),
    axis.title.x       = element_text(size = 15),
    axis.title.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    plot.title         = element_text(hjust = 0.5, size = 12, face = "bold")
  )

print(p)


abs_relative_bias <- matrix(NA, 7, 6)
for(i in 1:6){
  abs_relative_bias[, i] <- abs((my_est[, i] - my_est[1, i]) / my_est[1, i])
}

a <- cbind(abs_relative_bias[-1, ], apply(abs_relative_bias, MARGIN = 1, mean)[-1]) 
xtab <- xtable(a[c(2, 3, 1, 4, 5, 6), ] * 100, digits = 2)
print(xtab, include.rownames = FALSE)

