#install.packages("CompQuadForm")
library(CompQuadForm)

### PARAMETERS ###
# K:  number of phenotypes
# M:  number of sites/studies
# J:  number of populations 
# n:  total sample size

### INPUTS ###
# Y:     matrix (dimension n x K) of phenotypes 
# A:     matrix (dimension n x 1) of SNP/exposure of interest
# X:     matrix (dimension n x p+1) of confounders, including 1 column for intercept
# site:  vector (length n) of site ids, integers from 1 to M
# pop :  vector (length n) of population ids, integers from 1 to J
# model: choice of logistic or linear regression
# level: nominal FDR level

### OUTPUTS ###
# mean_effect:  selected set of non-null mean effects
# het_effect:   selected set of non-null heterogeneous effects
# p_mean:       (aggregated) p-values for mean effects
# p_het:        (aggregated) p-values for heterogeneous effects
# fitted_beta:  fitted coefficients (M matrices, each dimension J x K)
# se_beta:      standard error of coefficients (M matrices, each dimension J x K)

########################
### TESTING FUNCTION ###
########################

test_heterogeneity = function(Y, A, X, site, pop, family=c('binomial', 'gaussian'), level=0.1) {
  M = length(unique(site)) # number of sites/studies
  J = length(unique(pop)) # number of subpopulations
  K = ncol(Y) # number of phenotypes
  
  p_mu = rep(0, K) # p-values for mean effect
  p_T = rep(0, ncol=K) # p-values for heterogeneous effect
  
  beta_hat = list()
  sigma_hat = list()
  
  prev = list() # prevalence of outcome
  
  ### construct test statistics / p-values
  beta_hat = matrix(0, nrow=J, ncol=K)
  sigma_hat = matrix(0, nrow=J, ncol=K)
  prev = matrix(0, nrow=J, ncol=K)
  sig_T = list()
  
  # regressions for each phenotype (specify glm family)
  for (j in 1:J) {
    index = which(pop == j) # indexing for population/site
    
    prev[j,] = colMeans(Y[index,]) # prevalence for all outcomes per population
    
    # prev[j,] = mean(Y[index])
    
    for (k in 1:K) {
      fit = coef(summary(glm(Y[index,k] ~ 0 + A[index] + X[index,], family=family)))
      beta_hat[j,k] = fit[1,1]
      sigma_hat[j,k] = fit[1,2]
    }
  }
  
  ### test statistics / p-values
  for (k in 1:K) {
    # mean effect (inverse variance weighting)
    mu_hat = weighted.mean(beta_hat[,k], 1/sigma_hat[,k]^2) 
    
    p_mu[k] = 2*pnorm(-abs(mu_hat), 0, sqrt(1 / sum(1/sigma_hat[,k]^2)))
    
    # heterogeneity
    T_hat = sum((beta_hat[,k] - mean(beta_hat[,k]))^2)
    
    sig_T[[k]] = matrix(sum(sigma_hat[,k]^2/J^2), J, J) - 
      matrix(rep(sigma_hat[,k]^2/J, J), J) - 
      t(matrix(rep(sigma_hat[,k]^2/J, J), J)) + 
      diag(sigma_hat[,k]^2) # quadratic matrix
    lam_T = eigen(sig_T[[k]])$values # eigenvalues of sig_T
    
    q0 = liu(T_hat, lam_T)
    p_T[k] = q0
    
    #if (q0$Qq > 0) {
    #  p_T[k] = min(q0$Qq, 1)
    #} else { 
      # in case of algorithm errors
      # if p-value negative, add absolute error
      # if absolute error still not enough, just use absolute error (will be small)
     # p_T[k] = ifelse(q0$Qq + q0$abserr > 0, 
    #                  min(q0$Qq + q0$abserr, 1), min(q0$abserr, 1) )
    #}
  }
  
  ### FDR control 
  Sa1 = which(p.adjust(p_T, 'bonferroni') <= level/2) # intermediate non-null set
  
  pSm1 = p_mu; pSm1[Sa1] = 0
  Sm = which(p.adjust(pSm1, 'BH') <= level)
  
  pSa0 = p_T; pSa0[Sa1] = 0
  pSa = pSa0[Sm] # restrict to Sm
  Sa = Sm[which(p.adjust(pSa, 'BH') <= level)]
  
  ### results
  return(list(
    phecodes = colnames(Y), # phenotypes
    mean_effect = Sm, # index set of significant mean effect (index for columns of Y)
    het_effect = Sa, # index set of significant heterogeneous effect
    p_mean = p_mu, p_het = p_T, #  p-values
    fitted_beta = beta_hat, se_beta = sigma_hat, # regression estimates
    prevalences = prev, # prevalence of outcomes
    sigma_mats = sig_T # covariance matrices for heterogeneity
  ))
}


### INPUTS ###
# beta_hat: matrix (dimension J x K) of beta (for each population, each phenotype)
# sigma_hat: matrix (dimension J x K) of se.beta (for each population, each phenotype)

library(CompQuadForm)

test_heterogeneity_beta = function(beta_hat, sigma_hat, level=0.1) {

  J = nrow(beta_hat) # number of subpopulations
  K = ncol(beta_hat) # number of phenotypes
  
  p_mu = rep(0, K) # p-values for mean effect
  p_T = rep(0, ncol=K) # p-values for heterogeneous effect
  
  sig_T = list()
  ### test statistics / p-values
  for (k in 1:K) {
    # mean effect (inverse variance weighting)
    mu_hat = weighted.mean(beta_hat[,k], 1/sigma_hat[,k]^2) 
    
    p_mu[k] = 2*pnorm(-abs(mu_hat), 0, sqrt(1 / sum(1/sigma_hat[,k]^2)))
    
    # heterogeneity
    T_hat = sum((beta_hat[,k] - mean(beta_hat[,k]))^2)
    
    sig_T[[k]] = matrix(sum(sigma_hat[,k]^2/J^2), J, J) - 
      matrix(rep(sigma_hat[,k]^2/J, J), J) - 
      t(matrix(rep(sigma_hat[,k]^2/J, J), J)) + 
      diag(sigma_hat[,k]^2) # quadratic matrix
    lam_T = eigen(sig_T[[k]])$values # eigenvalues of sig_T
    
    q0 = liu(T_hat, lam_T)
    p_T[k] = q0
    
    ##### heterogeneit p-value for two sub-population:
    
    T_hat = beta_hat[1,k] - beta_hat[2,k]
    T_sd = sqrt(sigma_hat[1,k]^2 + sigma_hat[2,k]^2)
    p_T[k] = 2 * pnorm(-abs(T_hat), 0, T_sd)

  }
  
  ### FDR control 
  Sa1 = which(p.adjust(p_T, 'bonferroni') <= level/2) # intermediate non-null set
  
  pSm1 = p_mu; pSm1[Sa1] = 0
  Sm = which(p.adjust(pSm1, 'BH') <= level)
  
  pSa0 = p_T; pSa0[Sa1] = 0
  pSa = pSa0[Sm] # restrict to Sm
  Sa = Sm[which(p.adjust(pSa, 'BH') <= level)]
  
  ### results
  return(list(
    mean_effect = Sm, # index set of significant mean effect (index for columns of Y)
    het_effect = Sa, # index set of significant heterogeneous effect
    p_mean = p_mu, p_het = p_T #  p-values
  ))
}



get_sig_T = function(sigma_hat) {
  K = ncol(sigma_hat) # number of phenotypes
  J = nrow(sigma_hat) # number of populations
  
  sig_T = list()
  
  for (k in 1:K) {
    sig_T[[k]] = matrix(sum(sigma_hat[,k]^2/J^2), J, J) - 
      matrix(rep(sigma_hat[,k]^2/J, J), J) - 
      t(matrix(rep(sigma_hat[,k]^2/J, J), J)) + 
      diag(sigma_hat[,k]^2)
  }
  
  return(sig_T)
}


fisher.combine <- function(p_lst){
  stat <- sum(- 2 * log(p_lst))
  1 - pchisq(stat, 2 * length(p_lst))
}


logit <- function(a){
  log((a + 1e-50) / (1 - a + 1e-50))
}

expit <- function(a){
  exp(a) / (1 + exp(a))
}

adapt_test <- function(p_het, p_mean, tau = 0.001, level = 0.1){
  
  y <- ifelse(p_het < tau, 1, 0)
  x <- logit(p_mean)
  
  model.fit <- glm(y ~ x, family = binomial())  
  w.fit <- 1 - expit(- model.fit$linear.predictors) / (1 - tau)
  w.fit[w.fit<=0]=10^{-100}
  w.fit <- w.fit / (1 - w.fit)  
  w.fit <- w.fit / mean(w.fit)
  
  p_weight <- ifelse((p_het / w.fit) > 1, 1, p_het / w.fit)
  p.mean=p.adjust(p_mean, method = 'BH')
  p.heter=p.adjust(p_weight, method = 'BH')
  return(list(select = as.vector(which(p.adjust(p_weight, method = 'BH') < level)), weight.pvl = p_weight,
              p.mean=p.mean,p.heter=p.heter))
}


# ##########################
# ### EXAMPLE SIMULATION ###
# ##########################
# 
# # parameters
# set.seed(1000)
# K = 1000
# Smu = 50
# Salpha = 20
# n = 500 # sample size per site/population
# M = 3
# J = 2
# p = 3
# gamma = (1:(p+1)) / (p+1)
# 
# # coefficients
# mu = rbind(c(rep(0.25, Smu), rep(0, K - Smu)),
#            c(rep(0.25, Smu), rep(0, K - Smu)))
# alpha = rbind(c(rep(0.2, Salpha), rep(0, K - Salpha)),
#               c(rep(-0.2, Salpha), rep(0, K - Salpha)))
# beta = mu + alpha
# 
# # data / inputs
# site = rep(1:M, each=J*n) # site/study labels (=1, 2, 3)
# 
# pop = rep(rep(1:J, each=n), M) # population labels (=1, 2)
# 
# for (i in 1:10) {
#   print(i)
# A = rbinom(M*J*n, 2, 0.5) # SNP / exposure
# 
# X = cbind(rep(1, M*J*n), matrix(rnorm(M*J*n*p), ncol=p)) # confounders
# 
# Y = NULL
# for (m in 1:M) {
#   for (j in 1:J) {
#     Y0 = NULL
#     for (k in 1:K) {
#       index = intersect(which(pop == j), which(site == m))
#       Y0 = cbind(Y0, rbinom(n, 1, plogis(beta[j,k]*A[index] + X[index,] %*% gamma)))
#     }
#     Y = rbind(Y, Y0)
#   }
# }
# 
# # apply method
# test = test_heterogeneity(Y, A, X, site, pop, family='binomial', level=0.1)
# 
# sig_T0 = get_sig_T(test$se_beta)
# 
# names(test)
# 
# # evaluate simulation
# Smu_true = 1:Smu
# Smu_est = test$mean_effect
# 
# Salpha_true = 1:Salpha
# Salpha_est = test$het_effect
# 
# fdp_mu = length(setdiff(Smu_est, Smu_true)) / max(length(Smu_est), 1)
# tdp_mu = length(intersect(Smu_est, Smu_true)) / length(Smu_true)
# 
# fdp_alpha = length(setdiff(Salpha_est, Salpha_true)) / max(length(Salpha_est), 1)
# tdp_alpha = length(intersect(Salpha_est, Salpha_true)) / length(Salpha_true)
# 
# print(c(fdp_mu=fdp_mu, tdp_mu=tdp_mu, fdp_alpha=fdp_alpha, tdp_alpha=tdp_alpha))
# }
