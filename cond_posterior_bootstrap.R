## functions to generate random effects and covariance matrices from a Monte Carlo bootstrap


RE_bs_covar1 = function(model, x, t_obs, cont, mod_mat_long_f, W, gauss_rules, max.it, B = 200){
  
  
  re_bs_save = sim_par_save = score_save = log_lik_save =  omega1_save = NULL
  
  
  
  for(b in 1:B){
    
    ## generate simulated eta values
    
    eta_length = length(c(model$pars$beta, model$pars$gamma, model$pars$theta, model$pars$alpha))
    p = length(model$pars$beta)
    q = length(model$pars$gamma)
    m = length(model$pars$theta)
    r = length(model$pars$alpha)
    
    sim_par = MASS::mvrnorm(1, mu = c(model$pars$beta, model$pars$gamma, model$pars$theta, model$pars$alpha), 
                            Sigma = model$cov_H_RE[1:eta_length, 1:eta_length])
    
    sim_par_save = rbind(sim_par_save, sim_par)
    
    beta_bs = sim_par[1:p]
    gamma_bs = sim_par[(p+1):(p+q)]
    theta_bs = sim_par[(p+q+1):(p+q+m)]
    alpha_bs = sim_par[(p+q+m+1):(p+q+m+r)]
    
    # estimate RE up to time t
    
    
    w = ncol(model$pars$a_re_cols)
    
    quad.lambda = (max(t_obs) - 0)/2
    quad.mu = (max(t_obs) + 0)/2
    quad.y = t(as.matrix(quad.lambda) %*% rules[[gauss_rules]]$x + quad.mu)
    quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
    quad.phi.event = mod_mat_long_f(c(quad.y), tmax = max(t_obs), W_long =  c(sapply(W, rep, gauss_rules)))
    quad.w = rules[[gauss_rules]]$w
    
    
    h0_t_quad = quad.psi.event %*% theta_bs
    
    a_re_cols_t = matrix(0, nrow = 1, ncol = w)
    a_re_long_t = unlist(c(a_re_cols_t))
    
    a_re_cols_pad_t = matrix(rep(0, r), ncol = r)
    a_re_cols_pad_quad_t = matrix(rep(0, r*gauss_rules), ncol = r)
    
    phi.mu = mod_mat_long_f(t_obs, tmax = max(t_obs), W_long = rep(W, length(t_obs)))
    
    eXtB_r = eXtB = exp(as.matrix(x) %*% beta_bs)
    
    #step size for random effects
    #likelihood
    exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                              a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
    
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
    
    #right
    pl_r = - sum(eXtB_r * H0_t_r)
    
    #least squares
    alpha_re = alpha_bs
    alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
    z_t_ia = c(phi.mu %*% (alpha_re))
    
    pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
    
    #log likelihood
    log_lik = pl_r - pl_ls - 
      sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
    
    
    for(it in 1:max.it){
      
      alpha_re = alpha_bs
      alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
      z_t_ia = c(phi.mu %*% (alpha_re))
      
      exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
      
      
      B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
      Bphi_t_re_r = c(quad.lambda) * apply(data.frame(quad.w* B_t_re_quad), 2, sum)
      
      a1i_score_r = - as.numeric(gamma_bs) * as.numeric(eXtB_r) * (Bphi_t_re_r)
      
      a1i_score_s = a1i_score_r
      
      phi.mu.1 = phi.mu[,1:w]
      if(length(cont) > 1){
        a1i_score_ls = apply((c(cont-z_t_ia)*phi.mu.1), 2, sum)
      }else{
        a1i_score_ls = (c(cont-z_t_ia)*phi.mu.1)
      }
      
      a1i_score = a1i_score_s + (1/model$variance$sigma2_Et)*a1i_score_ls - a_re_cols_t/matrix(model$variance$sigma2_re, ncol = w, byrow = T)
      
      a_hess_i = t(rep(- as.numeric(gamma_bs^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
        t((1/model$variance$sigma2_Et) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
        diag(1/model$variance$sigma2_re)
      
      a_re_cols_old_t = a_re_cols_t
      a_re_cols_t = a_re_cols_old_t + (as.matrix(a1i_score)) %*% solve(-a_hess_i)
      
      a_re_long_old_t = a_re_long_t
      a_re_cols_pad_old_t = a_re_cols_pad_t
      a_re_cols_pad_quad_old_t = a_re_cols_pad_quad_t
      
      a_re_long_t = unlist(c(a_re_cols_t))
      a_re_cols_pad_t[,1:w] = as.matrix(a_re_cols_t)
      a_re_cols_pad_quad_t[,1:w] = matrix(sapply(a_re_cols_t, rep, gauss_rules), ncol = 2)
      
      ## step size
      
      #step size for random effects
      #likelihood
      exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
      
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
      
      #right
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #least squares
      alpha_re = alpha_bs
      alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
      z_t_ia = c(phi.mu %*% (alpha_re))
      
      pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_r - pl_ls - 
        sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
      
      #print(c("re", log_lik_old, log_lik))
      kappa = 2
      if(((log_lik < log_lik_old))){
        i = 0
        omega1 = 1/kappa
        
        while(log_lik<log_lik_old){
          
          a_re_cols_t = a_re_cols_old_t + omega1 * (as.matrix(a1i_score)) %*% solve(-a_hess_i)
          
          a_re_long_old_t = a_re_long_t
          a_re_cols_pad_old_t = a_re_cols_pad_t
          a_re_cols_pad_quad_old_t = a_re_cols_pad_quad_t
          
          a_re_long_t = unlist(c(a_re_cols_t))
          a_re_cols_pad_t[,1:w] = as.matrix(a_re_cols_t)
          a_re_cols_pad_quad_t[,1:w] = matrix(sapply(a_re_cols_t, rep, gauss_rules), ncol = 2)
          
          #step size for random effects
          #likelihood
          exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                    a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
          
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
          
          #right
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #least squares
          alpha_re = alpha_bs
          alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
          z_t_ia = c(phi.mu %*% (alpha_re))
          
          pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
          
          #log likelihood
          #log_lik_old = log_lik
          log_lik = pl_r - pl_ls - 
            sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
          
          omega1_save = c(omega1_save, omega1)
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>50){break}
          
          
        }
        #print(c(omega1))
        
      }
      if(all(abs(c(a_re_cols_t - a_re_cols_old_t)) < 1e-8)){
        break
      }
      
      
      #print(c(it, a_re_cols_t))
      
      
      
      
    }
    
    ##second derivative
    a_hess_i = t(rep(- as.numeric(model$pars$gamma^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
      t((1/model$variance$sigma2_Et) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
      diag(1/model$variance$sigma2_re)
    neg.hess.inv = solve(-a_hess_i)
    
    a_re_est = MASS::mvrnorm(1, mu = a_re_cols_t, Sigma = neg.hess.inv)
    
    #score_save = c(score_save, a1i_score)
    #log_lik_save = c(log_lik_save, log_lik)
    #cond_surv_bs = c(cond_surv_bs, exp(-H_t_u)/exp(-H_t_t))
    
    
    re_bs_save = rbind(re_bs_save, a_re_est)
    
    
    
    
    
  }
  
  
  
  M2_inv = matrix(0, nrow = (eta_length + w), ncol = (eta_length+w))
  M2_inv[1:eta_length, 1:eta_length] = model$M2_inv[1:eta_length, 1:eta_length]
  M2_inv[(eta_length+1):(eta_length+w), (eta_length+1):(eta_length+w)] = cov(re_bs_save)
  
  cov_H = matrix(0, nrow = (eta_length + w), ncol = (eta_length+w))
  cov_H[1:eta_length, 1:eta_length] = model$cov_H_RE[1:eta_length, 1:eta_length]
  cov_H[(eta_length+1):(eta_length+w), (eta_length+1):(eta_length+w)] = cov(re_bs_save)
  
  for(s in 1:w){
    
    for(par in 1:eta_length){
      M2_inv[(eta_length + s), par] = cov(re_bs_save[,s], sim_par_save[,par])
      cov_H[(eta_length + s), par] = cov(re_bs_save[,s], sim_par_save[,par])
      
      M2_inv[par, (eta_length + s)] = cov(re_bs_save[,s], sim_par_save[,par])
      cov_H[par, (eta_length + s)] = cov(re_bs_save[,s], sim_par_save[,par])
      
    }
    
    
  }
  
  cov_MC = cov(cbind(sim_par_save, re_bs_save))
  
  out = list(a_re_cols = apply(re_bs_save, 2, mean), 
             re_bs_save = re_bs_save, sim_par_save = sim_par_save, 
             M2_inv = M2_inv, cov_H = cov_H, cov_MC = cov_MC, score_save=score_save, 
             log_lik_save = log_lik_save, omega1_save = omega1_save)
  
  return(out)
  
  
}




RE_bs_covar2 = function(model, x, t_obs, cont, mod_mat_long_f, W, gauss_rules, max.it, B = 200){
  
  
  re_bs_save = sim_par_save = score_save = log_lik_save =  omega1_save = NULL
  
  
  
  for(b in 1:B){
    
    ## generate simulated eta values
    
    eta_length = length(c(model$pars$beta, model$pars$gamma, model$pars$theta, model$pars$alpha))
    p = length(model$pars$beta)
    q = length(model$pars$gamma)
    m = length(model$pars$theta)
    r = length(model$pars$alpha)
    
    sim_par = MASS::mvrnorm(1, mu = c(model$pars$beta, model$pars$gamma, model$pars$theta, model$pars$alpha), 
                            Sigma = model$cov_H_RE[1:eta_length, 1:eta_length])
    
    sim_par_save = rbind(sim_par_save, sim_par)
    
    beta_bs = sim_par[1:p]
    gamma_bs = sim_par[(p+1):(p+q)]
    theta_bs = sim_par[(p+q+1):(p+q+m)]
    alpha_bs = sim_par[(p+q+m+1):(p+q+m+r)]
    
    # estimate RE up to time t
    
    
    w = ncol(model$pars$a_re_cols)
    
    quad.lambda = (max(t_obs) - 0)/2
    quad.mu = (max(t_obs) + 0)/2
    quad.y = t(as.matrix(quad.lambda) %*% rules[[gauss_rules]]$x + quad.mu)
    quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
    quad.phi.event = mod_mat_long_f(c(quad.y), tmax = max(t_obs), W_long =  c(sapply(W, rep, gauss_rules)))
    quad.w = rules[[gauss_rules]]$w
    
    
    h0_t_quad = quad.psi.event %*% theta_bs
    
    a_re_cols_t = matrix(0, nrow = 1, ncol = w)
    a_re_long_t = unlist(c(a_re_cols_t))
    
    a_re_cols_pad_t = matrix(rep(0, r), ncol = r)
    a_re_cols_pad_quad_t = matrix(rep(0, r*gauss_rules), ncol = r)
    
    phi.mu = mod_mat_long_f(t_obs, tmax = max(t_obs), W_long = rep(W, length(t_obs)))
    
    eXtB_r = eXtB = exp(as.matrix(x) %*% beta_bs)
    
    #step size for random effects
    #likelihood
    exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                              a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
    
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
    
    #right
    pl_r = - sum(eXtB_r * H0_t_r)
    
    #least squares
    alpha_re = alpha_bs
    alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
    z_t_ia = c(phi.mu %*% (alpha_re))
    
    pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
    
    #log likelihood
    log_lik = pl_r - pl_ls - 
      sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
    
    
    for(it in 1:max.it){
      
      alpha_re = alpha_bs
      alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
      z_t_ia = c(phi.mu %*% (alpha_re))
      
      exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
      
      
      B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
      Bphi_t_re_r = c(quad.lambda) * apply(data.frame(quad.w* B_t_re_quad), 2, sum)
      
      a1i_score_r = - as.numeric(gamma_bs) * as.numeric(eXtB_r) * (Bphi_t_re_r)
      
      a1i_score_s = a1i_score_r
      
      phi.mu.1 = phi.mu[,1:w]
      if(length(cont) > 1){
        a1i_score_ls = apply((c(cont-z_t_ia)*phi.mu.1), 2, sum)
      }else{
        a1i_score_ls = (c(cont-z_t_ia)*phi.mu.1)
      }
      
      a1i_score = a1i_score_s + (1/model$variance$sigma2_Et)*a1i_score_ls - a_re_cols_t/matrix(model$variance$sigma2_re[[1]], ncol = w, byrow = T)
      
      a_hess_i = t(rep(- as.numeric(gamma_bs^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
        t((1/model$variance$sigma2_Et) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
        diag(1/model$variance$sigma2_re)
      
      a_re_cols_old_t = a_re_cols_t
      a_re_cols_t = a_re_cols_old_t + (as.matrix(a1i_score)) %*% solve(-a_hess_i)
      
      a_re_long_old_t = a_re_long_t
      a_re_cols_pad_old_t = a_re_cols_pad_t
      a_re_cols_pad_quad_old_t = a_re_cols_pad_quad_t
      
      a_re_long_t = unlist(c(a_re_cols_t))
      a_re_cols_pad_t[,1:w] = as.matrix(a_re_cols_t)
      a_re_cols_pad_quad_t[,1:w] = matrix(sapply(a_re_cols_t, rep, gauss_rules), ncol = 2)
      
      ## step size
      
      #step size for random effects
      #likelihood
      exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
      
      h0_t_star_quad = h0_t_quad * exp_zTg_quad
      H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
      
      #right
      pl_r = - sum(eXtB_r * H0_t_r)
      
      #least squares
      alpha_re = alpha_bs
      alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
      z_t_ia = c(phi.mu %*% (alpha_re))
      
      pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
      
      #log likelihood
      log_lik_old = log_lik
      log_lik = pl_r - pl_ls - 
        sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
      
      #print(c("re", log_lik_old, log_lik))
      kappa = 2
      if(((log_lik < log_lik_old))){
        i = 0
        omega1 = 1/kappa
        
        while(log_lik<log_lik_old){
          
          a_re_cols_t = a_re_cols_old_t + omega1 * (as.matrix(a1i_score)) %*% solve(-a_hess_i)
          
          a_re_long_old_t = a_re_long_t
          a_re_cols_pad_old_t = a_re_cols_pad_t
          a_re_cols_pad_quad_old_t = a_re_cols_pad_quad_t
          
          a_re_long_t = unlist(c(a_re_cols_t))
          a_re_cols_pad_t[,1:w] = as.matrix(a_re_cols_t)
          a_re_cols_pad_quad_t[,1:w] = matrix(sapply(a_re_cols_t, rep, gauss_rules), ncol = 2)
          
          #step size for random effects
          #likelihood
          exp_zTg_quad = exp(c(gamma_bs) * apply((matrix(rep(alpha_bs, gauss_rules), ncol = r, byrow = TRUE) + 
                                                    a_re_cols_pad_quad_t) * quad.phi.event, 1, sum))
          
          h0_t_star_quad = h0_t_quad * exp_zTg_quad
          H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
          
          #right
          pl_r = - sum(eXtB_r * H0_t_r)
          
          #least squares
          alpha_re = alpha_bs
          alpha_re[1:w] = alpha_re[1:w] + a_re_cols_t
          z_t_ia = c(phi.mu %*% (alpha_re))
          
          pl_ls = (1/(2*model$variance$sigma2_Et))*sum((cont - z_t_ia)^2)
          
          #log likelihood
          #log_lik_old = log_lik
          log_lik = pl_r - pl_ls - 
            sum((1/(2*model$variance$sigma2_re)) * apply(a_re_cols_t * a_re_cols_t, 2, sum))
          
          omega1_save = c(omega1_save, omega1)
          
          #update value of omega1
          if(omega1>=1e-2){
            omega1 = omega1/kappa
          }else if(omega1<1e-2 & omega1>=1e-5){
            omega1 = omega1*(5e-2)
          }else if(omega1<1e-5){
            omega1 = omega1*(1e-5)
          }
          i = i+1
          if(i>50){break}
          
          
        }
        #print(c(omega1))
        
      }
      if(all(abs(c(a_re_cols_t - a_re_cols_old_t)) < 1e-8)){
        break
      }
      
      
      #print(c(it, a_re_cols_t))
      
      
      
      
    }
    
    ##second derivative
    a_hess_i = t(rep(- as.numeric(model$pars$gamma^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
      t((1/model$variance$sigma2_Et) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
      diag(1/model$variance$sigma2_re)
    neg.hess.inv = solve(-a_hess_i)
    
    
    a_re_est = mvtnorm::rmvt(n=1, delta = a_re_cols_t, sigma = neg.hess.inv, df = 4)
    
    #score_save = c(score_save, a1i_score)
    #log_lik_save = c(log_lik_save, log_lik)
    #cond_surv_bs = c(cond_surv_bs, exp(-H_t_u)/exp(-H_t_t))
    
    
    re_bs_save = rbind(re_bs_save, a_re_est)
    
    
    
    
    
  }
  
  
  
  M2_inv = matrix(0, nrow = (eta_length + w), ncol = (eta_length+w))
  M2_inv[1:eta_length, 1:eta_length] = model$M2_inv[1:eta_length, 1:eta_length]
  M2_inv[(eta_length+1):(eta_length+w), (eta_length+1):(eta_length+w)] = cov(re_bs_save)
  
  cov_H = matrix(0, nrow = (eta_length + w), ncol = (eta_length+w))
  cov_H[1:eta_length, 1:eta_length] = model$cov_H_RE[1:eta_length, 1:eta_length]
  cov_H[(eta_length+1):(eta_length+w), (eta_length+1):(eta_length+w)] = cov(re_bs_save)
  
  for(s in 1:w){
    
    for(par in 1:eta_length){
      M2_inv[(eta_length + s), par] = cov(re_bs_save[,s], sim_par_save[,par])
      cov_H[(eta_length + s), par] = cov(re_bs_save[,s], sim_par_save[,par])
      
      M2_inv[par, (eta_length + s)] = cov(re_bs_save[,s], sim_par_save[,par])
      cov_H[par, (eta_length + s)] = cov(re_bs_save[,s], sim_par_save[,par])
      
    }
    
    
  }
  
  cov_MC = cov(cbind(sim_par_save, re_bs_save))
  
  out = list(a_re_cols = apply(re_bs_save, 2, mean), 
             re_bs_save = re_bs_save, sim_par_save = sim_par_save, 
             M2_inv = M2_inv, cov_H = cov_H, cov_MC = cov_MC, score_save=score_save, 
             log_lik_save = log_lik_save, omega1_save = omega1_save)
  
  return(out)
  
  
}

