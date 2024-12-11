## function to estimate random effects by optimising the conditional posterior using a Newton algorithm


estimate_RE_indep = function(model, x, t_obs, cont, mod_mat_long_f, W, gauss_rules, max.it){
  
  w = ncol(model$pars$a_re_cols[[1]])
  
  quad.lambda = (max(t_obs) - 0)/2
  quad.mu = (max(t_obs) + 0)/2
  quad.y = t(as.matrix(quad.lambda) %*% rules[[gauss_rules]]$x + quad.mu)
  quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                            Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event = mod_mat_long_f(c(quad.y), tmax = max(t_obs), W_long =  c(sapply(W, rep, gauss_rules)))
  quad.w = rules[[gauss_rules]]$w
  
  h0_t_quad = quad.psi.event %*% model$pars$theta
  
  a_re_cols = matrix(0, nrow = 1, ncol = w)
  a_re_long = unlist(c(a_re_cols))
  
  a_re_cols_pad = matrix(rep(0, length(model$pars$alpha[[1]])), ncol = length(model$pars$alpha[[1]]))
  a_re_cols_pad_quad = matrix(rep(0,length(model$pars$alpha[[1]])*gauss_rules), ncol = length(model$pars$alpha[[1]]))
  
  phi.mu = mod_mat_long_f(t_obs, tmax = max(t_obs), W_long = rep(W, length(t_obs)))
  
  eXtB_r = eXtB = exp(as.matrix(x) %*% model$pars$beta)
  
  #step size for random effects
  #likelihood
  exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                    a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  
  #exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  h0_t_star_quad = h0_t_quad * exp_zTg_quad
  #h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
  H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
  
  #right
  pl_r = - sum(eXtB_r * H0_t_r)
  
  #least squares
  alpha_re = model$pars$alpha[[1]]
  alpha_re[1:w] = alpha_re[1:w] + a_re_cols
  z_t_ia = c(phi.mu %*% (alpha_re))
  
  pl_ls = (1/(2*model$variance$sigma2_Et[[1]]))*sum((cont - z_t_ia)^2)
  
  #log likelihood
  log_lik = pl_r - pl_ls - 
    sum((1/(2*model$variance$sigma2_re[[1]])) * apply(a_re_cols * a_re_cols, 2, sum))
  
  
  for(it in 1:max.it){
    
    alpha_re = model$pars$alpha[[1]]
    alpha_re[1:w] = alpha_re[1:w] + a_re_cols
    z_t_ia = c(phi.mu %*% (alpha_re))
    
    exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                      a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    
    
    B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
    df = data.frame(quad.w * B_t_re_quad)
    Bphi_t_re_r = c(quad.lambda) * apply(data.frame(quad.w* B_t_re_quad), 2, sum)
    
    #a1i_score_e = as.numeric(gamma) * phi_e[,1:w] - as.numeric(gamma) * as.numeric(eXtB_e) * (Bphi_t_re_e)
    a1i_score_r = - as.numeric(model$pars$gamma) * as.numeric(eXtB_r) * (Bphi_t_re_r)
    
    a1i_score_s = a1i_score_r
    
    phi.mu.1 = phi.mu[,1:w]
    if(length(cont) > 1){
      a1i_score_ls = apply((c(cont-z_t_ia)*phi.mu.1), 2, sum)
    }else{
      a1i_score_ls = (c(cont-z_t_ia)*phi.mu.1)
    }
    
    
    a1i_score = a1i_score_s + (1/model$variance$sigma2_Et[[1]])*a1i_score_ls - a_re_cols/matrix(model$variance$sigma2_re[[1]], ncol = w, byrow = T)
    
    #a_cols_update = matrix(0, nrow = n, ncol = w)
    
    
    a_hess_i = t(rep(- as.numeric(model$pars$gamma^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
      t((1/model$variance$sigma2_Et[[1]]) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
      diag(1/model$variance$sigma2_re[[1]])
    
    a_re_cols_old = a_re_cols
    a_re_cols = a_re_cols_old + (as.matrix(a1i_score)) %*% solve(-a_hess_i)
    
    a_re_long_old = a_re_long
    a_re_cols_pad_old = a_re_cols_pad
    a_re_cols_pad_quad_old = a_re_cols_pad_quad
    
    a_re_long = unlist(c(a_re_cols))
    a_re_cols_pad[,1:w] = as.matrix(a_re_cols)
    a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules), ncol = 2)
    
    ## step size
    
    #step size for random effects
    #likelihood
    exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                      a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    
    #exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
    h0_t_star_quad = h0_t_quad * exp_zTg_quad
    #h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
    H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
    
    #right
    pl_r = - sum(eXtB_r * H0_t_r)
    
    #least squares
    alpha_re = model$pars$alpha[[1]]
    alpha_re[1:w] = alpha_re[1:w] + a_re_cols
    z_t_ia = c(phi.mu %*% (alpha_re))
    
    pl_ls = (1/(2*model$variance$sigma2_Et[[1]]))*sum((cont - z_t_ia)^2)
    
    #log likelihood
    log_lik_old = log_lik
    log_lik = pl_r - pl_ls - 
      sum((1/(2*model$variance$sigma2_re[[1]])) * apply(a_re_cols * a_re_cols, 2, sum))
    
    #print(c("re", log_lik_old, log_lik))
    kappa = 2
    if(((log_lik < log_lik_old))){
      i = 0
      omega1 = 1/kappa
      
      
      while(log_lik<log_lik_old){
        
        a_re_cols = a_re_cols_old + omega1 * (as.matrix(a1i_score)) %*% solve(-a_hess_i)
        
        a_re_long_old = a_re_long
        a_re_cols_pad_old = a_re_cols_pad
        a_re_cols_pad_quad_old = a_re_cols_pad_quad
        
        a_re_long = unlist(c(a_re_cols))
        a_re_cols_pad[,1:w] = as.matrix(a_re_cols)
        a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules), ncol = 2)
        
        #step size for random effects
        #likelihood
        exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                          a_re_cols_pad_quad) * quad.phi.event, 1, sum))
        
        #exp_zTg_quad = exp(c(gamma) * apply((matrix(rep(alpha, gauss_rules*n), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum))
        h0_t_star_quad = h0_t_quad * exp_zTg_quad
        #h0_t_star_quad = matrix(h0_t_star_quad, ncol = n, nrow = gauss_rules, byrow = FALSE)
        H0_t_r = quad.lambda * sum(c(quad.w) * h0_t_star_quad)
        
        #right
        pl_r = - sum(eXtB_r * H0_t_r)
        
        #least squares
        alpha_re = model$pars$alpha[[1]]
        alpha_re[1:w] = alpha_re[1:w] + a_re_cols
        z_t_ia = c(phi.mu %*% (alpha_re))
        
        pl_ls = (1/(2*model$variance$sigma2_Et[[1]]))*sum((cont - z_t_ia)^2)
        
        #log likelihood
        #log_lik_old = log_lik
        log_lik = pl_r - pl_ls - 
          sum((1/(2*model$variance$sigma2_re[[1]])) * apply(a_re_cols * a_re_cols, 2, sum))
        
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
    
    
    #print(c(it, a_re_cols))
    
    if(all(abs(c(a_re_cols - a_re_cols_old)) < 1e-8)){
      break
    }
    
    
    
  }
  
  
  #a1i_hess_ls = sum(phi.mu.1)
  
  ###new covariance matrix????
  
  #H0_t_r = matrix(mapply(h0_Quad_function, upper.time = max(t_obs), 
  #                       a0 = a0, a1 = a1, 
  #                       MoreArgs = list(g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                       event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)), ncol = 1)
  
  
  #Az_t_r = matrix(mapply(h0_zb_Quad_function, upper.time = max(t_obs), 
  #                       a0 = a0, a1 = a1, 
  #                       MoreArgs = list(g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                       event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)), ncol = 1)
  
  
  #PsiStar_r = t(mapply(vf_PsiStar_u, upper.time = max(t_obs), 
  #                     a0 = a0, a1 = a1, 
  #                     MoreArgs = list(u = c(1:(length(model$pars$theta))), g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                     event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)))
  
  #Bphi_t_r = t(mapply(vf_h0_phi_l, upper.time = max(t_obs), 
  #                    a0 = a0, a1 = a1, 
  #                    MoreArgs = list(l = 1:(length(model$pars$alpha)), g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                    event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)))
  
  
  #Bstar_t_r = t(mapply(vf_B_star, upper.time = max(t_obs), 
  #                     a0 = a0, a1 = a1, 
  #                     MoreArgs = list(l = 1:(length(model$pars$alpha)), g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                     event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)))
  
  #D_t_r = t(mapply(vf_D_ilu, upper.time = max(t_obs), 
  #                 a0 = a0, a1 = a1, 
  #                 MoreArgs = list(l = 2, u = 1:(length(model$pars$theta)), g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                 event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)))
  
  #phiBphi_t_r = t(mapply(vf_h0_phi_ls, upper.time = max(t_obs), 
  #                       a0 = a0, a1 = a1, 
  #                       MoreArgs = list(l = 1:(length(model$pars$alpha)), s=2, g=model$pars$gamma,th=model$pars$theta, a=model$pars$alpha,
  #                                       event.kn = model$knots$theta.knots, tvc.kn.mu = model$knots$eta.knots)))
  
  
  
  p = length(model$pars$beta)
  q = length(model$pars$gamma)
  m = length(model$pars$theta)
  r = length(model$pars$alpha[[1]])
  
  Hess = matrix(0, nrow = (p+q+m+r + w), ncol = (p+q+m+r + w))
  
  Hess[1:(p+q+m+r), 1:(p+q+m+r)] = model$H_full[1:(p+q+m+r), 1:(p+q+m+r)]
  
  #beta kappa
  exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                    a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  B_t_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event
  B_t = c(quad.lambda) * apply(c(quad.w) * B_t_quad, 2, sum)
  
  betaa1_hess = do.call(rbind, lapply(seq_len(w), function(X) -x)) *
    c(model$pars$gamma) * c(eXtB) * B_t[1:w] 
  Hess[1:p, (p+q+m+r+1):(p+q+m+r+w)] = t(as.matrix(betaa1_hess))
  Hess[(p+q+m+r+1):(p+q+m+r+w), 1:p] = t(Hess[1:p, (p+q+m+r+1):(p+q+m+r+w)])
  
  #gamma kappa
  E_t_quad = c(model$pars$gamma) *  c(h0_t_quad * exp_zTg_quad *apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = r, byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event, 1, sum)) * quad.phi.event + 
    c(h0_t_quad * exp_zTg_quad ) * quad.phi.event
  E_t = c(quad.lambda) * apply(c(quad.w) * E_t_quad, 2, sum)
  
  gammaa1_hess = c( - (c(eXtB) * E_t[1:w])) 
  Hess[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w)] = as.matrix(gammaa1_hess)
  Hess[(p+q+m+r+1):(p+q+m+r+w), (p+1):(p+q)] = t(Hess[(p+1):(p+q), (p+q+m+r+1):(p+q+m+r+w)])
  
  #theta a1
  psi_t_star_quad = c(exp_zTg_quad) * quad.psi.event
  psi_t_quad2 = c(quad.w) * c(sapply(quad.lambda, rep, gauss_rules)) * psi_t_star_quad 
  
  temp = c(quad.phi.event[,1:w]) * 
    do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(model$pars$gamma), rep, gauss_rules)) * psi_t_quad2))
  df = data.frame(re_no = c(sapply(1:w, rep, gauss_rules)), temp) 
  thetaa1_hess = as.matrix(aggregate( df[,2:(m+1)], list(df[,1]), FUN = sum )[,-c(1)]) ## CHECK THIS IS CORRECT
  Hess[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w)] = t(thetaa1_hess)
  Hess[(p+q+m+r+1):(p+q+m+r+w),(p+q+1):(p+q+m)] = t(Hess[(p+q+1):(p+q+m), (p+q+m+r+1):(p+q+m+r+w)])
  
  #alpha a1
  B_t_quad2 = c(quad.w) * c(sapply(quad.lambda, rep, gauss_rules)) * B_t_quad
  
  temp = c(quad.phi.event[,1:w]) * 
    do.call(rbind, lapply(seq_len(w), function(X) - c(sapply(c(eXtB) * c(model$pars$gamma^2), rep, gauss_rules)) * B_t_quad2))
  df = data.frame(re_no = c(sapply(1:w, rep, gauss_rules)), temp)
  alphaa1_hess = as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-c(1)]) ## CHECK THIS IS CORRECT
  Hess[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w)] = t(alphaa1_hess)
  Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+1):(p+q+m+r)] = t(Hess[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w)])
  
  #kappa kappa
  a_hess_i = t(rep(- as.numeric(model$pars$gamma^2) * as.numeric(eXtB), gauss_rules) * B_t_re_quad) %*%  quad.phi.event[,1:w] - 
    t((1/model$variance$sigma2_Et[[1]]) * matrix(phi.mu.1, ncol = w)) %*% matrix(phi.mu.1, ncol = w) - 
    diag(1/model$variance$sigma2_re[[1]])
  Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+r+1):(p+q+m+r+w)] = a_hess_i
  
  
  #least squares part
  #alpha kappa
  temp = c(phi.mu[,1:w]) * 
    do.call(rbind, lapply(seq_len(w), function(X) phi.mu))
  df = data.frame(re_no = c(sapply(1:w, rep, length(t_obs))), temp)
  alphaa1_hess_ls = - c(1/model$variance$sigma2_Et[[1]]) *  as.matrix(aggregate( df[,2:(r+1)], list(df[,1]), FUN = sum )[,-c(1)])
  
  Hess[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w)] = Hess[(p+q+m+1):(p+q+m+r), (p+q+m+r+1):(p+q+m+r+w)] + t(alphaa1_hess_ls)
  Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+1):(p+q+m+r)] = Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+1):(p+q+m+r)] + alphaa1_hess_ls
  
  #kappa kappa already included above
  #temp = c(phi.mu[,1:w]) * 
  #  do.call(rbind, lapply(seq_len(w), function(X) phi.mu[,1:w]))
  #df = data.frame(re_no = c(sapply(1:w, rep, length(t_obs))), temp)
  #a1a1_hess_ls = - c(1/model$variance$sigma2_Et) *  as.matrix(aggregate( df[,2:(w+1)], list(df[,1]), FUN = sum )[,-c(1)])
  
  #RE term
  #already in kappa kappa term above i think?
  #Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+r+1):(p+q+m+r+w)]  = Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+r+1):(p+q+m+r+w)] + diag(-c(1/c(model$variance$sigma2_re)))
  
  
  M2 = -Hess[(p+q+m+r+1):(p+q+m+r+w), (p+q+m+r+1):(p+q+m+r+w)]
  
  M2_inv = matrix(0, nrow(M2), nrow(M2))
  #diag(corr) = pos
  #corr[!pos,]=0
  
  M2_inv = solve(M2)
  
  #A_omega = corr %*% M2_inv %*% t(corr)
  
  #get rid of penalty part for the inside of the sandwich formula?
  #Hess[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] = Hess[(p+q+1):(p+q+m), (p+q+1):(p+q+m)] + model$variance$theta.lambda * model$penalty$R
  
  
  #cov_H = A_omega %*% (-(Hess)) %*% A_omega
  
  alpha_re = model$pars$alpha[[1]]
  alpha_re[1:w] = alpha_re[1:w] + a_re_cols
  z_t_ia = c(phi.mu %*% (alpha_re))
  
  out = list(a_re_cols = a_re_cols, z_t_ia = z_t_ia, M2_inv = M2_inv)
  
  return(out)
  
  
}

