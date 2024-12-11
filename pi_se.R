## functions to compute the standard error and confidence intervals around \pi(u|t) using either the direct conditional posterior or a Monte Carlo method


conditional_survival_prob_se5 = function(t, u, model, new.re, x, W, S_t_last_observed, S_t_max, gauss_rules, mod_mat_long_f){
  
  #gaussian quadrature stuff for t and u
  quad.lambda.t = (t - 0)/2
  quad.mu.t = (t + 0)/2
  quad.y.t = t(as.matrix(quad.lambda.t) %*% rules[[gauss_rules]]$x + quad.mu.t)
  quad.psi.event.t = t(sapply(quad.y.t, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event.t = mod_mat_long_f(c(quad.y.t), tmax = t, W_long =  c(sapply(W, rep, gauss_rules)))
  
  quad.lambda.u = (u - 0)/2
  quad.mu.u = (u + 0)/2
  quad.y.u = t(as.matrix(quad.lambda.u) %*% rules[[gauss_rules]]$x + quad.mu.u)
  quad.psi.event.u = t(sapply(quad.y.u, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event.u = mod_mat_long_f(c(quad.y.u), tmax = u, W_long =  c(sapply(W, rep, gauss_rules)))
  
  quad.w = rules[[gauss_rules]]$w
  
  a_re_cols = new.re$a_re_cols
  a_re_long = unlist(c(new.re$a_re_cols))
  w = length(new.re$a_re_cols)
  
  a_re_cols_pad = matrix(rep(0,length(model$pars$alpha)), ncol = length(model$pars$alpha))
  a_re_cols_pad_quad = matrix(rep(0,length(model$pars$alpha)*gauss_rules), ncol = length(model$pars$alpha))
  a_re_cols_pad[,1:w] = a_re_cols
  a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules))
  
  #cumulative baseline hazard function
  h0_t_quad.t = quad.psi.event.t %*% model$pars$theta
  exp_zTg_quad.t = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum))
  h0_t_star_quad.t = h0_t_quad.t * exp_zTg_quad.t
  H0_t_last_observed = quad.lambda.t * apply(quad.w *  h0_t_star_quad.t, 2, sum)
  
  h0_t_quad.u = quad.psi.event.u %*% model$pars$theta
  exp_zTg_quad.u = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum))
  h0_t_star_quad.u = h0_t_quad.u * exp_zTg_quad.u
  H0_t_max = quad.lambda.u * apply(quad.w *  h0_t_star_quad.u, 2, sum)
  
  #A function - gamma score
  
  A_t_quad.t = h0_t_quad.t * exp_zTg_quad.t * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum)
  Az_t_last = quad.lambda.t * apply(c(quad.w) * A_t_quad.t, 2, sum)
  
  A_t_quad.u = h0_t_quad.u * exp_zTg_quad.u * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum)
  Az_t_max = quad.lambda.u * apply(c(quad.w) * A_t_quad.u, 2, sum)
  
  
  
  psi_t_star_quad.t = c(exp_zTg_quad.t) * quad.psi.event.t
  PsiStar_last = c(quad.lambda.t) * apply(c(quad.w) * psi_t_star_quad.t, 2, sum)
  
  psi_t_star_quad.u = c(exp_zTg_quad.u) * quad.psi.event.u
  PsiStar_max = c(quad.lambda.u) * apply(c(quad.w) * psi_t_star_quad.u, 2, sum)
  
  B_t_quad.t = c(h0_t_quad.t * exp_zTg_quad.t) * quad.phi.event.t
  Bphi_t_last = c(quad.lambda.t) * apply(c(quad.w) * B_t_quad.t, 2, sum)
  
  B_t_quad.u = c(h0_t_quad.u * exp_zTg_quad.u) * quad.phi.event.u
  Bphi_t_max = c(quad.lambda.u) * apply(c(quad.w) * B_t_quad.u, 2, sum)
  
  
  
  eXtB = exp(as.matrix(x) %*%  model$pars$beta)
  H_t_last_observed = H0_t_last_observed * c(eXtB)
  H_t_max = H0_t_max * c(eXtB)
  
  d_logitSuSt_beta = x * S_t_last_observed * (H_t_last_observed - H_t_max)/(S_t_last_observed - S_t_max)
  
  d_logitSuSt_gamma = S_t_last_observed * (Az_t_last - Az_t_max)/(S_t_last_observed - S_t_max)
  
  d_logitSuSt_theta= c(S_t_last_observed/(S_t_last_observed - S_t_max)) * (PsiStar_last - PsiStar_max)
  
  d_logitSuSt_alpha= c(S_t_last_observed * model$pars$gamma/(S_t_last_observed - S_t_max)) * (Bphi_t_last - Bphi_t_max)
  
  d_logitSuSt_kappa = c(S_t_last_observed * model$pars$gamma/(S_t_last_observed - S_t_max)) * (Bphi_t_last[1:w] - Bphi_t_max[1:w])
  
  d_logitSuSt_dBGTA = as.vector(c(as.matrix(d_logitSuSt_beta), d_logitSuSt_gamma, d_logitSuSt_theta, d_logitSuSt_alpha))
  #d_logitSuSt_dAK = as.vector(c(d_logitSuSt_kappa))
  
  pqma = length(d_logitSuSt_dBGTA)
  
  conditional_pi_se = sqrt(t(d_logitSuSt_dBGTA) %*% model$M2_inv[1:pqma, 1:pqma] %*% (d_logitSuSt_dBGTA) + 
                             d_logitSuSt_kappa %*% new.re$M2_inv %*% (d_logitSuSt_kappa) )
  
  
  ll_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) - qnorm(0.975)*conditional_pi_se
  
  ll = exp(ll_logit)/(1 + exp(ll_logit))
  
  ul_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) + qnorm(0.975)*conditional_pi_se
  
  ul = exp(ul_logit)/(1 + exp(ul_logit))
  
  
  out = list(conditional_pi_se = conditional_pi_se, ll = ll, ul = ul)
  
  
  return(out)
  
}






conditional_survival_prob_se6 = function(t, u, model, new.re, x, W, S_t_last_observed, S_t_max, gauss_rules, mod_mat_long_f){
  
  #gaussian quadrature stuff for t and u
  quad.lambda.t = (t - 0)/2
  quad.mu.t = (t + 0)/2
  quad.y.t = t(as.matrix(quad.lambda.t) %*% rules[[gauss_rules]]$x + quad.mu.t)
  quad.psi.event.t = t(sapply(quad.y.t, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event.t = mod_mat_long_f(c(quad.y.t), tmax = t, W_long =  c(sapply(W, rep, gauss_rules)))
  
  quad.lambda.u = (u - 0)/2
  quad.mu.u = (u + 0)/2
  quad.y.u = t(as.matrix(quad.lambda.u) %*% rules[[gauss_rules]]$x + quad.mu.u)
  quad.psi.event.u = t(sapply(quad.y.u, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                              Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event.u = mod_mat_long_f(c(quad.y.u), tmax = u, W_long =  c(sapply(W, rep, gauss_rules)))
  
  quad.w = rules[[gauss_rules]]$w
  
  a_re_cols = new.re$a_re_cols
  a_re_long = unlist(c(new.re$a_re_cols))
  w = length(new.re$a_re_cols)
  
  a_re_cols_pad = matrix(rep(0,length(model$pars$alpha)), ncol = length(model$pars$alpha))
  a_re_cols_pad_quad = matrix(rep(0,length(model$pars$alpha)*gauss_rules), ncol = length(model$pars$alpha))
  a_re_cols_pad[,1:w] = a_re_cols
  a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules))
  
  #cumulative baseline hazard function
  h0_t_quad.t = quad.psi.event.t %*% model$pars$theta
  exp_zTg_quad.t = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum))
  h0_t_star_quad.t = h0_t_quad.t * exp_zTg_quad.t
  H0_t_last_observed = quad.lambda.t * apply(quad.w *  h0_t_star_quad.t, 2, sum)
  
  h0_t_quad.u = quad.psi.event.u %*% model$pars$theta
  exp_zTg_quad.u = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum))
  h0_t_star_quad.u = h0_t_quad.u * exp_zTg_quad.u
  H0_t_max = quad.lambda.u * apply(quad.w *  h0_t_star_quad.u, 2, sum)
  
  #A function - gamma score
  
  A_t_quad.t = h0_t_quad.t * exp_zTg_quad.t * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum)
  Az_t_last = quad.lambda.t * apply(c(quad.w) * A_t_quad.t, 2, sum)
  
  A_t_quad.u = h0_t_quad.u * exp_zTg_quad.u * apply((matrix(rep(model$pars$alpha, gauss_rules), ncol = length(model$pars$alpha), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum)
  Az_t_max = quad.lambda.u * apply(c(quad.w) * A_t_quad.u, 2, sum)
  
  
  
  psi_t_star_quad.t = c(exp_zTg_quad.t) * quad.psi.event.t
  PsiStar_last = c(quad.lambda.t) * apply(c(quad.w) * psi_t_star_quad.t, 2, sum)
  
  psi_t_star_quad.u = c(exp_zTg_quad.u) * quad.psi.event.u
  PsiStar_max = c(quad.lambda.u) * apply(c(quad.w) * psi_t_star_quad.u, 2, sum)
  
  B_t_quad.t = c(h0_t_quad.t * exp_zTg_quad.t) * quad.phi.event.t
  Bphi_t_last = c(quad.lambda.t) * apply(c(quad.w) * B_t_quad.t, 2, sum)
  
  B_t_quad.u = c(h0_t_quad.u * exp_zTg_quad.u) * quad.phi.event.u
  Bphi_t_max = c(quad.lambda.u) * apply(c(quad.w) * B_t_quad.u, 2, sum)
  
  
  
  eXtB = exp(as.matrix(x) %*%  model$pars$beta)
  H_t_last_observed = H0_t_last_observed * c(eXtB)
  H_t_max = H0_t_max * c(eXtB)
  
  d_logitSuSt_beta = x * S_t_last_observed * (H_t_last_observed - H_t_max)/(S_t_last_observed - S_t_max)
  
  d_logitSuSt_gamma = S_t_last_observed * (Az_t_last - Az_t_max)/(S_t_last_observed - S_t_max)
  
  d_logitSuSt_theta= c(S_t_last_observed/(S_t_last_observed - S_t_max)) * (PsiStar_last - PsiStar_max)
  
  d_logitSuSt_alpha= c(S_t_last_observed * model$pars$gamma/(S_t_last_observed - S_t_max)) * (Bphi_t_last - Bphi_t_max)
  
  d_logitSuSt_kappa = c(S_t_last_observed * model$pars$gamma/(S_t_last_observed - S_t_max)) * (Bphi_t_last[1:w] - Bphi_t_max[1:w])
  
  d_logitSuSt_dBGTA = as.vector(c(as.matrix(d_logitSuSt_beta), d_logitSuSt_gamma, d_logitSuSt_theta, d_logitSuSt_alpha))
  #d_logitSuSt_dAK = as.vector(c(d_logitSuSt_kappa))
  
  pqma = length(d_logitSuSt_dBGTA)
  
  conditional_pi_se = sqrt(d_logitSuSt_dBGTA %*% new.re$cov_MC[1:pqma, 1:pqma] %*% (d_logitSuSt_dBGTA) + 
                             d_logitSuSt_kappa %*% new.re$cov_MC[(pqma+1):(pqma+w), (pqma+1):(pqma+w)] %*% (d_logitSuSt_kappa) )
  
  
  ll_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) - qnorm(0.975)*conditional_pi_se
  
  ll = exp(ll_logit)/(1 + exp(ll_logit))
  
  ul_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) + qnorm(0.975)*conditional_pi_se
  
  ul = exp(ul_logit)/(1 + exp(ul_logit))
  
  
  out = list(conditional_pi_se = conditional_pi_se, ll = ll, ul = ul)
  
  
  return(out)
  
}


