## additional functions used in to compute predictive accuracy in interval censored simulations

#truth
h0_function_true = function(time, a_true, g_true, re_true, W, mod_mat_long_f, tmax){
  
  h0 = ((1/(sqrt(2*pi)*0.5*time+1e-4)) * exp(-0.5 * ((log(time+1e-4)+0.3)/0.5)^2))/(1 - pnorm((log(time+1e-4)+0.3)/0.5))

  #phi.tvc.mu = mod_mat_long_f(time, tmax, W)
  
  #a_re = as.matrix(a_true)
  #a_re[1:length(re_true)] = a_re[1:length(re_true)] + as.matrix(re_true)
  
  z_t_i = -0.6 + (1-(0.75)*exp(4*time)/(1 + exp(4*time)) ) + re_true[,1] + time*re_true[,2]
  
  out = h0 * as.numeric(exp(as.matrix(z_t_i)%*%g_true))
  out
}


h0_Quad_function_true = function(upper.time, a_true, g_true, re_true, W, mod_mat_long_f, tmax){
  
  legendre.quadrature(h0_function_true, rules[[15]], lower = 0, upper = upper.time,
                      weighted = TRUE, a_true = a_true, g_true = g_true, re_true = re_true,  W=W,tmax=tmax,
                      mod_mat_long_f=mod_mat_long_f)
  
  
}




h0_function = function(time, model, new.re, W, mod_mat_long_f, tmax){
  
  psi.events = mSpline(time, degree = 3, knots = model$knots$theta.knots$int, 
                       Boundary.knots = model$knots$theta.knots$bound, intercept = F)
  h0 = psi.events %*% model$pars$theta
  
  phi.tvc.mu = mod_mat_long_f(time, tmax, W)
  
  a_re = model$pars$alpha[[1]]
  a_re[1:length(new.re)] = a_re[1:length(new.re)] + new.re
  
  z_t_i = phi.tvc.mu %*% (a_re)
  
  out = h0 * as.numeric(exp(z_t_i%*%model$pars$gamma))
  out
}


h0_Quad_function = function(upper.time, model, new.re, W, mod_mat_long_f, tmax){
  
  legendre.quadrature(h0_function, rules[[15]], lower = 0, upper = upper.time,
                      weighted = TRUE, model=model, new.re = new.re,  W=W,tmax=tmax,
                      mod_mat_long_f=mod_mat_long_f)
  
  
}



lin_f = function(t, tmax, W_long){
  
  phi = cbind(1, t, W_long)
  return(phi)
  
}



dyn_pred3_sims_IC_old = function(t, u, i, new_baseline_data, model, new_long_data_i){
  
  save_results = matrix(0, nrow = 1, ncol = 16)
  
  eXtB_true = exp(as.matrix(new_baseline_data[i,7]) %*% as.matrix(c(-1)))
  
  H0_tt_true = h0_Quad_function_true(t, a_true = c(-0.1, -0.1, -0.3), g_true = -0.3, 
                                     re_true = new_baseline_data[i,10:11], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = lin_f, tmax = 1)
  S_t_true = exp(-eXtB_true * H0_tt_true)
  
  ##true S(t = t)
  H0_tu_true = h0_Quad_function_true(u, a_true = c(-0.1, -0.1, -0.3), g_true = -0.3, 
                                     re_true = new_baseline_data[i,10:11], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = lin_f, tmax = 1)
  
  S_u_true = exp(-eXtB_true * H0_tu_true)
  save_results[,1] = S_u_true/S_t_true
  
  
  #MPL method with "asymptotic" covariance
  new.re_t0.2 = estimate_RE_indep(model, x = new_baseline_data[i,7], 
                                  t_obs = new_long_data_i$visit.time[which(new_long_data_i$visit.time < t)], 
                                  cont = new_long_data_i$W[which(new_long_data_i$visit.time < t)], 
                                  mod_mat_long_f = lin_f, W = new_baseline_data$fixed_x3[i], 
                                  gauss_rules = 15, max.it = 1000)
  
  
  eXtB = exp(as.matrix(new_baseline_data[i,7]) %*% model$pars$beta)
  
  H0_tt_star = h0_Quad_function(t, model = model, new.re = new.re_t0.2$a_re_cols, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = lin_f, tmax = 1)
  
  S_t = exp(-eXtB * H0_tt_star)
  
  H0_tu_star = h0_Quad_function(u, model = model, new.re = new.re_t0.2$a_re_cols, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = lin_f, tmax = 1)
  
  S_u = exp(-eXtB * H0_tu_star)
  save_results[,2] = S_u/S_t
  
  mpl_ci = conditional_survival_prob_se5(t, u, model, new.re_t0.2, 
                                         new_baseline_data[i,7], new_baseline_data$fixed_x3[i], 
                                         S_t, S_u, gauss_rules = 15, lin_f)
  
  save_results[,3] = mpl_ci$ll
  save_results[,4] = mpl_ci$ul
  
  
  #MPL method with monte carlo covariance (t distribution 4 df distribution)
  new.re_t0.2_rs_t4df = RE_bs_covar2(model, x = new_baseline_data[i,7], 
                                     t_obs = new_long_data_i$visit.time[which(new_long_data_i$visit.time < t)], 
                                     cont = new_long_data_i$W[which(new_long_data_i$visit.time < t)], 
                                     mod_mat_long_f = lin_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_t4df = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_t_rs_t4df = exp(-eXtB * H0_tt_star_rs_t4df)
  
  H0_tu_star_rs_t4df = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_u_rs_t4df = exp(-eXtB * H0_tu_star_rs_t4df)
  save_results[,5] = S_u_rs_t4df/S_t_rs_t4df
  
  mpl_ci_rs_t4df = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_t4df, 
                                                 new_baseline_data[i,7], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_t4df, S_u_rs_t4df, gauss_rules = 15, lin_f)
  
  save_results[,6] = mpl_ci_rs_t4df$ll
  save_results[,7] = mpl_ci_rs_t4df$ul
  
  
  #MPL method with monte carlo covariance (t distribution 2 df distribution)
  new.re_t0.2_rs_t2df = RE_bs_covar3(model, x = new_baseline_data[i,7], 
                                     t_obs = new_long_data_i$visit.time[which(new_long_data_i$visit.time < t)], 
                                     cont = new_long_data_i$W[which(new_long_data_i$visit.time < t)], 
                                     mod_mat_long_f = lin_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_t2df = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_t2df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_t_rs_t2df = exp(-eXtB * H0_tt_star_rs_t2df)
  
  H0_tu_star_rs_t2df = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_t2df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_u_rs_t2df = exp(-eXtB * H0_tu_star_rs_t2df)
  save_results[,8] = S_u_rs_t2df/S_t_rs_t2df
  
  mpl_ci_rs_t2df = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_t2df, 
                                                 new_baseline_data[i,7], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_t2df, S_u_rs_t2df, gauss_rules = 15, lin_f)
  
  save_results[,9] = mpl_ci_rs_t2df$ll
  save_results[,10] = mpl_ci_rs_t2df$ul
  
  
  
  
  #MPL method with monte carlo covariance (t distribution 2 df distribution)
  new.re_t0.2_rs_t2df = RE_bs_covar3(model, x = new_baseline_data[i,7], 
                                     t_obs = new_long_data_i$visit.time[which(new_long_data_i$visit.time < t)], 
                                     cont = new_long_data_i$W[which(new_long_data_i$visit.time < t)], 
                                     mod_mat_long_f = lin_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_t2df = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_t2df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_t_rs_t2df = exp(-eXtB * H0_tt_star_rs_t2df)
  
  H0_tu_star_rs_t2df = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_t2df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_u_rs_t2df = exp(-eXtB * H0_tu_star_rs_t2df)
  save_results[,11] = S_u_rs_t2df/S_t_rs_t2df
  
  mpl_ci_rs_t2df = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_t2df, 
                                                 new_baseline_data[i,7], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_t2df, S_u_rs_t2df, gauss_rules = 15, lin_f)
  
  save_results[,12] = mpl_ci_rs_t2df$ll
  save_results[,13] = mpl_ci_rs_t2df$ul
  
  
  #MPL method with monte carlo covariance (normal distribution)
  new.re_t0.2_rs_norm = RE_bs_covar1(model, x = new_baseline_data[i,7], 
                                     t_obs = new_long_data_i$visit.time[which(new_long_data_i$visit.time < t)], 
                                     cont = new_long_data_i$W[which(new_long_data_i$visit.time < t)], 
                                     mod_mat_long_f = lin_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_norm = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_t_rs_norm = exp(-eXtB * H0_tt_star_rs_norm)
  
  H0_tu_star_rs_norm = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = lin_f, tmax = 1)
  
  S_u_rs_norm = exp(-eXtB * H0_tu_star_rs_norm)
  save_results[,14] = S_u_rs_norm/S_t_rs_norm
  
  mpl_ci_rs_norm = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_norm, 
                                                 new_baseline_data[i,7], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_norm, S_u_rs_norm, gauss_rules = 15, lin_f)
  
  save_results[,15] = mpl_ci_rs_norm$ll
  save_results[,16] = mpl_ci_rs_norm$ul
  
  return(save_results)
  
  
}

ss_cov_f = function(t, tmax, W_long){
  
  #kn = quantile(df[,1], c(0.25, 0.5))
  ss1 = t-0.5
  ss1[which(ss1 < 0)] = 0
  
  
  phi = cbind(1, t, t^2, t^3, ss1^3)
  #phi = cbind(1, t, t^2, t^3)
  
  #phi.re = cbind(1, df[,1])
  
  return(phi)
  
  #return(list(phi = phi, phi.re = phi.re))
  
}


dyn_pred3_sims_IC = function(t, u, i, new_baseline_data, model, new_long_data_i){
  
  save_results = matrix(0, nrow = 1, ncol = 16)
  
  eXtB_true = exp(as.matrix(new_baseline_data[i,7:8]) %*% as.matrix(c(0.2, -0.4)))
  
  H0_tt_true = h0_Quad_function_true(t, a_true = c(-0.1, -0.1, -0.3), g_true = 1, 
                                     re_true = new_baseline_data[i,c(10,15)], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = ss_cov_f, tmax = 1)
  S_t_true = exp(-eXtB_true * H0_tt_true)
  
  ##true S(t = t)
  H0_tu_true = h0_Quad_function_true(u, a_true = c(-0.1, -0.1, -0.3), g_true = 1, 
                                     re_true = new_baseline_data[i,c(10,15)], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u_true = exp(-eXtB_true * H0_tu_true)
  save_results[,1] = S_u_true/S_t_true
  
  
  #MPL method with "asymptotic" covariance
  new.re_t0.2 = estimate_RE_indep(model, x = new_baseline_data[i,7:8], 
                                  t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                                  cont = new_long_data_i$wt_samp[which(new_long_data_i$obs_time_samp < t)], 
                                  mod_mat_long_f = ss_cov_f, W = new_baseline_data$fixed_x3[i], 
                                  gauss_rules = 15, max.it = 1000)
  
  
  eXtB = exp(as.matrix(new_baseline_data[i,7:8]) %*% model$pars$beta)
  
  H0_tt_star = h0_Quad_function(t, model = model, new.re = c(new.re_t0.2$a_re_cols), 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_t = exp(-eXtB * H0_tt_star)
  
  H0_tu_star = h0_Quad_function(u, model = model, new.re = c(new.re_t0.2$a_re_cols), 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u = exp(-eXtB * H0_tu_star)
  save_results[,2] = S_u/S_t
  
  mpl_ci = conditional_survival_prob_se5(t, u, model, new.re_t0.2, 
                                         new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                         S_t, S_u, gauss_rules = 15, ss_cov_f)
  
  save_results[,3] = mpl_ci$ll
  save_results[,4] = mpl_ci$ul
  
  
  #MPL method with monte carlo covariance (t distribution 4 df distribution)
  new.re_t0.2_rs_t4df = RE_bs_covar2(model, x = new_baseline_data[i,7:8], 
                                     t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     cont = new_long_data_i$wt_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     mod_mat_long_f = ss_cov_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_t4df = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_t_rs_t4df = exp(-eXtB * H0_tt_star_rs_t4df)
  
  H0_tu_star_rs_t4df = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u_rs_t4df = exp(-eXtB * H0_tu_star_rs_t4df)
  save_results[,5] = S_u_rs_t4df/S_t_rs_t4df
  
  mpl_ci_rs_t4df = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_t4df, 
                                                 new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_t4df, S_u_rs_t4df, gauss_rules = 15, ss_cov_f)
  
  save_results[,6] = mpl_ci_rs_t4df$ll
  save_results[,7] = mpl_ci_rs_t4df$ul
  
  #MPL method with monte carlo covariance (normal distribution)
  new.re_t0.2_rs_norm = RE_bs_covar1(model, x = new_baseline_data[i,7:8], 
                                     t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     cont = new_long_data_i$wt_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     mod_mat_long_f = ss_cov_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_norm = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_t_rs_norm = exp(-eXtB * H0_tt_star_rs_norm)
  
  H0_tu_star_rs_norm = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u_rs_norm = exp(-eXtB * H0_tu_star_rs_norm)
  save_results[,14] = S_u_rs_norm/S_t_rs_norm
  
  mpl_ci_rs_norm = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_norm, 
                                                 new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_norm, S_u_rs_norm, gauss_rules = 15, ss_cov_f)
  
  save_results[,15] = mpl_ci_rs_norm$ll
  save_results[,16] = mpl_ci_rs_norm$ul
  
  return(save_results)
  
  
}


true_pi = function(t, u, new_baseline_data, i){
  
  
  eXtB_true = exp(as.matrix(new_baseline_data[i,7:8]) %*% as.matrix(c(0.2, -0.4)))

  H0_tt_true = h0_Quad_function_true(t, a_true = c(-0.1, -0.1, -0.3), g_true = 1, 
                                     re_true = new_baseline_data[i,c(10,15)], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = ss_cov_f, tmax = 1)
  S_t_true = exp(-eXtB_true * H0_tt_true)
  
  ##true S(t = t)
  H0_tu_true = h0_Quad_function_true(u, a_true = c(-0.1, -0.1, -0.3), g_true = 1, 
                                     re_true = new_baseline_data[i,c(10,15)], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u_true = exp(-eXtB_true * H0_tu_true)
  return(S_u_true/S_t_true)
  
}

within_sample_pi_est = function(t, u, model, new_baseline_data, i){
  
  new.re_t0.2 = model$pars$a_re_cols[[1]][i,]
  
  eXtB = exp(as.matrix(new_baseline_data[i,7:8]) %*% model$pars$beta)
  
  H0_tt_star = h0_Quad_function(t, model = model, new.re = new.re_t0.2, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_t = exp(-eXtB * H0_tt_star)
  
  H0_tu_star = h0_Quad_function(u, model = model, new.re = new.re_t0.2, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = ss_cov_f, tmax = 1)
  
  S_u = exp(-eXtB * H0_tu_star)
  return(list(S_t = S_t, S_u = S_u, pi = S_u/S_t))
  
  
}

## function for SE of pi using original model 
conditional_survival_prob_se7 = function(t, u, model, x, i, W, S_t_last_observed, S_t_max, 
                                         gauss_rules, mod_mat_long_f){
  
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
  
  a_re_cols = model$pars$a_re_cols[[1]][i,]
  a_re_long = unlist(c(model$pars$a_re_cols[[1]][i,]))
  w = length(model$pars$a_re_cols[[1]][i,])
  n = nrow(model$pars$a_re_cols[[1]])
  
  a_re_cols_pad = matrix(rep(0,length(model$pars$alpha[[1]])), ncol = length(model$pars$alpha[[1]]))
  a_re_cols_pad_quad = matrix(rep(0,length(model$pars$alpha[[1]])*gauss_rules), ncol = length(model$pars$alpha[[1]]))
  a_re_cols_pad[,1:w] = a_re_cols
  a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules))
  
  #cumulative baseline hazard function
  h0_t_quad.t = quad.psi.event.t %*% model$pars$theta
  exp_zTg_quad.t = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum))
  h0_t_star_quad.t = h0_t_quad.t * exp_zTg_quad.t
  H0_t_last_observed = quad.lambda.t * apply(quad.w *  h0_t_star_quad.t, 2, sum)
  
  h0_t_quad.u = quad.psi.event.u %*% model$pars$theta
  exp_zTg_quad.u = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum))
  h0_t_star_quad.u = h0_t_quad.u * exp_zTg_quad.u
  H0_t_max = quad.lambda.u * apply(quad.w *  h0_t_star_quad.u, 2, sum)
  
  #A function - gamma score
  
  A_t_quad.t = h0_t_quad.t * exp_zTg_quad.t * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.t, 1, sum)
  Az_t_last = quad.lambda.t * apply(c(quad.w) * A_t_quad.t, 2, sum)
  
  A_t_quad.u = h0_t_quad.u * exp_zTg_quad.u * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + a_re_cols_pad_quad) * quad.phi.event.u, 1, sum)
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
  
  conditional_pi_se = sqrt(d_logitSuSt_dBGTA %*% model$cov_H_RE[1:pqma, 1:pqma] %*% (d_logitSuSt_dBGTA) + 
                             d_logitSuSt_kappa %*% model$cov_H_RE[c(pqma+i, pqma+n+i), c(pqma+i, pqma+n+i)] %*% (d_logitSuSt_kappa) )
  
  
  ll_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) - qnorm(0.975)*conditional_pi_se
  
  ll = exp(ll_logit)/(1 + exp(ll_logit))
  
  ul_logit = log((S_t_max/S_t_last_observed)/(1-S_t_max/S_t_last_observed)) + qnorm(0.975)*conditional_pi_se
  
  ul = exp(ul_logit)/(1 + exp(ul_logit))
  
  
  out = list(conditional_pi_se = conditional_pi_se, ll = ll, ul = ul)
  
  
  return(out)
  
}


## conditional covariance of an in sample pi

cond_covar_insample = function(model, i, x, t_obs, cont, mod_mat_long_f, W, gauss_rules){
  
  w = ncol(model$pars$a_re_cols[[1]])
  a_re_cols = model$pars$a_re_cols[[1]][i,]
  a_re_long = unlist(c(a_re_cols))
  a_re_cols_pad = matrix(rep(0, length(model$pars$alpha[[1]])), ncol = length(model$pars$alpha[[1]]))
  a_re_cols_pad_quad = matrix(rep(0,length(model$pars$alpha[[1]])*gauss_rules), ncol = length(model$pars$alpha[[1]]))
  a_re_cols_pad[,1:w] = as.matrix(a_re_cols)
  a_re_cols_pad_quad[,1:w] = matrix(sapply(a_re_cols, rep, gauss_rules), ncol = 2)
  
  
  quad.lambda = (max(t_obs) - 0)/2
  quad.mu = (max(t_obs) + 0)/2
  quad.y = t(as.matrix(quad.lambda) %*% rules[[gauss_rules]]$x + quad.mu)
  quad.psi.event = t(sapply(quad.y, mSpline, degree = 3, knots = model$knots$theta.knots$int, 
                            Boundary.knots = model$knots$theta.knots$bound))
  quad.phi.event = mod_mat_long_f(c(quad.y), tmax = max(t_obs), W_long =  c(sapply(W, rep, gauss_rules)))
  quad.w = rules[[gauss_rules]]$w
  
  h0_t_quad = quad.psi.event %*% model$pars$theta
  
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
  
  alpha_re = model$pars$alpha[[1]]
  alpha_re[1:w] = alpha_re[1:w] + a_re_cols
  z_t_ia = c(phi.mu %*% (alpha_re))
  
  exp_zTg_quad = exp(c(model$pars$gamma) * apply((matrix(rep(model$pars$alpha[[1]], gauss_rules), ncol = length(model$pars$alpha[[1]]), byrow = TRUE) + 
                                                    a_re_cols_pad_quad) * quad.phi.event, 1, sum))
  
  
  B_t_re_quad = c(h0_t_quad * exp_zTg_quad) * quad.phi.event[,1:w]
  df = data.frame(quad.w * B_t_re_quad)
  Bphi_t_re_r = c(quad.lambda) * apply(data.frame(quad.w* B_t_re_quad), 2, sum)
  
  
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
  phi.mu.1 = phi.mu[,1:w]
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
  
  alpha_re = model$pars$alpha[[1]]
  alpha_re[1:w] = alpha_re[1:w] + a_re_cols
  z_t_ia = c(phi.mu %*% (alpha_re))
  
  out = list(a_re_cols = a_re_cols, z_t_ia = z_t_ia, M2_inv = M2_inv)
  
  return(out)
  
  
}


cond_survs_temp_withinsample = function(t, u, model, new_baseline_data){
  
  cond_surv = S_u_save = S_t_save = S_tL_save = S_tR_save = NULL
  n = nrow(new_baseline_data)
  t_L_fin = new_baseline_data$t_L
  t_L_fin[which(is.infinite(t_L_fin))] = 0
  
  
  for(i in 1:n){
    new.re_t0.2 = model$pars$a_re_cols[[1]][i,]
    
    eXtB = exp(as.matrix(new_baseline_data[i,7:8]) %*% model$pars$beta)
    
    H0_tt_star = h0_Quad_function(t, model = model, new.re = new.re_t0.2, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_t = exp(-eXtB * H0_tt_star)
    S_t_save = c(S_t_save, S_t)
    
    H0_tu_star = h0_Quad_function(u, model = model, new.re = new.re_t0.2, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_u = exp(-eXtB * H0_tu_star)
    S_u_save = c(S_u_save, S_u)
    
    H0_tL_star = h0_Quad_function(t_L_fin[i], model = model, new.re = new.re_t0.2, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_tL = exp(-eXtB * H0_tL_star)
    
    S_tL_save = c(S_tL_save, S_tL)
    
    if(!is.infinite(new_baseline_data$t_R[i])){
      H0_tR_star = h0_Quad_function(new_baseline_data$t_R[i], model = model, new.re = new.re_t0.2, 
                                    W = new_baseline_data$fixed_x3[i], 
                                    mod_mat_long_f = ss_cov_f, tmax = 1)
      
      S_tR = exp(-eXtB * H0_tR_star)
    }else{
      S_tR = NA
    }
    
   S_tR_save = c(S_tR_save, S_tR)
    
    cond_surv = c(cond_surv, S_u/S_t)
    
  }
  
  out = data.frame(id = 1:n, cond_surv, S_t_save, S_u_save, S_tL_save, S_tR_save)
  
  return(out)
}



cond_survs_temp_newsample = function(t, u, model, new_baseline_data, new_long_data){
  
  cond_surv = S_u_save = S_t_save = S_tL_save = S_tR_save = NULL
  n = nrow(new_baseline_data)
  t_L_fin = new_baseline_data$t_L
  t_L_fin[which(is.infinite(t_L_fin))] = 0
  
  
  for(i in 1:n){
    new_long_data_i = new_long_data[which(new_long_data$id_long == i),]
    
    new.re_t0.2 = estimate_RE_indep(model, x = new_baseline_data[i,7:8], 
                      t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                      cont = new_long_data_i$wt_samp[which(new_long_data_i$obs_time_samp < t)], 
                      mod_mat_long_f = ss_cov_f, W = new_baseline_data$fixed_x3[i], 
                      gauss_rules = 15, max.it = 1000)
    
    
    eXtB = exp(as.matrix(new_baseline_data[i,7:8]) %*% model$pars$beta)
    
    H0_tt_star = h0_Quad_function(t, model = model, new.re = new.re_t0.2$a_re_cols, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_t = exp(-eXtB * H0_tt_star)
    S_t_save = c(S_t_save, S_t)
    
    H0_tu_star = h0_Quad_function(u, model = model, new.re = new.re_t0.2$a_re_cols, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_u = exp(-eXtB * H0_tu_star)
    S_u_save = c(S_u_save, S_u)
    
    H0_tL_star = h0_Quad_function(t_L_fin[i], model = model, new.re = new.re_t0.2$a_re_cols, 
                                  W = new_baseline_data$fixed_x3[i], 
                                  mod_mat_long_f = ss_cov_f, tmax = 1)
    
    S_tL = exp(-eXtB * H0_tL_star)
    
    S_tL_save = c(S_tL_save, S_tL)
    
    if(!is.infinite(new_baseline_data$t_R[i])){
      H0_tR_star = h0_Quad_function(new_baseline_data$t_R[i], model = model, new.re = new.re_t0.2$a_re_cols, 
                                    W = new_baseline_data$fixed_x3[i], 
                                    mod_mat_long_f = ss_cov_f, tmax = 1)
      
      S_tR = exp(-eXtB * H0_tR_star)
    }else{
      S_tR = NA
    }
    
    S_tR_save = c(S_tR_save, S_tR)
    
    cond_surv = c(cond_surv, S_u/S_t)
    
  }
  
  out = data.frame(id = 1:n, cond_surv, S_t_save, S_u_save, S_tL_save, S_tR_save)
  
  return(out)
}


case_w_ij_auc = function(t, u, tLi, tRi, tLj, tRj, S_t_i, S_u_i, S_t_j, S_u_j, StL_i, StR_i, StL_j, StR_j){
  if(t < tLi & tLi <= tRi & tRi < u){
    if(tRj > u){
      case = 1
      w_ij = 1
    }else if(tLj < u & u < tRj){
      case = 2
      w_ij = (S_u_j - StR_j)/(StL_j - StR_j)
    }else{
      case = NA
      w_ij = NA
    }
  }else if(tLi < t & t < u & u < tRi){
    if(tRj > u){
      case = 3
      w_ij = (S_t_i - S_u_i)/(StL_i - StR_i)
    }else if(tLj < u & u < tRj){
      case = 4
      w_ij = ((S_t_i - S_u_i)/(StL_i - StR_i)) * (S_u_j - StR_j)/(StL_j - StR_j)
    }else{
      case = NA
      w_ij = NA
    }
  }else if(t < tLi & tLi < u & u < tRi){
    if(tRj > u){
      case = 5
      w_ij = (StL_i - S_u_i)/(StL_i - StR_i)
    }else if(tLj < u & u < tRj){
      case = 6
      w_ij = ((StL_i - S_u_i)/(StL_i - StR_i)) * ((S_u_i - StR_j)/(StL_j - StR_j))
    }else{
      case = NA
      w_ij = NA
    }
  }else if(tLi < t & t < tRi & tRi < u){
    if(tRj > u){
      case = 7
      w_ij = (S_t_i - StR_i)/(StL_i - StR_i)
    }else if(tLj < u & u < tRj){
      case = 8
      w_ij = ((S_t_i - StR_i)/(StL_i - StR_i)) * ((S_u_j - StR_j)/(StL_j - StR_j))
    }else{
      case = NA
      w_ij = NA
    }
  }else{
    case = NA
    w_ij = NA
  }
  return(cbind(case, w_ij))
}


sens_spec_auc_pe = function(p, t, u, dat.baseline, survs){
  
  w_it = p_it = rep(NA, nrow(dat.baseline))
  
  cond = rep(0, nrow(dat.baseline))
  
  cond[which(dat.baseline$event == 1 & t < dat.baseline$t_L & dat.baseline$t_L < u)] = 1
  cond[which(dat.baseline$event == 1 & t < dat.baseline$t_L & u < dat.baseline$t_L)] = 2
  
  cond[which(dat.baseline$right == 1 & dat.baseline$t_L < t)] = 3
  cond[which(dat.baseline$right == 1 & t < dat.baseline$t_L & dat.baseline$t_L < u)] = 4
  cond[which(dat.baseline$right == 1 & u < dat.baseline$t_L)] = 5
  
  cond[which((dat.baseline$left == 1 | dat.baseline$interval == 1) & t < dat.baseline$t_L & dat.baseline$t_R < u)] = 6
  cond[which((dat.baseline$left == 1 | dat.baseline$interval == 1) & t < dat.baseline$t_L & dat.baseline$t_L < u & u < dat.baseline$t_R)] = 7
  cond[which((dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_L < t & u < dat.baseline$t_R )] = 8
  cond[which((dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_L < t & t < dat.baseline$t_R & dat.baseline$t_R < u)] = 9
  
  w_it[which(cond == 1)] = 1
  w_it[which(cond == 2)] = 0
  
  w_it[which(cond == 3)] = ((survs$S_t_save-survs$S_u_save)/survs$S_tL_save)[which(cond == 3)]
  w_it[which(cond == 4)] = ((survs$S_tL_save-survs$S_u_save)/survs$S_tL_save)[which(cond == 4)]
  w_it[which(cond == 5)] = 0
  
  w_it[which(cond == 6)] = 1
  w_it[which(cond == 7)] = ((survs$S_tL_save-survs$S_u_save)/(survs$S_tL_save-survs$S_tR_save))[which(cond == 7)]
  w_it[which(cond == 8)] = ((survs$S_t_save-survs$S_u_save)/(survs$S_tL_save-survs$S_tR_save))[which(cond == 8)]
  w_it[which(cond == 9)] = ((survs$S_t_save-survs$S_tR_save)/(survs$S_tL_save-survs$S_tR_save))[which(cond == 9)]
  
  p_it[which(cond == 1)] = 1
  p_it[which(cond == 2)] = 1
  
  p_it[which(cond == 3)] = ((survs$S_t_save)/survs$S_tL_save)[which(cond == 3)]
  p_it[which(cond == 4)] = 1
  p_it[which(cond == 5)] = 1
  
  p_it[which(cond == 6)] = 1
  p_it[which(cond == 7)] = 1
  p_it[which(cond == 8)] = ((survs$S_t_save-survs$S_tR_save)/(survs$S_tL_save-survs$S_tR_save))[which(cond == 8)]
  p_it[which(cond == 9)] = ((survs$S_t_save-survs$S_tR_save)/(survs$S_tL_save-survs$S_tR_save))[which(cond == 9)]
  
  w_it[which((dat.baseline$event == 1 | dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_R < t)] = 0
  p_it[which((dat.baseline$event == 1 | dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_R < t)] = 0
  
  
  #sensitivity
  sens = sum(w_it*as.numeric(survs$cond_surv < p), na.rm = TRUE)/sum(w_it, na.rm = TRUE)
  
  #specificity
  spec = sum((p_it - w_it)*as.numeric(survs$cond_surv >= p), na.rm = TRUE)/sum((p_it - w_it), na.rm = TRUE)
  
  if(any((dat.baseline$event == 1 | dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_R < t)){
    id.dont.count = which((dat.baseline$event == 1 | dat.baseline$left == 1 | dat.baseline$interval == 1) & dat.baseline$t_R < t)
    
    #auc
    dat.baseline.inc = dat.baseline[-id.dont.count,]
    S_t_save_inc = survs$S_t_save[-id.dont.count]
    S_u_save_inc = survs$S_u_save[-id.dont.count]
    S_tL_save_inc = survs$S_tL_save[-id.dont.count]
    S_tR_save_inc = survs$S_tR_save[-id.dont.count]
    survs_include = survs[-id.dont.count,]
    
    
  }else{
    dat.baseline.inc = dat.baseline
    S_t_save_inc = survs$S_t_save
    S_u_save_inc = survs$S_u_save
    S_tL_save_inc = survs$S_tL_save
    S_tR_save_inc = survs$S_tR_save
    survs_include = survs
    
    
  }
  
  
  case_save = matrix(0, ncol = 4, nrow = sum(cumsum(1:(nrow(dat.baseline.inc)-1))))
  
  for(i in 1:(nrow(dat.baseline.inc) - 1)){
    for(j in ((i+1):nrow(dat.baseline.inc))){
      case_ij = case_w_ij_auc(t, u, dat.baseline.inc$t_L[i], dat.baseline.inc$t_R[i], 
                              dat.baseline.inc$t_L[j], dat.baseline.inc$t_R[j], 
                              S_t_save_inc[i], S_u_save_inc[i], S_t_save_inc[j], S_u_save_inc[j], 
                              S_tL_save_inc[i], S_tR_save_inc[i], S_tL_save_inc[j], S_tR_save_inc[j])
      case_save[(((i-1)*(nrow(dat.baseline.inc)-1))+(j-2)),1:4] = c(i,j,case_ij)
      #print(c((((i-1)*(nrow(dat.baseline.inc)-1))+(j-2))))
      #case_save = rbind(c(i, j, case_ij), case_save)
    }
    #print(c(i))
    
  }
  
  
  case_save_clean = case_save[-which(is.na(case_save[,4])),]
  case_save_clean = case_save_clean[-which(case_save_clean[,1]==0),]
  
  numer_save = denom_save = rep(0, nrow(case_save_clean))
  
  for(pair in 1:nrow(case_save_clean)){
    i = case_save_clean[pair,1]
    j = case_save_clean[pair,2]
    
    cond_surv_comp = as.numeric(survs_include$cond_surv[i] < survs_include$cond_surv[j])
    
    numer = cond_surv_comp * case_save_clean[pair,4]
    #denom = case_save_clean[pair,4]
    
    numer_save[pair] = numer
    #denom_save[pair] = denom
  }
  
  
  auc = sum(numer_save)/sum(case_save_clean[,4])
  
  ## prediction error - absolute prediction error??
  pe_abs = sum(p_it * ((p_it - w_it) * abs(1 - survs$cond_surv) + (1 - (p_it - w_it)) * abs(0 - survs$cond_surv)), na.rm = TRUE)/sum(p_it, na.rm = TRUE)
  
  ## prediction error - square prediction error??
  
  pe_sq = sum(p_it * ((p_it - w_it) * ((1 - survs$cond_surv)^2) + (1 - (p_it - w_it)) * ((0 - survs$cond_surv)^2)), na.rm = TRUE)/sum(p_it, na.rm = TRUE)
  
  
  out = list(sens = sens, spec = spec, auc = auc, pe_abs = pe_abs, pe_sq = pe_sq)
  return(out)
  
}







