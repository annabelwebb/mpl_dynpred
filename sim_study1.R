## function and loop used for simulaltion study 1


dyn_pred2_sims_RC = function(t, u, i, new_baseline_data, model, fj, new_long_data_i){
  
  save_results = matrix(0, nrow = 1, ncol = 19)
  
  eXtB_true = exp(as.matrix(new_baseline_data[i,7:8]) %*% as.matrix(c(-0.5, 0.5)))
  
  H0_tt_true = h0_Quad_function_true(t, a_true = c(0.5, 0.1, 0.2, -0.5, -0.8, 0.2), g_true = 0.25, 
                                     re_true = new_baseline_data[i,10:11], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = poly_temp_f, tmax = 1)
  S_t_true = exp(-eXtB_true * H0_tt_true)
  
  ##true S(t = t)
  H0_tu_true = h0_Quad_function_true(u, a_true = c(0.5, 0.1, 0.2, -0.5, -0.8, 0.2), g_true = 0.25, 
                                     re_true = new_baseline_data[i,10:11], W = new_baseline_data[i,9], 
                                     mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_u_true = exp(-eXtB_true * H0_tu_true)
  save_results[,1] = S_u_true/S_t_true
  
  
  #MPL method with "asymptotic" covariance
  new.re_t0.2 = estimate_RE_indep(model, x = new_baseline_data[i,7:8], 
                                  t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                                  cont = new_long_data_i$wt_i_samp[which(new_long_data_i$obs_time_samp < t)], 
                                  mod_mat_long_f = poly_temp_f, W = new_baseline_data$fixed_x3[i], 
                                  gauss_rules = 15, max.it = 1000)
  
  
  eXtB = exp(as.matrix(new_baseline_data[i,7:8]) %*% model$pars$beta)
  
  H0_tt_star = h0_Quad_function(t, model = model, new.re = new.re_t0.2$a_re_cols, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_t = exp(-eXtB * H0_tt_star)
  
  H0_tu_star = h0_Quad_function(u, model = model, new.re = new.re_t0.2$a_re_cols, 
                                W = new_baseline_data$fixed_x3[i], 
                                mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_u = exp(-eXtB * H0_tu_star)
  save_results[,2] = S_u/S_t
  
  mpl_ci = conditional_survival_prob_se5(t, u, model, new.re_t0.2, 
                                         new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                         S_t, S_u, gauss_rules = 15, poly_temp_f)
  
  save_results[,3] = mpl_ci$ll
  save_results[,4] = mpl_ci$ul
  

  #MPL method with monte carlo covariance (t distribution 4 df distribution)
  new.re_t0.2_rs_t4df = RE_bs_covar2(model, x = new_baseline_data[i,7:8], 
                               t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                               cont = new_long_data_i$wt_i_samp[which(new_long_data_i$obs_time_samp < t)], 
                               mod_mat_long_f = poly_temp_f, W = new_baseline_data$fixed_x3[i], 
                               gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_t4df = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                   W = new_baseline_data$fixed_x3[i], 
                                   mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_t_rs_t4df = exp(-eXtB * H0_tt_star_rs_t4df)

  H0_tu_star_rs_t4df = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_t4df$a_re_cols, 
                                   W = new_baseline_data$fixed_x3[i], 
                                   mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_u_rs_t4df = exp(-eXtB * H0_tu_star_rs_t4df)
  save_results[,5] = S_u_rs_t4df/S_t_rs_t4df
  
  mpl_ci_rs_t4df = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_t4df, 
                                            new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                            S_t_rs_t4df, S_u_rs_t4df, gauss_rules = 15, poly_temp_f)
  
  save_results[,6] = mpl_ci_rs_t4df$ll
  save_results[,7] = mpl_ci_rs_t4df$ul
  
  
  #MPL method with monte carlo covariance (normal distribution)
  new.re_t0.2_rs_norm = RE_bs_covar1(model, x = new_baseline_data[i,7:8], 
                                     t_obs = new_long_data_i$obs_time_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     cont = new_long_data_i$wt_i_samp[which(new_long_data_i$obs_time_samp < t)], 
                                     mod_mat_long_f = poly_temp_f, W = new_baseline_data$fixed_x3[i], 
                                     gauss_rules = 15, max.it = 1000)
  
  H0_tt_star_rs_norm = h0_Quad_function(t, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_t_rs_norm = exp(-eXtB * H0_tt_star_rs_norm)
  
  H0_tu_star_rs_norm = h0_Quad_function(u, model = model, new.re = new.re_t0.2_rs_norm$a_re_cols, 
                                        W = new_baseline_data$fixed_x3[i], 
                                        mod_mat_long_f = poly_temp_f, tmax = 1)
  
  S_u_rs_norm = exp(-eXtB * H0_tu_star_rs_norm)
  save_results[,14] = S_u_rs_norm/S_t_rs_norm
  
  mpl_ci_rs_norm = conditional_survival_prob_se6(t, u, model, new.re_t0.2_rs_norm, 
                                                 new_baseline_data[i,7:8], new_baseline_data$fixed_x3[i], 
                                                 S_t_rs_norm, S_u_rs_norm, gauss_rules = 15, poly_temp_f)
  
  save_results[,15] = mpl_ci_rs_norm$ll
  save_results[,16] = mpl_ci_rs_norm$ul
  
  
  
  
  #jm estimation
  jm.surv = survfitJM(fj, new_long_data_i[which(new_long_data_i$obs_time_samp<t),], idVar = "id_long", survTimes = u, last.time = t)
  
  save_results[,17] = unlist(jm.surv$summaries)[3]
  save_results[,18] = unlist(jm.surv$summaries)[4]
  save_results[,19] = unlist(jm.surv$summaries)[5]
  
  return(save_results)
  
  
}



save_prediction_performance = NULL
save_prediction_performance_kappa = NULL
save_prediction_performance_z_JM = NULL
save_fitted_model = matrix(0, nrow = 500, ncol = 42)

for(s in 1:500){
  
  #generate sample to fit the model to
  sample = gen_right_censoring_jointmodel(n=1000)
  dat.baseline = sample$dat.baseline
  dat.long = sample$dat.long
  
  #fit mpl
  baseline_model_mat = dat.baseline[,7:8]
  baseline_model_mat = as.matrix(baseline_model_mat)
  
  #fit mpl
  samp.lme = lmer(dat.long$wt_samp ~ dat.long$W_long * dat.long$obs_time_samp + 
                    I(dat.long$obs_time_samp^2) +
                    I(dat.long$obs_time_samp^3) + 
                    (1 + dat.long$obs_time_samp || dat.long$id_long))
  
  try.MLE = try(jm_fitting_right_profile_new(dat.baseline, dat.long, baseline_model_mat,
                                             max.iter = c(3, 2000), 
                                             n.knots.h0 = 3,  step_size=1, range = c(0.1, 0.85),
                                             mod_mat_long_f = poly_temp_f, re_ind = c(1,2),
                                             init_re = as.matrix(ranef(samp.lme)$`dat.long$id_long`, ncol = 2),
                                             init_sigma_eps = as.data.frame(VarCorr(samp.lme))$vcov[3], 
                                             init_sigma_re = as.data.frame(VarCorr(samp.lme))$vcov[1:2], gauss_rules = 15))

  #fit jm
  dat.baseline$status = as.numeric(!dat.baseline$right)
  #library(JM)
  samp.lme = lme(wt_i_samp ~ obs_time_samp*W_long + I(obs_time_samp^2) + I(obs_time_samp^3), 
                 random = ~ 1 + obs_time_samp | id_long, data = dat.long, 
                 control = lmeControl(opt = "optim"))
  samp.surv = coxph(Surv(mid_point, status) ~ fixed_x1 + fixed_x2, data = dat.baseline, x = TRUE)
  fit.joint = jointModel(samp.lme, samp.surv, timeVar = "obs_time_samp", method = "spline-PH-GH")
  
  
  
  
  #save regression parameters and standard error estimates
  save_fitted_model[s,1:2] = try.MLE$pars$beta
  save_fitted_model[s,3] = try.MLE$pars$gamma
  
  save_fitted_model[s,4:5] = sqrt(diag(try.MLE$cov_H_RE)[1:2])
  save_fitted_model[s,6] = sqrt(diag(try.MLE$cov_H_RE)[3])
  
  save_fitted_model[s,7:8] = fit.joint$coefficients$gammas
  save_fitted_model[s,9] = fit.joint$coefficients$alpha
  
  save_fitted_model[s,10:11] = sqrt(diag(vcov(fit.joint)[8:9,8:9]))
  save_fitted_model[s,12] = sqrt(vcov(fit.joint)[10,10])
  
  save_fitted_model[s,13] = try.MLE$variance$sigma2_re[1]
  save_fitted_model[s,14] = try.MLE$variance$sigma2_re[2]
  
  save_fitted_model[s,15] = fit.joint$coefficients$D[1]
  save_fitted_model[s,16] = fit.joint$coefficients$D[4]
  
  save_fitted_model[s,17] = sqrt(try.MLE$variance$sigma2_Et)
  save_fitted_model[s,18] = fit.joint$coefficients$sigma
  
  #alpha = try.MLE$par$alpha
  #covar_alpha_mpl = try.MLE$M2_inv_trunc[7:10, 7:10]
  save_fitted_model[s,19:24] = try.MLE$par$alpha
  save_fitted_model[s,25:30] = sqrt(diag(try.MLE$cov_H_RE[10:15, 10:15]))
  
  save_fitted_model[s,31:36] = fit.joint$coefficients$betas
  save_fitted_model[s,37:42] = sqrt(diag(vcov(fit.joint)[1:6, 1:6]))
  
  
  ##generate sample to make predictions about
  sample = gen_right_censoring_jointmodel(n=1)
  dat.baseline.new = sample$dat.baseline
  dat.long.new = sample$dat.long
  
  save_prediction_performance_i = matrix(0, ncol=23, nrow = 12)
  save_prediction_performance_i[,1] = s
  
  save_prediction_performance_i_kappa = matrix(0, ncol=13, nrow = 3)
  save_prediction_performance_i_kappa [,1] = s
  
  for(i in 1:1){
    
    ##save information about new individual
    save_prediction_performance_i[,2] = i
    save_prediction_performance_i_kappa[,2] = i
    
    
    dat.long.new_i = dat.long.new[which(dat.long.new$id_long == i),]
    dat.long.new_i$fixed_x1 = dat.baseline.new$fixed_x1[i]
    dat.long.new_i$fixed_x2 = dat.baseline.new$fixed_x2[i]
    
    ##
    save_prediction_performance_i[1:4,3] = 0.2
    save_prediction_performance_i[1:4,4] = c(0.2, 0.5, 0.8, 1)
    
    save_prediction_performance_i[1,5:23]=dyn_pred2_sims_RC(0.2, 0.4, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[2,5:23]=dyn_pred2_sims_RC(0.2, 0.7, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[3,5:23]=dyn_pred2_sims_RC(0.2, 1, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[4,5:23]=dyn_pred2_sims_RC(0.2, 1.2, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    
    ##
    save_prediction_performance_i[5:8,3] = 0.5
    save_prediction_performance_i[5:8,4] = c(0.2, 0.5, 0.8, 1)
    
    save_prediction_performance_i[5,5:23]=dyn_pred2_sims_RC(0.5, 0.7, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[6,5:23]=dyn_pred2_sims_RC(0.5, 1, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[7,5:23]=dyn_pred2_sims_RC(0.5, 1.3, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[8,5:23]=dyn_pred2_sims_RC(0.5, 1.5, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    
    ##
    save_prediction_performance_i[9:12,3] = 0.8
    save_prediction_performance_i[9:12,4] = c(0.2, 0.5, 0.8, 1)
    
    save_prediction_performance_i[9,5:23]=dyn_pred2_sims_RC(0.8, 1, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[10,5:23]=dyn_pred2_sims_RC(0.8, 1.3, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[11,5:23]=dyn_pred2_sims_RC(0.8, 1.6, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    save_prediction_performance_i[12,5:23]=dyn_pred2_sims_RC(0.8, 1.8, i, dat.baseline.new, try.MLE, fit.joint, dat.long.new_i)
    
    save_prediction_performance = rbind(save_prediction_performance, save_prediction_performance_i)
    
    
    ## estimate of kappa
    save_prediction_performance_i_kappa[,3] = dat.baseline.new$g0[i]
    save_prediction_performance_i_kappa[,4] = dat.baseline.new$g1[i]
    
    save_prediction_performance_i_kappa[1,5] = 0.2
    
    # MPL asymptotic covariance
    kappa_est = estimate_RE_indep(try.MLE, x = dat.baseline.new[i,7:8], 
                      t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                      cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                      mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                      gauss_rules = 15, max.it = 1000)
    
    save_prediction_performance_i_kappa[1,6:7] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 4 df
    
    kappa_est = RE_bs_covar2(try.MLE, x = dat.baseline.new[i,7:8], 
                       t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                       cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                       mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                       gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[1,8:9] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 2 df
    kappa_est = RE_bs_covar3(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[1,10:11] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; normal distribution
    kappa_est = RE_bs_covar1(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.2)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[1,12:13] = kappa_est$a_re_cols
    
    
    # JM 
    preds_02 = c(predict(fit.joint, dat.long.new_i[which(dat.long.new_i$obs_time_samp < 0.2),], type = c("Subject"), idVar = "id_long",
                      FtTimes = seq(0, 1.8, length.out = 500), M = 300))
    
    save_prediction_performance_z_JM = cbind(save_prediction_performance_z_JM, preds_02)
    
    
    
    save_prediction_performance_i_kappa[2,5] = 0.5 
    
    # MPL asymptotic covariance
    kappa_est = estimate_RE_indep(try.MLE, x = dat.baseline.new[i,7:8], 
                                  t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                                  cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                                  mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                                  gauss_rules = 15, max.it = 1000)
    
    save_prediction_performance_i_kappa[2,6:7] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 4 df
    
    kappa_est = RE_bs_covar2(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[2,8:9] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 2 df
    kappa_est = RE_bs_covar3(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[2,10:11] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; normal distribution
    kappa_est = RE_bs_covar1(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.5)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[2,12:13] = kappa_est$a_re_cols
    
    
    # JM 
    preds_05 = c(predict(fit.joint, dat.long.new_i[which(dat.long.new_i$obs_time_samp < 0.5),], type = c("Subject"), idVar = "id_long",
                         FtTimes = seq(0, 1.8, length.out = 500), M = 300))
    
    save_prediction_performance_z_JM = cbind(save_prediction_performance_z_JM, preds_05)
    
    
    save_prediction_performance_i_kappa[3,5] = 0.8
    
    # MPL asymptotic covariance
    kappa_est = estimate_RE_indep(try.MLE, x = dat.baseline.new[i,7:8], 
                                  t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                                  cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                                  mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                                  gauss_rules = 15, max.it = 1000)
    
    save_prediction_performance_i_kappa[3,6:7] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 4 df
    
    kappa_est = RE_bs_covar2(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[3,8:9] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; t dist 2 df
    kappa_est = RE_bs_covar3(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[3,10:11] = kappa_est$a_re_cols
    
    # MPL monte carlo covariance; normal distribution
    kappa_est = RE_bs_covar1(try.MLE, x = dat.baseline.new[i,7:8], 
                             t_obs = dat.long.new_i$obs_time_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             cont = dat.long.new_i$wt_i_samp[which(dat.long.new_i$obs_time_samp < 0.8)], 
                             mod_mat_long_f = poly_temp_f, W = dat.baseline.new$fixed_x3[i], 
                             gauss_rules = 15, max.it = 100, B=500)
    
    save_prediction_performance_i_kappa[3,12:13] = kappa_est$a_re_cols
    
    
    # JM 
    preds_08 = c(predict(fit.joint, dat.long.new_i[which(dat.long.new_i$obs_time_samp < 0.8),], type = c("Subject"), idVar = "id_long",
                         FtTimes = seq(0, 1.8, length.out = 500), M = 300))
    
    save_prediction_performance_z_JM = cbind(save_prediction_performance_z_JM, preds_08)
    
    save_prediction_performance_kappa = rbind(save_prediction_performance_kappa, save_prediction_performance_i_kappa)
    
    as.numeric(save_prediction_performance_i[1,7] < save_prediction_performance_i[1,5] & 
                 save_prediction_performance_i[1,5] < save_prediction_performance_i[1,8])
    
    print(c(s, save_prediction_performance_i[1,5:6]))
    
  }
  
}




