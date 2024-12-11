## loop used for simulation study 2b and 2c - predictive accuracy

within_sample_pred_acc = NULL
out_of_sample_pred_acc = NULL

for(s in 1:200){
  
  name3 = paste("model_obj/n1000_nobs5_right25/save_model_", s, ".rds", sep = "")
  try.MLE = readRDS(name3)
  
  ## part 1 predictive accuracy of a new sample
  
  sample = gen_interval_censoring_jointmodel(n=1000, nobs = 5)
  
  dat.baseline.new = sample$dat.baseline
  dat.long.new = sample$dat.long
  
  out_of_sample_pred_acc_i = matrix(0, nrow = 9, ncol = 12)
  
  out_of_sample_pred_acc_i[1:3,1] = 0.1
  out_of_sample_pred_acc_i[1:3,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_newsample(0.1, 0.4, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[1,3:7] = unlist(true_dynpred_insample(0.1, 0.4, 0.95, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[1,8:12] = unlist(sens_spec_auc_pe(0.95, 0.1, 0.4, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.1, 0.7, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[2,3:7] = unlist(true_dynpred_insample(0.1, 0.7, 0.7, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[2,8:12] = unlist(sens_spec_auc_pe(0.7, 0.1, 0.7, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.1, 1, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[3,3:7] = unlist(true_dynpred_insample(0.1, 1, 0.45, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[3,8:12] = unlist(sens_spec_auc_pe(0.45, 0.1, 1, dat.baseline.new, cs))
  
  out_of_sample_pred_acc_i[4:6,1] = 0.4
  out_of_sample_pred_acc_i[4:6,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_newsample(0.4, 0.7, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[4,3:7] = unlist(true_dynpred_insample(0.4, 0.7, 0.75, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[4,8:12] = unlist(sens_spec_auc_pe(0.75, 0.4, 0.7, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.4, 1, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[5,3:7] = unlist(true_dynpred_insample(0.4, 1, 0.5, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[5,8:12] = unlist(sens_spec_auc_pe(0.5, 0.4, 1, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.4, 1.3, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[6,3:7] = unlist(true_dynpred_insample(0.4, 1.3, 0.35, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[6,8:12] = unlist(sens_spec_auc_pe(0.35, 0.4, 1.3, dat.baseline.new, cs))
  
  out_of_sample_pred_acc_i[7:9,1] = 0.7
  out_of_sample_pred_acc_i[7:9,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_newsample(0.7, 1, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[7,3:7] = unlist(true_dynpred_insample(0.7, 1, 0.65, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[7,8:12] = unlist(sens_spec_auc_pe(0.65, 0.7, 1, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.7, 1.3, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[8,3:7] = unlist(true_dynpred_insample(0.7, 1.3, 0.45, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[8,8:12] = unlist(sens_spec_auc_pe(0.45, 0.7, 1.3, dat.baseline.new, cs))
  cs = cond_survs_temp_newsample(0.7, 1.6, try.MLE, dat.baseline.new, dat.long.new)
  out_of_sample_pred_acc_i[9,3:7] = unlist(true_dynpred_insample(0.7, 1.6, 0.3, dat.baseline.new, cs))
  out_of_sample_pred_acc_i[9,8:12] = unlist(sens_spec_auc_pe(0.3, 0.7, 1.6, dat.baseline.new, cs))
  
  out_of_sample_pred_acc = rbind(out_of_sample_pred_acc, out_of_sample_pred_acc_i)

  print(c(s, out_of_sample_pred_acc_i[1,5]-out_of_sample_pred_acc_i[1,10], 
          out_of_sample_pred_acc_i[7,5]-out_of_sample_pred_acc_i[7,10]))
  
  
  ### part 2 predictive accuracy of in-sample validation
  rand_i = sample(1:200, 1)
  name1 = paste("data/n1000_nobs5_right25/dat_baseline", s, ".csv", sep = "")
  dat.baseline = read.csv(name1)
  dat.baseline = dat.baseline[,-1]  
  name2 = paste("data/n1000_nobs5_right25/dat_long", s, ".csv", sep = "")
  dat.long = read.csv(name2)
  dat.long = dat.long[,-1]  
  
  within_sample_pred_acc_i = matrix(0, nrow = 9, ncol = 12)
  
  within_sample_pred_acc_i[1:3,1] = 0.1
  within_sample_pred_acc_i[1:3,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_withinsample(0.1, 0.4, try.MLE, dat.baseline)
  within_sample_pred_acc_i[1,3:7] = unlist(true_dynpred_insample(0.1, 0.4, 0.95, dat.baseline, cs))
  within_sample_pred_acc_i[1,8:12] = unlist(sens_spec_auc_pe(0.95, 0.1, 0.4, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.1, 0.7, try.MLE, dat.baseline)
  within_sample_pred_acc_i[2,3:7] = unlist(true_dynpred_insample(0.1, 0.7, 0.7, dat.baseline, cs))
  within_sample_pred_acc_i[2,8:12] = unlist(sens_spec_auc_pe(0.7, 0.1, 0.7, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.1, 1, try.MLE, dat.baseline)
  within_sample_pred_acc_i[3,3:7] = unlist(true_dynpred_insample(0.1, 1, 0.45, dat.baseline, cs))
  within_sample_pred_acc_i[3,8:12] = unlist(sens_spec_auc_pe(0.45, 0.1, 1, dat.baseline, cs))
  
  within_sample_pred_acc_i[4:6,1] = 0.4
  within_sample_pred_acc_i[4:6,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_withinsample(0.4, 0.7, try.MLE, dat.baseline)
  within_sample_pred_acc_i[4,3:7] = unlist(true_dynpred_insample(0.4, 0.7, 0.75, dat.baseline, cs))
  within_sample_pred_acc_i[4,8:12] = unlist(sens_spec_auc_pe(0.75, 0.4, 0.7, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.4, 1, try.MLE, dat.baseline)
  within_sample_pred_acc_i[5,3:7] = unlist(true_dynpred_insample(0.4, 1, 0.5, dat.baseline, cs))
  within_sample_pred_acc_i[5,8:12] = unlist(sens_spec_auc_pe(0.5, 0.4, 1, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.4, 1.3, try.MLE, dat.baseline)
  within_sample_pred_acc_i[6,3:7] = unlist(true_dynpred_insample(0.4, 1.3, 0.35, dat.baseline, cs))
  within_sample_pred_acc_i[6,8:12] = unlist(sens_spec_auc_pe(0.35, 0.4, 1.3, dat.baseline, cs))
  
  within_sample_pred_acc_i[7:9,1] = 0.7
  within_sample_pred_acc_i[7:9,2] = c(0.3, 0.6, 0.9)
  
  cs = cond_survs_temp_withinsample(0.7, 1, try.MLE, dat.baseline)
  within_sample_pred_acc_i[7,3:7] = unlist(true_dynpred_insample(0.7, 1, 0.65, dat.baseline, cs))
  within_sample_pred_acc_i[7,8:12] = unlist(sens_spec_auc_pe(0.65, 0.7, 1, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.7, 1.3, try.MLE, dat.baseline)
  within_sample_pred_acc_i[8,3:7] = unlist(true_dynpred_insample(0.7, 1.3, 0.45, dat.baseline, cs))
  within_sample_pred_acc_i[8,8:12] = unlist(sens_spec_auc_pe(0.45, 0.7, 1.3, dat.baseline, cs))
  cs = cond_survs_temp_withinsample(0.7, 1.6, try.MLE, dat.baseline)
  within_sample_pred_acc_i[9,3:7] = unlist(true_dynpred_insample(0.7, 1.6, 0.3, dat.baseline, cs))
  within_sample_pred_acc_i[9,8:12] = unlist(sens_spec_auc_pe(0.3, 0.7, 1.6, dat.baseline, cs))
  
  within_sample_pred_acc = rbind(within_sample_pred_acc, within_sample_pred_acc_i)

}







