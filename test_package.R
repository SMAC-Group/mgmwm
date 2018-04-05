

library(iterpc)
#### Check why it does not load this library


n1 = 1000000
n2 = 1000000
n3 = 1000000

model1 =  AR1(.995, sigma2 = 1e-6) + WN(.005) + RW (1e-7)
model2 = AR1() + WN() + RW ()
Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n3, model1)


mimu = make_wvar_mimu_obj(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))

start_time <- Sys.time()
test_mgmwm1 = mgmwm(mimu, model2, CI = FALSE, stationarity_test = FALSE, B_stationarity_test = NULL,
                    alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300)
end_time <- Sys.time()

end_time - start_time

test_mgmwm2 = mgmwm(test_mgmwm1, CI = FALSE, stationarity_test = FALSE, B_stationarity_test = NULL,
                    alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300)



plot(test.optim1, process.decomp = T)

# Resutlat du test
start_time <- Sys.time()
test_model_selection = model_selection(model,mimu, s_est = NULL, stationarity_test = FALSE,
                                       B_stationarity_test = NULL,
                                       alpha_near_test = NULL, paired_test = TRUE,
                                       alpha_paired_test = NULL)
end_time <- Sys.time()

end_time - start_time

plot(test_model_selection, process.decomp = TRUE)

ci_test_model_select = ci_mgmwm(test_model_selection)

ci_test_model_select$ci_low - ci_test_optim1$ci_low


