library(simts)
n1 = 100000
n2 = 500000
n3 = 1000000

model1 = AR1(.995, sigma2 = 1e-6) + WN(.005) + RW (1e-7)
model = 2*AR1() + RW() + DR()
Wt =  gen_gts(n3, model1)
Xt =  gen_gts(n1, model1)
Yt =  gen_gts(n2, model1)
Zt =  gen_gts(n3, model1)


mimu = make_mimu(Wt ,Xt, Yt, Zt, freq = 100, unit = "s",
         sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))

start_time <- Sys.time()
test_mgmwm1 = mgmwm(mimu, model, CI = FALSE, stationarity_test = FALSE, B_stationarity_test = NULL,
                    alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300)
end_time <- Sys.time()

end_time - start_time

plot(test_mgmwm1, decomp = T)

test_mgmwm2 = mgmwm(test_mgmwm1, model = 2*AR1(), CI = FALSE, stationarity_test = FALSE, B_stationarity_test = NULL,
                    alpha_ci = NULL, alpha_near_test = NULL, seed = 2710, n_boot_ci_max = 300)

plot(test_mgmwm2, decomp = T)


# Resutlat du test
start_time <- Sys.time()
test_model_selection = model_selection(mimu, model, s_est = NULL,
                                       alpha_paired_test = NULL, seed = 2710)
end_time <- Sys.time()

end_time - start_time

plot(test_model_selection)



ci_test_model_select = ci_mgmwm(test_model_selection)

ci_test_model_select$ci_low - ci_test_optim1$ci_low


