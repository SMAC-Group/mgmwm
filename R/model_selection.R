

#' @export
model_selection = function(mimu, model, s_est = NULL,
                           paired_test = FALSE,
                           alpha_paired_test = NULL, seed = 2710){

  # Level of confidence for Wilcoxon paired test
  if(is.null(alpha_paired_test)){
    alpha_paired_test = .05
  }else{
    alpha_paired_test = alpha_paired_test
  }

  # Number of replicates
  n_replicates = length(mimu)

  # Number of time series replicates for estimation
  if(is.null(s_est)){
    s_est = ceiling(n_replicates/2)
  }else{
    s_est = s_est
  }

  # Number of time series replicates for validation
  s_valid = n_replicates - s_est

  # Matrix of possible combination of Time series to estimate
  pair = t(combn(n_replicates,s_est))

  # Number of possible combination
  n_permutation = (dim(pair)[1])

  # Create model object with all possible nested model in model_max
  model_nested = model_combination(model_max = model)

  # List of model nested within model
  model_test = list()

  #Number of nested model
  n_models = length(model_nested)


  # Extract tau max
  length_tau = rep(NA,n_replicates)

  for (i in 1:n_replicates){
    length_tau[i] = length(mimu[[i]]$tau)
  }

  # Vector of maximum tau and scales
  tau_max_vec = mimu[[which.max(length_tau)]]$tau
  scales_max_vec = mimu[[which.max(length_tau)]]$scales

  cv_wvic = rep(NA,n_models)
  obj_out_sample = matrix(NA,n_permutation, n_models)
  model_complexity = rep(NA,n_models)


  set.seed(seed)

  pb <- progress_bar$new(
    format = "  Model :current of :total Models. Time remaining:  :eta",
    clear = FALSE, total = n_models, width = 100)

  for (i in 1:n_models){

    # Extract the type of model
    desc = model_nested[[i]]$desc

    # Number of latent pocesses in model
    np = model_nested[[i]]$plength

    # Name of parameters in model
    obj_desc = model_nested[[i]]$obj.desc

    # Value of parameters in model
    theta = model_nested[[i]]$theta

    N = rep(NA,n_replicates)
    param_starting = matrix(NA,n_replicates,np)

    for (k in 1:n_replicates){
      # Length of each error signal
      N[k] = length(mimu[[k]]$data)

      # Extract the data from mimu object
      data = mimu[[k]]$data

      # Set the number of boostrap replicates to compute covariance matrix
      if(N[k] > 10000){
        G = 1e6
      }else{
        G = 20000
      }

      uni_gmwm = .Call('gmwm_gmwm_master_cpp', PACKAGE = 'gmwm', data, theta, desc, obj = obj_desc,
                       model.type = 'imu' , starting = model_nested[[i]]$starting,
                       p = 0.05, compute_v = "fast", K = 1, H = 100, G = G,
                       robust=FALSE, eff = 1)[[1]]

      # Store the univariate gmwm
      param_starting[k,] = uni_gmwm
    }

    ## Select as starting value that minimize the mgmwm objective function.
    obj_value_starting_value = rep(NA, n_replicates)
    mgmwm_list = list()
    for (j in 1:n_replicates){
      starting_value = inv_param_transform(model_nested[[i]], param_starting[j,])
      mgmwm_list[[j]] = optim(starting_value, mgmwm_obj_function, model = model_nested[[i]], mimu = mimu)
      obj_value_starting_value[j] = mgmwm_list[[j]]$value
    }

    mgmwm_full_dataset = mgmwm_list[[which.min(obj_value_starting_value)]]


    for (d in 1:n_permutation){

      mimu_est = list()
      class(mimu_est) = "mimu"
      mimu_test = mimu[-pair[d,sequence(s_est)]]

      for (s in 1:s_est){
        mimu_est[[s]] = mimu[[pair[d,s]]]
      }

      out = optim(mgmwm_full_dataset$par, mgmwm_obj_function, model = model_nested[[i]], mimu = mimu_est)

      model_test[[i]] = model_nested[[i]]

      # Pass on the estimated paramters onto the model.
      model_test[[i]]$starting = FALSE
      model_test[[i]]$theta = out$par

      # Compute the WVIC on each permutation
      obj_out_sample[d,i] = mgmwm_obj_function(model_test[[i]]$theta, model_test[[i]], mimu_test)
    }

    # Pass on the estimated paramters onto the model.
    model_nested[[i]]$starting = FALSE
    model_nested[[i]]$theta =  param_transform(model_nested[[i]], mgmwm_full_dataset$par)

    # Pass obj_value to the model
    model_nested[[i]]$obj_value = mgmwm_full_dataset$value
    model_complexity[i] = model_nested[[i]]$plength


    # Compute the CV-WVIC
    cv_wvic[i] = mean(obj_out_sample[,i])
    # Update progress bar
    pb$tick(tokens = list(what = "foo   "))
    Sys.sleep(1 / n_models)
  }
  # Compute the average on all permutation of the out of sample WVIC
  mod_selected_cvwvic = which.min(cv_wvic)

  #Wilcoxon paired test on nested models
  wilcox_test = rep(FALSE,n_models)

  # Compute the Wilcoxon test
  for (i in 1:n_models){
    if(model_complexity[mod_selected_cvwvic] >= model_complexity[i]){
      #Compute the p-value for asll nested model
      wilcox_test[i] = wilcox.test(obj_out_sample[,i],obj_out_sample[,mod_selected_cvwvic],
                                   paired = T,alternative = "greater")$p.val
    }
  }
  # Decision rule on which model is equivalent
  test_wilcox_result = (wilcox_test > alpha_paired_test)

  # Index of model in "model_nested" which is equivalent to mod_selected_cvwvic
  index_select_wilcox_list = which(test_wilcox_result[1:n_models] == TRUE)

  # If index_select_wilcox_list has multiple model, select the smallest one
  if(length(index_select_wilcox_list) != 1){

    # Select equivalent models (0 means not equivalent)
    equivalent_models = model_complexity*test_wilcox_result
    index_equivalent_model  = which(equivalent_models != 0)
    equivalent_models[equivalent_models == 0] = NA

    # model complexity of the smallest one
    model_complexity_select_wilcox = model_complexity[which.min(equivalent_models)]

    # If more than one model of same size select, pick the one with the smallest cvwvic
    if(equivalent_models[1:n_models] == model_complexity_select_wilcox){

      equivalent_model_star = equivalent_models
      equivalent_model_star[equivalent_model_star != model_complexity_select_wilcox] = NA
      # select the one which min the cvwvic knowing the minimum size
      model_select_wilcox  = which.min(cv_wvic*equivalent_model_star)

    }else{
      model_select_wilcox = which.min(equivalent_models)
    }

  # Put the selected one in chosen model object
  model_hat = model_nested[[model_select_wilcox]]
  }else{
  model_hat = model_nested[[mod_selected_cvwvic]]
  }

   # Put the decision rule in model object
  for (i in 1:n_models){
    if(length(index_select_wilcox_list) != 1){
      if(i == model_select_wilcox){
        model_nested[[i]]$decision = "Model selected"
      }else if (i == index_equivalent_model & i != mod_selected_cvwvic){
        model_nested[[i]]$decision = "Bigger equivalent model"
      }else if (i == mod_selected_cvwvic){
        model_nested[[i]]$decision = "Model selected cvwvic"
      }else{
        model_nested[[i]]$decision = "Model not appropriate"
      }
    }else{
      if(i == mod_selected_cvwvic){
        model_nested[[i]]$decision = "Model selected cvwvic"
      }else{
        model_nested[[i]]$decision = "Model not appropriate"
      }
    }
  }

  ## Create the ouput for the selected model
  estimate = as.matrix(model_hat$theta)
  rownames(estimate) = model_hat$process.desc
  colnames(estimate) = "Estimates"

  # Selected model objective function value
  obj_value = model_hat$obj_value
  names(obj_value) = "Value Objective Function"


  # WV implied by the parameter
  wv_implied = wv_theo(model_hat, tau_max_vec)

  # Extact individual model for theoretical decomposition
  desc_decomp_theo = desc_decomp_theo_fun(model_hat, n_process)


  # Compute individual theoretical wv
  decomp_theo = list()
  for (i in 1:n_process){
    decomp_theo[[i]] =  wv_theo(desc_decomp_theo[[i]], tau_max_vec)
  }


}








