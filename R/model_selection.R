

#' @export
#' @import progress
model_selection = function(mimu, model, s_est = NULL,
                           alpha_paired_test = NULL, seed = 2710){

  # Level of confidence for Wilcoxon paired test
  if(is.null(alpha_paired_test)){
    alpha_paired_test = .05
  }

  # Number of replicates
  n_replicates = length(mimu)

  # Number of time series replicates for estimation
  if(is.null(s_est)){
    s_est = ceiling(n_replicates/2)
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
  test_wilcox_result = wilcox_test > alpha_paired_test

  # Which model are equivalent in the sense of Wilcoxon test
  equiv_mod = which(test_wilcox_result)

  # Which model is the smallest in terms of parameter
  smallest_equiv = (min(model_complexity[test_wilcox_result]) == model_complexity) & test_wilcox_result

  # If there is more than one model, pick the one with the samllest cvwvic
  if (sum(smallest_equiv) > 1){
    smallest_equiv = (min(cv_wvic[smallest_equiv]) == cv_wvic) & smallest_equiv
  }
  smallest_equiv = which(smallest_equiv)

  selection_decision = rep(NA,n_models)

   # Put the decision rule in model object
  for (i in 1:n_models){
    if(sum(equiv_mod) > 1){
      if(i == smallest_equiv){
        selection_decision[i] = "Model selected"
        model_hat =  model_nested[[i]]
      }else if (i == equiv_mod && i != mod_selected_cvwvic){
        selection_decision[i] = "Bigger equivalent model"
      }else if (i == mod_selected_cvwvic){
        selection_decision[i] = "Model selected cv-wvic"
      }else{
        selection_decision[i] = "Model not appropriate"
      }
    }else{
      if(i == mod_selected_cvwvic){
        selection_decision[i] = "Model selected"
        model_hat =  model_nested[[i]]
      }else{
        selection_decision[i] = "Model not appropriate"
      }
    }

    # WV implied by the parameter
    model_nested[[i]]$wv_implied = wv_theo(model_nested[[i]], tau_max_vec)

    # Extact individual model for theoretical decomposition
    model_nested[[i]]$desc_decomp_theo = desc_decomp_theo_fun(model_nested[[i]], length(model_nested[[i]]$desc))

    # Compute individual theoretical wv
    decomp_theo = list()
    for (j in 1:length(model_nested[[i]]$desc)){
      decomp_theo[[j]] =  wv_theo(model_nested[[i]]$desc_decomp_theo[[j]], tau_max_vec)
    }

    model_nested[[i]]$decomp_theo = decomp_theo
  }

  model_name = rep(NA,n_models)
  for (i in 1:n_models){
    model_name[i] = model_names(model_nested[[i]])
  }

  ## Create the ouput for the selected model
  estimate = as.matrix(model_hat$theta)
  rownames(estimate) = model_hat$process.desc
  colnames(estimate) = "Estimates"

  # Selected model objective function value
  obj_value = model_hat$obj_value
  names(obj_value) = "Value Objective Function"

  # Transform the parameter
  #model_hat$theta = param_transform(model_hat,model_hat$theta)

  # WV implied by the parameter
  model_hat$wv_implied = wv_theo(model_hat, tau_max_vec)

  # Extact individual model for theoretical decomposition
  model_hat$desc_decomp_theo = desc_decomp_theo_fun(model_hat, length(model_hat$desc))

  # Compute individual theoretical wv
  decomp_theo = list()
  for (i in 1:length(model_hat$desc)){
    decomp_theo[[i]] =  wv_theo(model_hat$desc_decomp_theo[[i]], tau_max_vec)
  }
  model_hat$decomp_theo = decomp_theo

  # Cancel previous seed
  set.seed(as.numeric(format(Sys.time(),"%s"))/10)

  out_model_selection = structure(list(estimate = estimate,
                                  obj_value = obj_value,
                                  model_hat = model_hat,
                                  model_nested = model_nested,
                                  selection_decision = selection_decision,
                                  cv_wvic = cv_wvic,
                                  model_name = model_name,
                                  scales_max_vec = scales_max_vec,
                                  obj_out_sample = obj_out_sample,
                                  mimu = mimu), class = "cvwvic")
  invisible(out_model_selection)

}


#' @export
plot.cvwvic = function(obj_list, decomp = TRUE, type = NULL, model = NULL,
                       add_legend_mgwmw = TRUE, legend_pos = NULL,
                       ylab_cvwvic = NULL, couleur_axis = FALSE){


  n_models = length(obj_list$model_nested)

  if (is.ts.model(model)){
    if (!is.null(type)){
      warning("type set to NULL.")
    }
    type = NULL
    model_name = model_names(model)

    if (sum(obj_list$model_name %in% model_name) == 0){
      stop("This model has not been estimated. Use this model has an input of the `model_selection` fuction.")
    }

    index_model = which(obj_list$model_name %in% model_name)
  }else{
    if (is.null(type)){
      type = "selected"
    }
  }

  if (is.null(ylab_cvwvic)){
    ylab = expression(paste("Wavelet Variance ", nu^2, sep = ""))
  }

  obj_plot = list()

  for (i in 1:n_models){
    obj_plot[[i]] = list(mimu = obj_list$mimu, model_hat = obj_list$model_nested[[i]],
                        scales_max_vec = obj_list$scales_max_vec)
  }

  # Extract index of equivalent model
  index_equiv = which(obj_list$selection_decision == "Model selected" | obj_list$selection_decision == "Model selected cv-wvic" | obj_list$selection_decision == "Bigger equivalent model")

  if (!is.null(type)){
    if(type == "selected"){
      plot.mgmwm(obj_plot[[which(obj_list$selection_decision == "Model selected")]], decomp = decomp,
                 add_legend_mgwmw = TRUE, legend_pos = NULL, ylab_mgmwm = NULL)
      title(main = "Model selected")
    }else if((type == "cvwvic")){
      if (sum((obj_list$selection_decision %in% "Model selected cv-wvic")) == 1){
        compare_to = "Model selected cv-wvic"
      }else{
        compare_to = "Model selected"
      }
      plot.mgmwm(obj_plot[[which(obj_list$selection_decision == compare_to)]],
                 decomp = decomp,add_legend_mgwmw = TRUE, legend_pos = NULL, ylab_mgmwm = NULL)
      title(main = compare_to)
    }else if(type == "equivalent"){
      decomp = FALSE

      index_equiv = which(obj_list$selection_decision == "Model selected" | obj_list$selection_decision == "Model selected cv-wvic" | obj_list$selection_decision == "Bigger equivalent model")

      plot(obj_plot[[1]]$mimu, add_legend = FALSE, transparency_wv = 0.4, transparency_ci = 0.05, ylab = ylab)

      U = length(index_equiv)
      col_wv = hcl(h = seq(100, 375, length = U + 1), l = 65, c = 200, alpha = 1)[1:U]
      legend_names = rep(NA,length(index_equiv))
      col_legend = rep(NA,length(index_equiv))
      for (i in 1:length(index_equiv)){
        # Plot implied WV
        lines(t(obj_list$scales_max_vec),obj_plot[[index_equiv[i]]]$model_hat$wv_implied, type = "l", lwd = 3, col = col_wv[[i]], pch = 1, cex = 1.5)
        lines(t(obj_list$scales_max_vec),obj_plot[[index_equiv[i]]]$model_hat$wv_implied, type = "p", lwd = 2, col = col_wv[[i]], pch = 1, cex = 1.5)

        legend_names[i] = obj_list$model_name[index_equiv[i]]
        col_legend[i] = col_wv[i]
        p_cex_legend = rep(c(1.5,NA),length(obj_list$wv_implied))

        if (is.null(legend_pos)){
          legend_pos = "bottomleft"
        }
        if (add_legend_mgwmw == TRUE){
          legend(legend_pos, legend_names, bty = "n", lwd = 1, pt.cex = 1.5, pch = p_cex_legend, col = col_legend)
        }
      }
    par(old.par)
    }else if (type == "compare"){
      model_WVIC_CI(obj_list, boot_ci = 500, alpha = 0.05, couleur_axis = FALSE)
    }else{
      stop("Please define a valid model")
    }
  }else{
    if (is.ts.model(model)){
      plot.mgmwm(obj_plot[[index_model]], decomp = decomp,
                 add_legend_mgwmw = TRUE, legend_pos = NULL,
                 ylab_mgmwm = NULL)
    }else{
      stop("Please define a valid model")
    }
  }
}



model_WVIC_CI = function(test_model_selection, boot_ci = 500, alpha = 0.05, couleur_axis = FALSE){

  dec = test_model_selection$selection_decision

  n_models = length(test_model_selection$model_name)
  n_permutation = dim(test_model_selection$obj_out_sample)[1]
  model_ord = order(test_model_selection$cv_wvic)

  hues = seq(15, 375, length = n_models + 1)
  couleur = hcl(h = hues, l = 65, c = 100, alpha = 1)[seq_len(n_models+1)]

  if (couleur_axis){
    col_desc = c("Black", couleur[3], couleur[5], couleur[7])
  }else{
    col_desc = rep("Black",4)
  }
  names(col_desc) = c("Model not appropriate", "Bigger selected model","Model selected", "Model selected cv-wvic")

  matrix_ci_cvwvic = matrix(NA, boot_ci, n_models)
  matrix_boot = matrix(NA,boot_ci, n_models)

  for (b in 1:boot_ci){
    matrix_boot[b,] =  apply(test_model_selection$obj_out_sample[sample(1:n_permutation, replace = TRUE),],2,mean)
  }

  sd_cvwvic = rep(NA,n_models)
  ci_low = rep(NA,n_models)
  ci_high = rep(NA,n_models)
  coef = qnorm(1-alpha/2)
  for (i in 1:n_models){
    sd_cvwvic[i] = sd(matrix_boot[,i])
    ci_low[i] = test_model_selection$cv_wvic[i] - coef*sd_cvwvic[i]
    ci_high[i] = test_model_selection$cv_wvic[i] + coef*sd_cvwvic[i]
  }




  #par(mar = c(0.7, 2,0,0), oma = c(4,6.2,1,1))

  xlab = "CV-WVIC"
  ylab = " "
  main = "CI for CV-WVIC"
  plot(NA, xlim = c(min(ci_low), max(ci_high)), ylim = c(1,n_models), xlab = xlab, ylab = ylab,
       xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE, log = "x")
  win_dim = par("usr")
  par(new = TRUE)
  plot(NA, xlim = c(min(ci_low), max(ci_high)), ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
       ylab = ylab, xlab = xlab, xaxt = 'n', yaxt = 'n', bty = "n", log = "x")
  mtext(xlab, side = 1, line = 2.5)
  mtext(ylab, side = 2, line = 0)
  win_dim = par("usr")

  # Add grid
  grid(NULL, NA, lty = 1, col = "grey95")
  abline(h = 1:n_models, lty = 1, col = "grey95")

  # Add title
  x_vec = 10^c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
  y_vec = c(win_dim[4], win_dim[4],
            win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
            win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
  polygon(x_vec, y_vec, col = "grey95", border = NA)
  text(x = 10^mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

  # Add axes and box
  lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = "grey50")
  box(col = "grey50")

  axis(1, padj = 0.3)

  y_axis = axis(2, labels = FALSE, tick = FALSE)
  y_axis = 1:n_models

  for (i in 1:n_models){
    axis(2, padj = 0.5, at = i, las = 1, labels = test_model_selection$model_name[model_ord[i]],
         col.axis = col_desc[test_model_selection$selection_decision[model_ord[i]]])
  }

  for (i in 1:n_models){
    lines(c(ci_low[model_ord[i]], ci_high[model_ord[i]]), c(i,i), col = couleur[i])
    points(c(ci_low[model_ord[i]], ci_high[model_ord[i]]), c(i,i), pch = "|", col = couleur[i])
    points(test_model_selection$cv_wvic[model_ord[i]], i, pch = 16, cex = 1.5, col = couleur[i])
  }

}


