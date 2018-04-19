# Copyright (C) 2014 - 2018  Gaetan Bakalli, Stephane Guerrier.
#
# This file is part of classimu R Methods Package
#
# The `mgmwm` R package is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# The `mgmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Multivariate Generalized Method of Wavelet Moments (MGMWM) for IMUs
#'
#' @export
#' @import progress
model_selection = function(mimu, model, s_est = NULL, alpha_ci_cvwvic = NULL,
                           b_ci_wvic = 500, seed = 2710){

  # Level of confidence for Wilcoxon paired test
  if(is.null(alpha_ci_cvwvic)){
    alpha_ci_cvwvic = .05
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
  model_nested = model_combination(model)

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
    format = "  Model :current of :total Models. Estimated remaining time:  :eta",
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
    options(warn=-1)
    for (j in 1:n_replicates){
      starting_value = inv_param_transform(model_nested[[i]], param_starting[j,])
      mgmwm_list[[j]] = optim(starting_value, mgmwm_obj_function, model = model_nested[[i]], mimu = mimu)
      obj_value_starting_value[j] = mgmwm_list[[j]]$value
    }
    options(warn=0)
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


  ############## Redefine equivalent models #####################



  # Create number of model for which to plot adjusted CI
  n_models_adj = n_models - 1

  # Extract vector of out-of-sample wvic for selected model
  selected_model_vector_wvic = as.vector(obj_out_sample[,mod_selected_cvwvic])

  # define the new matrix without selected model
  mat_criteria = obj_out_sample[,-mod_selected_cvwvic]

  # substract the vector of objective function from selected model
  mat_criteria_scaled = mat_criteria - matrix(rep(selected_model_vector_wvic, each = n_models_adj),
                                              n_permutation, n_models_adj, byrow = TRUE)

  # Initialize
  matrix_boot_wvic = matrix(NA,b_ci_wvic, n_models_adj)

  for (b in 1:b_ci_wvic){
    matrix_boot_wvic[b,] =  apply(mat_criteria_scaled[sample(1:n_permutation, replace = TRUE),],2,mean)
  }

  ci_low_wvic = rep(NA,n_models_adj)
  ci_high_wvic = rep(NA,n_models_adj)
  med_wvic = rep(NA,n_models_adj)

  for (i in 1:n_models_adj){
    med_wvic[i] = quantile(matrix_boot_wvic[,i], probs = 0.5)
    ci_low_wvic[i] = quantile(matrix_boot_wvic[,i], probs = alpha_ci_cvwvic/2)
    ci_high_wvic[i] = quantile(matrix_boot_wvic[,i], probs = 1 - alpha_ci_cvwvic/2)
  }

  ############## Redefine equivalent models

  selection_decision = rep(NA,n_models)
  selection_decision[mod_selected_cvwvic] = "Model selected"
  selection_decision_seq = which(is.na(selection_decision))
  model_hat =  model_nested[[i]]

  # Put the decision rule in model object
  for (i in 1:n_models_adj){
    if (ci_low_wvic[i] <=0){
      if(model_complexity[selection_decision_seq[i]] <= model_complexity[mod_selected_cvwvic]){
        selection_decision[selection_decision_seq[i]] = "Equivalent model"
      }else{
        selection_decision[selection_decision_seq[i]] = "Bigger equivalent model"
      }
    }else{
      selection_decision[selection_decision_seq[i]] = "Model not appropriate"
    }
  }

  for (i in 1:n_models){
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
  colnames(estimate) = "Estimates model selected"

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
  options(warn=-1)
  set.seed(as.numeric(format(Sys.time(),"%s"))/10)
  options(warn=0)

  out_model_selection = structure(list(estimate = estimate,
                                       obj_value = obj_value,
                                       model_hat = model_hat,
                                       model_nested = model_nested,
                                       selection_decision = selection_decision,
                                       cv_wvic = cv_wvic,
                                       med_wvic = med_wvic,
                                       ci_low_wvic = ci_low_wvic,
                                       ci_high_wvic = ci_high_wvic,
                                       model_name = model_name,
                                       scales_max_vec = scales_max_vec,
                                       obj_out_sample = obj_out_sample,
                                       mimu = mimu), class = "cvwvic")
  invisible(out_model_selection)

}


#' @export plot.cvwvic
plot.cvwvic = function(obj_list, decomp = TRUE, type = NULL, model = NULL,
                       add_legend_mgwmw = TRUE, legend_pos = NULL,
                       ylab = NULL, couleur_axis = FALSE){


  n_models = length(obj_list$model_nested)

  if (is.tsmodel(model)){
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

  if (is.null(ylab)){
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
                 add_legend_mgwmw = TRUE, legend_pos = NULL, ylab_mgmwm = ylab)
      title(main = "Model selected")
    }else if(type == "equivalent"){
      decomp = FALSE

      index_equiv = which(obj_list$selection_decision == "Model selected" | obj_list$selection_decision == "Equivalent model")

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
    }else if (type == "wvic_all"){
      model_WVIC_CI(obj_list, type = type)
    }else if (type == "wvic_equivalent"){
      model_WVIC_CI(obj_list, type = type)
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


#' @export model_WVIC_CI
model_WVIC_CI = function(obj_list, type = "wvic_all"){

  # Compute the number of model to plot
  n_models = length(obj_list$model_name)
  n_models_adj = n_models -1

  # Extract index of selectec model
  selected_model_index = which.min(obj_list$cv_wvic)

  # Order models with respect to the wvic
  model_ord = order(obj_list$cv_wvic[-selected_model_index])

  #Order model name
  model_names = obj_list$model_name[-selected_model_index]
  model_names = model_names[model_ord]

  # Order model description
  mod_des = obj_list$selection_decision[-selected_model_index]
  mod_des_ord = mod_des[model_ord]

  # Order median and ci
  med_wvic_ord = obj_list$med_wvic[model_ord]

  # Order median and ci
  ci_low_wvic_ord = obj_list$ci_low_wvic[model_ord]

  # Order median and ci
  ci_high_wvic_ord = obj_list$ci_high_wvic[model_ord]



  hues = seq(15, 375, length = n_models)
  couleur = hcl(h = hues, l = 65, c = 100, alpha = 1)[seq_len(n_models)]

  if(type == "wvic_all"){


    col_desc = c("Black", "Black", couleur[3])

    names(col_desc) = c("Model not appropriate", "Bigger equivalent model", "Equivalent model")


    ci_low_neg_index = which(ci_low_wvic_ord <= 0)
    ci_low_neg = ci_low_wvic_ord[ci_low_neg_index]
    ci_low_pos_index = which(ci_low_wvic_ord > 0)
    ci_low_pos =  ci_high_wvic_ord[ci_low_pos_index]
    left_bound = min(med_wvic_ord)


    xlab = paste("CV-WVIC (Model) - CV-WVIC (", obj_list$model_name[selected_model_index],")", sep="")
    ylab = " "
    main = "CI for CV-WVIC"
    par(oma = c(0.1,6.5,0,0))
    plot(NA, xlim = c(left_bound, max(ci_high_wvic_ord)), ylim = c(1,n_models_adj), ylab = ylab, xlab = NULL,
         xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE, log = "x")
    win_dim = par("usr")
    par(new = TRUE)
    plot(NA, xlim = c(left_bound , max(ci_high_wvic_ord)), ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
         ylab = ylab, xlab = " ", xaxt = 'n', yaxt = 'n', bty = "n", log = "x")
    graphics::box()
    mtext(xlab, side = 1, line = 2.5)
    mtext(ylab, side = 2, line = 0)
    win_dim = par("usr")

    # Add grid
    grid(NULL, NA, lty = 1, col = "grey95")
    abline(h = 1:n_models_adj, lty = 1, col = "grey95")

    # Add title
    x_vec = 10^c(win_dim[1] , win_dim[2], win_dim[2], win_dim[1])
    y_vec = c(win_dim[4], win_dim[4],
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
    polygon(x_vec, y_vec, col = "grey95", border = NA)
    text(x = 10^mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

    # Add axes and box
    lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = "grey50")

    axis(1, padj = 0.3)

    y_axis = axis(2, labels = FALSE, tick = FALSE)
    y_axis = 1:n_models_adj

    for (i in 1:n_models_adj){
      axis(2, padj = 0.5, at = i, las = 1, labels = model_names[i],
      col.axis = col_desc[mod_des_ord[i]], cex.axis = 0.8)
    }

    for (i in ci_low_neg_index){
      lines(c(10^(-30), ci_high_wvic_ord[i]), c(i,i), col = couleur[i], lty = 2)
      points(c(ci_low_wvic_ord[i], ci_high_wvic_ord[i]), c(i,i), pch = "|", col = couleur[i])
      points(med_wvic_ord[i], i, pch = 16, cex = 1.5, col = couleur[i])
    }

    for (i in ci_low_pos_index){
      lines(c(ci_low_wvic_ord[i], ci_high_wvic_ord[i]), c(i,i), col = couleur[i])
      points(c(ci_low_wvic_ord[i], ci_high_wvic_ord[i]), c(i,i), pch = "|", col = couleur[i])
      points(med_wvic_ord[i], i, pch = 16, cex = 1.5, col = couleur[i])
    }
  }else if(type == "wvic_equivalent"){

    mod_equivalent_ci_index = which(mod_des_ord == "Equivalent model")

    n_models_equiv = length(mod_equivalent_ci_index)

    med_equivalent_ci = med_wvic_ord[mod_equivalent_ci_index]
    ci_low_equivalent_ci = ci_low_wvic_ord[mod_equivalent_ci_index]
    ci_high_equivalent_ci = ci_high_wvic_ord[mod_equivalent_ci_index]
    model_names_equivalent_ci = model_names[mod_equivalent_ci_index]


    xlab = paste("CV-WVIC (Model) - CV-WVIC (", obj_list$model_name[selected_model_index],")", sep="")
    ylab = " "
    main = "CI for CV-WVIC"
    par(oma = c(0.1,6.5,0,0))
    plot(NA, xlim = c(min(ci_low_equivalent_ci), max(ci_high_equivalent_ci)), ylim = c(1,n_models_equiv), ylab = ylab, xlab = NULL,
         xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
    win_dim = par("usr")
    par(new = TRUE)
    plot(NA, xlim = c(min(ci_low_equivalent_ci) , max(ci_high_equivalent_ci)), ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
         ylab = ylab, xlab = " ", xaxt = 'n', yaxt = 'n', bty = "n")
    graphics::box()
    mtext(xlab, side = 1, line = 2.5)
    mtext(ylab, side = 2, line = 0)
    win_dim = par("usr")

    # Add grid
    grid(NULL, NA, lty = 1, col = "grey95")
    abline(h = 1:n_models_equiv, lty = 1, col = "grey95")

    # Add title
    x_vec = c(win_dim[1] , win_dim[2], win_dim[2], win_dim[1])
    y_vec = c(win_dim[4], win_dim[4],
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
    polygon(x_vec, y_vec, col = "grey95", border = NA)
    text(x = mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

    # Add axes and box
    lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = "grey50")

    axis(1, padj = 0.3)

    y_axis = axis(2, labels = FALSE, tick = FALSE)
    y_axis = 1:n_models_equiv

    for (i in 1:n_models_equiv){
      axis(2, padj = 0.5, at = i, las = 1, labels = model_names_equivalent_ci[i])
    }

    for (i in 1:n_models_equiv){
      lines(c(ci_low_equivalent_ci[i], ci_high_equivalent_ci[i]), c(i,i), col = couleur[i], lty = 2)
      points(c(ci_low_equivalent_ci[i], ci_high_equivalent_ci[i]), c(i,i), pch = "|", col = couleur[i])
      points(med_equivalent_ci[i], i, pch = 16, cex = 1.5, col = couleur[i])
    }
  }
}


#'@export summary.cvwvic
summary.cvwvic = function(object){

  out = round(object$cv_wvic, digits = 3)


  ms  = matrix(NA, length(out), 1)
  rownames(ms) = object$model_name
  colnames(ms) = c("CV-WVIC")
  ms[,1] =  out



  x = structure(list(estimates = object$estimate,
                     obj_value = object$obj_value,
                     cv_wvic_nested_model = ms))
  x
}


