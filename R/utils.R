#' @export
is.mimu = function(obj){class(obj) == "mimu"}

#' @export
is.mgmwm = function(obj){class(obj) == "mgmwm"}

#' @export
is.cvwvic = function(obj){class(obj) == "cvwvic"}

#' @export
is.tsmodel = function(obj){class(obj) == "ts.model"}

#' @export
# ----- phi
transform_phi = function(phi_R){
  phi = 2*(inv_logit(phi_R)-1/2)
  phi
}

#' @export
# ----- phi
inv_transform_phi = function(phi){
  phi_R = log(2/(1-phi)-1)
  phi_R
}

#' @export
# ----- weights
# use inv_logit function to transform weights (to map R -> [0,1])
inv_logit = function(x){
  exp(x)/(1 + exp(x))
}

#' @export
comb.mat = function(n){
  c = rep(list(1:0), n)
  expand.grid(c)
}

#' @export
param_transform = function(model, theta1){

  np = model$plength

  # Initialise counter
  counter = 1

  for (j in 1:np){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      theta1[counter] = exp(theta1[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      theta1[counter] = exp(theta1[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      theta1[counter] = exp(theta1[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      theta1[counter] = exp(theta1[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      theta1[counter] = transform_phi(theta1[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      theta1[counter] = exp(theta1[counter])
      counter = counter + 1
    }
  }
  theta1
}

#' @export
inv_param_transform = function(model,theta1){

  np = model$plength

  # Initialise counter
  counter = 1

  for (j in 1:np){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      theta1[counter] = log(theta1[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      theta1[counter] = log(theta1[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      theta1[counter] = log(theta1[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      theta1[counter] = log(theta1[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      theta1[counter] = inv_transform_phi(theta1[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      theta1[counter] = log(theta1[counter])
      counter = counter + 1
    }
  }
  theta1
}

#' @export
model_combination = function(model){

  # String description of model
  model_desc_max = model$desc

  # number of latent process in model max

  n_process_max = length(model_desc_max)


  # Build matrix of possible combination
  m = as.matrix(comb.mat(n_process_max))
  m = m[-nrow(m),]

  models_names =build_model_set(m,model_desc_max)

  n_models =length(models_names)

  all_model = list()


  for (i in 1:n_models){

    model_test = model

    model_test$desc = models_names[[i]]

    n_process = length(model_test$desc)

    model_test$starting = TRUE


    n_para = 0
    for (j in 1: n_process){
      if(model_test$desc[[j]] == "AR1"){
        n_para = n_para +  2
      }else{
        n_para = n_para + 1
      }
      model_test$plength =  n_para
    }

    process_desc_change = rep(NA,model_test$plength)
    theta_test = rep(NA,model_test$plength)
    obj_desc_list = list()

    counter = 1

    for (j in 1: n_process){
      if (model_test$desc[[j]] == "AR1"){
        process_desc_change[counter] = "AR1"
        process_desc_change[counter + 1] = "SIGMA2"
        theta_test[counter] = 0
        theta_test[counter + 1] = 1
        obj_desc_list[[j]] = c(1,1)

        counter = counter + 2
      }

      if (model_test$desc[[j]] == "WN"){
        process_desc_change[counter] = "WN"
        theta_test[counter] = 3
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "RW"){
        process_desc_change[counter] = "RW"
        theta_test[counter] = 4
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "QN"){
        process_desc_change[counter] = "QN"
        theta_test[counter] = 2
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }

      if (model_test$desc[[j]] == "DR"){
        process_desc_change[counter] = "DR"
        theta_test[counter] = 5
        obj_desc_list[[j]] = 1
        counter = counter + 1
      }
      model_test$process.desc =  process_desc_change
      model_test$theta =  theta_test
      model_test$obj.desc =  obj_desc_list
    }
    all_model[[i]] = model_test
  }
  all_model
}

#' @title Objective function Multivariate GMWM
#'
#' @description
#' Compute the objective function Multivariate GMWM, based on an \code{mimu} object and a \code{model}
#' @keywords internal
#' @param theta            A \code{vector} of parameter from a model.
#' @param model            A model \code{object} of class \code{ts.model}.
#' @param mimu             A mimu \code{object}.
#' @return \code{numeric} value of the objective function
#' @author Gaetan Bakalli
#' @keywords internal
#' @export
mgmwm_obj_function = function(theta, model, mimu){

  M = length(model$process.desc)

  # Initialise counter
  counter = 1

  for (j in 1:M){
    # is random walk?
    if (model$process.desc[j] == "RW"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){
      theta[counter] = transform_phi(theta[counter])
      counter = counter + 1
    }

    # is SIGMA2?
    if (model$process.desc[j] == "SIGMA2"){
      theta[counter] = exp(theta[counter])
      counter = counter + 1
    }
  }

  model$theta = theta

  # Step 1: compute theoretical WV
  tau = list()
  wv.theo = list()
  obj_len  = length(mimu)

  for (i in 1:obj_len) {
    tau[[i]] = 2^(1:length(mimu[[i]]$scales))
    # Compute theoretical wvar for each latent process
    wv.theo[[i]] = wv_theo(model, tau[[i]])
  }


  # Step 2: compute Omega
  Omega = list()
  for (i in 1:obj_len){
    Omega[[i]] = diag(1/(mimu[[i]]$ci_high - mimu[[i]]$ci_low)^2)
  }

  # Step 3: compute actual objective function
  out = 0
  for (i in 1:obj_len){
    dif_vect = wv.theo[[i]] - mimu[[i]]$variance
    out = out + t(dif_vect)%*%Omega[[i]]%*%dif_vect
  }
  out
}


#' @title Theoretical Wavelet Variance
#'
#' @description
#' Compute theroretical Wavelet Variance
#' @keywords internal
#' @param model            A \code{model} object.
#' @param tau              A \code{vector} of \eqn{2^1:number of scales}
#' @return \code{vector} of theoretical Wavelet Variance for a given model
#' @author Stephane Guerrier, Nathanael Claussen, and Justin Lee
#' @export
#' @examples
#' #' Xt = rnorm(n/4)
#' Yt = rnorm(n/2) + cumsum(rnorm(n/2, 0, 10^(-2)))
#' Zt = rnorm(n) + cumsum(rnorm(n, 0, 10^(-3)))
#' obj = make_wvar_mimu_obj(Xt, Yt, Zt, freq = 100, unit = "s",
#' sensor.name = "MTiG - Gyro. X", exp.name = c("today", "yesterday", "a few days ago"))
#' model = AR1() + WN()
#'
#' Extract the scales of first replicates of mimu \code{object} and compute tau
#' tau = 2^(1:length(sensor.name[[1]]$scales))
#' wv_theo = function(model, tau)



wv_theo = function(model, tau){
  if (!class(model) == "ts.model"){
    # Error
  }

  # Number of scales
  J = length(tau)

  # Extract process description
  desc = model$desc

  # Number of latent processes
  M = length(desc)

  # Extract parameters
  theta = model$theta

  # Initialise counter
  counter = 1

  # Compute theoretical wvar for each latent process
  wv = matrix(NA, M, J)
  for (i in 1:M){
    # is random walk?
    if (desc[i] == "RW"){
      wv[i, ] = wv::rw_to_wv(gamma2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is white noise?
    if (desc[i] == "WN"){
      wv[i, ] = wv::wn_to_wv(sigma2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is drift?
    if (desc[i] == "DR"){
      wv[i, ] = wv::dr_to_wv(omega = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is quantization noise?
    if (desc[i] == "QN"){
      wv[i, ] = wv::qn_to_wv(q2 = theta[counter], tau = tau)
      counter = counter + 1
    }

    # is AR1?
    if (desc[i] == "AR1"){
      wv[i, ] = wv::ar1_to_wv(phi = theta[counter], sigma2 = theta[counter + 1], tau = tau)
      counter = counter + 2
    }
  }
  # Summing them up and return
  apply(wv, 2, sum)
}

#' @export
#'
model_names = function(model){

  model_name = list()
  n_ar1 = sum(model$desc == "AR1")

  ar1_espression = bquote(paste(n_ar1,"*Wavelet Variance "))

  for (j in 1:length(model$desc)){
    counter = 1

    # is quantization noise?
    if (model$desc[j] == "QN"){
      model_name[[1]] = "QN"
      counter = counter + 1
    }

    # is white noise?
    if (model$desc[j] == "WN"){
      model_name[[2]] =  "WN"
      counter = counter + 1
    }

    # is AR1?
    if (model$desc[j] == "AR1"){
      if(n_ar1 == 1){
        model_name[[3]] = "AR1"
      }else{
        model_name[[3]] = paste0(c(n_ar1,"AR1"), collapse = '*')
      }
      counter = counter + 1

    }

    # is random walk?
    if (model$desc[j] == "RW"){
      model_name[[4]] = "RW"
      counter = counter + 1
    }

    # is drift?
    if (model$desc[j] == "DR"){
      model_name[[5]] =  "DR"
      counter = counter + 1
    }

  }
  x = unlist(model_name)
  paste0(unlist(model_name), collapse = '+')
}

