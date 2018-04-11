

param_name = function(model){

  param_name = list()

  counter = 1

  for (j in 1:length(model$desc)){

    # is quantization noise?
    if (model$process.desc[j] == "QN"){
      param_name[[counter]] = expression(paste(hat(Q)^2))
      counter = counter + 1
    }

    # is white noise?
    if (model$process.desc[j] == "WN"){
      param_name[[counter]] =  expression(paste(hat(sigma)^2))
      counter = counter + 1
    }

    # is AR1?
    if (model$process.desc[j] == "AR1"){

      param_name[[counter]] = expression(paste(hat(phi)))

      counter = counter + 1

    }

    # is random walk?
    if (model$process.desc[j] == "RW"){
      param_name[[counter]] = expression(paste(hat(gamma)^2))
      counter = counter + 1
    }

    # is random walk?
    if (model$process.desc[j] == "SIGMA2"){
      param_name[[counter]] =  expression(paste(hat(nu)^2))
      counter = counter + 1
    }

    # is drift?
    if (model$process.desc[j] == "DR"){
      param_name[[counter]] = expression(paste(hat(omega)^2))
      counter = counter + 1
    }

  }
  param_name
}


plot_CI = function(obj_list, units = NULL, xlab = NULL, ylab = NULL,
                     col_dens = NULL, col_ci = NULL, nb_ticks_x = NULL,
                     nb_ticks_y = NULL, ci_wv = NULL, point_cex = NULL,
                     point_pch = NULL,col_line = NULL, ...){

  n_param = obj_list$model_hat$plength


  windows_col = ceiling((obj_list$model_hat$plength)/2)
  windows_row = floor((obj_list$model_hat$plength)/2)

  #par(mfrow=c(windows_row,windows_row), mar = c(4,4,4,2))

  for (i in 1:n_param){

    density_param = (density(na.omit(obj_list$distrib_param[,i])))



    # Labels
    if (is.null(xlab)){
      xlab = obj_list$model_hat$process.desc[i]
    }else{
      xlab = xlab
    }
    ### to do
    if (is.null(ylab)){
      ylab = "Smoothed Density"
    }else{
      ylab = ylab
    }

    param_name = param_name(obj_list$model_hat)


    main = param_name[[i]]


    # Line and CI colors
    if (is.null(col_dens)){
      col_dens = hcl(h = 240, c = 65, l =70, alpha = 0.2, fixup = TRUE)
    }

    if (is.null(col_ci)){
      col_ci = hcl(h = 210, l = 65, c = 100, alpha = 1)
    }

    # Line and CI colors
    if (is.null(col_line)){
      col_line = "darkblue"
    }

    # Range
    x_range = c(obj_list$ci_low[[i]],obj_list$ci_high[[i]])
    x_low = obj_list$ci_low[[i]] - obj_list$ci_low[[i]]/10
    x_high = obj_list$ci_high[[i]] + obj_list$ci_high[[i]]/10

    y_range = range(density_param$y)
    y_low = min(y_range)
    y_high = max(y_range)

    x_ticks = obj_list$model_hat$theta[[i]]


    x_labels =  obj_list$model_hat$process.desc[i]

    # Main Plot
    plot(NA, xlim = x_range, ylim = c(y_low, y_high), xlab = xlab, ylab = ylab
         , yaxt = 'n' , bty = "n", ann = FALSE)
    win_dim = par("usr")

    par(new = TRUE)

    plot(NA, xlim = x_range, ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
         xlab = xlab, ylab = ylab, xaxt = 'n', yaxt = 'n', bty = "n")
    win_dim = par("usr")


    # Add Title
    x_vec = c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
    y_vec = c(win_dim[4], win_dim[4],
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
              win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
    polygon(x_vec, y_vec, col = "grey95", border = NA)
    text(x = mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

    # Add Axes and Box
    lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = 1)
    #y_ticks = y_ticks[(2^y_ticks) < 10^(win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))]
    box()
    axis(1, at = x_ticks, padj = 0.3)
    lines(density_param$x,density_param$y, type = "l", col = col_line)
    polygon(c(obj_list$ci_low[[i]],density_param$x[density_param$x>=obj_list$ci_low[[i]] & density_param$x<=obj_list$ci_high[[i]]],obj_list$ci_high[[i]])
            ,c(obj_list$ci_low[[i]],density_param$y[density_param$x>=obj_list$ci_low[[i]] & density_param$x<=obj_list$ci_high[[i]]],obj_list$ci_hig[[i]]),
            col=col_dens, border = F)
    lines(rep(obj_list$model_hat$theta[[i]],2), c(win_dim[3],
                                                  win_dim[4] - 0.09*(win_dim[4] - win_dim[3])), col = "red", lwd = 2)

  }

}




