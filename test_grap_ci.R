# Distribution of p1

dat1 <- with(density(na.omit(test.optim1$distrib_param[,1])))
quant11 = as.numeric(quantile(na.omit(test.optim1$distrib_param[,1]),0.975))
quant21 = as.numeric(quantile(na.omit(test.optim1$distrib_param[,1]),0.025))


plot_CI = function(x, units = NULL, xlab = NULL, ylab = NULL, main = NULL,
                     col_dens = NULL, col_ci = NULL, nb_ticks_x = NULL,
                     nb_ticks_y = NULL, ci_wv = NULL, point_cex = NULL,
                     point_pch = NULL, ...){

  density_param = (density(na.omit(test.optim1$distrib_param[,1])))
  quant11 = as.numeric(quantile(na.omit(test.optim1$distrib_param[,1]),0.975))
  quant21 = as.numeric(quantile(na.omit(test.optim1$distrib_param[,1]),0.025))


  # Labels
  if (is.null(xlab)){
    xlab = test.optim1$model_hat$process.desc[1]
  }else{
    xlab = xlab
  }
  ### to do
  if (is.null(ylab)){
    ylab = "Smoothed Density"
  }else{
    ylab = ylab
  }

  # Main Title
  if (is.null(main)){
    main = "CI"
  }

  # Line and CI colors
  if (is.null(col_dens)){
    col_dens = "darkblue"
  }

  if (is.null(col_ci)){
    col_ci = hcl(h = 210, l = 65, c = 100, alpha = 0.2)
  }

  # Line and CI colors
  if (is.null(col_line)){
    col_line = "darkblue"
  }

  # Range
  x_range = range(density_param$x)
  x_low = min(x_range)
  x_high = 1

  y_range = range(density_param$y)
  y_low = 0
  y_high = max(y_range)

  x_ticks = test.optim1$model_hat$theta[[1]]


  # Main Plot
  plot(NA, xlim = x_range, ylim = c(y_low, y_high), xlab = xlab, ylab = ylab
       , xaxt = 'n', bty = "n", ann = FALSE)
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
  axis(1, at = x_ticks, labels = x_labels, padj = 0.3)
  lines(density_param$x,density_param$y, type = "l", col = col_line)
  polygon(c(quant21,density_param$x[density_param$x>=quant21 & density_param$x<=quant11],quant11)
          ,c(quant21,density_param$y[density_param$x>=quant21 & density_param$x<=quant11],quant11), col="lightgrey")
  lines(rep(test.optim1$model_hat$theta[[1]],2), c(win_dim[3],
                                                   win_dim[4] - 0.09*(win_dim[4] - win_dim[3])))

}








#QN, WN, AR1, RW DR



