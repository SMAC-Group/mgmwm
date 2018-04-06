
# CI for cvwvic
boot_ci = 500

matrix_ci_cvwvic = matrix(NA, boot_ci, n_models)


n_permutation = 6

matrix_boot = matrix(NA,boot_ci, n_models)

for (b in 1:boot_ci){
  matrix_boot[b,] =  apply(test_model_selection$obj_out_sample[sample(1:n_permutation, replace = TRUE),],2,mean)
}

sd_cvwvic = rep(NA,n_models)
ci_low = rep(NA,n_models)
ci_high = rep(NA,n_models)
for (i in 1:n_models){
  sd_cvwvic[i] = sd(matrix_boot[,i])
  ci_low[i] = test_model_selection$cv_wvic[i] - sd_cvwvic[i]
  ci_high[i] = test_model_selection$cv_wvic[i] + sd_cvwvic[i]
}



range(test_model_selection$cv_wvic)


par(mar = c(0.7, 2,0,0), oma = c(4,6.2,1,1))

xlab = "CV-WVIC"
ylab = "Models"
main = "CI for CV-WVIC"
plot(NA, xlim = range(test_model_selection$cv_wvic), ylim = c(1,n_models), xlab = xlab, ylab = ylab,
     xaxt = 'n', yaxt = 'n', bty = "n", ann = FALSE)
win_dim = par("usr")
par(new = TRUE)
plot(NA, xlim = range(test_model_selection$cv_wvic), ylim = c(win_dim[3], win_dim[4] + 0.09*(win_dim[4] - win_dim[3])),
     ylab = ylab, xlab = xlab, xaxt = 'n', yaxt = 'n', bty = "n")
mtext(xlab, side = 1, line = 2.5)
mtext(ylab, side = 2, line = 0)
win_dim = par("usr")

# Add grid
grid(NULL, 11, lty = 1, col = "grey95")

# Add title
x_vec = c(win_dim[1], win_dim[2], win_dim[2], win_dim[1])
y_vec = c(win_dim[4], win_dim[4],
          win_dim[4] - 0.09*(win_dim[4] - win_dim[3]),
          win_dim[4] - 0.09*(win_dim[4] - win_dim[3]))
polygon(x_vec, y_vec, col = "grey95", border = NA)
text(x = mean(c(win_dim[1], win_dim[2])), y = (win_dim[4] - 0.09/2*(win_dim[4] - win_dim[3])), main)

# Add axes and box
lines(x_vec[1:2], rep((win_dim[4] - 0.09*(win_dim[4] - win_dim[3])),2), col = "grey50")
box(col = "grey50")

axis(1, padj = 0.3)

y_axis = axis(2, labels = FALSE, tick = FALSE)
y_axis = 1:n_models
axis(2, padj = -0.2, at = y_axis, las = 1, labels = test_model_selection$model_name)

hues = seq(15, 375, length = n_models + 1)
couleur = hcl(h = hues, l = 65, c = 100, alpha = 1)[seq_len(n_models+1)]


for (i in 1:n_models){
  lines(c(ci_low[i], ci_high[i]), c(i,i), col = couleur[i])
  points(c(ci_low[i], ci_high[i]), c(i,i), pch = "|", col = couleur[i])
  points(test_model_selection$cv_wvic[i], i, pch = 16, cex = 1.5, col = couleur[i])
}
