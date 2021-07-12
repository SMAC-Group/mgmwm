# bug steffan

library(imudata)
library(simts)

axisImuData <- list()

data(adis_16405_imu1)
axisImuData[length(axisImuData)+1] <- list(as.vector(adis_16405_imu1[,1]))

data(adis_16405_imu2)
axisImuData[length(axisImuData)+1] <- list(as.vector(adis_16405_imu2[,1]))

titleLabel = 'ADIS 16405 IMU1 and IMU2'
xLabel = 'Time Scale (s)'
yLabel = 'Variance ((deg/s)^2)'
axisLabel = 'Gyro_X'
expLabels = c('IMU1', 'IMU2')

axisImuObjs <- do.call(make_mimu, c(axisImuData, list(freq=100, unit='s', sensor.name=paste(titleLabel, ', ', axisLabel, sep=''), exp.name=paste(expLabels, sep =''))))

plot(axisImuObjs, xlab=xLabel, ylab=yLabel)

model <- AR1() + RW()
estimModel = mgmwm(axisImuObjs, model, CI=TRUE)




# test without specifying class
class(estimModel)
summary(estimModel)
plot(estimModel)

# test with specifying class
summary.mgmwm(estimModel)
plot.mgmwm(estimModel)
plot.mgmwm(estimModel, decomp=TRUE, ylab_mgmwm=yLabel)

