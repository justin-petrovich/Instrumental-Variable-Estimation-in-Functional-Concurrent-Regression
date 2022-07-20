## Read in results
setwd("./Simulation Results/Scenario 1")
scen1.1.N100 <- get(load("./Scheme 1/mse-N100-s1scheme1.RData"))
scen1.1.N200 <- get(load("./Scheme 1/mse-N200-s1scheme1.RData"))
scen1.1.N400 <- get(load("./Scheme 1/mse-N400-s1scheme1.RData"))
scen1.1.N800 <- get(load("./Scheme 1/mse-N800-s1scheme1.RData"))

scen1.2.N100 <- get(load("./Scheme 2/mse-N100-s1scheme2.RData"))
scen1.2.N200 <- get(load("./Scheme 2/mse-N200-s1scheme2.RData"))
scen1.2.N400 <- get(load("./Scheme 2/mse-N400-s1scheme2.RData"))
scen1.2.N800 <- get(load("./Scheme 2/mse-N800-s1scheme2.RData"))

setwd("../Scenario 2")
scen2.1.N100 <- get(load("./Scheme 1/mse-N100-s2scheme1.RData"))
scen2.1.N200 <- get(load("./Scheme 1/mse-N200-s2scheme1.RData"))
scen2.1.N400 <- get(load("./Scheme 1/mse-N400-s2scheme1.RData"))
scen2.1.N800 <- get(load("./Scheme 1/mse-N800-s2scheme1.RData"))

scen2.2.N100 <- get(load("./Scheme 2/mse-N100-s2scheme2.RData"))
scen2.2.N200 <- get(load("./Scheme 2/mse-N200-s2scheme2.RData"))
scen2.2.N400 <- get(load("./Scheme 2/mse-N400-s2scheme2.RData"))
scen2.2.N800 <- get(load("./Scheme 2/mse-N800-s2scheme2.RData"))

setwd("../Scenario 3")
scen3.1.N100 <- get(load("./Scheme 1/mse-N100-s3scheme1.RData"))
scen3.1.N200 <- get(load("./Scheme 1/mse-N200-s3scheme1.RData"))
scen3.1.N400 <- get(load("./Scheme 1/mse-N400-s3scheme1.RData"))
scen3.1.N800 <- get(load("./Scheme 1/mse-N800-s3scheme1.RData"))

scen3.2.N100 <- get(load("./Scheme 2/mse-N100-s3scheme2.RData"))
scen3.2.N200 <- get(load("./Scheme 2/mse-N200-s3scheme2.RData"))
scen3.2.N400 <- get(load("./Scheme 2/mse-N400-s3scheme2.RData"))
scen3.2.N800 <- get(load("./Scheme 2/mse-N800-s3scheme2.RData"))

rm(results)

## Create matrices that form the table: Median ISE
scen1.1mat <- round(t(cbind(apply(scen1.1.N100$msedf,2,median)[1:4],
                            apply(scen1.1.N200$msedf,2,median)[1:4],
                            apply(scen1.1.N400$msedf,2,median)[1:4],
                            apply(scen1.1.N800$msedf,2,median)[1:4])),2)

scen1.2mat <- round(t(cbind(apply(scen1.2.N100$msedf,2,median)[1:4],
                            apply(scen1.2.N200$msedf,2,median)[1:4],
                            apply(scen1.2.N400$msedf,2,median)[1:4],
                            apply(scen1.2.N800$msedf,2,median)[1:4])),2)

scen2.1mat <- round(t(cbind(apply(scen2.1.N100$msedf,2,median)[1:4],
                            apply(scen2.1.N200$msedf,2,median)[1:4],
                            apply(scen2.1.N400$msedf,2,median)[1:4],
                            apply(scen2.1.N800$msedf,2,median)[1:4])),2)

scen2.2mat <- round(t(cbind(apply(scen2.2.N100$msedf,2,median)[1:4],
                            apply(scen2.2.N200$msedf,2,median)[1:4],
                            apply(scen2.2.N400$msedf,2,median)[1:4],
                            apply(scen2.2.N800$msedf,2,median)[1:4])),2)

scen3.1mat <- round(t(cbind(apply(scen3.1.N100$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.1.N200$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.1.N400$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.1.N800$msedf,2,median)[c(1,2,4,5)])),2)

scen3.2mat <- round(t(cbind(apply(scen3.2.N100$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.2.N200$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.2.N400$msedf,2,median)[c(1,2,4,5)],
                            apply(scen3.2.N800$msedf,2,median)[c(1,2,4,5)])),2)




## Create matrices that form the table: Mean ISE
scen1.1mat <- round(t(cbind(apply(scen1.1.N100$msedf,2,mean)[1:4],
                            apply(scen1.1.N200$msedf,2,mean)[1:4],
                            apply(scen1.1.N400$msedf,2,mean)[1:4],
                            apply(scen1.1.N800$msedf,2,mean)[1:4])),2)

scen1.2mat <- round(t(cbind(apply(scen1.2.N100$msedf,2,mean)[1:4],
                            apply(scen1.2.N200$msedf,2,mean)[1:4],
                            apply(scen1.2.N400$msedf,2,mean)[1:4],
                            apply(scen1.2.N800$msedf,2,mean)[1:4])),2)

scen2.1mat <- round(t(cbind(apply(scen2.1.N100$msedf,2,mean)[1:4],
                            apply(scen2.1.N200$msedf,2,mean)[1:4],
                            apply(scen2.1.N400$msedf,2,mean)[1:4],
                            apply(scen2.1.N800$msedf,2,mean)[1:4])),2)

scen2.2mat <- round(t(cbind(apply(scen2.2.N100$msedf,2,mean)[1:4],
                            apply(scen2.2.N200$msedf,2,mean)[1:4],
                            apply(scen2.2.N400$msedf,2,mean)[1:4],
                            apply(scen2.2.N800$msedf,2,mean)[1:4])),2)

scen3.1mat <- round(t(cbind(apply(scen3.1.N100$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.1.N200$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.1.N400$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.1.N800$msedf,2,mean)[c(1,2,4,5)])),2)

scen3.2mat <- round(t(cbind(apply(scen3.2.N100$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.2.N200$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.2.N400$msedf,2,mean)[c(1,2,4,5)],
                            apply(scen3.2.N800$msedf,2,mean)[c(1,2,4,5)])),2)

