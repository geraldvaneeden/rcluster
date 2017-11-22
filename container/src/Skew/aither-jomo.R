#!/usr/bin/env Rscript
# install.packages('MASS', repos = 'http://cran.mirror.ac.za')
# install.packages('sn', repos = 'http://cran.mirror.ac.za')
library(MASS)
library(sn)

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

#Covariance matrix from an actual luminex dataset
Sigma <- read.table(args[1])
Sigma <- abs(Sigma)
Sigma <- Sigma[-c(17),-c(17)]

#Select dependent variable
colnum <- 18
depvar <- noquote(names(Sigma[colnum]))

source('functionsSkew.R')

#CONTINUOUS DATA
#BIG_TABLE
set.seed(8397)
bt <- genBigTable(100000, 4, Sigma, T, "right")

#generate indexes sampled from bt for small(100), medium(1000) and large(10000) databases 1000 times
indexesList <- indexes(bt, 100, 1000, 10000, 1000)

table <- c()
rsqDiffTable <- c()
time <- c()
rmseTable <- c()
coefficientsDiff <- c()

# Set missingness model formula to be used for all analysis
print(mc5names)
formula <- as.formula(Ferritin ~ IP.10+SAA+Fibrinogen+TPA)

#sC5
name = "sC5"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc5_small)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 100 obs, 5 highly corr vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- resultsDiff(results[[3]], name, "jomo")
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#sC0
formula <- as.formula(Ferritin ~ TGF.a+IFN.g+TNF.a+CRP+PCT+IP.10+SAA+Fibrinogen+TPA)
name = "sC0"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc10_small)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 100 obs, 10 highly corr vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lC5
formula <- as.formula(Ferritin ~ IP.10+SAA+Fibrinogen+TPA)
name = "lC5"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc5_medium)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1 000 obs, 5 highly corr vars. All variables are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lC0
formula <- as.formula(Ferritin ~ TGF.a+IFN.g+TNF.a+CRP+PCT+IP.10+SAA+Fibrinogen+TPA)
name = "lC0"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc10_medium)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1 000 obs, 10 highly corr vars. All variables are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#vC5
formula <- as.formula(Ferritin ~ IP.10+SAA+Fibrinogen+TPA)
name = "vC5"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc5_large)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 10 000 obs, 5 highly corr vars. All variables are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#vC0
formula <- as.formula(Ferritin ~ TGF.a+IFN.g+TNF.a+CRP+PCT+IP.10+SAA+Fibrinogen+TPA)
name = "vC0"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$mc10_large)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 10 000 obs, 10 highly corr vars. All variables are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5c1
formula <- as.formula(Ferritin ~ IP.10+SAA+Fibrinogen+TPA)
name = "lM5c1"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 1)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 1 Binary variable.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5c2
name = "lM5c2"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 2)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 2 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5c3
name = "lM5c3"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 3 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

# Set missingness model formula to be used for all analysis
print(mc10names)
formula <- as.formula(Ferritin ~ TGF.a+IFN.g+TNF.a+CRP+PCT+IP.10+SAA+Fibrinogen+TPA)

#lM0c1
name = "lM0c1"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 1)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 1 Binary variable.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0c4
name = "lM0c4"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 4)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 4 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0c7
name = "lM0c7"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 7)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 7 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

write.csv(as.data.frame(table), "jomo-continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "jomo-rsqDiffTable.csv")
table <- c()
residDiffTable <- c()

# Set missingness model formula to be used for all analysis
print(mc5names)
formula <- as.formula(Ferritin ~ IP.10+SAA+Fibrinogen+TPA)

#Mixed data
#lM5d0
name = "lM5d0"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 0)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dep var is binary, rest are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- resultsDiff(results[[3]], name, "jomo")
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5d2
name = "lM5d2"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 2)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = binary. 2 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5d3
name = "lM5d3"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = binary. 3 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

# Set missingness model formula to be used for all analysis
print(mc10names)
formula <- as.formula(Ferritin ~ TGF.a+IFN.g+TNF.a+CRP+PCT+IP.10+SAA+Fibrinogen+TPA)

#lM0d0
name = "lM0d0"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 0)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dep var is binary, rest are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0d4
name = "lM0d4"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 4)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = binary. 4 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0d7
name = "lM0d7"
start.time <- Sys.time()
results <- jomoMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 7)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = binary. 7 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)
names(time) <- c("DB", "Time taken", "Start time", "End time")

write.csv(as.data.frame(table), "jomo-mixed-results.csv")
write.csv(as.data.frame(residDiffTable), "jomo-residDiffTable.csv")
write.csv(as.data.frame(time), "jomo-time.csv")
write.csv(as.data.frame(rmseTable), "jomo-rmseTable.csv")
write.csv(as.data.frame(coefficientsDiff), "coefficientsDiff_jomo.csv")