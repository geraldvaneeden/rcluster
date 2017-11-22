#!/usr/bin/env Rscript
# install.packages('MASS', repos = 'http://cran.mirror.ac.za')
# install.packages('sn', repos = 'http://cran.mirror.ac.za')
library(MASS)
library(sn)
library(snow)
library(doMC)

rm(list=ls())

#Covariance matrix from an actual luminex dataset
Sigma <- read.table('covMatrix.txt')
Sigma <- abs(Sigma)
Sigma <- Sigma[-c(17),-c(17)]

#Select dependent variable
colnum <- 18
depvar <- noquote(names(Sigma[colnum]))

source('functionsSkew.R')

clus <- makeCluster("10.2.0.99", "10.2.0.100", "10.2.0.101", "10.2.0.102", "10.2.0.103", "10.2.0.104", "10.2.0.105", "10.2.0.106", master='10.2.2.122', type="SOCK")
clusterExport(clus, c("st", "MCAR", "MAR", "MNAR", "regAnalysis", "complete", "calcP", "mice", "mice.impute.pmm", "mice.impute.norm", "logAnalysis", "resultsDiff", "rmse", "registerDoMC", "graphme", "resultsTable", "checkMethod"))

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

#sC5
name = "sC5"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc5_small, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 100 obs, 5 highly corr vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- resultsDiff(results[[3]], name, "knn")
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#sC0
name = "sC0"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc10_small, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 100 obs, 10 highly corr vars. All vars are continuous.")
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
name = "lC5"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc5_medium, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1 000 obs, 5 highly corr vars. All variables are continuous.")
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
name = "lC0"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc10_medium, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1 000 obs, 10 highly corr vars. All variables are continuous.")
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
name = "vC5"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc5_large, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 10 000 obs, 5 highly corr vars. All variables are continuous.")
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
name = "vC0"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$mc10_large, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 10 000 obs, 10 highly corr vars. All variables are continuous.")
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
name = "lM5c1"
start.time <- Sys.time()
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 1, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 1 Binary variable.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 2, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 2 Binary variables.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 3, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = continuous. 3 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0c1
name = "lM0c1"
start.time <- Sys.time()
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 1, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 1 Binary variable.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 4, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 4 Binary variables.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, F, 7, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = continuous. 7 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

write.csv(as.data.frame(table), "knn-continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "knn-rsqDiffTable.csv")
table <- c()
residDiffTable <- c()

#Mixed data
#lM5d0
name = "lM5d0"
start.time <- Sys.time()
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 0, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dep var is binary, rest are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- resultsDiff(results[[3]], name, "knn")
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM5d2
name = "lM5d2"
start.time <- Sys.time()
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 2, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = binary. 2 Binary variables.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 3, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 5 highly corr vars. Dependent var = binary. 3 Binary variables.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
residDiffTable <- cbind(residDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)

#lM0d0
name = "lM0d0"
start.time <- Sys.time()
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 0, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dep var is binary, rest are continuous.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 4, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = binary. 4 Binary variables.")
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
results <- knnMixedAnalysis(bt, 1000, indexesList$mc10_medium, T, 7, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTableMixed(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "residual deviance", "Based on 1000 samples. 1000 obs, 10 highly corr vars. Dependent var = binary. 7 Binary variables.")
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

write.csv(as.data.frame(table), "knn-mixed-results.csv")
write.csv(as.data.frame(residDiffTable), "knn-residDiffTable.csv")
write.csv(as.data.frame(time), "knn-time.csv")
write.csv(as.data.frame(rmseTable), "knn-rmseTable.csv")
write.csv(as.data.frame(coefficientsDiff), "coefficientsDiff_knn.csv")