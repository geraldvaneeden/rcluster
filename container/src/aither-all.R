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
print(allnames)
formula <- as.formula(Ferritin ~ IL.1RA+TGF.a+IP.10+TNF.a+IFN.a2+IFN.g+VEGF+MMP.2+MMP.9+Apo.AI+Apo.CIII+Transthyretin+Comp.FH+A2M+Haptoglobin+CRP+PCT+TPA+Fibrinogen+SAA)

#knn
name = "knn"
start.time <- Sys.time()
results <- knnAnalysis(bt, 1000, indexesList$all_medium, 3)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "knn", "r-squared value", "Based on 1000 samples. 1000 obs, 21 vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- resultsDiff(results[[3]], name, NA)
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)
write.csv(as.data.frame(table), "continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "rsqDiffTable.csv")

#mice
name = "mice"
start.time <- Sys.time()
results <- miceAnalysis(bt, 1000, indexesList$all_medium)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "mice", "r-squared value", "Based on 1000 samples. 1000 obs, 21 vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)
write.csv(as.data.frame(table), "continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "rsqDiffTable.csv")

#jomo
name = "jomo"
start.time <- Sys.time()
results <- jomoAnalysis(bt, 1000, indexesList$all_medium)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "jomo", "r-squared value", "Based on 1000 samples. 1000 obs, 21 vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)
write.csv(as.data.frame(table), "continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "rsqDiffTable.csv")

#rf
name = "rf"
start.time <- Sys.time()
results <- rfAnalysis(bt, 1000, indexesList$all_medium)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "rf", "r-squared value", "Based on 1000 samples. 1000 obs, 21 vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[3]], name, NA))
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsTemp <- cbind.data.frame(results[[5]], "SD" = results[[6]][1:3])
coefficientsDiff <- rbind(coefficientsDiff, coefficientsTemp)
write.csv(as.data.frame(table), "continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "rsqDiffTable.csv")
write.csv(as.data.frame(rmseTable), "all-rmseTable.csv")
write.csv(as.data.frame(coefficientsDiff), "coefficientsDiff_all.csv")