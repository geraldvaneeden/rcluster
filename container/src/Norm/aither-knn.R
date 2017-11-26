#!/usr/bin/env Rscript
# install.packages('MASS', repos = 'http://cran.mirror.ac.za')
# install.packages('sn', repos = 'http://cran.mirror.ac.za')
library(MASS)
library(sn)
library(doMC)

rm(list=ls())

#Covariance matrix from an actual luminex dataset
Sigma <- read.table('covMatrix.txt')
Sigma <- abs(Sigma)
Sigma <- Sigma[-c(17),-c(17)]

#Select dependent variable
colnum <- 18
depvar <- noquote(names(Sigma[colnum]))

source('functionsNorm.R')

clus <- makeCluster(c("10.2.0.47", "10.2.3.20", "10.2.4.25", "10.2.5.27", "10.2.7.32", "10.2.1.21", "10.2.2.25", "10.2.6.23", "10.2.9.24"), master='10.2.0.47', type="SOCK")
clusterExport(clus, c("st", "MCAR", "MAR", "MNAR", "regAnalysis", "complete", "calcP", "neighbour", "kNN", "knnAnalysis", "knnMixedAnalysis", "mice", "mice.impute.pmm", "mice.impute.norm", "logAnalysis", "resultsDiff", "rmse", "registerDoMC", "graphme", "resultsTable", "checkMethod", "genMixedData", "createFrames", "writeTables"))

#CONTINUOUS DATA
#BIG_TABLE
set.seed(8397)
bt <- genBigTable(100000, 4, Sigma, F, "right")

#generate indexes sampled from bt for small(100), medium(1000) and large(10000) databases 1000 times
indexesList <- indexes(bt, 100, 1000, 10000, 1000)

table <- c()
rsqDiffTable <- c()
rmseTable <- c()
time <- c()
coefficientsDiff <- c()

registerDoSNOW(clus)

#Set parallel iteration scheme
n = 1000
tbl <- matrix(1:n, 125, 8)
k <- 3
myNames <- c("sC5", "sC0", "lC5", "lC0", "vC5", "vC0")
myIndexes <- list(indexesList$mc5_small, indexesList$mc10_small, indexesList$mc5_medium, indexesList$mc10_medium, indexesList$mc5_large, indexesList$mc10_large)

for(i in 1:6){
  name = myNames[i]
  start.time <- Sys.time()
  results <- knnAnalysis(bt, tbl, myIndexes[[i]], k)
  end.time <- Sys.time()
  time.taken <- cbind(name, end.time-start.time, start.time, end.time)
  time <- rbind(time, time.taken)
  writeTables(name, results)
  write.csv(as.data.frame(results[[1]]), paste(name, "rmseTable.csv", sep = "_"))
}

myNames <- c("lM5c1", "lM5c2", "lM5c3")
binVars <- c(1, 2, 3)

for(i in 1:3){
  name = myNames[i]
  start.time <- Sys.time()
  results <- knnMixedAnalysis(bt, tbl, indexesList$mc5_medium, F, binVars[i], k)
  end.time <- Sys.time()
  time.taken <- cbind(name, end.time-start.time, start.time, end.time)
  time <- rbind(time, time.taken)
  writeTables(name, results)
  write.csv(as.data.frame(results[[1]]), paste(name, "rmseTable.csv", sep = "_"))
}

myNames <- c("lM0c1", "lM0c4", "lM0c7")
binVars <- c(1, 4, 7)

for(i in 1:3){
  name = myNames[i]
  start.time <- Sys.time()
  results <- knnMixedAnalysis(bt, tbl, indexesList$mc10_medium, F, binVars[i], k)
  end.time <- Sys.time()
  time.taken <- cbind(name, end.time-start.time, start.time, end.time)
  time <- rbind(time, time.taken)
  writeTables(name, results)
  write.csv(as.data.frame(results[[1]]), paste(name, "rmseTable.csv", sep = "_"))
}

myNames <- c("lM5d0", "lM5d2", "lM5d3")
binVars <- c(1, 2, 3)

for(i in 1:3){
  name = myNames[i]
  mixed = T
  start.time <- Sys.time()
  results <- knnMixedAnalysis(bt, tbl, indexesList$mc5_medium, T, binVars[i], k)
  end.time <- Sys.time()
  time.taken <- cbind(name, end.time-start.time, start.time, end.time)
  time <- rbind(time, time.taken)
  writeTables(name, results, mixed)
  write.csv(as.data.frame(results[[1]]), paste(name, "rmseTable.csv", sep = "_"))
}

myNames <- c("lM0d0", "lM0d4", "lM0d7")
binVars <- c(1, 4, 7)

for(i in 1:3){
  name = myNames[i]
  start.time <- Sys.time()
  results <- knnMixedAnalysis(bt, tbl, indexesList$mc10_medium, T, binVars[i], k)
  end.time <- Sys.time()
  time.taken <- cbind(name, end.time-start.time, start.time, end.time)
  time <- rbind(time, time.taken)
  writeTables(name, results)
  write.csv(as.data.frame(results[[1]]), paste(name, "rmseTable.csv", sep = "_"))
}

