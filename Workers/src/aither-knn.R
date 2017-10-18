#!/usr/bin/env Rscript
#install.packages('MASS', repos = 'http://cran.mirror.ac.za')
#install.packages('sn', repos = 'http://cran.mirror.ac.za')
library(MASS)
library(sn)
library(snow)

rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

#Covariance matrix from an actual luminex dataset
Sigma <- read.table('covMatrix.txt')
Sigma <- abs(Sigma)
Sigma <- Sigma[-c(17),-c(17)]

#Select dependent variable
colnum <- 18
depvar <- noquote(names(Sigma[colnum]))

source('functionsNorm.R')

clus <- makeCluster(c("172.18.0.2", "172.18.0.3", "172.18.0.4"), type = "SOCK")

clusterExport(clus, "miceAnalysis")

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

#sC5
name = "sC5"
start.time <- Sys.time()
registerDoSNOW(clus)
results <- parRapply(clus, indexesList, miceAnalysis(bt, 1000, indexesList$mc5_small))
stopCluster(clus)
end.time <- Sys.time()
time.taken <- cbind(name, end.time-start.time, start.time, end.time)
table <- resultsTable(data.frame(results[[1]]), c(name))
graphme(setorder(results[[2]], ID), name, "mice", "r-squared value", "Based on 1000 samples. 100 obs, 5 highly corr vars. All vars are continuous.")
rmseResults <- aggregate(results[[4]][,2] ~ as.character(results[[4]][,1]), results[4], mean)
for(i in 1:3){rmseResults[i,1] <- paste(name, rmseResults[i,1])}
rmseResults <- data.frame("DB"=rmseResults[,1], "rmse"=rmseResults[,2])
rmseTable <- rbind(rmseTable, rmseResults)
rsqDiffTable <- resultsDiff(results[[3]], name, "mice")
time <- rbind(time, time.taken)
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

write.csv(as.data.frame(table), "mice-continuous-results.csv")
write.csv(as.data.frame(rsqDiffTable), "mice-rsqDiffTable.csv")
