#!/usr/bin/env Rscript
# install.packages('MASS', repos = 'http://cran.mirror.ac.za')
# install.packages('sn', repos = 'http://cran.mirror.ac.za')
library(MASS)
library(sn)

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
registerDoMC(4)

#CONTINUOUS DATA
#BIG_TABLE
set.seed(8397)
bt <- genBigTable(100000, 4, Sigma, F, "right")

#generate indexes sampled from bt for small(100), medium(1000) and large(10000) databases 1000 times
indexesList <- indexes(bt, 100, 1000, 10000, 1000)

table <- c()
table.aov <- c()
resultPlots <- c()
rsqDiffTable <- c()
coefficientsDiff <- c()

#sC5
name = "sC5"
results <- analysis(bt, 1000, indexesList$mc5_small)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- data.frame(results[[3]][,1:2])
rsqDiffTable <- resultsDiff(results[[4]], name, "")
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#sC_5
name = "sC_5"
results <- analysis(bt, 1000, indexesList$lc5_small)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])


#sC0
name = "sC0"
results <- analysis(bt, 1000, indexesList$mc10_small)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#sCa
name = "sCa"
results <- analysis(bt, 1000, indexesList$all_small)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

sCPlots <- resultPlots
resultPlots <- c()

#lC5
name = "lC5"
results <- analysis(bt, 1000, indexesList$mc5_medium)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lC_5
name = "lC_5"
results <- analysis(bt, 1000, indexesList$lc5_medium)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lC0
name = "lC0"
results <- analysis(bt, 1000, indexesList$mc10_medium)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lCa
name = "lCa"
results <- analysis(bt, 1000, indexesList$all_medium)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

lCPlots <- resultPlots
resultPlots <- c()

#vC5
name = "vC5"
results <- analysis(bt, 1000, indexesList$mc5_large)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#vC_5
name = "vC_5"
results <- analysis(bt, 1000, indexesList$lc5_large)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#vC0
name = "vC0"
results <- analysis(bt, 1000, indexesList$mc10_large)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#vCa
name = "vCa"
results <- analysis(bt, 1000, indexesList$all_large)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
resultPlots <- rbind(resultPlots, results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

vCPlots <- resultPlots

write.csv(as.data.frame(table), "table_continuous.csv")
write.csv(as.data.frame(table.aov), "table_aov_continuous.csv")
write.csv(as.data.frame(rsqDiffTable), "norm-rsqDiffTable.csv")

#Density curves depicting the effect of missingness on the linear regression models of continuous data
name <- c("sC5", "sC_5", "lC5", "lC_5", "vC5", "vC_5")
filename <- c("sC", "lC", "vC")
resultPlots <- list(sCPlots, lCPlots, vCPlots)
corr <- c("highly", "weakly")
x = 1
for(i in 1:3){
  png(filename = paste(filename[i], ".png"), width = 720, height = 480, units = "px")
  plots <- list()
  for(j in 1:2){
    plotdata <- resultPlots[[i]][(1+4000*(j-1)):(4000*j),]
    sumdat <- aggregate(value ~ ID, plotdata, mean, na.action = na.pass)
    p <- ggplot(plotdata, aes(x=value, color=ID)) +
      geom_density(size = 0.8) + 
      scale_colour_manual(values=c("#000000", "#000099", "#990000", "#006600")) + 
      geom_vline(data= sumdat, aes(xintercept=value,  colour=ID), linetype="dashed", size= 0.5) + 
      labs(title= paste("Density plot of the r-squared values for five", corr[j], "correlated variables and", 10**(i+1), "observations"), colour = "Dataset", x = "R-squared Value", caption = paste("based on 1000 samples of the dataset in question -", name[x]))
    plots <- list(plots, p)
    x = x + 1
  }
  grid.arrange(plots[[1]][[2]], plots[[2]], nrow = 2)
  dev.off() 
}

resultPlots <- list(sCPlots, lCPlots, vCPlots)
png(filename = paste("SizeCompare", ".png"), width = 720, height = 480, units = "px")
plotdata <- resultPlots[[2]][1:4000,]
sumdat <- aggregate(value ~ ID, plotdata, mean, na.action = na.pass)
right <- round(max(density(resultPlots[[2]][1:4000,2], na.rm = T)$x), 2)
left <- round(min(density(resultPlots[[2]][1:4000,2], na.rm = T)$x), 2)
p <- ggplot(plotdata, aes(x=value, color=ID)) + 
  geom_density(size = 0.8) + 
  scale_x_continuous(limits = c(left, right)) +
  scale_colour_manual(values=c("#000000", "#000099", "#990000", "#006600")) + 
  geom_vline(data= sumdat, aes(xintercept=value,  colour=ID), linetype="dashed", size= 0.5) + 
  labs(title= paste("Density plot of the r-squared values for five most correlated variables and 1000 observations"), colour = "Dataset", x = "R-squared Value", caption = paste("based on 1000 samples of the dataset in question -", name[x]))
plotdata <- resultPlots[[3]][1:4000,]
sumdat <- aggregate(value ~ ID, plotdata, mean, na.action = na.pass)
q <- ggplot(plotdata, aes(x=value, color=ID)) + 
  geom_density(size = 0.8) + 
  scale_x_continuous(limits = c(left, right)) +
  scale_colour_manual(values=c("#000000", "#000099", "#990000", "#006600")) + 
  geom_vline(data= sumdat, aes(xintercept=value,  colour=ID), linetype="dashed", size= 0.5) + 
  labs(title= paste("Density plot of the r-squared values for five most correlated variables and 10 000 observations"), colour = "Dataset", x = "R-squared Value", caption = paste("based on 1000 samples of the dataset in question -", name[x]))
grid.arrange(p, q, nrow = 2)
dev.off()


#Density curve depicting the effect of missingness on the data
COMcurve <- st(indexesList$mc5_large, 1)
MCARcurve <- na.omit(MCAR(st(indexesList$mc5_large, 1), 0.05))
MARcurve <-  na.omit(MAR(st(indexesList$mc5_large, 1), 0.75))
MNARcurve <-  na.omit(MNAR(st(indexesList$mc5_large, 1), 0.75))
dataCurve <- list(COMcurve[2], MCARcurve[2], MARcurve[2], MNARcurve[2])
png(filename = paste("dataCurve", ".png"), width = 720, height = 480, units = "px")
par(mar = c(4,6,4,2))
x <- max(max(density(dataCurve[[1]][,1])$y), max(density(dataCurve[[2]][,1])$y), max(density(dataCurve[[3]][,1])$y), max(density(dataCurve[[4]][,1])$y))
plot(density(dataCurve[[1]][,1]), col = 1, lwd=3, ylim = c(0, x), main = "Density Curve Representing the Effect of Different Kinds of Missingness on Data")
lines(density(dataCurve[[2]][,1]), col = 2, lwd=3)
lines(density(dataCurve[[3]][,1]), col = 3, lwd=3)
lines(density(dataCurve[[4]][,1]), col = 4, lwd=3)
abline(v = mean(dataCurve[[1]][,1]), col = 1, lty = 2, lwd=2)
abline(v = mean(dataCurve[[2]][,1]), col = 2, lty = 2, lwd=2)
abline(v = mean(dataCurve[[3]][,1]), col = 3, lty = 2, lwd=2)
abline(v = mean(dataCurve[[4]][,1]), col = 4, lty = 2, lwd=2)
legend("topleft",legend=c("COM","MCAR","MAR", "MNAR"), col=(1:4), lwd=2, lty = 1, seg.len = 1)
dev.off()  


table <- c()
table.aov <- c()
resultPlots <- c()
residDiffTable <- c()
results <- c()

#Mixed data
#lM5d0
name = "lM5d0"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 0)
table <- resultsTableMixed(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
d0Plots <- data.frame(results[[3]][,1:2])
residDiffTable <-resultsDiff(results[[4]], name, "")
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lM5d2
name = "lM5d2"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 2)
table <- resultsTableMixed(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
d2Plots <- data.frame(results[[3]][,1:2])
residDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lM5d3
name = "lM5d3"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, T, 3)
table <- resultsTableMixed(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
d3Plots <- data.frame(results[[3]][,1:2])
residDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

write.csv(as.data.frame(table), "table_mixed_d.csv")
write.csv(as.data.frame(table.aov), "table_aov_mixed_d.csv")
write.csv(as.data.frame(residDiffTable), "norm-residDiffTable.csv")

table <- c()
table.aov <- c()

#lM5c1
name = "lM5c1"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 1)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
c1Plots <- data.frame(results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])


#lM5c2
name = "lM5c2"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 2)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
c2Plots <- data.frame(results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

#lM5c3
name = "lM5c3"
results <- mixedAnalysis(bt, 1000, indexesList$mc5_medium, F, 3)
table <- resultsTable(data.frame(results[1]), c(name))
table.aov <- resultsTable.aov(data.frame(results[2]), c(name))
c3Plots <- data.frame(results[[3]][,1:2])
rsqDiffTable <- cbind(rsqDiffTable, resultsDiff(results[[4]], name, ""))
rownames(results[[5]]) <- c(paste(name, rownames(results[[5]])[1]), paste(name, rownames(results[[5]])[2]), paste(name, rownames(results[[5]])[3]))
coefficientsDiff <- rbind(coefficientsDiff, results[[5]])

write.csv(as.data.frame(table), "table_mixed_c.csv")
write.csv(as.data.frame(table.aov), "table_aov_mixed_c.csv")
write.csv(as.data.frame(rsqDiffTable), "norm-rsqDiffTable.csv")
write.csv(as.data.frame(coefficientsDiff), "norm_coefficientsDiff_gen_data.csv")

#Density curves depicting the effect of missingness on the logistic regression models of mixed data
name <- c("lM5d0", "lM5d2", "lM5d3")
resultPlots <- list(d0Plots, d2Plots, d3Plots)
for(i in 1:3){
  plotdata <- resultPlots[[i]]
  sumdat <- aggregate(value ~ ID, plotdata, mean, na.action = na.pass)
  ggplot(plotdata, aes(x=value, color=ID)) + 
    geom_density(size = 0.8) + 
    scale_colour_manual(values=c("#000000", "#000099", "#990000", "#006600")) + 
    geom_vline(data= sumdat, aes(xintercept=value,  colour=ID), linetype="dashed", size= 0.5) + 
    labs(title= paste(name[i], " - Density plot for 5 highly correlated variables and 1000 observations"), colour = "Dataset", x = "Residual Deviance")
  ggsave(paste(name[i], ".png"), width = 19.05, height = 12.7, units = "cm")
}

#Density curves depicting the effect of missingness on the linear regression models of mixed data
name <- c("lM5c1", "lM5c2", "lM5c3")
resultPlots <- list(c1Plots, c2Plots, c3Plots)
for(i in 1:3){
  plotdata <- resultPlots[[i]]
  sumdat <- aggregate(value ~ ID, plotdata, mean, na.action = na.pass)
  ggplot(plotdata, aes(x=value, color=ID)) + 
    geom_density(size = 0.8) + 
    scale_colour_manual(values=c("#000000", "#000099", "#990000", "#006600")) + 
    geom_vline(data= sumdat, aes(xintercept=value,  colour=ID), linetype="dashed", size= 0.5) + 
    labs(title= paste(name[i], " - Density plot for 5 highly correlated variables and 1000 observations"), colour = "Dataset", x = "R-squared Value")
  ggsave(paste(name[i], ".png"), width = 19.05, height = 12.7, units = "cm")
}
