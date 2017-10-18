#install.packages('Metrics', repos = 'http://cran.mirror.ac.za')
#install.packages('reshape', repos = 'http://cran.mirror.ac.za')
#install.packages('VIM', repos = 'http://cran.mirror.ac.za')
#install.packages('jomo', repos = 'http://cran.mirror.ac.za')
#install.packages('mice', repos = 'http://cran.mirror.ac.za')
#install.packages('missForest', repos = 'http://cran.mirror.ac.za')
#install.packages('ggplot2', repos = 'http://cran.mirror.ac.za')
#install.packages('gridExtra', repos = 'http://cran.mirror.ac.za')
#install.packages('doRNG', repos = 'http://cran.mirror.ac.za')
#install.packages('doMC', repos = 'http://cran.mirror.ac.za')
library(Metrics)
library(reshape)
library(VIM)
library(jomo)
library(mice)
library(missForest)
library(ggplot2)
library(gridExtra)
library(data.table)
#Libraries necessary for reproducible, parallel computation.
library(doRNG)
library(doSNOW)

source('editedJomofunction.R')

#Drawing a random smaller complete dataset from the large dataset
inx <- function(db, obs, vars, n){
  db2 <- db[,vars]
  inx <- foreach(i=1:n) %dorng% {
    arrayInd(sample(nrow(db2), obs,  replace = F), dim(db2))
  }
  
  return(list(inx, db2))
}

st <- function(inx, dbNum) {
  db2 <- inx[[2]]
  return(data.frame(db2[inx[[1]][[dbNum]][,1], ]))
}

#Generating the indexes of 1000 draws of the various data types from the big table.
# The indexes are generated before ALL the analysis so that the different imputation techniques are applied to the same data,
# since each imputation technique would affect the seed differently and thus produce differing draws if the data was generated before each analysis.
indexes <- function(db, small, medium, large, n){
  mylist = list()
  mylist[["lc5_small"]] <- inx(db, small, lc5names, n)
  mylist[["lc5_medium"]] <- inx(db, medium, lc5names, n)
  mylist[["lc5_large"]] <- inx(db, large, lc5names, n)
  mylist[["mc5_small"]] <- inx(db, small, mc5names, n)
  mylist[["mc5_medium"]] <- inx(db, medium, mc5names, n)
  mylist[["mc5_large"]] <- inx(db, large, mc5names, n)
  mylist[["mc10_small"]] <- inx(db, small, mc10names, n)
  mylist[["mc10_medium"]] <- inx(db, medium, mc10names, n)
  mylist[["mc10_large"]] <- inx(db, large, mc10names, n)
  mylist[["all_small"]] <- inx(db, small, allnames, n)
  mylist[["all_medium"]] <- inx(db, medium, allnames, n)
  mylist[["all_large"]] <- inx(db, large, allnames, n)
  return(mylist)
}

#Missingness
#MCAR - randomly induces missingness in the dataset
MCAR <- function(db, prob){
  x = db
  prop.m = prob
  for(i in 2:ncol(x)) {
    y = x[,i]
    mcar   = runif(nrow(x), min=0, max=1)
    y.mcar = ifelse(mcar<prop.m, NA, y)
    x[,i] <- y.mcar
  }
  return(x)
}

#MAR - induces missingness based on whether or not the value to the right, of the value
#under consideration, is more than 2 sd's away from the mean.
MAR <- function(db, prob) {
  x = db
  db2 = x
  prop.m = prob
  for(i in ncol(x):ncol(x)) {
    y = x[,2]
    mcar   = runif(nrow(x), min=0, max=1)
    mean = mean(x[,2])
    sd = sqrt(var(x[,2]))
    mar   = y
    y.mar = ifelse((mar<(mean-2*sd)) & (mcar<prop.m), NA, y)
    x[,i] <- y.mar
  }
  for(i in 2:(ncol(x)-1)) {
    y = db2[,i+1]
    mcar   = runif(nrow(db2), min=0, max=1)
    mean = mean(db2[,i+1])
    sd = sqrt(var(db2[,i+1]))
    mar   = y
    y.mar = ifelse((mar<(mean-2*sd)) & (mcar<prop.m), NA, y)
    x[,i] <- y.mar
  }
  
  return(x)
}

#MNAR - induces missingness based on whether or not the value under consideration is more 
#than 2 sd's away from the mean.
MNAR <- function(db, prob) {
  x = db
  prop.m = prob
  for(i in 2:ncol(x)) {
    y = x[,i]
    mcar   = runif(nrow(x), min=0, max=1)
    mean = mean(x[,i])
    sd = sqrt(var(x[,i]))
    mnar   = y
    y.mnar = ifelse((mnar<(mean-2*sd)) & (mcar<prop.m), NA, y)
    x[,i] <- y.mnar
  }
  return(x)
}

#Analysis
#Calculates the p-value associated with an f statistic of a given linear regression
calcP <- function(fstat, dataReg) {
  pf(as.numeric(fstat), dataReg[[1]][1,3], dataReg[[1]][1,4], lower.tail=F)
}

genBigTable <- function(size, mu, Sigma, skew, direction){
  rawdb <- mvrnorm(size, rep(mu,length(Sigma)), Sigma)
  if(skew == T){
    pvars <- pnorm(rawdb, mean = mu)
    if(direction == "right"){
      skwbt <-  matrix(qsn(pvars, mu, 2, 5), nrow = nrow(rawdb), ncol = ncol(rawdb))
    }
    else if(direction == "left"){
      skwbt <-  matrix(qsn(pvars, mu, -2, 5), nrow = nrow(rawdb), ncol = ncol(rawdb))
    }
    else {
      stop("Error! Skew direction must be 'left' or 'right'!")
    }
    colnames(skwbt) <- colnames(rawdb)
    return(skwbt)
  }
  else{
    return(rawdb)
  }
}

#MIXED DATA
#Using MASS to generate binomial vars from the bt table generated from Sigma
genMixedData <- function(baseData, dep_var_discr, no_of_discr_vars) {
  bt2 <- c()
  if(dep_var_discr == T){
    #Dicrete Var generation
    pvars <- pnorm(baseData[,1], mean = 4)
    binomvar <- qbinom(1-pvars, 1, 0.39)
    dep_var <- binomvar
  }
  else {
    dep_var <- baseData[,1]
  }
  x <- baseData[,-c(1)]
  if(no_of_discr_vars != 0){
    for(i in 1:no_of_discr_vars){
      #Dicrete Var generation
      pvars <- pnorm(x[,i], mean = 4)
      binomvar <- qbinom(1-pvars, 1, 0.39)
      discr_var <- binomvar
      #Array of discrete vars
      bt2 <- cbind(bt2, discr_var)
    }
    bt2 <- cbind(dep_var, bt2, baseData[,(no_of_discr_vars+2):ncol(baseData)])
    colnames(bt2) <- names(baseData[1:ncol(baseData)])
  }
  else{
    bt2 <- cbind(dep_var, x)
    colnames(bt2) <- names(baseData[1:ncol(baseData)])
  }
  return(data.frame(bt2))
}

#Requirement by the 'jomo' package that all discrete vars be presented as factors
factorize <- function(db){
  db2 <- db
  for(i in 1:ncol(db2)){
    if(length(unique(db2[,i])) <= 3){
      db2[,i] <- factor(db2[,i])
    }
  }
  return(db2)
}

checkMethod <- function(db){
  db2 <- db
  meth <- c()
  for(i in 1:ncol(db2)){
    if(length(unique(db2[,i])) <= 3){
      x <- "logreg"
    }
    else{
      x <- "pmm"
    }
    meth <- c(meth, x)
  }
  return(meth)
}

#Anova between estimates of the parameters of a linear regression
coeffANOVA <- function(db) {
  y <- c()
  for(i in 3:ncol(db)) {
    x <- data.frame(unclass(summary(aov(data = db, formula = db[,i] ~ db[,2]))))
    rownames(x) <- c(names(db[i]), paste(names(db[i]), "resid", sep = "-"))
    y <- rbind(y, x)
  }
  return(y)
}

#Function for the linear regression model. Returns the Rsquared value, f statistic and coefficients
regAnalysis <- function(db){
  params <- c()
  fit <- lm(data = db[,-1], db[,1]~.)
  x <-  summary(fit)$r.squared
  f <- summary(fit)$fstatistic
  c <- summary(fit)$coefficients
  colnames(c) <- c("Estimate", "Std Error", "t-value", "Pr")
  params <-  cbind(x, f[1], f[2], f[3])
  return(list(params, c))
}

#Mixed Analysis
#Function for the generalised linear model. Returns the null deviance, residual deviance, null df and residual df
logAnalysis <- function(db){
  params <- c()
  fit <- glm(data = db[,-1], db[,1]~., family = binomial)
  c <- summary(fit)$coefficients
  colnames(c) <- c("Estimate", "Std Error", "z-value", "Pr")
  params <-  cbind(fit$null.deviance, fit$df.null, fit$deviance, fit$df.residual)
  return(list(params, c))
}

#Calls regAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions
#It does this n times and averages the values returned for each respective version of the complete tabel
#Also returns the ANOVA between the complete table and its MCAR, MAR and MNAR versions and the ANOVA
#for the coefficients (see "coeffANOVA"). Also returns values for density plot
analysis <- function(db, n, inx){
  all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
    d <- st(inx, i)
    regCOM <- regAnalysis(d)
    regMCAR <- regAnalysis(MCAR(d, 0.05))
    regMAR <- regAnalysis(MAR(d, 0.75))
    regMNAR <- regAnalysis(MNAR(d, 0.75))
    x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2])
    y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "MCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMCAR[[2]][2:nrow(regMCAR[[2]]),1:3]), "MAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMAR[[2]][2:nrow(regMAR[[2]]),1:3]), "MNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMNAR[[2]][2:nrow(regMNAR[[2]]),1:3]))
    cbind(x,y)
  }
  vars.rs <- rbind(matrix(rep("COMrs", 1000), ncol = 1), matrix(rep("MCARrs", 1000), ncol = 1), matrix(rep("MARrs", 1000), ncol = 1), matrix(rep("MNARrs", 1000), ncol = 1))
  tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7)]), c(1,3,5,7)])[,2])
  vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1))
  tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8)]), c(2,4,6,8)])[,2])
  rsq.aov <- data.frame(unclass(summary(aov(data = tab1.rs, formula = tab1.rs$value ~ tab1.rs$ID))))
  stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
  sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
  stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
  reg <- regAnalysis(st(inx, 1))
  pvalue <- c()
  rsquared <- c()
  sd.rs <- c()
  for(i in 1:4) {
    p <- calcP(stats.f[i,2], reg)
    pvalue <- c(pvalue, p)
    r <- round(stats.rs[i,2], 5)
    rsquared <- c(rsquared, r)
    s <- round(sd.stats.rs[i,2],5)
    sd.rs <- c(sd.rs, s)
  }
  #Rsq(COM)-Rsq(Imp)
  rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[2001:3000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[1001:2000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[3001:4000,2])
  rsq <- matrix(c(rbind("COM", "MAR", "MCAR", "MNAR"), rsquared, sd.rs, pvalue), nrow = 4, ncol = 4)
  lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,10] ~ all[,9], all[9:12], mean)[,1], "MCAR" = aggregate(all[,10] ~ all[,9], all[9:12], mean)[,2], "MAR" = aggregate(all[,11] ~ all[,9], all[9:12], mean)[,2], "MNAR" = aggregate(all[,12] ~ all[,9], all[9:12], mean)[,2])[,2:4])
  colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
  return(list(rsq, rsq.aov, tab1.rs, rsqDiff, lin.coefficients))
}

#Calls logAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions
#It does this n times and averages the values returned for each respective version of the complete tabel
#Also returns the ANOVA between the complete table and its MCAR, MAR and MNAR versions. Also returns values for density plot
mixedAnalysis <- function(db, n, inx, dep_var_discr, no_of_discr_vars){
  if(dep_var_discr == T){
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      logCOM <- logAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      logMCAR <- logAnalysis(genMixedData(MCAR(d, 0.05), dep_var_discr, no_of_discr_vars))
      logMAR <- logAnalysis(genMixedData(MAR(d, 0.75), dep_var_discr, no_of_discr_vars))
      logMNAR <- logAnalysis(genMixedData(MNAR(d, 0.75), dep_var_discr, no_of_discr_vars))
      x <- cbind.data.frame(logCOM[[1]][1,1],logCOM[[1]][1,2], logCOM[[1]][1,3],logCOM[[1]][1,4], logMCAR[[1]][1,1],logMCAR[[1]][1,2], logMCAR[[1]][1,3],logMCAR[[1]][1,4], logMAR[[1]][1,1],logMAR[[1]][1,2], logMAR[[1]][1,3],logMAR[[1]][1,4], logMNAR[[1]][1,1],logMNAR[[1]][1,2], logMNAR[[1]][1,3],logMNAR[[1]][1,4])
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "z-value difference"), "MCAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logMCAR[[2]][2:nrow(logMCAR[[2]]),1:3]), "MAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logMAR[[2]][2:nrow(logMAR[[2]]),1:3]), "MNAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logMNAR[[2]][2:nrow(logMNAR[[2]]),1:3]))
      cbind(x,y)
    }
    vars.dev <- rbind(matrix(rep("COMdev", 1000), ncol = 1), matrix(rep("MCARdev", 1000), ncol = 1), matrix(rep("MARdev", 1000), ncol = 1), matrix(rep("MNARdev", 1000), ncol = 1))
    tab1.dev <- data.frame("ID" = vars.dev, "value" = melt(all[!duplicated(all[,c(3,7,11,15)]), c(3,7,11,15)])[,2])
    tempTab <- data.frame("COMdevdif" = (all[,1]-all[,3]), "MCARdevdif" = (all[,5]-all[,7]), "MARdevdif" = (all[,9]-all[,11]), "MNARdevdif" = (all[,13]-all[,15]))
    vars.p <- rbind(matrix(rep("COMp", 1000), ncol = 1), matrix(rep("MCARp", 1000), ncol = 1), matrix(rep("MARp", 1000), ncol = 1), matrix(rep("MNARp", 1000), ncol = 1))
    tab1.p <- data.frame("ID" = vars.p, "value" = melt(tempTab[!duplicated(tempTab), ])[,2])
    resid.dev.aov <- data.frame(unclass(summary(aov(data = tab1.dev, formula = tab1.dev$value ~ tab1.dev$ID))))
    mean.resid.dev <- aggregate(value ~ ID, tab1.dev, mean, na.action = na.pass)
    sd.resid.dev <- aggregate(value ~ ID, tab1.dev, sd, na.action = na.pass)
    stats.p <- aggregate(value ~ ID, tab1.p, mean, na.action = na.pass)
    logAna <- logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,2] - logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,4]
    pvalue <- c()
    mean.resid <- c()
    sd.resid <- c()
    for(i in 1:4) {
      p <- 1-pchisq(stats.p[i,2], logAna)
      pvalue <- c(pvalue, p)
      r <- round(mean.resid.dev[i,2], 5)
      mean.resid <- c(mean.resid, r)
      s <- round(sd.resid.dev[i,2], 5)
      sd.resid <- c(sd.resid, s)
    }
    #Rsq(COM)-Rsq(Imp)
    residDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[2001:3000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[1001:2000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[3001:4000,2])
    residual.deviance <- matrix(c(rbind("COM", "MAR", "MCAR", "MNAR"), mean.resid, sd.resid, pvalue), nrow = 4, ncol = 4)
    log.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,18] ~ all[,17], all[17:20], mean)[,1], "MCAR" = aggregate(all[,18] ~ all[,17], all[17:20], mean)[,2], "MAR" = aggregate(all[,19] ~ all[,17], all[17:20], mean)[,2], "MNAR" = aggregate(all[,20] ~ all[,17], all[17:20], mean)[,2])[,2:4])
    colnames(log.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    return(list(residual.deviance, resid.dev.aov, tab1.dev, residDiff, log.coefficients))
  }
  else{
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      regCOM <- regAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      regMCAR <- regAnalysis(genMixedData(MCAR(d, 0.05), dep_var_discr, no_of_discr_vars))
      regMAR <- regAnalysis(genMixedData(MAR(d, 0.75), dep_var_discr, no_of_discr_vars))
      regMNAR <- regAnalysis(genMixedData(MNAR(d, 0.75), dep_var_discr, no_of_discr_vars))
      x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2])
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "MCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMCAR[[2]][2:nrow(regMCAR[[2]]),1:3]), "MAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMAR[[2]][2:nrow(regMAR[[2]]),1:3]), "MNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regMNAR[[2]][2:nrow(regMNAR[[2]]),1:3]))
      cbind(x,y)
    }
    vars.rs <- rbind(matrix(rep("COMrs", 1000), ncol = 1), matrix(rep("MCARrs", 1000), ncol = 1), matrix(rep("MARrs", 1000), ncol = 1), matrix(rep("MNARrs", 1000), ncol = 1))
    tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7)]), c(1,3,5,7)])[,2])
    vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1))
    tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8)]), c(2,4,6,8)])[,2])
    rsq.aov <- data.frame(unclass(summary(aov(data = tab1.rs, formula = tab1.rs$value ~ tab1.rs$ID))))
    stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
    sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
    stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
    reg <- regAnalysis(st(inx, 1))
    pvalue <- c()
    rsquared <- c()
    sd.rs <- c()
    for(i in 1:4) {
      p <- calcP(stats.f[i,2], reg)
      pvalue <- c(pvalue, p)
      r <- round(stats.rs[i,2], 5)
      rsquared <- c(rsquared, r)
      s <- round(sd.stats.rs[i,2],5)
      sd.rs <- c(sd.rs, s)
    }
    #Rsq(COM)-Rsq(Imp)
    rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[2001:3000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[1001:2000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[3001:4000,2])
    rsq <- matrix(c(rbind("COM", "MAR", "MCAR", "MNAR"), rsquared, sd.rs, pvalue), nrow = 4, ncol = 4)
    lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,10] ~ all[,9], all[9:12], mean)[,1], "MCAR" = aggregate(all[,10] ~ all[,9], all[9:12], mean)[,2], "MAR" = aggregate(all[,11] ~ all[,9], all[9:12], mean)[,2], "MNAR" = aggregate(all[,12] ~ all[,9], all[9:12], mean)[,2])[,2:4])
    colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
    return(list(rsq, rsq.aov, tab1.rs, rsqDiff, lin.coefficients))
  }
}

neighbour <- function(db, k){
  imputed <- db[1]
  for(i in 2:ncol(db)){
    x <- kNN(db, variable = names(db[i]), k = k, dist_var = names(db[-i]))[i]
    imputed <- cbind(imputed, x)
  }
  return(imputed)
}

#Calls regAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#knn-imputation upon the newly generated datasets with missing values and produces imputed tables for each type. It does 
#this n times and averages the estimates calculated for each respective version of the tabel. 
knnAnalysis <- function(db, n, inx, k){
  all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
    d <- st(inx, i)
    mcar <- MCAR(d, 0.05)
    mar <- MAR(d, 0.75)
    mnar <- MNAR(d, 0.75)
    regCOM <- regAnalysis(d)
    regMCAR <- regAnalysis(mcar)
    regMAR <- regAnalysis(mar)
    regMNAR <- regAnalysis(mnar)
    ImpMCAR <- neighbour(mcar, k)
    ImpMAR <- neighbour(mar, k)
    ImpMNAR <- neighbour(mnar, k)
    regImpMCAR <- regAnalysis(ImpMCAR)
    regImpMAR <- regAnalysis(ImpMAR)
    regImpMNAR <- regAnalysis(ImpMNAR)
    rmseMCAR <- rmse(d,ImpMCAR)
    rmseMAR <- rmse(d,ImpMAR)
    rmseMNAR <- rmse(d,ImpMNAR)
    x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
    y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
    cbind(x, y)
  }
  vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
  tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
  vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
  tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
  stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
  sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
  stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
  reg <- regAnalysis(st(inx, 1))
  pvalue <- c()
  rsquared <- c()
  sd.rs <- c()
  for(i in 1:7) {
    p <- calcP(stats.f[i,2], reg)
    pvalue <- c(pvalue, p)
    r <- round(stats.rs[i,2], 5)
    rsquared <- c(rsquared, r)
    s <- round(sd.stats.rs[i,2],5)
    sd.rs <- c(sd.rs, s)
  }
  #Rsq(COM)-Rsq(Imp)
  rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
  rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
  lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
  colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
  rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
  rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
  return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
}

#Calls logAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#knn-imputation upon the newly generated datasets with missing values and produces imputed tables for each type. It does 
#this n times and averages the estimates calculated for each respective version of the tabel. 
knnMixedAnalysis <- function(db, n, inx, dep_var_discr, no_of_discr_vars, k){
  if(dep_var_discr == T){
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      ImpMCAR <- neighbour(genMixedData(mcar, dep_var_discr, no_of_discr_vars), k)
      ImpMAR <- neighbour(genMixedData(mar, dep_var_discr, no_of_discr_vars), k)
      ImpMNAR <- neighbour(genMixedData(mnar, dep_var_discr, no_of_discr_vars), k)
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      logCOM <- logAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      logMCAR <- logAnalysis(genMixedData(mcar, dep_var_discr, no_of_discr_vars))
      logMAR <- logAnalysis(genMixedData(mar, dep_var_discr, no_of_discr_vars))
      logMNAR <- logAnalysis(genMixedData(mnar, dep_var_discr, no_of_discr_vars))
      logImpMCAR <- logAnalysis(ImpMCAR)
      logImpMAR <- logAnalysis(ImpMAR)
      logImpMNAR <- logAnalysis(ImpMNAR)
      x <- cbind.data.frame(logCOM[[1]][1,1],logCOM[[1]][1,2], logCOM[[1]][1,3],logCOM[[1]][1,4], logMCAR[[1]][1,1],logMCAR[[1]][1,2], logMCAR[[1]][1,3],logMCAR[[1]][1,4], logMAR[[1]][1,1],logMAR[[1]][1,2], logMAR[[1]][1,3],logMAR[[1]][1,4], logMNAR[[1]][1,1],logMNAR[[1]][1,2], logMNAR[[1]][1,3],logMNAR[[1]][1,4], logImpMCAR[[1]][1,1],logImpMCAR[[1]][1,2], logImpMCAR[[1]][1,3],logImpMCAR[[1]][1,4], logImpMAR[[1]][1,1],logImpMAR[[1]][1,2],logImpMAR[[1]][1,3],logImpMAR[[1]][1,4],logImpMNAR[[1]][1,1],logImpMNAR[[1]][1,2],logImpMNAR[[1]][1,3],logImpMNAR[[1]][1,4], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMCAR[[2]][2:nrow(logImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMAR[[2]][2:nrow(logImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMNAR[[2]][2:nrow(logImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.dev <- rbind(matrix(rep("COMdev", 1000), ncol = 1), matrix(rep("MCARdev", 1000), ncol = 1), matrix(rep("MARdev", 1000), ncol = 1), matrix(rep("MNARdev", 1000)), matrix(rep("impMCARdev", 1000), ncol = 1), matrix(rep("impMARdev", 1000), ncol = 1), matrix(rep("impMNARdev", 1000), ncol = 1))
    tab1.dev <- data.frame("ID" = vars.dev, "value" = melt(all[!duplicated(all[,c(3,7,11,15,19,23,27)]), c(3,7,11,15,19,23,27)])[,2])
    tempTab <- data.frame("COMdevdif" = (all[,1]-all[,3]), "MCARdevdif" = (all[,5]-all[,7]), "MARdevdif" = (all[,9]-all[,11]), "MNARdevdif" = (all[,13]-all[,15]), "impMCARdevdif" = (all[,17]-all[,19]), "impMARdevdif" = (all[,21]-all[,23]), "impMNARdevdif" = (all[,25]-all[,27]))
    vars.p <- rbind(matrix(rep("COMp", 1000), ncol = 1), matrix(rep("MCARp", 1000), ncol = 1), matrix(rep("MARp", 1000), ncol = 1), matrix(rep("MNARp", 1000), ncol = 1), matrix(rep("impMCARp", 1000), ncol = 1), matrix(rep("impMARp", 1000), ncol = 1), matrix(rep("impMNARp", 1000), ncol = 1))
    tab1.p <- data.frame("ID" = vars.p, "value" = melt(tempTab[!duplicated(tempTab), ])[,2])
    mean.resid.dev <- aggregate(value ~ ID, tab1.dev, mean, na.action = na.pass)
    sd.resid.dev <- aggregate(value ~ ID, tab1.dev, sd, na.action = na.pass)
    stats.p <- aggregate(value ~ ID, tab1.p, mean, na.action = na.pass)
    logAna <- logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,2] - logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,4]
    pvalue <- c()
    mean.resid <- c()
    sd.resid <- c()
    for(i in 1:7) {
      p <- 1-pchisq(stats.p[i,2], logAna)
      pvalue <- c(pvalue, p)
      r <- round(mean.resid.dev[i,2], 5)
      mean.resid <- c(mean.resid, r)
      s <- round(sd.resid.dev[i,2], 5)
      sd.resid <- c(sd.resid, s)
    }
    #Rsq(COM)-Rsq(Imp)
    residDiff <- data.frame("diffMAR" = tab1.dev[1:1000,2]-tab1.dev[5001:6000,2], "diffMCAR" = tab1.dev[1:1000,2]-tab1.dev[4001:5000,2], "diffMNAR" = tab1.dev[1:1000,2]-tab1.dev[6001:7000,2])
    residual.deviance <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), mean.resid, sd.resid, pvalue), nrow = 7, ncol = 4)
    log.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,1], "ImpMCAR" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,2], "ImpMAR" = aggregate(all[,34] ~ all[,32], all[32:35], mean)[,2], "ImpMNAR" = aggregate(all[,35] ~ all[,32], all[32:35], mean)[,2])[,2:4])
    colnames(log.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(28,29,30)]), c(28,29,30)])[,2])
    return(list(residual.deviance, tab1.dev, residDiff, rmseAll, log.coefficients))
  }
  else{
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      ImpMCAR <- neighbour(genMixedData(mcar, dep_var_discr, no_of_discr_vars), k)
      ImpMAR <- neighbour(genMixedData(mar, dep_var_discr, no_of_discr_vars), k)
      ImpMNAR <- neighbour(genMixedData(mnar, dep_var_discr, no_of_discr_vars), k)
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      regCOM <- regAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      regMCAR <- regAnalysis(genMixedData(mcar, dep_var_discr, no_of_discr_vars))
      regMAR <- regAnalysis(genMixedData(mar, dep_var_discr, no_of_discr_vars))
      regMNAR <- regAnalysis(genMixedData(mnar, dep_var_discr, no_of_discr_vars))
      regImpMCAR <- regAnalysis(ImpMCAR)
      regImpMAR <- regAnalysis(ImpMAR)
      regImpMNAR <- regAnalysis(ImpMNAR)
      x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
    tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
    vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
    tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
    stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
    sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
    stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
    reg <- regAnalysis(st(inx, 1))
    pvalue <- c()
    rsquared <- c()
    sd.rs <- c()
    for(i in 1:7) {
      p <- calcP(stats.f[i,2], reg)
      pvalue <- c(pvalue, p)
      r <- round(stats.rs[i,2], 5)
      rsquared <- c(rsquared, r)
      s <- round(sd.stats.rs[i,2],5)
      sd.rs <- c(sd.rs, s)
    }
    #Rsq(COM)-Rsq(Imp)
    rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
    rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
    lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
    colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
    return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
  }
}

jointModel.lm <- function(db){
  attach(db)
  imputed <- jomo.lm.edited(formula = formula, data = db, output = 0, out.iter = 0, nburn = 50, nbetween = 250)
  x <- c()
  for(i in 1:nrow(imputed)){
    if(imputed[i,]$Imputation == 5){
      x <- rbind(x,imputed[i,1:ncol(db)])
    }  
  }
  imputed <- x
  return(imputed)
}

jointModel.glm <- function(db){
  attach(db)
  imputed <- jomo.glm.edited(formula = formula, data = db, output = 0, out.iter = 0, nburn = 50, nbetween = 50, family = "binomial")
  x <- c()
  for(i in 1:nrow(imputed)){
    if(imputed[i,]$Imputation == 5){
      x <- rbind(x,imputed[i,1:ncol(db)])
    }  
  }
  imputed <- x
  imputed[,1] <- imputed[,1]-1
  return(imputed)
}

#Calls regAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#knn-imputation upon the newly generated datasets with missing values and produces imputed tables for each type. It does 
#this n times and averages the estimates calculated for each respective version of the tabel. 
jomoAnalysis <- function(db, n, inx){
  all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
    d <- st(inx, i)
    mcar <- MCAR(d, 0.05)
    mar <- MAR(d, 0.75)
    mnar <- MNAR(d, 0.75)
    regCOM <- regAnalysis(d)
    regMCAR <- regAnalysis(mcar)
    regMAR <- regAnalysis(mar)
    regMNAR <- regAnalysis(mnar)
    ImpMCAR <- jointModel.lm(mcar)
    ImpMAR <- jointModel.lm(mar)
    ImpMNAR <- jointModel.lm(mnar)
    regImpMCAR <- regAnalysis(ImpMCAR)
    regImpMAR <- regAnalysis(ImpMAR)
    regImpMNAR <- regAnalysis(ImpMNAR)
    rmseMCAR <- rmse(d,ImpMCAR)
    rmseMAR <- rmse(d,ImpMAR)
    rmseMNAR <- rmse(d,ImpMNAR)
    x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
    y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
    cbind(x, y)
  }
  vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
  tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
  vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
  tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
  stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
  sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
  stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
  reg <- regAnalysis(st(inx, 1))
  pvalue <- c()
  rsquared <- c()
  sd.rs <- c()
  for(i in 1:7) {
    p <- calcP(stats.f[i,2], reg)
    pvalue <- c(pvalue, p)
    r <- round(stats.rs[i,2], 5)
    rsquared <- c(rsquared, r)
    s <- round(sd.stats.rs[i,2],5)
    sd.rs <- c(sd.rs, s)
  }
  #Rsq(COM)-Rsq(Imp)
  rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
  rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
  lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
  colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
  rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
  rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
  return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
}

#Calls logAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#knn-imputation upon the newly generated datasets with missing values and produces imputed tables for each type. It does 
#this n times and averages the estimates calculated for each respective version of the tabel. 
jomoMixedAnalysis <- function(db, n, inx, dep_var_discr, no_of_discr_vars){
  if(dep_var_discr == T){
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      d <- genMixedData(d, dep_var_discr, no_of_discr_vars)
      ImpMCAR <- jointModel.glm(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      ImpMCAR <- ImpMCAR[,c(colnames(d))]
      ImpMAR <- jointModel.glm(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      ImpMAR <- ImpMAR[,c(colnames(d))]
      ImpMNAR <- jointModel.glm(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      ImpMNAR <- ImpMNAR[,c(colnames(d))]
      logCOM <- logAnalysis(factorize(d))
      logMCAR <- logAnalysis(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      logMAR <- logAnalysis(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      logMNAR <- logAnalysis(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      logImpMCAR <- logAnalysis(ImpMCAR)
      logImpMAR <- logAnalysis(ImpMAR)
      logImpMNAR <- logAnalysis(ImpMNAR)
      for(i in 1:ncol(d)){
        ImpMCAR[,i]=as.numeric(ImpMCAR[,i])
        ImpMAR[,i]=as.numeric(ImpMAR[,i])
        ImpMNAR[,i]=as.numeric(ImpMNAR[,i])
      }
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      x <- cbind.data.frame(logCOM[[1]][1,1],logCOM[[1]][1,2], logCOM[[1]][1,3],logCOM[[1]][1,4], logMCAR[[1]][1,1],logMCAR[[1]][1,2], logMCAR[[1]][1,3],logMCAR[[1]][1,4], logMAR[[1]][1,1],logMAR[[1]][1,2], logMAR[[1]][1,3],logMAR[[1]][1,4], logMNAR[[1]][1,1],logMNAR[[1]][1,2], logMNAR[[1]][1,3],logMNAR[[1]][1,4], logImpMCAR[[1]][1,1],logImpMCAR[[1]][1,2], logImpMCAR[[1]][1,3],logImpMCAR[[1]][1,4], logImpMAR[[1]][1,1],logImpMAR[[1]][1,2],logImpMAR[[1]][1,3],logImpMAR[[1]][1,4],logImpMNAR[[1]][1,1],logImpMNAR[[1]][1,2],logImpMNAR[[1]][1,3],logImpMNAR[[1]][1,4], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "z-value difference"), "ImpMCAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMCAR[[2]][2:nrow(logImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMAR[[2]][2:nrow(logImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMNAR[[2]][2:nrow(logImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.dev <- rbind(matrix(rep("COMdev", 1000), ncol = 1), matrix(rep("MCARdev", 1000), ncol = 1), matrix(rep("MARdev", 1000), ncol = 1), matrix(rep("MNARdev", 1000)), matrix(rep("impMCARdev", 1000), ncol = 1), matrix(rep("impMARdev", 1000), ncol = 1), matrix(rep("impMNARdev", 1000), ncol = 1))
    tab1.dev <- data.frame("ID" = vars.dev, "value" = melt(all[!duplicated(all[,c(3,7,11,15,19,23,27)]), c(3,7,11,15,19,23,27)])[,2])
    tempTab <- data.frame("COMdevdif" = (all[,1]-all[,3]), "MCARdevdif" = (all[,5]-all[,7]), "MARdevdif" = (all[,9]-all[,11]), "MNARdevdif" = (all[,13]-all[,15]), "impMCARdevdif" = (all[,17]-all[,19]), "impMARdevdif" = (all[,21]-all[,23]), "impMNARdevdif" = (all[,25]-all[,27]))
    vars.p <- rbind(matrix(rep("COMp", 1000), ncol = 1), matrix(rep("MCARp", 1000), ncol = 1), matrix(rep("MARp", 1000), ncol = 1), matrix(rep("MNARp", 1000), ncol = 1), matrix(rep("impMCARp", 1000), ncol = 1), matrix(rep("impMARp", 1000), ncol = 1), matrix(rep("impMNARp", 1000), ncol = 1))
    tab1.p <- data.frame("ID" = vars.p, "value" = melt(tempTab[!duplicated(tempTab), ])[,2])
    mean.resid.dev <- aggregate(value ~ ID, tab1.dev, mean, na.action = na.pass)
    sd.resid.dev <- aggregate(value ~ ID, tab1.dev, sd, na.action = na.pass)
    stats.p <- aggregate(value ~ ID, tab1.p, mean, na.action = na.pass)
    logAna <- logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,2] - logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,4]
    pvalue <- c()
    mean.resid <- c()
    sd.resid <- c()
    for(i in 1:7) {
      p <- 1-pchisq(stats.p[i,2], logAna)
      pvalue <- c(pvalue, p)
      r <- round(mean.resid.dev[i,2], 5)
      mean.resid <- c(mean.resid, r)
      s <- round(sd.resid.dev[i,2], 5)
      sd.resid <- c(sd.resid, s)
    }
    #Resid(COM)-Resid(Imp)
    residDiff <- data.frame("diffMAR" = tab1.dev[1:1000,2]-tab1.dev[5001:6000,2], "diffMCAR" = tab1.dev[1:1000,2]-tab1.dev[4001:5000,2], "diffMNAR" = tab1.dev[1:1000,2]-tab1.dev[6001:7000,2])
    residual.deviance <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), mean.resid, sd.resid, pvalue), nrow = 7, ncol = 4)
    log.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,1], "ImpMCAR" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,2], "ImpMAR" = aggregate(all[,34] ~ all[,32], all[32:35], mean)[,2], "ImpMNAR" = aggregate(all[,35] ~ all[,32], all[32:35], mean)[,2])[,2:4])
    colnames(log.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(28,29,30)]), c(28,29,30)])[,2])
    return(list(residual.deviance, tab1.dev, residDiff, rmseAll, log.coefficients))
  }
  else{
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      d <- genMixedData(d, dep_var_discr, no_of_discr_vars)
      regCOM <- regAnalysis(factorize(d))
      regMCAR <- regAnalysis(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      regMAR <- regAnalysis(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      regMNAR <- regAnalysis(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      ImpMCAR <- jointModel.lm(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      ImpMCAR <- ImpMCAR[,c(colnames(d))]
      ImpMAR <- jointModel.lm(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      ImpMAR <- ImpMAR[,c(colnames(d))]
      ImpMNAR <- jointModel.lm(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      ImpMNAR <- ImpMNAR[,c(colnames(d))]
      regImpMCAR <- regAnalysis(ImpMCAR)
      regImpMAR <- regAnalysis(ImpMAR)
      regImpMNAR <- regAnalysis(ImpMNAR)
      for(i in 1:ncol(d)){
        ImpMCAR[,i]=as.numeric(ImpMCAR[,i])
        ImpMAR[,i]=as.numeric(ImpMAR[,i])
        ImpMNAR[,i]=as.numeric(ImpMNAR[,i])
      }
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
    tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
    vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
    tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
    stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
    sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
    stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
    reg <- regAnalysis(st(inx, 1))
    pvalue <- c()
    rsquared <- c()
    sd.rs <- c()
    for(i in 1:7) {
      p <- calcP(stats.f[i,2], reg)
      pvalue <- c(pvalue, p)
      r <- round(stats.rs[i,2], 5)
      rsquared <- c(rsquared, r)
      s <- round(sd.stats.rs[i,2],5)
      sd.rs <- c(sd.rs, s)
    }
    #Rsq(COM)-Rsq(Imp)
    rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
    rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
    lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
    colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
    return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
  }
}

#Calls regAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#mice-imputation upon the newly generated datasets with missing values and produces imputed tables for each type.
miceAnalysis <- function(db, n, inx){
  all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
    d <- st(inx, i)
    mcar <- MCAR(d, 0.05)
    mar <- MAR(d, 0.75)
    mnar <- MNAR(d, 0.75)
    regCOM <- regAnalysis(d)
    regMCAR <- regAnalysis(mcar)
    regMAR <- regAnalysis(mar)
    regMNAR <- regAnalysis(mnar)
    ImpMCAR <- complete(mice(mcar, printFlag = F))
    regImpMCAR <- regAnalysis(ImpMCAR)
    rmseMCAR <- rmse(d,ImpMCAR)
    if(sum(is.na(mar)) == 0){
      regImpMAR <- regMAR
      rmseMAR <- rmse(d,mar)
    }
    else{
      ImpMAR <- complete(mice(mar, printFlag = F))
      regImpMAR <- regAnalysis(ImpMAR)
      rmseMAR <- rmse(d,ImpMAR)
    }
    if(sum(is.na(mnar)) == 0){
      regImpMNAR <- regMNAR
      rmseMNAR <- rmse(d,mnar)
    }
    else{
      ImpMNAR <- complete(mice(mnar, printFlag = F))
      regImpMNAR <- regAnalysis(ImpMNAR)
      rmseMNAR <- rmse(d,ImpMNAR)
    }
    x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
    y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
    cbind(x,y)
  }
  vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
  tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
  vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
  tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
  stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
  sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
  stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
  reg <- regAnalysis(st(inx, 1))
  pvalue <- c()
  rsquared <- c()
  sd.rs <- c()
  for(i in 1:7) {
    p <- calcP(stats.f[i,2], reg)
    pvalue <- c(pvalue, p)
    r <- round(stats.rs[i,2], 5)
    rsquared <- c(rsquared, r)
    s <- round(sd.stats.rs[i,2],5)
    sd.rs <- c(sd.rs, s)
  }
  #Rsq(COM)-Rsq(Imp)
  rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
  rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
  lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
  colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
  rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
  rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
  return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
}

#Calls logAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#mice-imputation upon the newly generated datasets with missing values and produces imputed tables for each type.
miceMixedAnalysis <- function(db, n, inx, dep_var_discr, no_of_discr_vars){
  if(dep_var_discr == T){
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      method <- checkMethod(genMixedData(d, dep_var_discr, no_of_discr_vars))
      logCOM <- logAnalysis(factorize(genMixedData(d, dep_var_discr, no_of_discr_vars)))
      logMCAR <- logAnalysis(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      logMAR <- logAnalysis(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      logMNAR <- logAnalysis(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      ImpMCAR <- complete(mice(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
      logImpMCAR <- logAnalysis(ImpMCAR)
      for(i in 1:ncol(d)){ImpMCAR[,i]=as.numeric(ImpMCAR[,i])}
      rmseMCAR <- rmse(d,ImpMCAR)
      {if(sum(is.na(mar)) == 0){
        logImpMAR <- logMAR
        rmseMAR <- rmse(d,mar)
      }
      else{
        ImpMAR <- complete(mice(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
        logImpMAR <- logAnalysis(ImpMAR)
        for(i in 1:ncol(d)){ImpMAR[,i]=as.numeric(ImpMAR[,i])}
        rmseMAR <- rmse(d,ImpMAR)
      }}
      {if(sum(is.na(mnar)) == 0){
        logImpMNAR <- logMNAR
        rmseMNAR <- rmse(d,mnar)
      }
      else{
        ImpMNAR <- complete(mice(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
        logImpMNAR <- logAnalysis(ImpMNAR)
        for(i in 1:ncol(d)){ImpMNAR[,i]=as.numeric(ImpMNAR[,i])}
        rmseMNAR <- rmse(d,ImpMNAR)
        }}
      x <- cbind.data.frame(logCOM[[1]][1,1],logCOM[[1]][1,2], logCOM[[1]][1,3],logCOM[[1]][1,4], logMCAR[[1]][1,1],logMCAR[[1]][1,2], logMCAR[[1]][1,3],logMCAR[[1]][1,4], logMAR[[1]][1,1],logMAR[[1]][1,2], logMAR[[1]][1,3],logMAR[[1]][1,4], logMNAR[[1]][1,1],logMNAR[[1]][1,2], logMNAR[[1]][1,3],logMNAR[[1]][1,4], logImpMCAR[[1]][1,1],logImpMCAR[[1]][1,2], logImpMCAR[[1]][1,3],logImpMCAR[[1]][1,4], logImpMAR[[1]][1,1],logImpMAR[[1]][1,2],logImpMAR[[1]][1,3],logImpMAR[[1]][1,4],logImpMNAR[[1]][1,1],logImpMNAR[[1]][1,2],logImpMNAR[[1]][1,3],logImpMNAR[[1]][1,4], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMCAR[[2]][2:nrow(logImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMAR[[2]][2:nrow(logImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMNAR[[2]][2:nrow(logImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.dev <- rbind(matrix(rep("COMdev", 1000), ncol = 1), matrix(rep("MCARdev", 1000), ncol = 1), matrix(rep("MARdev", 1000), ncol = 1), matrix(rep("MNARdev", 1000)), matrix(rep("impMCARdev", 1000), ncol = 1), matrix(rep("impMARdev", 1000), ncol = 1), matrix(rep("impMNARdev", 1000), ncol = 1))
    tab1.dev <- data.frame("ID" = vars.dev, "value" = melt(all[!duplicated(all[,c(3,7,11,15,19,23,27)]), c(3,7,11,15,19,23,27)])[,2])
    tempTab <- data.frame("COMdevdif" = (all[,1]-all[,3]), "MCARdevdif" = (all[,5]-all[,7]), "MARdevdif" = (all[,9]-all[,11]), "MNARdevdif" = (all[,13]-all[,15]), "impMCARdevdif" = (all[,17]-all[,19]), "impMARdevdif" = (all[,21]-all[,23]), "impMNARdevdif" = (all[,25]-all[,27]))
    vars.p <- rbind(matrix(rep("COMp", 1000), ncol = 1), matrix(rep("MCARp", 1000), ncol = 1), matrix(rep("MARp", 1000), ncol = 1), matrix(rep("MNARp", 1000), ncol = 1), matrix(rep("impMCARp", 1000), ncol = 1), matrix(rep("impMARp", 1000), ncol = 1), matrix(rep("impMNARp", 1000), ncol = 1))
    tab1.p <- data.frame("ID" = vars.p, "value" = melt(tempTab[!duplicated(tempTab), ])[,2])
    mean.resid.dev <- aggregate(value ~ ID, tab1.dev, mean, na.action = na.pass)
    sd.resid.dev <- aggregate(value ~ ID, tab1.dev, sd, na.action = na.pass)
    stats.p <- aggregate(value ~ ID, tab1.p, mean, na.action = na.pass)
    logAna <- logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,2] - logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,4]
    pvalue <- c()
    mean.resid <- c()
    sd.resid <- c()
    for(i in 1:7) {
      p <- 1-pchisq(stats.p[i,2], logAna)
      pvalue <- c(pvalue, p)
      r <- round(mean.resid.dev[i,2], 5)
      mean.resid <- c(mean.resid, r)
      s <- round(sd.resid.dev[i,2], 5)
      sd.resid <- c(sd.resid, s)
    }
    #Resid(COM)-Resid(Imp)
    residDiff <- data.frame("diffMAR" = tab1.dev[1:1000,2]-tab1.dev[5001:6000,2], "diffMCAR" = tab1.dev[1:1000,2]-tab1.dev[4001:5000,2], "diffMNAR" = tab1.dev[1:1000,2]-tab1.dev[6001:7000,2])
    residual.deviance <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), mean.resid, sd.resid, pvalue), nrow = 7, ncol = 4)
    log.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,1], "ImpMCAR" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,2], "ImpMAR" = aggregate(all[,34] ~ all[,32], all[32:35], mean)[,2], "ImpMNAR" = aggregate(all[,35] ~ all[,32], all[32:35], mean)[,2])[,2:4])
    colnames(log.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(28,29,30)]), c(28,29,30)])[,2])
    return(list(residual.deviance, tab1.dev, residDiff, rmseAll, log.coefficients))
  }
  else{
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      method <- checkMethod(genMixedData(d, dep_var_discr, no_of_discr_vars))
      regCOM <- regAnalysis(factorize(genMixedData(d, dep_var_discr, no_of_discr_vars)))
      regMCAR <- regAnalysis(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)))
      regMAR <- regAnalysis(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)))
      regMNAR <- regAnalysis(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)))
      ImpMCAR <- complete(mice(factorize(genMixedData(mcar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
      regImpMCAR <- regAnalysis(ImpMCAR)
      rmseMCAR <- rmse(d,ImpMCAR)
      if(sum(is.na(mar)) == 0){
        regImpMAR <- regMAR
        rmseMAR <- rmse(d,mar)
      }
      else{
        ImpMAR <- complete(mice(factorize(genMixedData(mar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
        regImpMAR <- regAnalysis(ImpMAR)
        for(i in 1:ncol(d)){ImpMAR[,i]=as.numeric(ImpMAR[,i])}
        rmseMAR <- rmse(d,ImpMAR)
        }
      if(sum(is.na(mnar)) == 0){
        regImpMNAR <- regMNAR
        rmseMNAR <- rmse(d,mnar)
      }
      else{
        ImpMNAR <- complete(mice(factorize(genMixedData(mnar, dep_var_discr, no_of_discr_vars)), method = method, printFlag = F))
        regImpMNAR <- regAnalysis(ImpMNAR)
        for(i in 1:ncol(d)){ImpMNAR[,i]=as.numeric(ImpMNAR[,i])}
        rmseMNAR <- rmse(d,ImpMNAR)
        }
      x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
    tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
    vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
    tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
    stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
    sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
    stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
    reg <- regAnalysis(st(inx, 1))
    pvalue <- c()
    rsquared <- c()
    sd.rs <- c()
    for(i in 1:7) {
      p <- calcP(stats.f[i,2], reg)
      pvalue <- c(pvalue, p)
      r <- round(stats.rs[i,2], 5)
      rsquared <- c(rsquared, r)
      s <- round(sd.stats.rs[i,2],5)
      sd.rs <- c(sd.rs, s)
    }
    #Rsq(COM)-Rsq(Imp)
    rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
    rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
    lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
    colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
    return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
  }
}

rfAnalysis <- function(db, n, inx){
  all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
    d <- st(inx, i)
    mcar <- MCAR(d, 0.05)
    mar <- MAR(d, 0.75)
    mnar <- MNAR(d, 0.75)
    regCOM <- regAnalysis(d)
    regMCAR <- regAnalysis(mcar)
    regMAR <- regAnalysis(mar)
    regMNAR <- regAnalysis(mnar)
    ImpMCAR <- missForest(mcar)$ximp
    ImpMAR <- missForest(mar)$ximp
    ImpMNAR <- missForest(mnar)$ximp
    regImpMCAR <- regAnalysis(ImpMCAR)
    regImpMAR <- regAnalysis(ImpMAR)
    regImpMNAR <- regAnalysis(ImpMNAR)
    rmseMCAR <- rmse(d,ImpMCAR)
    rmseMAR <- rmse(d,ImpMAR)
    rmseMNAR <- rmse(d,ImpMNAR)
    x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
    y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
    cbind(x, y)
  }
  vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
  tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
  vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
  tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
  stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
  sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
  stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
  reg <- regAnalysis(st(inx, 1))
  pvalue <- c()
  rsquared <- c()
  sd.rs <- c()
  for(i in 1:7) {
    p <- calcP(stats.f[i,2], reg)
    pvalue <- c(pvalue, p)
    r <- round(stats.rs[i,2], 5)
    rsquared <- c(rsquared, r)
    s <- round(sd.stats.rs[i,2],5)
    sd.rs <- c(sd.rs, s)
  }
  #Rsq(COM)-Rsq(Imp)
  rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
  rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
  lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
  colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
  rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
  rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
  return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
}

#Calls logAnalysis on a certain complete table as well as its MCAR, MAR and MNAR versions and the function also performs 
#random forrest imputation upon the newly generated datasets with missing values and produces imputed tables for each type. It does 
#this n times and averages the estimates calculated for each respective version of the tabel. 
rfMixedAnalysis <- function(db, n, inx, dep_var_discr, no_of_discr_vars, k){
  if(dep_var_discr == T){
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      ImpMCAR <- missForest(genMixedData(mcar, dep_var_discr, no_of_discr_vars))$ximp
      ImpMAR <- missForest(genMixedData(mar, dep_var_discr, no_of_discr_vars))$ximp
      ImpMNAR <- missForest(genMixedData(mnar, dep_var_discr, no_of_discr_vars))$ximp
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      logCOM <- logAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      logMCAR <- logAnalysis(genMixedData(mcar, dep_var_discr, no_of_discr_vars))
      logMAR <- logAnalysis(genMixedData(mar, dep_var_discr, no_of_discr_vars))
      logMNAR <- logAnalysis(genMixedData(mnar, dep_var_discr, no_of_discr_vars))
      logImpMCAR <- logAnalysis(ImpMCAR)
      logImpMAR <- logAnalysis(ImpMAR)
      logImpMNAR <- logAnalysis(ImpMNAR)
      x <- cbind.data.frame(logCOM[[1]][1,1],logCOM[[1]][1,2], logCOM[[1]][1,3],logCOM[[1]][1,4], logMCAR[[1]][1,1],logMCAR[[1]][1,2], logMCAR[[1]][1,3],logMCAR[[1]][1,4], logMAR[[1]][1,1],logMAR[[1]][1,2], logMAR[[1]][1,3],logMAR[[1]][1,4], logMNAR[[1]][1,1],logMNAR[[1]][1,2], logMNAR[[1]][1,3],logMNAR[[1]][1,4], logImpMCAR[[1]][1,1],logImpMCAR[[1]][1,2], logImpMCAR[[1]][1,3],logImpMCAR[[1]][1,4], logImpMAR[[1]][1,1],logImpMAR[[1]][1,2],logImpMAR[[1]][1,3],logImpMAR[[1]][1,4],logImpMNAR[[1]][1,1],logImpMNAR[[1]][1,2],logImpMNAR[[1]][1,3],logImpMNAR[[1]][1,4])
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "z-value difference"), "ImpMCAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMCAR[[2]][2:nrow(logImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMAR[[2]][2:nrow(logImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(logCOM[[2]][2:nrow(logCOM[[2]]),1:3] - logImpMNAR[[2]][2:nrow(logImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.dev <- rbind(matrix(rep("COMdev", 1000), ncol = 1), matrix(rep("MCARdev", 1000), ncol = 1), matrix(rep("MARdev", 1000), ncol = 1), matrix(rep("MNARdev", 1000)), matrix(rep("impMCARdev", 1000), ncol = 1), matrix(rep("impMARdev", 1000), ncol = 1), matrix(rep("impMNARdev", 1000), ncol = 1))
    tab1.dev <- data.frame("ID" = vars.dev, "value" = melt(all[!duplicated(all[,c(3,7,11,15,19,23,27)]), c(3,7,11,15,19,23,27)])[,2])
    tempTab <- data.frame("COMdevdif" = (all[,1]-all[,3]), "MCARdevdif" = (all[,5]-all[,7]), "MARdevdif" = (all[,9]-all[,11]), "MNARdevdif" = (all[,13]-all[,15]), "impMCARdevdif" = (all[,17]-all[,19]), "impMARdevdif" = (all[,21]-all[,23]), "impMNARdevdif" = (all[,25]-all[,27]))
    vars.p <- rbind(matrix(rep("COMp", 1000), ncol = 1), matrix(rep("MCARp", 1000), ncol = 1), matrix(rep("MARp", 1000), ncol = 1), matrix(rep("MNARp", 1000), ncol = 1), matrix(rep("impMCARp", 1000), ncol = 1), matrix(rep("impMARp", 1000), ncol = 1), matrix(rep("impMNARp", 1000), ncol = 1))
    tab1.p <- data.frame("ID" = vars.p, "value" = melt(tempTab[!duplicated(tempTab), ])[,2])
    mean.resid.dev <- aggregate(value ~ ID, tab1.dev, mean, na.action = na.pass)
    sd.resid.dev <- aggregate(value ~ ID, tab1.dev, sd, na.action = na.pass)
    stats.p <- aggregate(value ~ ID, tab1.p, mean, na.action = na.pass)
    logAna <- logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,2] - logAnalysis(genMixedData(st(inx, 1), dep_var_discr, no_of_discr_vars))[[1]][,4]
    pvalue <- c()
    mean.resid <- c()
    sd.resid <- c()
    for(i in 1:7) {
      p <- 1-pchisq(stats.p[i,2], logAna)
      pvalue <- c(pvalue, p)
      r <- round(mean.resid.dev[i,2], 5)
      mean.resid <- c(mean.resid, r)
      s <- round(sd.resid.dev[i,2], 5)
      sd.resid <- c(sd.resid, s)
    }
    #Rsq(COM)-Rsq(Imp)
    residDiff <- data.frame("diffMAR" = tab1.dev[1:1000,2]-tab1.dev[5001:6000,2], "diffMCAR" = tab1.dev[1:1000,2]-tab1.dev[4001:5000,2], "diffMNAR" = tab1.dev[1:1000,2]-tab1.dev[6001:7000,2])
    residual.deviance <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), mean.resid, sd.resid, pvalue), nrow = 7, ncol = 4)
    log.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,1], "ImpMCAR" = aggregate(all[,33] ~ all[,32], all[32:35], mean)[,2], "ImpMAR" = aggregate(all[,34] ~ all[,32], all[32:35], mean)[,2], "ImpMNAR" = aggregate(all[,35] ~ all[,32], all[32:35], mean)[,2])[,2:4])
    colnames(log.coefficients) <- c("Estimate difference", "Std Err difference", "z-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(28,29,30)]), c(28,29,30)])[,2])
    return(list(residual.deviance, tab1.dev, residDiff, rmseAll, log.coefficients))
  }
  else{
    all <- foreach(i=1:n, .combine = rbind.data.frame) %dorng% {
      d <- st(inx, i)
      mcar <- MCAR(d, 0.05)
      mar <- MAR(d, 0.75)
      mnar <- MNAR(d, 0.75)
      regCOM <- regAnalysis(genMixedData(d, dep_var_discr, no_of_discr_vars))
      regMCAR <- regAnalysis(genMixedData(mcar, dep_var_discr, no_of_discr_vars))
      regMAR <- regAnalysis(genMixedData(mar, dep_var_discr, no_of_discr_vars))
      regMNAR <- regAnalysis(genMixedData(mnar, dep_var_discr, no_of_discr_vars))
      ImpMCAR <- missForest(genMixedData(mcar, dep_var_discr, no_of_discr_vars))$ximp
      ImpMAR <- missForest(genMixedData(mar, dep_var_discr, no_of_discr_vars))$ximp
      ImpMNAR <- missForest(genMixedData(mnar, dep_var_discr, no_of_discr_vars))$ximp
      regImpMCAR <- regAnalysis(ImpMCAR)
      regImpMAR <- regAnalysis(ImpMAR)
      regImpMNAR <- regAnalysis(ImpMNAR)
      rmseMCAR <- rmse(d,ImpMCAR)
      rmseMAR <- rmse(d,ImpMAR)
      rmseMNAR <- rmse(d,ImpMNAR)
      x <- cbind.data.frame(regCOM[[1]][1,1],regCOM[[1]][1,2], regMCAR[[1]][1,1],regMCAR[[1]][1,2], regMAR[[1]][1,1],regMAR[[1]][1,2], regMNAR[[1]][1,1],regMNAR[[1]][1,2], regImpMCAR[[1]][1,1],regImpMCAR[[1]][1,2], regImpMAR[[1]][1,1],regImpMAR[[1]][1,2], regImpMNAR[[1]][1,1],regImpMNAR[[1]][1,2], rmseMCAR, rmseMAR, rmseMNAR)
      y <- cbind.data.frame("coefficients"=c("Estimates difference", "Std Err difference", "t-value difference"), "ImpMCAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMCAR[[2]][2:nrow(regImpMCAR[[2]]),1:3]), "ImpMAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMAR[[2]][2:nrow(regImpMAR[[2]]),1:3]), "ImpMNAR" = colSums(regCOM[[2]][2:nrow(regCOM[[2]]),1:3] - regImpMNAR[[2]][2:nrow(regImpMNAR[[2]]),1:3]))
      cbind(x, y)
    }
    vars.rs <- rbind(matrix(rep("COM", 1000), ncol = 1), matrix(rep("MCAR", 1000), ncol = 1), matrix(rep("MAR", 1000), ncol = 1), matrix(rep("MNAR", 1000), ncol = 1), matrix(rep("impMCAR", 1000), ncol = 1), matrix(rep("impMAR", 1000), ncol = 1), matrix(rep("impMNAR", 1000), ncol = 1))
    tab1.rs <- data.frame("ID" = vars.rs, "value" = melt(all[!duplicated(all[,c(1,3,5,7,9,11,13)]), c(1,3,5,7,9,11,13)])[,2])
    vars.f <- rbind(matrix(rep("COMf", 1000), ncol = 1), matrix(rep("MCARf", 1000), ncol = 1), matrix(rep("MARf", 1000), ncol = 1), matrix(rep("MNARf", 1000), ncol = 1), matrix(rep("impMCARf", 1000), ncol = 1), matrix(rep("impMARf", 1000), ncol = 1), matrix(rep("impMNARf", 1000), ncol = 1))
    tab1.f <- data.frame("ID" = vars.f, "value" = melt(all[!duplicated(all[,c(2,4,6,8,10,12,14)]), c(2,4,6,8,10,12,14)])[,2])
    stats.rs <- aggregate(value ~ ID, tab1.rs, mean, na.action = na.pass)
    sd.stats.rs <- aggregate(value ~ ID, tab1.rs, sd, na.action = na.pass)
    stats.f <- aggregate(value ~ ID, tab1.f, mean, na.action = na.pass)
    reg <- regAnalysis(st(inx, 1))
    pvalue <- c()
    rsquared <- c()
    sd.rs <- c()
    for(i in 1:7) {
      p <- calcP(stats.f[i,2], reg)
      pvalue <- c(pvalue, p)
      r <- round(stats.rs[i,2], 5)
      rsquared <- c(rsquared, r)
      s <- round(sd.stats.rs[i,2],5)
      sd.rs <- c(sd.rs, s)
    }
    #Rsq(COM)-Rsq(Imp)
    rsqDiff <- data.frame("diffMAR" = tab1.rs[1:1000,2]-tab1.rs[5001:6000,2], "diffMCAR" = tab1.rs[1:1000,2]-tab1.rs[4001:5000,2], "diffMNAR" = tab1.rs[1:1000,2]-tab1.rs[6001:7000,2])
    rsq <- matrix(c(rbind("COM", "impMAR", "impMCAR", "impNMAR", "MAR", "MCAR", "NMAR"), rsquared, sd.rs, pvalue), nrow = 7, ncol = 4)
    lin.coefficients <- t(cbind.data.frame("Coefficients" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,1], "ImpMCAR" = aggregate(all[,19] ~ all[,18], all[18:21], mean)[,2], "ImpMAR" = aggregate(all[,20] ~ all[,18], all[18:21], mean)[,2], "ImpMNAR" = aggregate(all[,21] ~ all[,18], all[18:21], mean)[,2])[,2:4])
    colnames(lin.coefficients) <- c("Estimate difference", "Std Err difference", "t-value difference")
    rmse.vars <- rbind(matrix(rep("rmseMCAR", 1000), ncol = 1), matrix(rep("rmseMAR", 1000), ncol = 1), matrix(rep("rmseMNAR", 1000), ncol = 1))
    rmseAll <- data.frame("ID" = rmse.vars, "value" = melt(all[!duplicated(all[,c(15,16,17)]), c(15,16,17)])[,2])
    return(list(rsq, tab1.rs, rsqDiff, rmseAll, lin.coefficients))
  }
}

#Produces table with results for the linear regression
resultsTable <- function(param, rowname) {
  for(i in 1:nrow(param)) {
    x <- data.frame("Dependent Var" = c(names(Sigma[colnum])), "Mean Rsquared" = param[i,2], "SD Rsquared" = param[i,3], "p-value" = param[i,4])
    rownames(x) <- paste(rowname, param[i,1], sep = "")
    table <- rbind(table, x)
  }
  return(table)
}

#Produces table with results for the generalised linear model
resultsTableMixed <- function(param, rowname) {
  for(i in 1:nrow(param)) {
    x <- data.frame("Dependent Var" = c(names(Sigma[colnum])), "Mean Residual Deviance" = param[i,2], "SD of Residual Deviance" = param[i,3], "p-value" = param[i,4])
    rownames(x) <- paste(rowname, param[i,1], sep = "")
    table <- rbind(table, x)
  }
  return(table)
}

#Produces table with results for the ANOVA between the complete table and its MCAR, MAR and MNAR versions.
resultsTable.aov <- function(param, rowname) {
  x <- cbind("table" = rep(rowname, nrow(param)), param)
  table.aov <- rbind(table.aov, x)
}

#Produces table of differences between complete and imputed estimates
resultsDiff <- function(db, name, method){
  names(db) <- c(paste(name, "-MAR"), paste(name, "-MCAR"), paste(name, "-MNAR"))
  if(is.na(method)){
    db2 <- data.frame(db)
  }
  else{
    db2 <- data.frame("DB" = rep(method, 1000), db)
  }
  return(db2)
}

#Produces 3 separate density curves for each kind of missingness together with the imputed curve for each and a curve of the complete data.
graphme <- function(db, name, method, predictors, caption){
  top <- max(density(db[1:1000,2])$y)*1.1
  right <- round(max(density(db[1:7000,2], na.rm = T)$x), 2)
  left <- round(min(density(db[1:7000,2], na.rm = T)$x), 2)
  summean <- aggregate(value ~ ID, db, mean, na.action = na.pass)
  sumsd <- aggregate(value ~ ID, db, sd, na.action = na.pass)
  sumdat <- cbind(summean, sumsd[,2], c(top*1.12, rep(top*1.08,3),rep(top*1.04,3)))
  colnames(sumdat) <- c("ID", "mean", "sd", "y")
  filename <- c(paste(method, "ImpMAR"), paste(method, "ImpMCAR"), paste(method, "ImpMNAR"))
  mycolours <- list(c("#000000", "#000099", "#0066FF"), c("#000000", "#990000", "#FF0033"), c("#000000", "#006600", "#00CC66"))
  titles <- c("MAR", "MCAR", "MNAR")
  for(i in 1:3){
    ggplot(db[c(1:1000, (1000*(i)+1):(1000*(i+1)), (1000*(i+3)+1):(1000*(i+4))),], aes(x=value, color=ID)) +
      geom_density(size = 0.8) +
      scale_x_continuous(limits = c(left, right)) +
      scale_colour_manual(values=mycolours[[i]]) +
      geom_vline(data= sumdat[c(1, i+1, i+4),], aes(xintercept=mean,  colour=ID), linetype="dashed", size= 0.5) +
      labs(title= paste(titles[i], " - Density plot of the", predictors), colour = "Dataset", x = predictors, caption = caption) +
      geom_segment(data= sumdat[c(1, i+1, i+4),], aes(y=y, yend=y, x=mean-sd, xend=mean+sd), lwd=4) +
      annotate("text", label = formatC(sumsd[1,2], format = "e", digits = 3), x = sumdat[1,2]+sumdat[1,3]*2, y = top*1.12, size = 4) +
      annotate("text", label = formatC(sumsd[i+1,2], format = "e", digits = 3), x = sumdat[i+1,2]+sumdat[i+1,3]*2, y = top*1.08, size = 4) +
      annotate("text", label = formatC(sumsd[i+4,2], format = "e", digits = 3), x = sumdat[i+4,2]-sumdat[i+4,3]*2, y = top*1.04, size = 4)
    ggsave(paste(name, " - ", filename[i], ".png"), width = 19.05, height = 12.7, units = "cm")
  }
}

x <- Sigma[colnum]
z <- Sigma[-c(colnum)]
#Finding least correlated 4 variables to dependent var 
lc5 <- head(sort(x[,1]), 4)
lc5names <- c()
for(i in lc5){
  y <- which(x==i, arr.ind = T)
  lc5names <- c(lc5names, row.names(y))
}
lc5names <- c(names(x), lc5names)
#Finding most correlated 4 variables to dependent var 
mc5 <- tail(sort(x[,1]), 5)
mc5 <- mc5[1:4]
mc5names <- c()
for(i in mc5){
  y <- which(x==i, arr.ind = T)
  mc5names <- c(mc5names, row.names(y))
}
mc5names <- c(names(x), mc5names)
#Finding most correlated 9 variables to dependent var 
mc10 <- tail(sort(x[,1]), 10)
mc10 <- mc10[1:9]
mc10names <- c()
for(i in mc10){
  y <- which(x==i, arr.ind = T)
  mc10names <- c(mc10names, row.names(y))
}
mc10names <- c(names(x), mc10names)
#all columns
allnames <- c(names(x), names(z))
