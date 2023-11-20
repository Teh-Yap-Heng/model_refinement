##### Load libraries #####
library(readxl)
library(writexl)
library(dplyr)
library(rstudioapi)
library(zoo)
library(lmtest)
library(car)
library(olsrr)
library(tseries)
library(fUnitRoots)
library(urca)
library(DescTools)

##### Clear workspace #####
rm(list = ls()) 

##### Import excel files #####
merged_ODR <- read_excel("mergedData.xlsx", sheet = "Treasury" , col_names = TRUE) 
MEV_ori_fixed <- read_excel("mergedData.xlsx", sheet = "Treasury" , col_names = TRUE) 
Relationship <- read_excel("Relationship.xlsx", sheet = "Input for R" , col_names = TRUE)

##### Extract data #####
ODR_logit <- as.matrix(merged_ODR[ncol(merged_ODR)])
MEV_ori_fixed <- as.matrix(MEV_ori_fixed[,2:119])

##### Regression options ####
MEVs <- MEV_ori_fixed[,c(1:118)]
logitODR <- ODR_logit

##### Find all combos ######
combo2 <- t(combn(t(colnames(MEVs)), 2, FUN = NULL, simplify = TRUE))  # 2 MEV combos

#####Find Correlation####
MatCorr <- matrix(,nrow = ncol(MEVs), ncol=2)
colnames(MatCorr) <- c("MEV", "Correlation with ODR")
for(ia in 1: ncol(MEVs)){
  MEVCorr <- MEVs[,ia]
  MEVName <- c(colnames(MEV_ori_fixed))
  Corr <- cor.test(logitODR,MEVCorr)
  MatCorr[ia,] <- c(MEVName[ia],as.numeric(Corr$estimate))
}

##### Multivariate regression #####
## Preallocate memory for results matrix
rMat1 <- matrix(, nrow = ncol(MEV_ori_fixed)+nrow(combo2), ncol = 2)
colnames(rMat1) <- c("MEV1", "MEV2")
rMat2 <- matrix(, nrow = ncol(MEV_ori_fixed)+nrow(combo2), ncol = 7)
colnames(rMat2) <- c("Intercept", "MEV1 Coeff", "MEV2 Coeff", "AdjRsq","Intercept P-Value","MEV1 P-Value","MEV2 P-Value")
rMat3 <- matrix(, nrow = ncol(MEV_ori_fixed)+nrow(combo2), ncol = 10)
colnames(rMat3) <- c("MSE","T-test","Correlation","Stationarity (MEV1) P-Value","Stationarity (MEV2) P-Value","Stationarity (Residuals) P-Value","VIF","A-D test P-Value","B-P test P-Value","D-W test Statistics")
# Note: We have to use 2 matrices, one for the characters & one for the numbers - they're combined at the end.
# Update: Now we use 3, third one is for the test results

## Create progress bar for loop
progBar <- winProgressBar(title = "Multivariate Regression Progress", min = 0, max = ncol(MEV_ori_fixed)+nrow(combo2), width = 300)

## Start looping through single MEVs
for (iw in 1:ncol(MEV_ori_fixed)){
  
  ## select required MEVs
  MEVName <- c(colnames(MEV_ori_fixed))
  mvMEV1 <- MEVs[,iw]
  
  ## Run MV regression
  regFits <- lm(logitODR ~ mvMEV1) 
  
  ##### Extract Performance Data #####
  coeffs <- t(summary(regFits)$coefficients[,1]) # extract coeffs
  adjrsq <- summary(regFits)$adj.r.squared       # extract rsq
  pvalue <- t(summary(regFits)$coefficients[,4]) # extract pvalue
  if (length(coeffs)==2) {                       # Correct for all-zero MEVs
    coeffs <- c(coeffs,0)
    pvalue <- c(pvalue,0)
  }
  
  ##### Run Stats Tests #####
  #MAPEOutput <- MAPE(regFits)
  MSEOutput <- mean(regFits$residuals ^ 2)
  ttestOutput <-  t.test(regFits$residuals)
  adfOutput1 <- adfTest(mvMEV1,type = "nc")
  adfOutput <- adfTest(summary(regFits)$residuals,type = "nc")
  
  normOutput <- ols_test_normality(regFits)
  bpOutput <- bptest(regFits)
  dwOutput <- dwtest(regFits)
  #MEV1Cor <- cor.test(logitODR,mvMEV1)
  
  # Compile results matrices
  rMat1[iw,] <- c(MEVName[iw],"-")                 # collect results
  rMat2[iw,] <- c(coeffs,adjrsq,pvalue)            # collect results
  #rMat3[iw,] <- c(MAPEOutput,MSEOutput,ttestOutput$p.value,adfOutput@test$p.value, 0, normOutput$anderson$p.value, bpOutput$p.value, dwOutput$p.value, MEV1Cor$estimate,1) 
  rMat3[iw,] <- c(MSEOutput, ttestOutput$p.value, 1, adfOutput1@test$p.value, 0, adfOutput@test$p.value, 0, normOutput$anderson$p.value, bpOutput$p.value, dwOutput$statistic) 
  
  ## Update progress bar
  setWinProgressBar(progBar, iw, title=paste(round(iw/(ncol(MEV_ori_fixed)+nrow(combo2))*100, 0),"% done"))
}



## Start looping through all the combos
for (ix in 1:nrow(combo2)){
  
  ## select required MEVs
  mvMEV1 <- MEVs[,combo2[ix,1]]
  mvMEV2 <- MEVs[,combo2[ix,2]]
  
  ## Run MV regression
  regFits <- lm(logitODR ~ mvMEV1 + mvMEV2) 
  
  ##### Extract Performance Data #####
  coeffs <- t(summary(regFits)$coefficients[,1]) # extract coeffs
  adjrsq <- summary(regFits)$adj.r.squared       # extract rsq
  pvalue <- t(summary(regFits)$coefficients[,4]) # extract pvalue
  if (length(coeffs)==2) {                       # Correct for all-zero MEVs
    coeffs <- c(coeffs,0)
    pvalue <- c(pvalue,0)
  }
  
  ##### Run Stats Tests #####
  #MAPEOutput <- MAPE(regFits)
  MSEOutput <- mean(regFits$residuals ^ 2)
  ttestOutput <-  t.test(regFits$residuals)
  Corr <- cor(data.frame(logitODR,mvMEV1,mvMEV2),use="all.obs",method="pearson")
  adfOutput1 <- adfTest(mvMEV1,type = "nc")
  adfOutput2 <- adfTest(mvMEV2,type = "nc")
  adfOutput <- adfTest(summary(regFits)$residuals,type = "nc")
  vifOutput <- ols_vif_tol(regFits)$VIF[1]
  normOutput <- ols_test_normality(regFits)
  bpOutput <- bptest(regFits)
  dwOutput <- dwtest(regFits)
  #MEV1Cor <- cor.test(logitODR,mvMEV1)
  #MEV2Cor <- cor.test(logitODR,mvMEV2)
  
  # Compile results matrices
  rMat1[ix+ncol(MEV_ori_fixed),] <- cbind(combo2[ix,1] , combo2[ix,2]) # collect results
  rMat2[ix+ncol(MEV_ori_fixed),] <- c(coeffs,adjrsq,pvalue)            # collect results
  #rMat3[ix+ncol(MEV_ori_fixed),] <- c(MAPEOutput,MSEOutput,ttestOutput$p.value,adfOutput@test$p.value ,vifOutput[1], normOutput$anderson$p.value, bpOutput$p.value, dwOutput$p.value, MEV1Cor$estimate, MEV2Cor$estimate) 
  #rMat3[ix+ncol(MEV_ori_fixed),] <- c(MAPEOutput,MSEOutput,ttestOutput$p.value,adfOutput@test$p.value, normOutput$anderson$p.value, bpOutput$p.value, dwOutput$p.value, MEV1Cor$estimate, MEV2Cor$estimate)
  rMat3[ix+ncol(MEV_ori_fixed),] <- c(MSEOutput, ttestOutput$p.value, Corr[2,3], adfOutput1@test$p.value, adfOutput2@test$p.value, adfOutput@test$p.value, vifOutput[1], normOutput$anderson$p.value, bpOutput$p.value, dwOutput$statistic) 
  
  ## Update progress bar
  setWinProgressBar(progBar, ix+ncol(MEV_ori_fixed), title=paste(round((ix+ncol(MEV_ori_fixed))/(ncol(MEV_ori_fixed)+nrow(combo2))*100, 0),"% done"))
}
close(progBar)

##### Sort data by AdjRsq value #####
rMat1 <- rMat1[order(rMat2[, "AdjRsq"]), decreasing = TRUE]
rMat1 <- rMat1[(nrow(rMat1)):1,]  # flip it because it sorts ascending for some reason
rMat3 <- rMat3[order(rMat2[, "AdjRsq"]), decreasing = TRUE] 
rMat3 <- rMat3[(nrow(rMat3)):1,]  # flip it because it sorts ascending for some reason
rMat2 <- rMat2[order(rMat2[, "AdjRsq"]), decreasing = TRUE] 
rMat2 <- rMat2[(nrow(rMat2)):1,]  # flip it because it sorts ascending for some reason

##### Business intuition section #####
## Add in relationship file data
rMat1 <- left_join(as.data.frame(rMat1), Relationship, by = c("MEV1" = "VarName"), copy = FALSE)
rMat1 <- left_join(as.data.frame(rMat1), Relationship, by = c("MEV2" = "VarName"), copy = FALSE)
colnames(rMat1) <- c("MEV1", "MEV2","Intuition1","Raw1","Intuition2","Raw2")

##### Keep unsorted/unfiltered results data #####
rawResults <- cbind(rMat1,rMat2,rMat3)  # Merge results matrices 

##### Filter data by P-Value #####
pValLimit <- 0.05 # this can be changed
deleteIDs <- vector()
# Find out which IDs
for (iy in 1:nrow(rMat2)){
  if (abs(rMat2[iy,"MEV1 P-Value"]) > pValLimit || abs(rMat2[iy,"MEV2 P-Value"]) > pValLimit){
    deleteIDs <- c(deleteIDs,iy)
  }
}
rMat1 <- rMat1[-deleteIDs,]
rMat2 <- rMat2[-deleteIDs,]
rMat3 <- rMat3[-deleteIDs,]

## Loop through solutions and remove unsuitable options
deleteIDs <- vector()
existingCombo <- vector()
for (iz in 1:nrow(rMat1)){
  ## Create combination of labels
  testCombo <- c(rMat1[iz,"Raw1"],rMat1[iz,"Raw2"])
  testCombo <- testCombo[order(testCombo)]
  testCombo <- paste(testCombo[1],testCombo[2], sep = "", collapse = NULL)
  
  ## Remove if MEV1 = MEV2
  if (rMat1[iz,"Raw1"] == rMat1[iz,"Raw2"] && !is.na(rMat1[iz,"Raw2"])){
    deleteIDs <- c(deleteIDs,iz)
    
    ## Remove if trend doesnt match +- sign
  }else if ((sign(rMat2[iz,"MEV1 Coeff"]) != sign(rMat1[iz,"Intuition1"]) || sign(rMat2[iz,"MEV2 Coeff"]) != sign(rMat1[iz,"Intuition2"])) && !is.na(rMat1[iz,"Raw2"])){
    #deleteIDs <- c(deleteIDs,iz) ##Business intuition checking is turned off, will be done in the result excel file
    
    ## Remove if this pairing of MEVs is already in
  }else if ((testCombo %in% existingCombo) == TRUE){
    deleteIDs <- c(deleteIDs,iz)
    
    ## Add new pair to list of MEV pairs
  }else {
    existingCombo <- c(existingCombo,testCombo)
  }
}
rMat1 <- rMat1[-deleteIDs,]
rMat2 <- rMat2[-deleteIDs,]
rMat3 <- rMat3[-deleteIDs,]

##### Combine results matrices to output #####
results <- cbind(rMat1,rMat2,rMat3)  # Merge results matrices
results <- subset(results,select = -c(Intuition1,Raw1,Intuition2,Raw2))  # Drop extra columns

View(rawResults)
View(results)

write_xlsx(as.data.frame(rawResults) ,path='rawResults-Modified.xlsx',col_names=TRUE,format_headers = TRUE)
write_xlsx(as.data.frame(MatCorr), path='Correlation.xlsx',col_names=TRUE,format_headers = TRUE)
write_xlsx(as.data.frame(results) ,path='results-Modified.xlsx',col_names=TRUE,format_headers = TRUE)
