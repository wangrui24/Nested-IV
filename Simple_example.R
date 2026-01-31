#########################
######################### We provide a simple synthetic data example in this file
#########################

library(SuperLearner)
library(caret)
library(MASS)

######################### 1. Define the functions

# Functions


logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}


# A function for generating data:
# parameters for effect size : beta
# parameters for compliance rate ; alpha


dat_gen <- function(n=1000,  alpha = c(0.5,0.2,0.05), beta = c(1,1,1)){
  # n = 1000
  mu <- c(0,1,-0.5)
  var <- matrix(data = c(1,0.2,-0.3,0.2,1,0.1,-0.3,0.1,1), ncol = 3, nrow = 3)
  x <- mvrnorm(n, mu, var)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  
  x1[x1>4] <- 4
  x1[x1<-4] <- -4
  
  x2[x2>4] <- 4
  x2[x2<-4] <- -4
  
  x3[x3>4] <- 4
  x3[x3<-4] <- -4
  
  x4 <- rbinom(n, 1, 0.5)
  x7 <- runif(n, min = -3, max = 3)
  # 
  # 
  x8 <- rbinom(n, 4, 0.5)
  x5 <- rbinom(n, 1, 0.5)
  x6 <- rbinom(n, 1, 0.5)
  
  # Generate G
  pg <- expit(1+0.2*x1-0.1*x2+0.3*x3)
  g <- rbinom(n,1,pg)
  pa <- expit(1+0.5*x1-x2+0.7*x3)
  pb <- expit(0.5+0.6*x1+0.3*x2+0.4*x3)
  
  # Generate Z
  z_h <- rbinom(n,1,g*pa+(1-g)*pb)
  
  z <- c()
  z[g==0 & z_h==0] <- "a0"
  z[g==0 & z_h==1] <- "a1"
  z[g==1 & z_h==0] <- "b0"
  z[g==1 & z_h==1] <- "b1"
  
  z_a0 <- c()
  z_a0[z=="a0"] <- 1
  z_a0[z!="a0"] <- 0
  
  z_a1 <- c()
  z_a1[z=="a1"] <- 1
  z_a1[z!="a1"] <- 0
  
  z_b0 <- c()
  z_b0[z=="b0"] <- 1
  z_b0[z!="b0"] <- 0
  
  z_b1 <- c()
  z_b1[z=="b1"] <- 1
  z_b1[z!="b1"] <- 0
  
  # Generate unmeasured confounder
  u <- rnorm(n, 0, 0.6)
  
  # Generate principal stratum
  par_ant <- c(0,-1,0.7)
  par_swi_1 <- c(1,2,-1)
  par_swi_2 <- c(-1,-2,0)
  par_atnt <- c(0.5,1,0.5)
  par_aco <- c(0.8, -2, 2)
  par_ntat <- c(-0.5,-1,-0.5)
  par_aat <- c(2,0,2)
  
  # set the parameter
  
  
  num_ant <- c(exp(1+x%*%par_ant+0.3*u))
  num_swi_1 <- c(exp(alpha[1]+alpha[2]*x%*%par_swi_1+alpha[3]*u))
  num_swi_2 <- c(exp(alpha[1]+alpha[2]*x%*%par_swi_2+alpha[3]*u))
  num_atnt <- c(exp(1+x%*%par_atnt+0.5*u))
  num_aco <- c(exp(1+x%*%par_aco+0.5*u))
  num_ntat <- c(exp(1+x%*%par_ntat-0.5*u))
  num_aat <- c(exp(1+x%*%par_aat-u))
  
  
  dem <- num_ant + num_swi_1 + num_swi_2 + num_atnt + num_aco + num_ntat + num_aat
  p_ant <- num_ant/dem
  p_swi_1 <- num_swi_1/dem
  p_swi_2 <- num_swi_2/dem
  p_atnt <- num_atnt/dem
  p_aco <- num_aco/dem
  p_ntat <- num_ntat/dem
  p_aat <- num_aat/dem
  
  p_matrix <- cbind(p_ant, p_swi_1, p_swi_2, p_atnt, p_aco, p_ntat, p_aat)
  
  s <- c()
  strata <- c("ANT", "SWI1", "SWI2", "ATNT", "ACO", "NTAT", "AAT")
  for (i in c(1:n)) {
    p_i <- c(p_matrix[i,])
    dat_i <- c(rmultinom(1,1,p_i))
    s[i] <- strata[which(dat_i==1)]
  }
  
  d <- c()
  d[(s=="ANT")]<-0
  d[(s=="AAT")]<-1
  d[(s=="ATNT")&(g==0)]<- 1
  d[(s=="ATNT")&(g==1)]<- 0
  
  d[(s=="NTAT")&(g==0)]<- 0
  d[(s=="NTAT")&(g==1)]<- 1
  
  d[(s=="ACO")&(z=="a0")]<- 0
  d[(s=="ACO")&(z=="a1")]<- 1
  d[(s=="ACO")&(z=="b0")]<- 0
  d[(s=="ACO")&(z=="b1")]<- 1 
  
  
  d[(s=="SWI1")&(z=="a0")]<- 0
  d[(s=="SWI1")&(z=="a1")]<- 0
  d[(s=="SWI1")&(z=="b0")]<- 0
  d[(s=="SWI1")&(z=="b1")]<- 1 
  
  d[(s=="SWI2")&(z=="a0")]<- 1
  d[(s=="SWI2")&(z=="a1")]<- 1
  d[(s=="SWI2")&(z=="b0")]<- 0
  d[(s=="SWI2")&(z=="b1")]<- 1
  
  # Generate potential outcomes of Y
  y_0 <- 1*(s=="ANT"|s=="AAT")*(1+x1+x2+x3+x4+u)+
    1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+x3+x4+u)+
    1*(s=="SWI1"|s=="SWI2")*(1+x1+x2+x3+x4+u)+
    1*(s=="ACO")*(1+x1+x2+x3+x4+u)
  
  
  
  y_1 <- 1*(s=="ANT"|s=="AAT")*(1+x1+2*x2+2*x3+x4+u)+
    1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+2*x3+x4+u)+
    1*(s=="SWI1"|s=="SWI2")*(beta[1]+beta[2]*x1+2*x2+beta[3]*x3+x4+u)+
    1*(s=="ACO")*(1+x1+2*x2+0.2*x2^2+x3+x4+u)
  
  y <- (1-d)*y_0+(d)*y_1
  
  
  dat <- data.frame(x1,x2,x3,x4,
                    x5,x6,x7,x8,
                    g,z_h,z,z_a0,z_a1,z_b0,z_b1,d,y,y_0, y_1,s)
  return(dat)
}


# A function for generating the point estimates
npest<- function(dataset,V=2)
{
  ####
  # dataset <- read.csv("dat_true2.csv")
  # V = 2
  ####
  index <- createFolds(1:nrow(dataset),V)
  influ_swi_vec <- NULL
  os_swi_vec <- NULL
  ee_swi_vec <- NULL
  for(v in 1:V)
  {
    ####
    #v=1
    ####
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    
    # Model for P(Z=a1|G=a, X=x)
    trainingset_a <- trainingset[which(trainingset$g==0),]
    sl_z_a1 <- SuperLearner(Y = trainingset_a$z_h,
                            X = trainingset_a[,c("x1","x2","x3", "x4",
                                                 "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    # Model for P(Z=b1|G=b, X=x)
    trainingset_b <- trainingset[which(trainingset$g==1),]
    sl_z_b1 <- SuperLearner(Y = trainingset_b$z_h,
                            X = trainingset_b[,c("x1","x2","x3", "x4",
                                                 "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(G=1|X=x)
    sl_z_g1 <- SuperLearner(Y = trainingset$g,
                            X = trainingset[,c("x1","x2","x3", "x4",
                                               "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # prediction for P(Z=a1|G=a,X)
    predictset_a1 <- predictset
    predictset_a1$z_h <- 1
    p_z_a1_a <- predict(sl_z_a1, newdata = predictset_a1, onlySL = T)$pred
    # prediction for P(Z=a0|G=a,X)
    p_z_a0_a <- 1 - p_z_a1_a
    
    # prediction for P(Z=b1|G=b,X)
    predictset_b1 <- predictset
    predictset_b1$z_h <- 1   
    
    p_z_b1_b <- predict(sl_z_b1, newdata = predictset_b1, onlySL = T)$pred
    
    # prediction for P(Z=b0|G=b,X)
    p_z_b0_b <- 1 - p_z_b1_b
    
    # prediction for P(G=a|X)
    p_g1 <- predict(sl_z_g1, newdata = predictset, onlySL = T)$pred
    p_g0 <- 1 - p_g1
    
    z_a0_pred <- p_g0*p_z_a0_a
    z_a1_pred <- p_g0*p_z_a1_a
    z_b0_pred <- p_g1*p_z_b0_b
    z_b1_pred <- p_g1*p_z_b1_b
    
    # z_a0_pred + z_a1_pred + z_b0_pred + z_b1_pred
    
    # Model and prediction for E(D|Z,X)
    # sl_d <- SuperLearner(Y = trainingset$d,
    #                         X = trainingset[,c("z_a1","z_b0","z_b1",
    #                                            "x1","x2","x3")],
    #                         family = binomial(),
    #                         SL.library = c("SL.glm", "SL.randomForest"),
    #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    # 
    # d_1_pred <- predict(sl_d, newdata = predictset, onlySL = T)$pred
    # # d_pred <- predictset$d*d_1_pred+(1-predictset$d)*(1-d_1_pred)
    # 
    # # Prediction for E[D|Z=a0,X]
    # pred_d_a0 <- predictset
    # pred_d_a0$z_a1 <- 0
    # pred_d_a0$z_b0 <- 0
    # pred_d_a0$z_b1 <- 0
    # d_1_pred_a0 <- predict(sl_d,newdata = pred_d_a0,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[D|Z=a1,X]
    # pred_d_a1 <- predictset
    # pred_d_a1$z_a1 <- 1
    # pred_d_a1$z_b0 <- 0
    # pred_d_a1$z_b1 <- 0
    # d_1_pred_a1 <- predict(sl_d,newdata = pred_d_a1,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[D|Z=b0,X]
    # pred_d_b0 <- predictset
    # pred_d_b0$z_a1 <- 0
    # pred_d_b0$z_b0 <- 1
    # pred_d_b0$z_b1 <- 0
    # d_1_pred_b0 <- predict(sl_d,newdata = pred_d_b0,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[D|Z=b1,X]
    # pred_d_b1 <- predictset
    # pred_d_b1$z_a1 <- 0
    # pred_d_b1$z_b0 <- 0
    # pred_d_b1$z_b1 <- 1
    # d_1_pred_b1 <- predict(sl_d,newdata = pred_d_b1,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z,X)
    # sl_y <- SuperLearner(Y = trainingset$y,
    #                      X = trainingset[,c("z_a1","z_b0","z_b1",
    #                                         "x1","x2","x3")],
    #                      family = gaussian(),
    #                      SL.library = c("SL.glm", "SL.randomForest"),
    #                      cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    # 
    # y_pred <- predict(sl_y,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(D|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_d_a0 <- SuperLearner(Y = train_a0$d,
                            X = train_a0[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a0 <- predict(sl_d_a0,newdata = predictset,onlySL = T)$pred
    
    # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
    
    # Model and prediction for E(D|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_d_a1 <- SuperLearner(Y = train_a1$d,
                            X = train_a1[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a1 <- predict(sl_d_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(D|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_d_b0 <- SuperLearner(Y = train_b0$d,
                            X = train_b0[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b0 <- predict(sl_d_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(D|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_d_b1 <- SuperLearner(Y = train_b1$d,
                            X = train_b1[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b1 <- predict(sl_d_b1,newdata = predictset,onlySL = T)$pred
    
    
    #####################################
    
    # Model and prediction for E(Y|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_y_a0 <- SuperLearner(Y = train_a0$y,
                            X = train_a0[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a0 <- predict(sl_y_a0,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_y_a1 <- SuperLearner(Y = train_a1$y,
                            X = train_a1[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a1 <- predict(sl_y_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_y_b0 <- SuperLearner(Y = train_b0$y,
                            X = train_b0[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b0 <- predict(sl_y_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_y_b1 <- SuperLearner(Y = train_b1$y,
                            X = train_b1[,c("x1","x2","x3", "x4",
                                            "x5","x6","x7", "x8")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b1 <- predict(sl_y_b1,newdata = predictset,onlySL = T)$pred
    
    
    # # Prediction for E[Y|Z=a0,X]
    # pred_y_a0 <- predictset
    # pred_y_a0$z_a1 <- 0
    # pred_y_a0$z_b0 <- 0
    # pred_y_a0$z_b1 <- 0
    # y_pred_a0 <- predict(sl_y,newdata = pred_y_a0,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[Y|Z=a1,X]
    # pred_y_a1 <- predictset
    # pred_y_a1$z_a1 <- 1
    # pred_y_a1$z_b0 <- 0
    # pred_y_a1$z_b1 <- 0
    # y_pred_a1 <- predict(sl_y,newdata = pred_y_a1,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[Y|Z=b0,X]
    # pred_y_b0 <- predictset
    # pred_y_b0$z_a1 <- 0
    # pred_y_b0$z_b0 <- 1
    # pred_y_b0$z_b1 <- 0
    # y_pred_b0 <- predict(sl_y,newdata = pred_y_b0,onlySL = T)$pred
    # 
    # 
    # # Prediction for E[Y|Z=b1,X]
    # pred_y_b1 <- predictset
    # pred_y_b1$z_a1 <- 0
    # pred_y_b1$z_b0 <- 0
    # pred_y_b1$z_b1 <- 1
    # y_pred_b1 <- predict(sl_y,newdata = pred_y_b1,onlySL = T)$pred
    
    
    # Calculation for influence function of swi
    # E[eta_b-eta_a]
    coeff_y <- 1/(mean(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0))
    var_y_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$y-y_pred_b1)
    var_y_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$y-y_pred_b0)
    var_y_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$y-y_pred_a1)
    var_y_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$y-y_pred_a0)
    var_y_5 <- y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0
    
    # mean(var_y_1)
    # mean(var_y_2)
    # mean(var_y_3)
    # mean(var_y_4)
    # 
    # sd(var_y_1)
    # sd(var_y_2)
    # sd(var_y_3)
    # sd(var_y_4)
    
    psi_p <- mean(y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0)/
      mean(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0)
    
    coeff_d <- 1/(mean(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0))
    var_d_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$d-d_pred_b1)
    var_d_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$d-d_pred_b0)
    var_d_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$d-d_pred_a1)
    var_d_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$d-d_pred_a0)
    var_d_5 <- d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0
    
    # mean(var_d_1)
    # mean(var_d_2)
    # mean(var_d_3)
    # mean(var_d_4)
    # 
    # sd(var_d_1)
    # sd(var_d_2)
    # sd(var_d_3)
    # sd(var_d_4)
    # 
    # 
    # mean(var_d_1-var_d_2-var_d_3+var_d_4)
    # Calculation for influence function of aco
    
    
    influ_swi <- coeff_y*(var_y_1-var_y_2-var_y_3+var_y_4+var_y_5)-psi_p*coeff_d*
      (var_d_1-var_d_2-var_d_3+var_d_4+var_d_5)
    
    os_swi <- psi_p + mean(influ_swi)
    
    ee_swi <- mean(var_y_1-var_y_2-var_y_3+var_y_4+var_y_5)/mean(var_d_1-var_d_2-var_d_3+var_d_4+var_d_5)
    os_swi_vec[v] <- os_swi
    ee_swi_vec[v] <- ee_swi
    influ_swi_vec <- c(influ_swi_vec, influ_swi)
  }
  os_swi_est <- mean(os_swi_vec)
  ee_swi_est <- mean(ee_swi_vec)
  # se_swi <- sd(os_swi_vec)/sqrt(length(influ_swi_vec))
  influ_swi <- mean(influ_swi_vec)
  se_influ <- sd(influ_swi_vec)/sqrt(length(influ_swi_vec))
  return(c(os_swi_est,ee_swi_est,se_influ))
}

# A function for generating the testing results
nptest <- function(dataset,V=2)
{
  ####
  # V = 2
  ####
  index <- createFolds(1:nrow(dataset),V)
  influ_beta_vec <- NULL
  beta_vec <- NULL
  influ_beta_vec_2 <- NULL
  beta_vec_2 <- NULL
  influ_beta_vec_3 <- NULL
  beta_vec_3 <- NULL
  
  
  ee_beta_vec <- NULL
  ee_beta_vec_2 <- NULL
  ee_beta_vec_3 <- NULL
  for(v in 1:V)
  {
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    
    # Model for P(Z=a1|G=a, X=x)
    trainingset_a <- trainingset[which(trainingset$g==0),]
    sl_z_a1 <- SuperLearner(Y = trainingset_a$z_h,
                            X = trainingset_a[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    # Model for P(Z=b1|G=b, X=x)
    trainingset_b <- trainingset[which(trainingset$g==1),]
    sl_z_b1 <- SuperLearner(Y = trainingset_b$z_h,
                            X = trainingset_b[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(G=1|X=x)
    sl_z_g1 <- SuperLearner(Y = trainingset$g,
                            X = trainingset[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # prediction for P(Z=a1|G=a,X)
    predictset_a1 <- predictset
    predictset_a1$z_h <- 1
    p_z_a1_a <- predict(sl_z_a1, newdata = predictset_a1, onlySL = T)$pred
    # prediction for P(Z=a0|G=a,X)
    p_z_a0_a <- 1 - p_z_a1_a
    
    # prediction for P(Z=b1|G=b,X)
    predictset_b1 <- predictset
    predictset_b1$z_h <- 1   
    
    p_z_b1_b <- predict(sl_z_b1, newdata = predictset_b1, onlySL = T)$pred
    
    # prediction for P(Z=b0|G=b,X)
    p_z_b0_b <- 1 - p_z_b1_b
    
    # prediction for P(G=a|X)
    p_g1 <- predict(sl_z_g1, newdata = predictset, onlySL = T)$pred
    p_g0 <- 1 - p_g1
    
    z_a0_pred <- p_g0*p_z_a0_a
    z_a1_pred <- p_g0*p_z_a1_a
    z_b0_pred <- p_g1*p_z_b0_b
    z_b1_pred <- p_g1*p_z_b1_b
    
    
    # Model and prediction for E(D|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_d_a0 <- SuperLearner(Y = train_a0$d,
                            X = train_a0[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a0 <- predict(sl_d_a0,newdata = predictset,onlySL = T)$pred
    
    # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
    
    # Model and prediction for E(D|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_d_a1 <- SuperLearner(Y = train_a1$d,
                            X = train_a1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a1 <- predict(sl_d_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(D|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_d_b0 <- SuperLearner(Y = train_b0$d,
                            X = train_b0[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b0 <- predict(sl_d_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(D|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_d_b1 <- SuperLearner(Y = train_b1$d,
                            X = train_b1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b1 <- predict(sl_d_b1,newdata = predictset,onlySL = T)$pred
    
    
    #####################################
    
    # Model and prediction for E(Y|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_y_a0 <- SuperLearner(Y = train_a0$y,
                            X = train_a0[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a0 <- predict(sl_y_a0,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_y_a1 <- SuperLearner(Y = train_a1$y,
                            X = train_a1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a1 <- predict(sl_y_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_y_b0 <- SuperLearner(Y = train_b0$y,
                            X = train_b0[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b0 <- predict(sl_y_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_y_b1 <- SuperLearner(Y = train_b1$y,
                            X = train_b1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b1 <- predict(sl_y_b1,newdata = predictset,onlySL = T)$pred
    
    # E[eta_b-eta_a]
    # coeff_y <- 1/(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0)
    var_y_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$y-y_pred_b1)
    var_y_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$y-y_pred_b0)
    var_y_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$y-y_pred_a1)
    var_y_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$y-y_pred_a0)
    # var_y_5 <- y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0
    
    # denominator of swate
    dem_swate <- (d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0) #+0.01*sign(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0)
    
    
    # denominator of acoate
    dem_acoate <- (d_pred_a1-d_pred_a0)#+0.01*sign(d_pred_a1-d_pred_a0)
    dem_coate <- (d_pred_b1-d_pred_b0)
    
    swate_x <- (y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0)/
      dem_swate
    acoate_x <- (y_pred_a1 - y_pred_a0)/
      dem_acoate
    
    coate_x <- (y_pred_b1 - y_pred_b0)/
      dem_coate
    
    # coeff_d <- 1/dem_swate
    var_d_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$d-d_pred_b1)
    var_d_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$d-d_pred_b0)
    var_d_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$d-d_pred_a1)
    var_d_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$d-d_pred_a0)
    # var_d_5 <- d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0
    
    
    # compute d2
    d_2 <- (var_y_1-var_y_2-var_y_3+var_y_4)/dem_swate-swate_x/dem_swate*
      (var_d_1-var_d_2-var_d_3+var_d_4)
    
    # compute d1
    d_1 <- (var_y_3-var_y_4)/dem_acoate-acoate_x/dem_acoate*
      (var_d_3-var_d_4)
    
    
    # compute d3
    d_3 <- (var_y_1-var_y_2)/dem_coate-coate_x/dem_coate*
      (var_d_1-var_d_2)
    
    # compute theta_p
    theta_p_x <- acoate_x-swate_x
    
    theta_p_x_2 <- acoate_x-coate_x
    
    theta_p_x_3 <- swate_x-coate_x
    
    
    predictset$outcome_h  <- theta_p_x
    predictset$outcome_h_2  <- theta_p_x_2
    predictset$outcome_h_3 <- theta_p_x_3
    
    model_fit <- lm(outcome_h~x1+x2, data = predictset)
    model_fit_2 <- lm(outcome_h_2~x1+x2, data = predictset)
    model_fit_3 <- lm(outcome_h_3~x1+x2, data = predictset)
    
    # beta_h is the plug-in estimator
    beta_h <- coef(model_fit)
    beta_h_2 <- coef(model_fit_2)
    beta_h_3 <- coef(model_fit_3)
    
    predictset$predict_h <- predict(model_fit, data = predictset)
    predictset$predict_h_2 <- predict(model_fit_2, data = predictset)
    predictset$predict_h_3 <- predict(model_fit_3, data = predictset)
    
    # matrix for covariates
    x_m_pred <- cbind(rep(1,nrow(predictset)),predictset$x1,predictset$x2)
    
    j_hat <- t(x_m_pred)%*%x_m_pred/nrow(predictset)
    
    
    # calculate the \varphi in theorem 5
    eff_score <- t(x_m_pred)%*%diag(c(predictset$outcome_h+d_1-d_2-predictset$predict_h))
    eff_score_2 <- t(x_m_pred)%*%diag(c(predictset$outcome_h_2+d_1-d_3-predictset$predict_h_2))
    eff_score_3 <- t(x_m_pred)%*%diag(c(predictset$outcome_h_3+d_2-d_3-predictset$predict_h_3))
    
    
    # -solve(j_hat) is the normalizing constant
    influ <- solve(j_hat)%*%eff_score
    influ_2 <- solve(j_hat)%*%eff_score_2
    influ_3 <- solve(j_hat)%*%eff_score_3
    
    # one-step correction
    beta_os <- beta_h+apply(influ,1,mean)
    beta_os_2 <- beta_h_2+apply(influ_2,1,mean)
    beta_os_3 <- beta_h_3+apply(influ_3,1,mean)
    
    
    beta_vec <- cbind(beta_vec, beta_os)
    influ_beta_vec <- cbind(influ_beta_vec,influ)
    
    beta_vec_2 <- cbind(beta_vec_2, beta_os_2)
    influ_beta_vec_2 <- cbind(influ_beta_vec_2,influ_2)
    
    beta_vec_3 <- cbind(beta_vec_3, beta_os_3)
    influ_beta_vec_3 <- cbind(influ_beta_vec_3,influ_3)
    
    
    # Estimating equation based estimator
    predictset$outcome_ee_h  <- d_1-d_2+theta_p_x
    predictset$outcome_ee_h_2  <- d_1-d_3+theta_p_x_2
    predictset$outcome_ee_h_3 <- d_2-d_3+theta_p_x_3
    
    ee_model_fit <- lm(outcome_ee_h~x1+x2, data = predictset)
    ee_model_fit_2 <- lm(outcome_ee_h_2~x1+x2, data = predictset)
    ee_model_fit_3 <- lm(outcome_ee_h_3~x1+x2, data = predictset)
    
    ee_beta_vec <- cbind(ee_beta_vec, coef(ee_model_fit))
    ee_beta_vec_2 <- cbind(ee_beta_vec_2, coef(ee_model_fit_2))
    ee_beta_vec_3 <- cbind(ee_beta_vec_3, coef(ee_model_fit_3))
  }
  beta_est <- apply(beta_vec,1,mean)
  beta_est_2 <- apply(beta_vec_2,1,mean)
  beta_est_3 <- apply(beta_vec_3,1,mean)
  
  
  ee_beta_est <- apply(ee_beta_vec,1,mean)
  ee_beta_est_2 <- apply(ee_beta_vec_2,1,mean)
  ee_beta_est_3 <- apply(ee_beta_vec_3,1,mean)
  # covariance matrix estimator
  beta_cov_hat <- cov(t(influ_beta_vec))/nrow(dataset)
  beta_cov_hat_2 <- cov(t(influ_beta_vec_2))/nrow(dataset)
  beta_cov_hat_3 <- cov(t(influ_beta_vec_3))/nrow(dataset)
  
  return(list(list(beta_est,ee_beta_est,beta_cov_hat),
              list(beta_est_2,ee_beta_est_2,beta_cov_hat_2),
              list(beta_est_3,ee_beta_est_3,beta_cov_hat_3)))
}



####### A function for summarizing the projection test results

test_summary <- function(est,cov,n){
  # est <-1
  # cov <- 1
  # n <- 1
  test_error <- try(W <- t(est)%*%solve(cov)%*%est, silent=TRUE)
  
  if (inherits(test_error, 'try-error')) {
    return(c(rep(NA,(21+2*length(est))),"error"))
  } 
  t <- 1*(c(W)>qchisq(seq(0,1,0.05),length(est)))
  est_out <- est
  se_out <- sqrt(diag(cov))
  
  return(c(t,est_out,se_out,"no error"))
}


####### Function for generating the nonparametric test results
kstest <- function(dataset,V=2)
{
  ####
  # V = 2
  ####
  
  N <- nrow(dataset)
  
  index <- createFolds(1:nrow(dataset),V)
  d_1_star_pre_vec <- NULL
  d_2_star_pre_vec <- NULL
  d_3_star_pre_vec <- NULL
  
  x1_vec <- NULL
  x2_vec <- NULL
  
  for(v in 1:V)
  {
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    
    # Model for P(Z=a1|G=a, X=x)
    trainingset_a <- trainingset[which(trainingset$g==0),]
    sl_z_a1 <- SuperLearner(Y = trainingset_a$z_h,
                            X = trainingset_a[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    # Model for P(Z=b1|G=b, X=x)
    trainingset_b <- trainingset[which(trainingset$g==1),]
    sl_z_b1 <- SuperLearner(Y = trainingset_b$z_h,
                            X = trainingset_b[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(G=1|X=x)
    sl_z_g1 <- SuperLearner(Y = trainingset$g,
                            X = trainingset[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # prediction for P(Z=a1|G=a,X)
    predictset_a1 <- predictset
    predictset_a1$z_h <- 1
    p_z_a1_a <- predict(sl_z_a1, newdata = predictset_a1, onlySL = T)$pred
    # prediction for P(Z=a0|G=a,X)
    p_z_a0_a <- 1 - p_z_a1_a
    
    # prediction for P(Z=b1|G=b,X)
    predictset_b1 <- predictset
    predictset_b1$z_h <- 1   
    
    p_z_b1_b <- predict(sl_z_b1, newdata = predictset_b1, onlySL = T)$pred
    
    # prediction for P(Z=b0|G=b,X)
    p_z_b0_b <- 1 - p_z_b1_b
    
    # prediction for P(G=a|X)
    p_g1 <- predict(sl_z_g1, newdata = predictset, onlySL = T)$pred
    p_g0 <- 1 - p_g1
    
    z_a0_pred <- p_g0*p_z_a0_a
    z_a1_pred <- p_g0*p_z_a1_a
    z_b0_pred <- p_g1*p_z_b0_b
    z_b1_pred <- p_g1*p_z_b1_b
    
    
    # Model and prediction for E(D|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_d_a0 <- SuperLearner(Y = train_a0$d,
                            X = train_a0[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a0 <- predict(sl_d_a0,newdata = predictset,onlySL = T)$pred
    
    # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
    
    # Model and prediction for E(D|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_d_a1 <- SuperLearner(Y = train_a1$d,
                            X = train_a1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a1 <- predict(sl_d_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(D|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_d_b0 <- SuperLearner(Y = train_b0$d,
                            X = train_b0[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b0 <- predict(sl_d_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(D|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_d_b1 <- SuperLearner(Y = train_b1$d,
                            X = train_b1[,c("x1","x2")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b1 <- predict(sl_d_b1,newdata = predictset,onlySL = T)$pred
    
    
    #####################################
    
    # Model and prediction for E(Y|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_y_a0 <- SuperLearner(Y = train_a0$y,
                            X = train_a0[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a0 <- predict(sl_y_a0,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_y_a1 <- SuperLearner(Y = train_a1$y,
                            X = train_a1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a1 <- predict(sl_y_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_y_b0 <- SuperLearner(Y = train_b0$y,
                            X = train_b0[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b0 <- predict(sl_y_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_y_b1 <- SuperLearner(Y = train_b1$y,
                            X = train_b1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b1 <- predict(sl_y_b1,newdata = predictset,onlySL = T)$pred
    
    # E[eta_b-eta_a]
    # coeff_y <- 1/(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0)
    var_y_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$y-y_pred_b1)
    var_y_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$y-y_pred_b0)
    var_y_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$y-y_pred_a1)
    var_y_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$y-y_pred_a0)
    # var_y_5 <- y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0
    
    # denominator of swate
    dem_swate <- (d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0) 
    
    
    # denominator of acoate
    dem_acoate <- (d_pred_a1-d_pred_a0)
    dem_coate <- (d_pred_b1-d_pred_b0)
    
    swate_x <- (y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0)/
      dem_swate
    acoate_x <- (y_pred_a1 - y_pred_a0)/
      dem_acoate
    
    coate_x <- (y_pred_b1 - y_pred_b0)/
      dem_coate
    
    # coeff_d <- 1/dem_swate
    var_d_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$d-d_pred_b1)
    var_d_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$d-d_pred_b0)
    var_d_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$d-d_pred_a1)
    var_d_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$d-d_pred_a0)
    # var_d_5 <- d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0
    
    
    # compute d2
    d_2 <- (var_y_1-var_y_2-var_y_3+var_y_4)/dem_swate-swate_x/dem_swate*
      (var_d_1-var_d_2-var_d_3+var_d_4)
    
    # compute d1
    d_1 <- (var_y_3-var_y_4)/dem_acoate-acoate_x/dem_acoate*
      (var_d_3-var_d_4)
    
    
    # compute d3
    d_3 <- (var_y_1-var_y_2)/dem_coate-coate_x/dem_coate*
      (var_d_1-var_d_2)
    
    # compute theta_p
    theta_p_x <- acoate_x-swate_x
    
    theta_p_x_2 <- acoate_x-coate_x
    
    theta_p_x_3 <- swate_x-coate_x
    
    d_1_star_pre <- d_1-d_2+theta_p_x
    d_2_star_pre <- d_1-d_3+theta_p_x_2
    d_3_star_pre <- d_2-d_3+theta_p_x_3
    
    x1_vec <- c(x1_vec, predictset$x1)
    x2_vec <- c(x2_vec, predictset$x2)
    
    
    d_1_star_pre_vec <- c(d_1_star_pre_vec, d_1_star_pre)
    d_2_star_pre_vec <- c(d_2_star_pre_vec, d_2_star_pre)
    d_3_star_pre_vec <- c(d_3_star_pre_vec, d_3_star_pre)
    
    
  }
  
  d_1_star_c <- function(c){
    indi <- 1*(x1_vec<=c[1])*(x2_vec<=c[2])
    out1 <- d_1_star_pre_vec*indi-mean(d_1_star_pre_vec*indi)
    out2 <- d_1_star_pre_vec*indi
    return(list(out1,out2))
  }
  
  d_2_star_c <- function(c){
    indi <- 1*(x1_vec<=c[1])*(x2_vec<=c[2])
    out1 <- d_2_star_pre_vec*indi-mean(d_2_star_pre_vec*indi)
    out2 <- d_2_star_pre_vec*indi
    
    return(list(out1,out2))
  }
  
  d_3_star_c <- function(c){
    indi <- 1*(x1_vec<=c[1])*(x2_vec<=c[2])
    out1 <- d_3_star_pre_vec*indi-mean(d_3_star_pre_vec*indi)
    out2 <- d_3_star_pre_vec*indi
    
    return(list(out1, out2))
  }
  
  ###### sup Omega
  ###### bootstrap to simulate the supremum of Gaussian process
  ###### Here we suppose to calculate a large covariance matrix (n times n)
  
  sup_omega_1 <- 0
  sup_omega_2 <- 0
  sup_omega_3 <- 0
  
  l2_omega_1 <- 0
  l2_omega_2 <- 0
  l2_omega_3 <- 0
  
  
  A1_matrix <- NULL
  A2_matrix <- NULL
  A3_matrix <- NULL
  
  for (i in c(1:N)) {
    col_1_i <- d_1_star_c(c(x1_vec[i],x2_vec[i]))
    col_2_i <- d_2_star_c(c(x1_vec[i],x2_vec[i]))
    col_3_i <- d_3_star_c(c(x1_vec[i],x2_vec[i]))
    
    sup_omega_1 <- max(sup_omega_1, abs(mean(col_1_i[[2]])))
    sup_omega_2 <- max(sup_omega_2, abs(mean(col_2_i[[2]])))
    sup_omega_3 <- max(sup_omega_3, abs(mean(col_3_i[[2]])))
    
    
    l2_omega_1 <- l2_omega_1 +(mean(col_1_i[[2]]))^2
    l2_omega_2 <- l2_omega_2 +(mean(col_2_i[[2]]))^2
    l2_omega_3 <- l2_omega_3 +(mean(col_3_i[[2]]))^2
    
    
    A1_matrix <- cbind(A1_matrix, col_1_i[[1]])
    A2_matrix <- cbind(A2_matrix, col_2_i[[1]])
    A3_matrix <- cbind(A3_matrix, col_3_i[[1]])
  }
  
  cov_1 <- 1/N*t(A1_matrix)%*%A1_matrix
  cov_2 <- 1/N*t(A2_matrix)%*%A2_matrix
  cov_3 <- 1/N*t(A3_matrix)%*%A3_matrix
  
  
  l2_omega_1 <- sqrt(l2_omega_1)
  l2_omega_2 <- sqrt(l2_omega_2)
  l2_omega_3 <- sqrt(l2_omega_3)
  
  
  ###### Simulate the Gaussian distribution
  
  boot.samp <- 1000
  
  paths_1 <- mvrnorm(boot.samp,mu=rep(0,N),Sigma=cov_1)
  paths_2 <- mvrnorm(boot.samp,mu=rep(0,N),Sigma=cov_2)
  paths_3 <- mvrnorm(boot.samp,mu=rep(0,N),Sigma=cov_3)
  
  stats_1 <- apply(abs(paths_1), 1, max)
  stats_2 <- apply(abs(paths_2), 1, max)
  stats_3 <- apply(abs(paths_3), 1, max)
  
  l2_stats_1 <- apply(abs(paths_1), 1, function(x) sqrt(sum(x^2)))
  l2_stats_2 <- apply(abs(paths_2), 1, function(x) sqrt(sum(x^2)))
  l2_stats_3 <- apply(abs(paths_3), 1, function(x) sqrt(sum(x^2)))
  
  reject_1 <- 1*((sup_omega_1*N^{1/2})>quantile(stats_1,probs = seq(0,1,0.05)))
  reject_2 <- 1*((sup_omega_2*N^{1/2})>quantile(stats_2,probs = seq(0,1,0.05)))
  reject_3 <- 1*((sup_omega_3*N^{1/2})>quantile(stats_3,probs = seq(0,1,0.05)))
  
  l2_reject_1 <- 1*((l2_omega_1*N^{1/2})>quantile(l2_stats_1,probs = seq(0,1,0.05)))
  l2_reject_2 <- 1*((l2_omega_2*N^{1/2})>quantile(l2_stats_2,probs = seq(0,1,0.05)))
  l2_reject_3 <- 1*((l2_omega_3*N^{1/2})>quantile(l2_stats_3,probs = seq(0,1,0.05)))
  
  
  return(list(reject_1,reject_2,reject_3,l2_reject_1,l2_reject_2,l2_reject_3))
}





####
######################### 2. A synthetic data analysis example

set.seed(1234)
# generat data
dat <- dat_gen(n = 2000, alpha = c(1,1,1), beta = c(4,4,4))
# estimation results:
est_results <- npest(dat)

# projection test
###### We have 3 tests based on 3 different null hypothesis, we used the result from the second test

project_test <- nptest(dat)[[2]]
project_test_res <- test_summary(est = project_test[[1]], cov = project_test[[3]], n = 1000)

# Output the result for a size 0.05 test
print(project_test_res[2]) # 1 means reject


# NP test
np_test_res <- kstest(dat)

# We have 6 tests in total (3 null hypothesis * 2 nonparametric test)
# We print the test result from the KS test under the second null hypothesis here
# Output the result for a size 0.05 test
print((np_test_res[[1]])[2])
