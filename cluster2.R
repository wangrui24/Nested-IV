# Cluster TASK_ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(touse+2000)

library(MASS)

expit <- function(x){
  y <- 1/(1+exp(-x))
  return(y)
}

expit(1)


expit <- function(x){
  y <- 1/(1+exp(-x))
  return(y)
}

expit(1)

# parameters for effect size : beta
# parameters for compliance rate ; alpha


dat_gen <- function(n=1000, alpha = c(0.5,0.2,0.05), beta = c(0,2,-3)){
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
  x5 <- runif(n, min = -3, max = 3)
  # 
  # 
  x6 <- rbinom(n, 4, 0.5)
  x7 <- rbinom(n, 1, 0.5)
  x8 <- rbinom(n, 1, 0.5)
  
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
  p_y_1 <- expit(1*(s=="ANT"|s=="AAT")*(1+x1+x2+x3+u)+
                   1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+x3+u)+
                   1*(s=="SWI1"|s=="SWI2")*(1+x1+x2+x3+u)+
                   1*(s=="ACO")*(1+x1+x2+x3+u))
  
  y_1 <- rbinom(n,1,p_y_1)
  
  
  p_y_0 <- expit(1*(s=="ANT"|s=="AAT")*(1+x1+2*x2+2*x3+x4+x5-x6-x7+x8+u)+
                   1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+2*x3+x4+x5-x6-x7+x8+u)+
                   1*(s=="SWI1"|s=="SWI2")*(beta[1]+beta[2]*x1+2*x2+x3+beta[3]*x6-x7+2*x8+u)+
                   1*(s=="ACO")*(1+x1+2*x2+0.2*x2^2+x3+x4+x5-x6-x7+x8+u))
  
  y_0 <- rbinom(n,1,p_y_0)
  
  y <- (1-d)*y_0+(d)*y_1
  
  
  dat <- data.frame(x1,x2,x3,x4,
                    x5,x6,x7,x8,
                    g,z_h,z,z_a0,z_a1,z_b0,z_b1,d,y,y_0, y_1,s)
  return(dat)
}




# binay outcome case
library(SuperLearner)
library(caret)
library(MASS)


logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}

npest<- function(dataset,V=2)
{
  ####
  # dataset <- read.csv("dat_true2.csv")
  # V = 2
  ####
  index <- createFolds(1:nrow(dataset),V)
  influ_swi_vec <- NULL
  ee_swi_vec <- NULL
  os_swi_vec <- NULL
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
                            X = trainingset_a[,c("x1","x2","x3","x4",
                                                 "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    # Model for P(Z=b1|G=b, X=x)
    trainingset_b <- trainingset[which(trainingset$g==1),]
    sl_z_b1 <- SuperLearner(Y = trainingset_b$z_h,
                            X = trainingset_b[,c("x1","x2","x3","x4",
                                                 "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(G=1|X=x)
    sl_z_g1 <- SuperLearner(Y = trainingset$g,
                            X = trainingset[,c("x1","x2","x3","x4",
                                               "x5","x6","x7","x8")],
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
                            X = train_a0[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a0 <- predict(sl_d_a0,newdata = predictset,onlySL = T)$pred
    
    # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
    
    # Model and prediction for E(D|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_d_a1 <- SuperLearner(Y = train_a1$d,
                            X = train_a1[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a1 <- predict(sl_d_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(D|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_d_b0 <- SuperLearner(Y = train_b0$d,
                            X = train_b0[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b0 <- predict(sl_d_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(D|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_d_b1 <- SuperLearner(Y = train_b1$d,
                            X = train_b1[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b1 <- predict(sl_d_b1,newdata = predictset,onlySL = T)$pred
    
    
    #####################################
    
    # Model and prediction for E(Y|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    sl_y_a0 <- SuperLearner(Y = train_a0$y,
                            X = train_a0[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a0 <- predict(sl_y_a0,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_y_a1 <- SuperLearner(Y = train_a1$y,
                            X = train_a1[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a1 <- predict(sl_y_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_y_b0 <- SuperLearner(Y = train_b0$y,
                            X = train_b0[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b0 <- predict(sl_y_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_y_b1 <- SuperLearner(Y = train_b1$y,
                            X = train_b1[,c("x1","x2","x3","x4",
                                            "x5","x6","x7","x8")],
                            family = binomial(),
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
    os_swi_vec[v] <- os_swi
    
    ee_swi <- mean(var_y_1-var_y_2-var_y_3+var_y_4+var_y_5)/mean(var_d_1-var_d_2-var_d_3+var_d_4+var_d_5)
    ee_swi_vec[v] <- ee_swi
    
    influ_swi_vec <- c(influ_swi_vec, influ_swi)
  }
  os_swi_est <- mean(os_swi_vec)
  ee_swi_est <- mean(ee_swi_vec)
  se_swi <- sd(os_swi_vec)/sqrt(length(influ_swi_vec))
  influ_swi <- mean(influ_swi_vec)
  se_influ <- sd(influ_swi_vec)/sqrt(length(influ_swi_vec))
  return(c(os_swi_est,ee_swi_est,se_influ))
}

beta_list <- list(c(0,1,-1),c(0,2,-3))
alpha_list <- list(c(-0.2,0.1,0.005),c(0.5,0.2,0.05),c(0.3,0.5,0.1),c(1,1,1))


simu_res = NULL
for (n in c(1000, 2000, 5000, 10000)){
  for (i in c(1:2)) {
    beta_i <- beta_list[[i]]
    for (j in c(1:4)) {
      alpha_j <- alpha_list[[j]]
      dt = dat_gen(n,alpha = alpha_j,beta = beta_i)
      # calculating results
      res_temp = c(npest(dt), n, i, j)
      simu_res = rbind(simu_res, res_temp)
    }
  }
}
FileSave <- paste0("/home/bzhang3/Wang_Rui/nested_iv/output2/sim_",touse,".csv")

write.table(simu_res, file=FileSave, row.names = FALSE)
