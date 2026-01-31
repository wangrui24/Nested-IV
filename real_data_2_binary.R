library(SuperLearner)
library(caret)
library(MASS)
set.seed(1234)
######################
######################  baseline estimated charateristics
######################

dat_2 <- read.csv("dat_2.csv")
dat_hf <- read.csv("dat_hf.csv")

dat_2$age_2 <-(dat_2$age)^2
dat_hf$age_2 <-(dat_hf$age)^2
dat_2_1 <- dat_2[(dat_2$z %in% c("a0","a1","b0","b1")),]

dat_2_2 <- dat_2[(dat_2$z %in% c("a0","a1","c0","c1")),]
dat_2_2$z[dat_2_2$z=="c0"] <- "b0"
dat_2_2$z[dat_2_2$z=="c1"] <- "b1"
dat_2_2$g[dat_2_2$g==2] <- 1

base_characteristic <- function(data){
  dat_train <- data[,c("sex_m","age","race_m", "educ_m1",
                       "educ_m2","bmi_m","cig_m1","cig_m2","d","z")]
  dat_train_a1 <- dat_train[dat_train$z=="a1",]
  m_d_a1  <- SuperLearner(Y = dat_train_a1$d,
                          X = dat_train_a1[,c("sex_m","age","race_m", "educ_m1",
                                              "educ_m2","bmi_m","cig_m1","cig_m2")],
                          family = binomial(),
                          SL.library = c("SL.glm", "SL.randomForest"),
                          cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  
  d_a1_pred <- predict(m_d_a1, newdata = dat_train, onlySL = T)$pred
  
  
  dat_train_b1 <- dat_train[dat_train$z=="b1",]
  m_d_b1  <- SuperLearner(Y = dat_train_b1$d,
                          X = dat_train_b1[,c("sex_m","age","race_m", "educ_m1",
                                              "educ_m2","bmi_m","cig_m1","cig_m2")],
                          family = binomial(),
                          SL.library = c("SL.glm", "SL.randomForest"),
                          cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  
  d_b1_pred <- predict(m_d_b1, newdata = dat_train, onlySL = T)$pred
  # mean(d_b1_pred-d_a1_pred)
  # mean(d_b1_pred)
  # mean(d_a1_pred)
  
  sw_ratio <- (d_b1_pred-d_a1_pred)/(mean(d_b1_pred-d_a1_pred))
  ac_ratio <- d_a1_pred/mean(d_a1_pred)
  
  hist(sw_ratio)
  hist(ac_ratio)
  # mean(d_a1_pred)
  # mean(d_b1_pred)
  # hist(sw_ratio)
  # hist(ac_ratio)
  # 
  res <- rbind(c(mean(dat_train$age*sw_ratio),mean(dat_train$age*ac_ratio)),
               c((1-mean(dat_train$sex_m*sw_ratio)),(1-mean(dat_train$sex_m*ac_ratio))),
               c(mean(dat_train$sex_m*sw_ratio),mean(dat_train$sex_m*ac_ratio)),
               c((1-mean(dat_train$race_m*sw_ratio)),(1-mean(dat_train$race_m*ac_ratio))),
               c(mean(dat_train$race_m*sw_ratio),mean(dat_train$race_m*ac_ratio)),
               c((1-mean(dat_train$educ_m1*sw_ratio)-mean(dat_train$educ_m2*sw_ratio)),
                 (1-mean(dat_train$educ_m1*ac_ratio)-mean(dat_train$educ_m2*ac_ratio))),
               c(mean(dat_train$educ_m1*sw_ratio),mean(dat_train$educ_m1*ac_ratio)),
               c(mean(dat_train$educ_m2*sw_ratio),mean(dat_train$educ_m2*ac_ratio)),
               c((1-mean(dat_train$cig_m1*sw_ratio)-mean(dat_train$cig_m2*sw_ratio)),
                 (1-mean(dat_train$cig_m1*ac_ratio)-mean(dat_train$cig_m2*ac_ratio))),
               c(mean(dat_train$cig_m1*sw_ratio),mean(dat_train$cig_m1*ac_ratio)),
               c(mean(dat_train$cig_m2*sw_ratio),mean(dat_train$cig_m2*ac_ratio)),
               c((1-mean(dat_train$bmi_m*sw_ratio)),(1-mean(dat_train$bmi_m*ac_ratio))),
               c(mean(dat_train$bmi_m*sw_ratio),mean(dat_train$bmi_m*ac_ratio)))
  row.names(res) <- c("age","female","male","white","minority","No hight school",
                      "High school","College or above","Non-smoker","Current smoker",
                      "Former smoker", "<25", ">25")
  return(res)
}


baseline_char_est1 <- base_characteristic(dat_2_1)
baseline_char_est2 <- base_characteristic(dat_2_2)
baseline_char_est_hf <- base_characteristic(dat_hf)
write.csv(baseline_char_est1,"results/real_data/baseline_char_est1.csv")
write.csv(baseline_char_est2,"results/real_data/baseline_char_est2.csv")
write.csv(baseline_char_est_hf,"results/real_data/baseline_char_est_hf.csv")



cov_reg <- c("sex_m","age","race_m", 
             "educ_m2","bmi_m","cig_m1")


#### Consider the outcome up to 15 years
dat_hf$colo_cancer[which((dat_hf$colo_cancer==1)&(dat_hf$colo_year>15))] <- 0
dat_2_1$colo_cancer[which((dat_2_1$colo_cancer==1)&(dat_2_1$colo_year>15))] <- 0
dat_2_2$colo_cancer[which((dat_2_2$colo_cancer==1)&(dat_2_2$colo_year>15))] <- 0

######################
######################  Estimating SWATE, ACOATE using EIF
######################


est_eif <- function(dataset, cov_reg, cov_test){
  V = 2
  
  index <- createFolds(1:nrow(dataset),V)
  
  # For estimation purpose
  influ_swi_vec <- NULL
  influ_aco_vec <- NULL
  os_swi_vec <- NULL
  os_aco_vec <- NULL
  
  
  ee_swi_vec <- NULL
  ee_aco_vec <- NULL
  
  # For testing purpose
  influ_beta_vec <- NULL
  beta_vec <- NULL
  influ_beta_vec_2 <- NULL
  beta_vec_2 <- NULL
  influ_beta_vec_3 <- NULL
  beta_vec_3 <- NULL
  for(v in 1:V)
  {
    ####
    #v=1
    ####
    trainingset <- dataset[-index[[v]],]
    predictset <- dataset[index[[v]],]
    
    
    # Model for P(Z=a1|G=a, X=x)
    trainingset_a <- trainingset[which(trainingset$g==0),]
    
    formula_reg <- as.formula(paste0(c("z_h~"),paste0(cov_reg, collapse = "+")))
    
    sl_z_a1 <- SuperLearner(Y = trainingset_a$z_h,
                            X = trainingset_a[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    # Model for P(Z=b1|G=b, X=x)
    trainingset_b <- trainingset[which(trainingset$g==1),]
    sl_z_b1 <- SuperLearner(Y = trainingset_b$z_h,
                            X = trainingset_b[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    # Model for P(G=1|X=x)
    sl_z_g1 <- SuperLearner(Y = trainingset$g,
                            X = trainingset[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    formula_reg <- as.formula(paste0(c("g~"),paste0(cov_reg, collapse = "+")))
    mm <- glm(formula_reg, data = trainingset)
    # summary(sl_z_g1)
    min(predict(mm, type = "response"))
    # prediction for P(Z=a1|G=a,X)
    predictset_a1 <- predictset
    predictset_a1$z_h <- 1
    p_z_a1_a <- predict(sl_z_a1, newdata = predictset_a1[,cov_reg], type = "response")$pred
    # prediction for P(Z=a0|G=a,X)
    p_z_a0_a <- 1 - p_z_a1_a
    
    # prediction for P(Z=b1|G=b,X)
    predictset_b1 <- predictset
    predictset_b1$z_h <- 1   
    
    p_z_b1_b <- predict(sl_z_b1, newdata = predictset_b1[,cov_reg], type = "response")$pred
    
    # prediction for P(Z=b0|G=b,X)
    p_z_b0_b <- 1 - p_z_b1_b
    
    # prediction for P(G=a|X)
    p_g1 <- predict(sl_z_g1, newdata = predictset[,cov_reg], type = "response")$pred
    p_g0 <- 1 - p_g1
    
    z_a0_pred <- p_g0*p_z_a0_a
    z_a1_pred <- p_g0*p_z_a1_a
    z_b0_pred <- p_g1*p_z_b0_b
    z_b1_pred <- p_g1*p_z_b1_b
    
    
    # Model and prediction for E(D|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    d_pred_a0 <- rep(0, nrow(predictset))
    
    # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
    
    # Model and prediction for E(D|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_d_a1 <- SuperLearner(Y = train_a1$d,
                            X = train_a1[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_a1 <- predict(sl_d_a1,newdata = predictset[,cov_reg],onlySL = T)$pred
    
    
    # Model and prediction for E(D|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    d_pred_b0 <- rep(0, nrow(predictset))
    
    
    
    
    # Model and prediction for E(D|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_d_b1 <- SuperLearner(Y = train_b1$d,
                            X = train_b1[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    d_pred_b1 <- predict(sl_d_b1,newdata = predictset[,cov_reg],onlySL = T)$pred
    
    
    #####################################
    predictset_y <- predictset[,cov_reg]
    predictset_y$colo_year <- 1000
    # Model and prediction for E(Y|Z=a0,X)
    train_a0 <- trainingset[trainingset$z=="a0", ]
    formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
    
    sl_y_a0 <- SuperLearner(Y = train_a0$colo_cancer,
                            X = train_a0[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))

    
    y_pred_a0 <- predict(sl_y_a0, newdata = predictset_y,onlySL = T)$pred
    # print(c(summary(y_pred_a0),"y_pred_a0"))
    # mean(predict(sl_y_a0,newdata = train_a0, type = "response"))
    # 
    # mean(y_pred_a0)
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
    
    sl_y_a1 <- SuperLearner(Y = train_a1$colo_cancer,
                            X = train_a1[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    y_pred_a1 <- predict(sl_y_a1, newdata = predictset_y,onlySL = T)$pred
    # print(c(summary(y_pred_a1),"y_pred_a1"))
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
    
    sl_y_b0 <- SuperLearner(Y = train_b0$colo_cancer,
                            X = train_b0[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    
    y_pred_b0 <- predict(sl_y_b0, newdata = predictset_y,onlySL = T)$pred
    # print(c(summary(y_pred_b0),"y_pred_b0"))
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
    
    sl_y_b1 <- SuperLearner(Y = train_b1$colo_cancer,
                            X = train_b1[,cov_reg],
                            family = binomial(),
                            SL.library = c("SL.glm", "SL.randomForest"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    # summary(glm(formula_reg, data = train_b0, family = binomial()))
    # summary(glm(formula_reg, data = train_a1, family = binomial()))
    
    y_pred_b1 <- predict(sl_y_b1, newdata = predictset_y,onlySL = T)$pred
    
    # print(c(summary(y_pred_b1),"y_pred_b1"))
    
    
    
    # Calculation for influence function of swi
    # E[eta_b-eta_a]
    coeff_y <- 1/(mean(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0))
    var_y_1 <- 1*(predictset$z_b1==1)/(z_b1_pred)*(predictset$colo_cancer-y_pred_b1)
    var_y_2 <- 1*(predictset$z_b0==1)/(z_b0_pred)*(predictset$colo_cancer-y_pred_b0)
    var_y_3 <- 1*(predictset$z_a1==1)/(z_a1_pred)*(predictset$colo_cancer-y_pred_a1)
    var_y_4 <- 1*(predictset$z_a0==1)/(z_a0_pred)*(predictset$colo_cancer-y_pred_a0)
    var_y_5 <- y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0
    
    
    
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
    # Calculation for influence function of sw
    
    
    influ_swi <- coeff_y*(var_y_1-var_y_2-var_y_3+var_y_4+var_y_5)-psi_p*coeff_d*
      (var_d_1-var_d_2-var_d_3+var_d_4+var_d_5)
    
    #summary(influ_swi)
    
    #hist(influ_swi)
    # Calculation for influence function of aco
    influ_aco <- (var_y_3-var_y_4+(y_pred_a1 -  y_pred_a0))/(mean(d_pred_a1-d_pred_a0))-
      mean(y_pred_a1 -  y_pred_a0)/
      (mean(d_pred_a1-d_pred_a0))^2*
      (var_d_3-var_d_4+(d_pred_a1 -  d_pred_a0))
    
    os_swi <- psi_p + mean(influ_swi)
    os_aco <- mean(y_pred_a1 - y_pred_a0)/mean(d_pred_a1-d_pred_a0) + mean(influ_aco)
    
    
    ee_swi <- mean(var_y_1-var_y_2-var_y_3+var_y_4+var_y_5)/mean(var_d_1-var_d_2-var_d_3+var_d_4+var_d_5)
    ee_aco <- mean(var_y_3-var_y_4+(y_pred_a1 -  y_pred_a0))/mean(var_d_3-var_d_4+(d_pred_a1 -  d_pred_a0))
    
    
    print(mean(y_pred_a1 - y_pred_a0)/mean(d_pred_a1-d_pred_a0))
    print(mean(influ_aco))
    os_swi_vec[v] <- os_swi
    os_aco_vec[v] <- os_aco
    
    ee_swi_vec[v] <- ee_swi
    ee_aco_vec[v] <- ee_aco
    
    influ_swi_vec <- c(influ_swi_vec, influ_swi)
    influ_aco_vec <- c(influ_aco_vec, influ_aco)
    
    
    
    ############################################## Testing part
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
    
    formula_reg <- as.formula(paste0(c("outcome_h~"),paste0(cov_test, collapse = "+")))
    model_fit <- lm(formula_reg, data = predictset)
    
    formula_reg <- as.formula(paste0(c("outcome_h_2~"),paste0(cov_test, collapse = "+")))
    model_fit_2 <- lm(formula_reg, data = predictset)
    
    formula_reg <- as.formula(paste0(c("outcome_h_3~"),paste0(cov_test, collapse = "+")))
    model_fit_3 <- lm(formula_reg, data = predictset)
    
    # beta_h is the plug-in estimator
    beta_h <- coef(model_fit)
    beta_h_2 <- coef(model_fit_2)
    beta_h_3 <- coef(model_fit_3)
    
    predictset$predict_h <- predict(model_fit, data = predictset)
    predictset$predict_h_2 <- predict(model_fit_2, data = predictset)
    predictset$predict_h_3 <- predict(model_fit_3, data = predictset)
    
    # matrix for covariates
    x_m_pred <- cbind(rep(1,nrow(predictset)),predictset[,cov_test])
    x_m_pred <- as.matrix(x_m_pred)
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
    
  }
  
  
  #### Results estimation
  os_swi_est <- round(mean(os_swi_vec),digits = 3)
  ee_swi_est <- round(mean(ee_swi_vec),digits = 3)
  # se_swi <- sd(os_swi_vec)/sqrt(length(influ_swi_vec))
  # influ_swi <- mean(influ_swi_vec)
  se_swi_influ <- round(sd(influ_swi_vec)/sqrt(length(influ_swi_vec)),digits = 3)
  
  os_aco_est <- round(mean(os_aco_vec),digits = 3)
  ee_aco_est <- round(mean(ee_aco_vec),digits = 3)
  # se_swi <- sd(os_swi_vec)/sqrt(length(influ_swi_vec))
  # influ_swi <- mean(influ_swi_vec)
  se_aco_influ <- round(sd(influ_aco_vec)/sqrt(length(influ_aco_vec)),digits = 3)
  
  
  
  os_swi_CI  <- paste("[",os_swi_est-1.96*se_swi_influ,",",os_swi_est+1.96*se_swi_influ,"]")
  os_aco_CI  <- paste("[",os_aco_est-1.96*se_aco_influ,",",os_aco_est+1.96*se_aco_influ,"]")
  
  ee_swi_CI  <- paste("[",ee_swi_est-1.96*se_swi_influ,",",ee_swi_est+1.96*se_swi_influ,"]")
  ee_aco_CI  <- paste("[",ee_aco_est-1.96*se_aco_influ,",",ee_aco_est+1.96*se_aco_influ,"]")
  
  est_res <- rbind(c(os_swi_est,os_swi_CI,ee_swi_est,ee_swi_CI),
                   c(os_aco_est,os_aco_CI,ee_aco_est,ee_aco_CI))
  
  
  
  ##########
  
  
  
  #### Results testing
  
  beta_est <- apply(beta_vec,1,mean)
  beta_est_2 <- apply(beta_vec_2,1,mean)
  beta_est_3 <- apply(beta_vec_3,1,mean)
  # covariance matrix estimator
  beta_cov_hat <- cov(t(influ_beta_vec))/nrow(dataset)
  beta_cov_hat_2 <- cov(t(influ_beta_vec_2))/nrow(dataset)
  beta_cov_hat_3 <- cov(t(influ_beta_vec_3))/nrow(dataset)
  
  
  W1 <- t(beta_est)%*%solve(beta_cov_hat)%*%beta_est
  W2 <- t(beta_est_2)%*%solve(beta_cov_hat_2)%*%beta_est_2
  W3 <- t(beta_est_3)%*%solve(beta_cov_hat_3)%*%beta_est_3
  
  return(list(est = est_res, test = c(W1,W2,W3,qchisq(0.95,length(beta_est)))))
}






######################
######################  Estimating SWATE, ACOATE using WALD
######################
est_para <- function(dataset,cov_reg){
  
  # Model for P(Z=a1|G=a, X=x)
  dataset_a <- dataset[which(dataset$g==0),]
  formula_reg <- as.formula(paste0(c("z_h~"),paste0(cov_reg, collapse = "+")))
  
  sl_z_a1 <- glm(formula_reg,
                 family = binomial, data = dataset_a)
  
  # summary(sl_z_a1)
  
  # Model for P(Z=b1|G=b, X=x)
  dataset_b <- dataset[which(dataset$g==1),]
  
  formula_reg <- as.formula(paste0(c("z_h~"),paste0(cov_reg, collapse = "+")))
  
  sl_z_b1 <- glm(formula_reg,
                 family = binomial, data = dataset_b)
  
  # summary(sl_z_b1)
  
  # Model for P(G=1|X=x)
  formula_reg <- as.formula(paste0(c("g~"),paste0(cov_reg, collapse = "+")))
  
  sl_z_g1 <- glm(formula_reg,
                 family = binomial, data = dataset)
  
  # summary(sl_z_g1)
  
  # prediction for P(Z=a1|G=a,X)
  dataset_a1 <- dataset
  dataset_a1$z_h <- 1
  p_z_a1_a <- predict(sl_z_a1, newdata = dataset_a1, type = "response")
  # prediction for P(Z=a0|G=a,X)
  p_z_a0_a <- 1 - p_z_a1_a
  
  # prediction for P(Z=b1|G=b,X)
  dataset_b1 <- dataset
  dataset_b1$z_h <- 1   
  
  p_z_b1_b <- predict(sl_z_b1, newdata = dataset_b1, type = "response")
  
  # prediction for P(Z=b0|G=b,X)
  p_z_b0_b <- 1 - p_z_b1_b
  
  # prediction for P(G=a|X)
  p_g1 <- predict(sl_z_g1, newdata = dataset, type = "response")
  p_g0 <- 1 - p_g1
  
  z_a0_pred <- p_g0*p_z_a0_a
  z_a1_pred <- p_g0*p_z_a1_a
  z_b0_pred <- p_g1*p_z_b0_b
  z_b1_pred <- p_g1*p_z_b1_b
  
  
  
  # Model and prediction for E(D|Z=a0,X)
  dataset_a0 <- dataset[dataset$z=="a0", ]
  # sl_d_a0 <- SuperLearner(Y = train_a0$d,
  #                         X = train_a0[,c("x1","x2","x3", "x4",
  #                                         "x5","x6","x7", "x8")],
  #                         family = binomial(),
  #                         SL.library = c("SL.glm", "SL.randomForest"),
  #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  
  # sl_d_a0 <- glm(d~sex_m+age+race_m+educ_m1+educ_m2+bmi_m+cig_m1+cig_m2,
  #                family = binomial, data = train_a0)
  # 
  # summary(sl_d_a0)
  
  d_pred_a0 <- rep(0, nrow(dataset))
  
  # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
  
  # Model and prediction for E(D|Z=a1,X)
  dataset_a1 <- dataset[dataset$z=="a1", ]
  # sl_d_a1 <- SuperLearner(Y = train_a1$d,
  #                         X = train_a1[,c("x1","x2","x3", "x4",
  #                                         "x5","x6","x7", "x8")],
  #                         family = binomial(),
  #                         SL.library = c("SL.glm", "SL.randomForest"),
  #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  formula_reg <- as.formula(paste0(c("d~"),paste0(cov_reg, collapse = "+")))
  sl_d_a1 <- glm(formula_reg,
                 family = binomial, data = dataset_a1)
  
  d_pred_a1 <- predict(sl_d_a1,newdata = dataset, type = "response")
  
  
  # Model and prediction for E(D|Z=b0,X), due to one sided non-compliance, they should be zero
  dataset_b0 <- dataset[dataset$z=="b0", ]
  # sl_d_b0 <- SuperLearner(Y = train_b0$d,
  #                         X = train_b0[,c("x1","x2","x3", "x4",
  #                                         "x5","x6","x7", "x8")],
  #                         family = binomial(),
  #                         SL.library = c("SL.glm", "SL.randomForest"),
  #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  
  d_pred_b0 <- rep(0, nrow(dataset))
  
  
  
  # Model and prediction for E(D|Z=b1,X)
  dataset_b1 <- dataset[dataset$z=="b1", ]
  # sl_d_b1 <- SuperLearner(Y = train_b1$d,
  #                         X = train_b1[,c("x1","x2","x3", "x4",
  #                                         "x5","x6","x7", "x8")],
  #                         family = binomial(),
  #                         SL.library = c("SL.glm", "SL.randomForest"),
  #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  formula_reg <- as.formula(paste0(c("d~"),paste0(cov_reg, collapse = "+")))
  sl_d_b1 <- glm(formula_reg,
                 family = binomial, data = dataset_b1)
  
  d_pred_b1 <- predict(sl_d_b1,newdata = dataset, type = "response")
  
  
  #####################################
  dataset_y <- dataset
  dataset_y$colo_year <- 1000
  # Model and prediction for E(Y|Z=a0,X)
  dataset_a0 <- dataset[dataset$z=="a0", ]
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  
  sl_y_a0 <- glm(formula_reg,family = binomial(), data = dataset_a0)
  
  y_pred_a0 <- predict(sl_y_a0,newdata = dataset_y, type = "response")
  
  # mean(predict(sl_y_a0,newdata = train_a0, type = "response"))
  # 
  # mean(y_pred_a0)
  # Model and prediction for E(Y|Z=a1,X)
  dataset_a1 <- dataset[dataset$z=="a1", ]
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  
  sl_y_a1 <- glm(formula_reg,family = binomial(), data = dataset_a1)
  
  y_pred_a1 <- predict(sl_y_a1,newdata = dataset_y, type = "response")
  
  
  # Model and prediction for E(Y|Z=b0,X)
  dataset_b0 <- dataset[dataset$z=="b0", ]
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  
  sl_y_b0 <- glm(formula_reg, family = binomial(), data = dataset_b0)
  
  y_pred_b0 <- predict(sl_y_b0,newdata = dataset_y, type = "response")
  
  
  
  # Model and prediction for E(Y|Z=b1,X)
  dataset_b1 <- dataset[dataset$z=="b1", ]
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  
  sl_y_b1 <- glm(formula_reg,family = binomial(), data = dataset_b1)
  
  y_pred_b1 <- predict(sl_y_b1,newdata = dataset_y, type = "response")
  
  
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
  var_y_1 <- 1*(dataset$z_b1==1)/(z_b1_pred)*(dataset$colo_cancer-y_pred_b1)
  var_y_2 <- 1*(dataset$z_b0==1)/(z_b0_pred)*(dataset$colo_cancer-y_pred_b0)
  var_y_3 <- 1*(dataset$z_a1==1)/(z_a1_pred)*(dataset$colo_cancer-y_pred_a1)
  var_y_4 <- 1*(dataset$z_a0==1)/(z_a0_pred)*(dataset$colo_cancer-y_pred_a0)
  var_y_5 <- y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0
  
  
  
  swate <- mean(y_pred_b1 -  y_pred_b0 - y_pred_a1 + y_pred_a0)/
    mean(d_pred_b1-d_pred_b0-d_pred_a1+d_pred_a0)
  
  acoate <- mean(y_pred_a1 - y_pred_a0)/
    mean(d_pred_a1-d_pred_a0)
  
  
  #################################
  ################################# Usual IV analysis
  #################################
  
  # Model and prediction for E(D|Z=a0,X)
  dataset_0 <- dataset[dataset$z_h==0, ]
  
  d_pred_0 <- rep(0, nrow(dataset))
  
  # mean(predictset$d[predictset$z=="a0"]-d_pred_a0[predictset$z=="a0"])
  
  # Model and prediction for E(D|Z=a1,X)
  dataset_1 <- dataset[dataset$z_h==1, ]
  # sl_d_a1 <- SuperLearner(Y = train_a1$d,
  #                         X = train_a1[,c("x1","x2","x3", "x4",
  #                                         "x5","x6","x7", "x8")],
  #                         family = binomial(),
  #                         SL.library = c("SL.glm", "SL.randomForest"),
  #                         cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
  formula_reg <- as.formula(paste0(c("d~"),paste0(cov_reg, collapse = "+")))
  
  sl_d_1 <- glm(formula_reg,
                family = binomial, data = dataset_1)
  
  d_pred_1 <- predict(sl_d_1,newdata = dataset, type = "response")
  
  
  
  
  #####################################
  dataset_y <- dataset
  dataset_y$colo_year <- 1000
  # Model and prediction for E(Y|Z=a0,X)
  dataset_0 <- dataset[dataset$z_h==0, ]
  
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  sl_y_0 <- glm(formula_reg,family = binomial(), data = dataset_0)
  
  y_pred_0 <- predict(sl_y_0,newdata = dataset_y, type = "response")
  
  # mean(predict(sl_y_a0,newdata = train_a0, type = "response"))
  # 
  # mean(y_pred_a0)
  # Model and prediction for E(Y|Z=a1,X)
  dataset_1 <- dataset[dataset$z_h==1, ]
  
  formula_reg <- as.formula(paste0(c("colo_cancer~"),paste0(cov_reg, collapse = "+")))
  sl_y_1 <- glm(formula_reg,family = binomial(), data = dataset_1)
  
  y_pred_1 <- predict(sl_y_1,newdata = dataset_y, type = "response")
  
  
  
  cate <- mean(y_pred_1 -  y_pred_0)/
    mean(d_pred_1-d_pred_0)
  
  return(c(swate,acoate,cate))
}


cov_reg <- c("sex_m","age","race_m", "educ_m1",
             "educ_m2","bmi_m","cig_m1","cig_m2")
cov_test <- c("sex_m","age","race_m", "educ_m1",
              "educ_m2","bmi_m","cig_m1","cig_m2")
# cov_reg <- c("sex_m")
# cov_test <- c("sex_m")
set.seed(1234)
est_eif_1 <- est_eif(dat_2_1,cov_reg,cov_test)
est_eif_2 <- est_eif(dat_2_2,cov_reg,cov_test)
est_eif_hf <- est_eif(dat_hf,cov_reg,cov_test)
write.csv(est_eif_1$est,"results/real_data/est_eif_1.csv")
write.csv(est_eif_2$est,"results/real_data/est_eif_2.csv")
write.csv(est_eif_hf$est,"results/real_data/est_eif_hf.csv")


write.csv(est_eif_1$test,"results/real_data/test_eif_1.csv")
write.csv(est_eif_2$test,"results/real_data/test_eif_2.csv")
write.csv(est_eif_hf$test,"results/real_data/test_eif_hf.csv")


set.seed(1234)
wald_est_1 <- est_para(dat_2_1,cov_reg)
wald_est_2 <- est_para(dat_2_2,cov_reg)
wald_est_3 <- est_para(dat_hf,cov_reg)
#### Bootstrap



swate_boots <- c()
acoate_boots <- c()
cate_boots <- c()

iter <- 1000
for (i in c(1:iter)) {
  boots_sample <- sample(c(1:nrow(dat_hf)), replace = T, size = nrow(dat_hf))
  dat_boots <- dat_hf[boots_sample,]
  res_i <- est_para(dat_boots,cov_reg)
  swate_boots[i] <- res_i[1]
  acoate_boots[i] <- res_i[2]
  cate_boots[i] <- res_i[3]
  print(i)
}
sd_hf <- c(sd(swate_boots),sd(acoate_boots),sd(cate_boots))
CI_low_hf <- wald_est_3 -1.96*sd_hf
CI_high_hf <- wald_est_3 +1.96*sd_hf

CI_hf <- paste0("[",CI_low_hf,",",CI_high_hf,"]")







swate_boots <- c()
acoate_boots <- c()
cate_boots <- c()

iter <- 1000
for (i in c(1:iter)) {
  boots_sample <- sample(c(1:nrow(dat_2_1)), replace = T, size = nrow(dat_2_1))
  dat_boots <- dat_2_1[boots_sample,]
  res_i <- est_para(dat_boots,cov_reg)
  swate_boots[i] <- res_i[1]
  acoate_boots[i] <- res_i[2]
  cate_boots[i] <- res_i[3]
  print(i)
}
sd_1 <- c(sd(swate_boots),sd(acoate_boots),sd(cate_boots))
CI_low_1 <- wald_est_1 -1.96*sd_1
CI_high_1 <- wald_est_1 +1.96*sd_1

CI_1 <- paste0("[",CI_low_1,",",CI_high_1,"]")




swate_boots <- c()
acoate_boots <- c()
cate_boots <- c()

iter <- 1000
for (i in c(1:iter)) {
  boots_sample <- sample(c(1:nrow(dat_2_2)), replace = T, size = nrow(dat_2_2))
  dat_boots <- dat_2_2[boots_sample,]
  res_i <- est_para(dat_boots,cov_reg)
  swate_boots[i] <- res_i[1]
  acoate_boots[i] <- res_i[2]
  cate_boots[i] <- res_i[3]
  # print(i)
}
sd_2 <- c(sd(swate_boots),sd(acoate_boots),sd(cate_boots))
CI_low_2 <- wald_est_2 -1.96*sd_2
CI_high_2 <- wald_est_2 +1.96*sd_2

CI_2 <- paste0("[",CI_low_2,",",CI_high_2,"]")


res_wald <- rbind(c(wald_est_3, CI_hf),c(wald_est_1, CI_1),c(wald_est_2, CI_2))
write.csv(res_wald,"results/real_data/res_wald.csv")