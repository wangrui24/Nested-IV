# Cluster TASK_ID
touse<-as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(touse+2000)

library(SuperLearner)
library(caret)
library(MASS)


logit <- function(x){log(x/(1-x))}
expit <- function(x){exp(x)/(1+exp(x))}


# parameters for effect size : beta
# parameters for compliance rate ; alpha


dat_gen <- function(n=1000,alpha = 0.2, beta = c(1,2,2)){
  # n = 1000
  mu <- c(0,0)
  var <- matrix(data = c(1,0,0,1), ncol = 2, nrow = 2)
  x <- mvrnorm(n, mu, var)
  x1 <- x[,1]
  x2 <- x[,2]
  
  x_star <- 1*(x>0)
  
  
  # Generate G
  pg <- expit(0+0.1*(x1>0)-0.1*(x2>0))
  g <- rbinom(n,1,pg)
  # pa <- expit(1+(x1>0)-2*(x2>0))
  # pb <- expit(1+2*(x1>0)-1*(x2>0))
  
  
  # pg <- 0.5
  # g <- rbinom(n,1,pg)
  # pa <- 0.5
  # pb <- 0.5
  # 
  # # Generate Z
  # z_h <- rbinom(n,1,g*pa+(1-g)*pb)
  
  
  #  pb <- expit(1+(x1>0)-2*(x2>0))
  #  pa <- expit(1+2*(x1>0)-1*(x2>0))
  # 
  # 
  # pg <- 0.5
  # g <- rbinom(n,1,pg)
  pa <- 0.5+0.1*1*(x1>0)-0.1*1*(x2>0)
  pb <- 0.5-0.1*1*(x1>0)+0.1*1*(x2>0)
  # 
  # # Generate Z
  z_h <- rbinom(n,1,g*pb+(1-g)*pa)
  
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
  u <- rnorm(n, -0.3, 0.3)
  
  # Generate principal stratum
  par_ant <- c(0,-1)
  par_swi_1 <- c(1,2)
  par_swi_2 <- c(1,2)
  par_atnt <- c(0.5,1)
  par_aco <- c(1, 1)
  par_ntat <- c(-0.5,-1)
  par_aat <- c(2,0)
  
  # set the parameter
  
  
  num_ant <- c(exp(1+x_star%*%par_ant+0.3*(u>0)))
  num_swi_1 <- c(exp(3.5+0.5*x_star%*%par_swi_1+0.1*(u>0)))
  num_swi_2 <- c(exp(3.5+0.5*x_star%*%par_swi_2+0.1*(u>0)))
  num_atnt <- c(exp(1+x_star%*%par_atnt+0.5*(u>0)))
  num_aco <- c(exp(2.5+x_star%*%par_aco+0.1*(u>0)))
  num_ntat <- c(exp(1+x_star%*%par_ntat-0.5*(u>0)))
  num_aat <- c(exp(1+x_star%*%par_aat-(u>0)))
  
  
  
  # num_ant <- 0
  # num_swi_1 <- c(exp(alpha[1]+alpha[2]*x%*%par_swi_1+alpha[3]*u))
  # num_swi_2 <- 0
  # num_atnt <- 0
  # num_aco <- c(exp(3+x%*%par_aco+0.5*u))
  # num_ntat <- 0
  # num_aat <- 0
  
  dem <- num_ant + num_swi_1 + num_swi_2 + num_atnt + num_aco + num_ntat + num_aat
  p_ant <- num_ant/dem
  p_swi_1 <- (num_swi_1+num_swi_2+num_aco)/dem*(1-alpha)/2
  p_swi_2 <- (num_swi_1+num_swi_2+num_aco)/dem*(1-alpha)/2
  p_atnt <- num_atnt/dem
  p_aco <- (num_swi_1+num_swi_2+num_aco)/dem*alpha
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
  epsilon <- rnorm(n,0,1)
  # Generate potential outcomes of Y
  y_0 <- 1*(s=="ANT"|s=="AAT")*(1+x1+x2+u)+
    1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+u)+
    1*(s=="SWI1"|s=="SWI2")*(1+x1+x2+u)+
    1*(s=="ACO")*(1+x1+x2+u)+epsilon
  
  
  
  
  # y_1 <- 1*(s=="ANT"|s=="AAT")*(1+x1+2*x2+2*x3+u)+
  #   1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+2*x3+u)+
  #   1*(s=="SWI1"|s=="SWI12")*(beta[1]+beta[2]*x1+2*x2+x3+beta[3]*x3^2+u)+
  #   1*(s=="ACO")*(1+x1+2*x2+0.2*x2^2+x3+u)
  
  y_1 <- 1*(s=="ANT"|s=="AAT")*(1+x1+x2+u)+
    1*(s=="ATNT"|s=="NTAT")*(1+x1+x2+u)+
    1*(s=="SWI1"|s=="SWI2")*(beta[1]+2*x1+2*x2+beta[2]*cos(x1)+beta[3]*expit(x2)+u)+
    1*(s=="ACO")*(1+2*x1+2*x2+2*cos(x1)+2*expit(x2)+u)+epsilon
  
  y <- (1-d)*y_0+(d)*y_1
  
  
  dat <- data.frame(x1,x2,
                    g,z_h,z,z_a0,z_a1,z_b0,z_b1,d,y,y_0, y_1,s,u)
  return(dat)
}

beta_list <- list(c(1,2,2),c(1.5,2.5,2.5))
alpha_list <- c(0.1,0.2,0.4,0.6,0.8,0.9)
matrix_prop <- NULL
for (i in c(1:length(alpha_list))) {
  dataset <- dat_gen(10000,alpha = alpha_list[i])
  matrix_prop <- rbind(matrix_prop,c(mean(dataset$s=="SWI1"|dataset$s=="SWI2"),mean(dataset$s=="ACO")))
}
matrix_prop
# dataset <- dat_gen(5000,alpha = c(2,0.5,0.1))
# 
# dat_swi <- dataset[(dataset$s=="SWI1"|dataset$s=="SWI2"),]
# nrow(dat_swi)/5000
# mean(dataset$s=="ACO")
# dat_swi$swate <- dat_swi$y_1-dat_swi$y_0
# 
# dat_aco <- dataset[(dataset$s=="ACO"),]
# dat_aco$acoate <- dat_aco$y_1-dat_aco$y_0
# 
# dat_co <- dataset[(dataset$s=="ACO"|dataset$s=="SWI1"|dataset$s=="SWI2"),]
# dat_co$coate <- dat_co$y_1-dat_co$y_0
# 
# lm(swate~x1+x2, dat = dat_swi)
# lm(acoate~x1+x2, dat = dat_aco)
# lm(coate~x1+x2, dat = dat_co)



nptest <- function(dataset,V=2)
{
  ####
  # V = 2
  ####
  #dataset$sinx1 <- cos(dataset$x1)
  #dataset$expitx2 <- expit(dataset$x2)
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
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a0 <- predict(sl_y_a0,newdata = predictset,onlySL = T)$pred
    
    # Model and prediction for E(Y|Z=a1,X)
    train_a1 <- trainingset[trainingset$z=="a1", ]
    sl_y_a1 <- SuperLearner(Y = train_a1$y,
                            X = train_a1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_a1 <- predict(sl_y_a1,newdata = predictset,onlySL = T)$pred
    
    
    # Model and prediction for E(Y|Z=b0,X)
    train_b0 <- trainingset[trainingset$z=="b0", ]
    sl_y_b0 <- SuperLearner(Y = train_b0$y,
                            X = train_b0[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
                            cvControl  = list(V = 5L,shuffle = TRUE, validRows = NULL))
    
    y_pred_b0 <- predict(sl_y_b0,newdata = predictset,onlySL = T)$pred
    
    
    
    # Model and prediction for E(Y|Z=b1,X)
    train_b1 <- trainingset[trainingset$z=="b1", ]
    sl_y_b1 <- SuperLearner(Y = train_b1$y,
                            X = train_b1[,c("x1","x2")],
                            family = gaussian(),
                            SL.library = c("SL.glm", "SL.randomForest","SL.gam"),
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
# length(test_summary(c(1,1),diag(c(1,1)),2))
# length(test_summary(c(0,1),diag(c(0,0)),2))
# 
# # This function is used for calculating the standard error and if reject
# test_summary <- function(est,cov,n){
#   W <- t(est)%*%solve(cov)%*%est
#   t <- 1*(c(W)>qchisq(seq(0,1,0.05),length(est)))
#   est_out <- est
#   se_out <- sqrt(diag(cov))
#   return(c(t,est_out,se_out))
# }


beta_list <- list(c(1,2,2),c(1.5,2.5,2.5), c(2,3,3))
alpha_list <- c(0.1,0.2,0.4,0.6,0.8,0.9)


simu_res = NULL
for (n in c(2000, 5000, 10000, 15000, 20000)){
  for (i in c(1:length(beta_list))) {
    for (j in c(1:length(alpha_list))) {
    dt = dat_gen(n, beta = beta_list[[i]],alpha = alpha_list[j])
    res_list = nptest(dt)
    res_test_1 <- res_list[[1]]
    res_test_2 <- res_list[[2]]
    res_test_3 <- res_list[[3]]
    simu_res = rbind(simu_res,c(test_summary(res_test_1[[1]],res_test_1[[3]],n),n,i,j,"test 1","os"))
    simu_res = rbind(simu_res,c(test_summary(res_test_2[[1]],res_test_2[[3]],n),n,i,j,"test 2","os"))
    simu_res = rbind(simu_res,c(test_summary(res_test_3[[1]],res_test_3[[3]],n),n,i,j,"test 3","os"))
    
    simu_res = rbind(simu_res,c(test_summary(res_test_1[[2]],res_test_1[[3]],n),n,i,j,"test 1","ee"))
    simu_res = rbind(simu_res,c(test_summary(res_test_2[[2]],res_test_2[[3]],n),n,i,j,"test 2","ee"))
    simu_res = rbind(simu_res,c(test_summary(res_test_3[[2]],res_test_3[[3]],n),n,i,j,"test 3","ee"))
    
    }
  }
}


FileSave <- paste0("/home/bzhang3/Wang_Rui/nested_iv/output_test_wald_nonlinear/sim_",touse,".csv")

write.table(simu_res, file=FileSave, row.names = FALSE)
