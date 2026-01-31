set.seed(1234)

# setwd("E:/research/projects/Nested IV/Nested IV")
library(tidyverse)
library(tableone)
library(SuperLearner)

dat_full <- read.csv("colo_data_mar22_d032222.csv")

dat <- dat_full %>% select(arm, center, dual,bq_returned, fsg_result0,colo_eligible_bq,
                           sex,age, agelevel,race7,educat,cig_stat,bmi_curc,
                           marital, occupat, colo_fh, rndyear,colo_cancer,colo_exitdays)


dat_full_1 <- dat_full[dat_full$arm==1,]
sum(dat_full_1$colo_cancer)
median(dat_full_1$colo_cancer_diagdays, na.rm = T)
median(dat_full_1$colo_exitdays)/365


# dat <- dat %>% 
#   mutate(D = case_when(
#     fsg_result0 == "C" ~ 0,
#     fsg_result0 == "1" ~ 1,
#     fsg_result0 == "2" ~ 1,
#     fsg_result0 == "3" ~ 1,
#     fsg_result0 == "4" ~ 1,
#     fsg_result0 == "8" ~ 0,
#     fsg_result0 == "9" ~ 0,
#   )) %>% 
#   mutate(Z_h = case_when(
#     arm == 1 ~ 1,
#     arm == 2 ~ 0
#   ))
# 
# dat_treat <- dat %>% filter(arm == 1) %>% 
#   filter(colo_eligible_bq == 1) %>% 
#   mutate(compliance = case_when(
#     Z_h == D ~ 1,
#     Z_h != D ~ 0
#   )) %>% 
#   drop_na()
# 
# 
# 
# 
# covars <- c("age", "arm", "sex", "bq_returned","D",
#             "center", "dual", "agelevel","race7","educat","cig_stat",
#             "bmi_curc",
#             "marital", "occupat")
# 
# cat_covars <- c("arm", "sex", "bq_returned","D",
#                 "center", "dual", "agelevel","race7","educat","cig_stat",
#                 "bmi_curc",
#                 "marital", "occupat")
# 
# tab1 <- CreateTableOne(vars = covars, strata = "compliance",
#                        data = dat_treat, factorVars = cat_covars)
# tab1_out <- print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
# 
# write.csv(tab1_out, "tab1_out.csv")

dat_hf_pre <- dat %>% filter(center == 4) %>% filter(colo_eligible_bq == 1)
table(dat_hf_pre$rndyear)
table(dat_hf_pre$rndyear[(dat_hf_pre$dual==1)&(dat_hf_pre$arm==1)])
table(dat_hf_pre$rndyear[(dat_hf_pre$dual==0)&(dat_hf_pre$arm==1)])

############## Data for analysis, Henry Ford Health system
dat_hf <- dat %>% filter(center == 4) %>% filter(colo_eligible_bq == 1) %>% 
  filter(rndyear != 1997) %>% 
  mutate(z_h = case_when(
    arm == 1 ~ 1,
    arm == 2 ~ 0
  )) %>% 
  mutate(d = case_when(
    fsg_result0 == "C" ~ 0,
    fsg_result0 == "1" ~ 1,
    fsg_result0 == "2" ~ 1,
    fsg_result0 == "3" ~ 1,
    fsg_result0 == "4" ~ 1,
    fsg_result0 == "8" ~ 0,
    fsg_result0 == "9" ~ 0,
    z_h == 0 ~ 0
  )) %>% mutate(g = case_when(
    rndyear<1997 ~ 0,
    rndyear>1997 ~ 1
  )) %>% 
  mutate(z = case_when(
    (z_h == 0)&(g == 0) ~ "a0",
    (z_h == 1)&(g == 0) ~ "a1",
    (z_h == 0)&(g == 1) ~ "b0",
    (z_h == 1)&(g == 1) ~ "b1",
  )) %>% 
  mutate(z_a0 = case_when(
    z == "a0" ~ 1,
    z != "a0" ~ 0
  )) %>% 
  mutate(z_a1 = case_when(
    z == "a1" ~ 1,
    z != "a1" ~ 0
  )) %>% 
  mutate(z_b0 = case_when(
    z == "b0" ~ 1,
    z != "b0" ~ 0
  )) %>% 
  mutate(z_b1 = case_when(
    z == "b1" ~ 1,
    z != "b1" ~ 0
  )) %>% mutate(colo_year = colo_exitdays/365) %>%drop_na("d", "age", "sex","race7","educat","cig_stat",
                 "bmi_curc",
                 "marital", "occupat", "colo_cancer","colo_exitdays") %>% 
  mutate(sex_m = 2-sex) %>% 
  mutate(race_m = case_when(
    race7 == 1 ~ 0,
    race7 != 1 ~ 1
  )) %>% 
  mutate(educ_m1 = case_when(
    (educat==3) |(educat==4) ~ 1,
    .default = 0
  )) %>% 
  mutate(educ_m2 =case_when(
    educat >4 ~ 1,
    .default = 0
  )) %>% 
  mutate(bmi_m =case_when(
    bmi_curc <=2 ~ 0,
    bmi_curc >2 ~ 1
  )) %>% 
  mutate(cig_m1 =case_when(
    cig_stat ==1 ~ 1,
    .default = 0
  ))%>% 
  mutate(cig_m2 =case_when(
    cig_stat ==2 ~ 1,
    .default = 0
  )) %>% 
  mutate(educ_m = case_when(
    educat <=2 ~ 0,
    (educat==3) |(educat==4) ~ 1,
    educat >4 ~ 2
  )) %>% 
  mutate(cig_m = case_when(
    cig_stat ==1 ~ 1,
    cig_stat ==2 ~ 2,
    .default = 0
  ))


table(dat_hf$rndyear[dat_hf$dual==1])
table(dat_hf$rndyear[dat_hf$dual==0])

table(dat_hf$rndyear[(dat_hf$dual==0)&(dat_hf$arm==1)])
table(dat_hf$rndyear[(dat_hf$dual==1)&(dat_hf$arm==1)])

# Descriptive statistics

# covars <- c("d", "age", "sex","race7","educat","cig_stat",
#             "bmi_curc",
#             "marital", "occupat")
# 
# cat_covars <- c("d", "sex","race7","educat","cig_stat",
#                 "bmi_curc", "marital", "occupat")
# 
# tab1 <- CreateTableOne(vars = covars, strata = "z",
#                        data = dat_hf, factorVars = cat_covars)
# des_hf <- print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))
# 
# write.csv(des_hf, "des_hf.csv")



covars <- c("d", "age", "sex_m","race_m","educ_m","cig_m",
            "bmi_m", "colo_cancer","colo_year")

cat_covars <- c("d", "sex_m","race_m","educ_m","cig_m",
                "bmi_m","colo_cancer")

tab1 <- CreateTableOne(vars = covars, strata = "z",
                       data = dat_hf, factorVars = cat_covars)
des_hf <- print(tab1, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

write.csv(des_hf, "des_hf.csv")
write.csv(dat_hf, "dat_hf.csv")


mean(dat_hf$d[dat_hf$z=="a1"])-mean(dat_hf$d[dat_hf$z=="a0"])
mean(dat_hf$d[dat_hf$z=="b1"])-mean(dat_hf$d[dat_hf$z=="b0"])


# ITT effect
m1 <- glm(colo_cancer~d+g+d:g+sex_m+age+race_m+educ_m1+educ_m2+bmi_m+cig_m1+cig_m2,
    offset = colo_year, family = poisson, data = dat_hf)
summary(m1)
confint(m1)
vcov(m1)
sd_b <- sqrt(0.03947+0.07563-2*0.03946)
exp(-0.453+0.09647+1.96*sd_b)
exp(-0.453+0.09647-1.96*sd_b)

dat$colo_year <- dat$colo_exitdays/365
dat$d <- 0
dat$d[dat$arm==1] <- 1

m1 <- glm(colo_cancer~d, offset = colo_year, family = poisson, data = dat)
summary(m1)

# produce the summary statistics for age, sex and race
dat_train <- dat_hf[,c("sex_m","age","race_m", "educ_m1",
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
mean(d_b1_pred-d_a1_pred)
mean(d_b1_pred)
mean(d_a1_pred)

sw_ratio <- (d_b1_pred-d_a1_pred)/(mean(d_b1_pred-d_a1_pred))
ac_ratio <- d_a1_pred/mean(d_a1_pred)
mean(d_a1_pred)
mean(d_b1_pred)
hist(sw_ratio)
hist(ac_ratio)

mean(dat_train$sex_m*sw_ratio)
mean(dat_train$sex_m*ac_ratio)


mean(dat_train$age*sw_ratio)
mean(dat_train$age*ac_ratio)

mean(dat_train$race_m*sw_ratio)
mean(dat_train$race_m*ac_ratio)

mean(dat_train$educ_m1*sw_ratio)
mean(dat_train$educ_m1*ac_ratio)


mean(dat_train$educ_m2*sw_ratio)
mean(dat_train$educ_m2*ac_ratio)


mean(dat_train$bmi_m*sw_ratio)
mean(dat_train$bmi_m*ac_ratio)

mean(dat_train$cig_m1*sw_ratio)
mean(dat_train$cig_m1*ac_ratio)

mean(dat_train$cig_m2*sw_ratio)
mean(dat_train$cig_m2*ac_ratio)



mean(d_x[dat_hf$z=="a1"])-mean(d_x[dat_hf$z=="a0"])




############## Data for analysis, different centers


dat_2 <- dat %>% filter(colo_eligible_bq == 1) %>% 
  filter((center == 1)|(center == 4)|(center == 11)) %>% 
  filter(!(rndyear >= 1997 & center == 4)) %>% 
  mutate(z_h = case_when(
    arm == 1 ~ 1,
    arm == 2 ~ 0
  )) %>% 
  mutate(d = case_when(
    fsg_result0 == "C" ~ 0,
    fsg_result0 == "1" ~ 1,
    fsg_result0 == "2" ~ 1,
    fsg_result0 == "3" ~ 1,
    fsg_result0 == "4" ~ 1,
    fsg_result0 == "8" ~ 0,
    fsg_result0 == "9" ~ 0,
    z_h == 0 ~ 0
  )) %>% mutate(g = case_when(
    center== 4~ 0,
    center== 1~ 1,
    center== 11~ 2,
  )) %>% 
  mutate(z = case_when(
    (z_h == 0)&(g == 0) ~ "a0",
    (z_h == 1)&(g == 0) ~ "a1",
    (z_h == 0)&(g == 1) ~ "b0",
    (z_h == 1)&(g == 1) ~ "b1",
    (z_h == 0)&(g == 2) ~ "c0",
    (z_h == 1)&(g == 2) ~ "c1",
  )) %>% 
  mutate(z_a0 = case_when(
    z == "a0" ~ 1,
    z != "a0" ~ 0
  )) %>% 
  mutate(z_a1 = case_when(
    z == "a1" ~ 1,
    z != "a1" ~ 0
  )) %>% 
  mutate(z_b0 = case_when(
    z == "b0" ~ 1,
    z != "b0" ~ 0
  )) %>% 
  mutate(z_b1 = case_when(
    z == "b1" ~ 1,
    z != "b1" ~ 0
  )) %>% mutate(colo_year = colo_exitdays/365) %>% drop_na("d", "age", "sex","race7","educat","cig_stat",
                                                          "bmi_curc",
                                                          "marital", "occupat", "colo_cancer","colo_exitdays") %>% 
  mutate(sex_m = 2-sex) %>% 
  mutate(race_m = case_when(
    race7 == 1 ~ 0,
    race7 != 1 ~ 1
  )) %>% 
  mutate(educ_m1 = case_when(
    (educat==3) |(educat==4) ~ 1,
    .default = 0
  )) %>% 
  mutate(educ_m2 =case_when(
    educat >4 ~ 1,
    .default = 0
  )) %>% 
  mutate(bmi_m =case_when(
    bmi_curc <=2 ~ 0,
    bmi_curc >2 ~ 1
  )) %>% 
  mutate(cig_m1 =case_when(
    cig_stat ==1 ~ 1,
    .default = 0
  ))%>% 
  mutate(cig_m2 =case_when(
    cig_stat ==2 ~ 1,
    .default = 0
  )) %>% 
  mutate(educ_m = case_when(
    educat <=2 ~ 0,
    (educat==3) |(educat==4) ~ 1,
    educat >4 ~ 2
  )) %>% 
  mutate(cig_m = case_when(
    cig_stat ==1 ~ 1,
    cig_stat ==2 ~ 2,
    .default = 0
  )) 


tab2 <- CreateTableOne(vars = covars, strata = "z",
                       data = dat_2, factorVars = cat_covars)
des_2 <- print(tab2, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

write.csv(des_2, "des_2.csv")
write.csv(dat_2, "dat_2.csv")
