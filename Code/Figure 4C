# how the persistence of an invader species to the end of the 4th year 
library(openxlsx)
library(car)
library(lme4)
library(lmerTest)
library(ggplot2)
library(glmmTMB)
library(tibble)
library(dplyr)

###### 对于每盆中先加入以及后引入植物总体生物量的影响
pot_cwm_corr_data = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "pot_mass", rowNames = F, colNames = T)
colnames(pot_cwm_corr_data)
pot_cwm_corr_data$yr = as.character(pot_cwm_corr_data$Time)

pot_cwm_corr_data = subset(pot_cwm_corr_data, yr %in% c(2024, 2025))
dim(pot_cwm_corr_data)
shapiro.test(log10(pot_cwm_corr_data$E_NAG))
hist(log10(pot_cwm_corr_data$E_NAG))

colnames(pot_cwm_corr_data)
var_select <- c("E_BG", "E_NAG", "E_AP", "soilN", "X2_F_hgt","X2_F_LA","X2_F_LDMC","X2_F_germT","X2_F_flowT","X2_F_seeds","X2_F_SLA")
pot_cwm_corr_data[var_select] <- scale(pot_cwm_corr_data[var_select])

pot_cwm_corr_E_BG = subset(pot_cwm_corr_data, X2_F_LA != "NA" & E_BG != "NA")

E_BG_model <- lmer(E_BG ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                (1|yr) + 
                (1|Origin2:SR2:pot) + 
                (1|res_sp_list) + 
                (1|add_sp_list), 
              data = pot_cwm_corr_E_BG, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(E_BG_model)

#anova(model)
shapiro.test(resid(E_BG_model))
car::Anova(E_BG_model)
effectsize::effectsize(E_BG_model)

options(na.action = "na.fail")
library(MuMIn)
E_BG_dd12 <- dredge(E_BG_model, trace = 2, rank = "AICc")

de6 <- model.avg(E_BG_dd12, subset = delta < 4, fit = TRUE)
sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)

# get the coefficient and the standard error
E_BG_resultModel <- summary(object = MuMIn::model.avg(object = E_BG_dd12, 
                                                 subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

E_BG_resultModel$Parameter = rownames(E_BG_resultModel)
E_BG_resultModel = E_BG_resultModel %>% left_join(importance_df)
print(E_BG_resultModel)

E_BG_final <- lmer(E_BG ~ X2_F_hgt + X2_F_LA + X2_F_germT + X2_F_seeds + 
                     (1|yr) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_E_BG, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(E_BG_final)

MuMIn::r.squaredGLMM(E_BG_final)

################################################################################
pot_cwm_corr_E_NAG = subset(pot_cwm_corr_data, X2_F_LA != "NA" & E_NAG != "NA")

E_NAG_model <- lmer(E_NAG ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                     (1|yr) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_E_NAG, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(E_NAG_model)
#anova(model)
shapiro.test(resid(E_NAG_model))
car::Anova(E_NAG_model)
effectsize::effectsize(E_NAG_model)

options(na.action = "na.fail")
library(MuMIn)

E_NAG_dd12 <- dredge(E_NAG_model, trace = 2, rank = "AICc")

de6 <- model.avg(E_NAG_dd12, subset = delta < 4, fit = TRUE)
sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)

# get the coefficient and the standard error
E_NAG_resultModel <- summary(object = MuMIn::model.avg(object = E_NAG_dd12, 
                                                      subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

E_NAG_resultModel$Parameter = rownames(E_NAG_resultModel)
E_NAG_resultModel = E_NAG_resultModel %>% left_join(importance_df)
print(E_NAG_resultModel)

E_NAG_final <- lmer(E_NAG ~ X2_F_LA + 
                      (1|yr) + 
                      (1|Origin2:SR2:pot) + 
                      (1|res_sp_list) + 
                      (1|add_sp_list), 
                    data = pot_cwm_corr_E_NAG, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(E_NAG_final)
MuMIn::r.squaredGLMM(E_NAG_final)

################################################################################
pot_cwm_corr_E_AP = subset(pot_cwm_corr_data, X2_F_LA != "NA" & E_AP != "NA")

E_AP_model <- lmer(E_AP ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                      (1|yr) + 
                      (1|Origin2:SR2:pot) + 
                      (1|res_sp_list) + 
                      (1|add_sp_list), 
                    data = pot_cwm_corr_E_AP, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(E_AP_model)
#anova(model)
shapiro.test(resid(E_AP_model))
car::Anova(E_AP_model)
effectsize::effectsize(E_AP_model)

options(na.action = "na.fail")
library(MuMIn)

E_AP_dd12 <- dredge(E_AP_model, trace = 2, rank = "AICc")

de6 <- model.avg(E_AP_dd12, subset = delta < 4, fit = TRUE)
sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)

# get the coefficient and the standard error
E_AP_resultModel <- summary(object = MuMIn::model.avg(object = E_AP_dd12, 
                                                       subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

E_AP_resultModel$Parameter = rownames(E_AP_resultModel)
E_AP_resultModel = E_AP_resultModel %>% left_join(importance_df)
print(E_AP_resultModel)

E_AP_final <- lmer(E_AP ~ X2_F_seeds + 
                     (1|yr) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_E_AP, control = lmerControl(optimizer = "bobyqa"), REML = T)
MuMIn::r.squaredGLMM(E_AP_final)

################################################################################
colnames(pot_cwm_corr_data)
pot_cwm_corr_soilN = subset(pot_cwm_corr_data, X2_F_LA != "NA" & soilN != "NA")

soilN_model <- lmer(soilN ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                     #(1|yr) + 
                     #(1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_soilN, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(soilN_model)
#anova(model)
shapiro.test(resid(soilN_model))
car::Anova(soilN_model)
effectsize::effectsize(soilN_model)

options(na.action = "na.fail")
library(MuMIn)
soilN_dd12 <- dredge(soilN_model, trace = 2, rank = "AICc")

de6 <- model.avg(soilN_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)

# get the coefficient and the standard error
soilN_resultModel <- summary(object = MuMIn::model.avg(object = soilN_dd12, 
                                                      subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

soilN_resultModel$Parameter = rownames(soilN_resultModel)
soilN_resultModel = soilN_resultModel %>% left_join(importance_df)
print(soilN_resultModel)

soilN_model <- lmer(soilN ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                      #(1|yr) + 
                      #(1|Origin2:SR2:pot) + 
                      (1|res_sp_list) + 
                      (1|add_sp_list), 
                    data = pot_cwm_corr_soilN, control = lmerControl(optimizer = "bobyqa"), REML = T)
vif(soilN_model)

################################################################################
# Effects on litter loss estimated in situ. ((invaded vs. control)) ##
pot_mass_Linsitu = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "Linsitu", rowNames = F, colNames = T)
pot_mass_Linsitu$Origin2 = as.factor(pot_mass_Linsitu$Origin2)
pot_mass_Linsitu$SR2 = as.factor(pot_mass_Linsitu$SR2)
pot_mass_Linsitu$pot = as.factor(pot_mass_Linsitu$pot)
pot_mass_Linsitu$add_comm = as.factor(pot_mass_Linsitu$add_comm)
pot_mass_Linsitu$Res_comm = as.factor(pot_mass_Linsitu$Res_comm)

# 添加2025年群落加权值
CWM_traits_2025 = subset(pot_cwm_corr_data, yr == "2025")

pot_mass_Linsitu_add = pot_mass_Linsitu %>% left_join(CWM_traits_2025[, c("pot", predict_list)]) 

colnames(pot_mass_Linsitu_add)
var_select <- c("Mloss", "X2_F_hgt","X2_F_LA","X2_F_LDMC","X2_F_germT","X2_F_flowT","X2_F_seeds","X2_F_SLA")
pot_mass_Linsitu_add[var_select] <- scale(pot_mass_Linsitu_add[var_select])

pot_cwm_corr_Linsitu = subset(pot_mass_Linsitu_add, X2_F_LA != "NA" & Mloss != "NA")
cor.test(pot_cwm_corr_Linsitu$X2_F_hgt, pot_cwm_corr_Linsitu$X2_F_LA)

Linsitu_model <- lmer(Mloss ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_Linsitu, control = lmerControl(optimizer = "bobyqa"))
vif(Linsitu_model)
#anova(model)
shapiro.test(resid(Linsitu_model))
car::Anova(Linsitu_model)
effectsize::effectsize(Linsitu_model)

options(na.action = "na.fail")
library(MuMIn)
Linsitu_dd12 <- dredge(Linsitu_model, trace = 2, rank = "AICc")

de6 <- model.avg(Linsitu_dd12, subset = delta < 4, fit = TRUE)
sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)


# get the coefficient and the standard error
Linsitu_resultModel <- summary(object = MuMIn::model.avg(object = Linsitu_dd12, 
                                                      subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

Linsitu_resultModel$Parameter = rownames(Linsitu_resultModel)
Linsitu_resultModel = Linsitu_resultModel %>% left_join(importance_df)
print(Linsitu_resultModel)

Linsitu_final <- lmer(Mloss ~ X2_F_flowT + X2_F_seeds + 
                        (1|Origin2:SR2:pot) + 
                        (1|res_sp_list) + 
                        (1|add_sp_list), 
                      data = pot_cwm_corr_Linsitu, control = lmerControl(optimizer = "bobyqa"))
vif(Linsitu_final)
MuMIn::r.squaredGLMM(Linsitu_final)

###### 对于每盆引入、基底、整体群落生物量、病虫害水平影响
pot_cwm_corr_data = read.xlsx("field experiment cover biomass survival 260409.xlsx", sheet = "pot_mass", rowNames = F, colNames = T)
colnames(pot_cwm_corr_data)
pot_cwm_corr_data$SQRTMass = sqrt(pot_cwm_corr_data$mass)
pot_cwm_corr_data$CWM_F_Disease = sqrt(pot_cwm_corr_data$CWM_F_Disease+1)
pot_cwm_corr_data$CWM_F_Pest = sqrt(pot_cwm_corr_data$CWM_F_Pest+1)
pot_cwm_corr_data$yr = as.character(pot_cwm_corr_data$Time)

cor.test(pot_cwm_corr_data$X2_F_hgt, pot_cwm_corr_data$X2_F_LDMC)

var_select <- c("SQRTX2mass", "SQRTRmass", "SQRTMass", "CWM_F_Disease", "CWM_F_Pest", "X2_F_hgt","X2_F_LA","X2_F_LDMC","X2_F_germT","X2_F_flowT","X2_F_seeds","X2_F_SLA")
pot_cwm_corr_data[var_select] <- scale(pot_cwm_corr_data[var_select])

pot_cwm_corr_SQRTX2mass = subset(pot_cwm_corr_data, X2_F_LA != "NA" & SQRTX2mass != "NA")

SQRTX2mass_model <- lmer(SQRTX2mass ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                     (1|yr) + 
                     (1|Origin2:SR2:pot) + 
                     (1|res_sp_list) + 
                     (1|add_sp_list), 
                   data = pot_cwm_corr_SQRTX2mass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTX2mass_model)
#anova(model)
shapiro.test(resid(SQRTX2mass_model))
car::Anova(SQRTX2mass_model)
effectsize::effectsize(SQRTX2mass_model)

options(na.action = "na.fail")
library(MuMIn)
SQRTX2mass_dd12 <- dredge(SQRTX2mass_model, trace = 2, rank = "AICc")

Final_model <- get.models(SQRTX2mass_dd12,1)[[1]] 
summary(Final_model)
vif(Final_model)

de6 <- model.avg(SQRTX2mass_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)

# get the coefficient and the standard error
SQRTX2mass_resultModel <- summary(object = MuMIn::model.avg(object = SQRTX2mass_dd12, 
                                                      subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

SQRTX2mass_resultModel$Parameter = rownames(SQRTX2mass_resultModel)
SQRTX2mass_resultModel = SQRTX2mass_resultModel %>% left_join(importance_df)
print(SQRTX2mass_resultModel)


SQRTX2mass_final <- lmer(SQRTX2mass ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_flowT + X2_F_seeds + 
                           (1|yr) + 
                           (1|Origin2:SR2:pot) + 
                           (1|res_sp_list) + 
                           (1|add_sp_list), 
                         data = pot_cwm_corr_SQRTX2mass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTX2mass_final)
MuMIn::r.squaredGLMM(SQRTX2mass_final)

################################################################################
pot_cwm_corr_SQRTRmass = subset(pot_cwm_corr_data, X2_F_LA != "NA" & SQRTRmass != "NA")

SQRTRmass_model <- lmer(SQRTRmass ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                           (1|yr) + 
                           (1|Origin2:SR2:pot) + 
                           (1|res_sp_list) + 
                           (1|add_sp_list), 
                         data = pot_cwm_corr_SQRTRmass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTRmass_model)
MuMIn::r.squaredGLMM(SQRTRmass_model)
#anova(model)
shapiro.test(resid(SQRTRmass_model))
car::Anova(SQRTRmass_model)
summary(SQRTRmass_model)
effectsize::effectsize(SQRTRmass_model)

options(na.action = "na.fail")
library(MuMIn)
SQRTRmass_dd12 <- dredge(SQRTRmass_model, trace = 2, rank = "AICc")

Final_model <- get.models(SQRTRmass_dd12,1)[[1]] 
summary(Final_model)
vif(Final_model)

de6 <- model.avg(SQRTRmass_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)


# get the coefficient and the standard error
SQRTRmass_resultModel <- summary(object = MuMIn::model.avg(object = SQRTRmass_dd12, 
                                                            subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

SQRTRmass_resultModel$Parameter = rownames(SQRTRmass_resultModel)
SQRTRmass_resultModel = SQRTRmass_resultModel %>% left_join(importance_df)
print(SQRTRmass_resultModel)

SQRTRmass_final <- lmer(SQRTRmass ~ X2_F_hgt + 
                          (1|yr) + 
                          (1|Origin2:SR2:pot) + 
                          (1|res_sp_list) + 
                          (1|add_sp_list), 
                        data = pot_cwm_corr_SQRTRmass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTRmass_model)
MuMIn::r.squaredGLMM(SQRTRmass_model)

################################################################################
pot_cwm_corr_SQRTMass = subset(pot_cwm_corr_data, X2_F_LA != "NA" & SQRTMass != "NA")

SQRTMass_model <- lmer(SQRTMass ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                          (1|yr) + 
                          (1|Origin2:SR2:pot) + 
                          (1|res_sp_list) + 
                          (1|add_sp_list), 
                        data = pot_cwm_corr_SQRTMass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTMass_model)
MuMIn::r.squaredGLMM(SQRTMass_model)
#anova(model)
shapiro.test(resid(SQRTMass_model))
car::Anova(SQRTMass_model)
summary(SQRTMass_model)
effectsize::effectsize(SQRTMass_model)

options(na.action = "na.fail")
library(MuMIn)
SQRTMass_dd12 <- dredge(SQRTMass_model, trace = 2, rank = "AICc")

Final_model <- get.models(SQRTMass_dd12,1)[[1]] 
summary(Final_model)
vif(Final_model)

de6 <- model.avg(SQRTMass_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)


# get the coefficient and the standard error
SQRTMass_resultModel <- summary(object = MuMIn::model.avg(object = SQRTMass_dd12, 
                                                           subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

SQRTMass_resultModel$Parameter = rownames(SQRTMass_resultModel)
SQRTMass_resultModel = SQRTMass_resultModel %>% left_join(importance_df)
print(SQRTMass_resultModel)


SQRTMass_final <- lmer(SQRTMass ~ X2_F_LDMC + X2_F_flowT + X2_F_seeds + 
                         (1|yr) + 
                         (1|Origin2:SR2:pot) + 
                         (1|res_sp_list) + 
                         (1|add_sp_list), 
                       data = pot_cwm_corr_SQRTMass, control = lmerControl(optimizer = "bobyqa"))
vif(SQRTMass_final)
MuMIn::r.squaredGLMM(SQRTMass_final)


################################################################################
pot_cwm_corr_CWM_F_Disease = subset(pot_cwm_corr_data, X2_F_LA != "NA" & CWM_F_Disease != "NA")

CWM_F_Disease_model <- lmer(CWM_F_Disease ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                         (1|yr) + 
                         (1|Origin2:SR2:pot) + 
                         (1|res_sp_list) + 
                         (1|add_sp_list), 
                       data = pot_cwm_corr_CWM_F_Disease, control = lmerControl(optimizer = "bobyqa"))
vif(CWM_F_Disease_model)
#anova(model)
shapiro.test(resid(CWM_F_Disease_model))
car::Anova(CWM_F_Disease_model)
summary(CWM_F_Disease_model)
effectsize::effectsize(CWM_F_Disease_model)

options(na.action = "na.fail")
library(MuMIn)
CWM_F_Disease_dd12 <- dredge(CWM_F_Disease_model, trace = 2, rank = "AICc")

Final_model <- get.models(CWM_F_Disease_dd12,1)[[1]] 
summary(Final_model)
vif(Final_model)

de6 <- model.avg(CWM_F_Disease_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)


# get the coefficient and the standard error
CWM_F_Disease_resultModel <- summary(object = MuMIn::model.avg(object = CWM_F_Disease_dd12, 
                                                          subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

CWM_F_Disease_resultModel$Parameter = rownames(CWM_F_Disease_resultModel)
CWM_F_Disease_resultModel = CWM_F_Disease_resultModel %>% left_join(importance_df)
print(CWM_F_Disease_resultModel)

CWM_F_Disease_final <- lmer(CWM_F_Disease ~ X2_F_flowT + 
                              (1|yr) + 
                              (1|Origin2:SR2:pot) + 
                              (1|res_sp_list) + 
                              (1|add_sp_list), 
                            data = pot_cwm_corr_CWM_F_Disease, control = lmerControl(optimizer = "bobyqa"))
MuMIn::r.squaredGLMM(CWM_F_Disease_final)

################################################################################
pot_cwm_corr_CWM_F_Pest = subset(pot_cwm_corr_data, X2_F_LA != "NA" & CWM_F_Pest != "NA")

CWM_F_Pest_model <- lmer(CWM_F_Pest ~ X2_F_hgt + X2_F_LA + X2_F_LDMC + X2_F_germT + X2_F_flowT + X2_F_seeds + 
                              (1|yr) + 
                              (1|Origin2:SR2:pot) + 
                              (1|res_sp_list) + 
                              (1|add_sp_list), 
                            data = pot_cwm_corr_CWM_F_Pest, control = lmerControl(optimizer = "bobyqa"))
vif(CWM_F_Pest_model)
#anova(model)
shapiro.test(resid(CWM_F_Pest_model))
car::Anova(CWM_F_Pest_model)
summary(CWM_F_Pest_model)
effectsize::effectsize(CWM_F_Pest_model)

options(na.action = "na.fail")
library(MuMIn)
CWM_F_Pest_dd12 <- dredge(CWM_F_Pest_model, trace = 2, rank = "AICc")

Final_model <- get.models(CWM_F_Pest_dd12,1)[[1]] 
summary(Final_model)
vif(Final_model)

de6 <- model.avg(CWM_F_Pest_dd12, subset = delta < 4, fit = TRUE)
de6$msTable
print(de6$msTable[1,])

sw_unclass <- unclass(sw(de6))
weights <- sw_unclass
n_models <- attr(sw_unclass, "n.models")
variables <- names(sw_unclass)
importance_df <- data.frame(
  Parameter = variables,
  Sum_of_weights = as.numeric(weights),
  N_containing_models = as.numeric(n_models))
print(importance_df)


# get the coefficient and the standard error
CWM_F_Pest_resultModel <- summary(object = MuMIn::model.avg(object = CWM_F_Pest_dd12, 
                                                               subset = delta < 4))$coefmat.subset %>% 
  as.data.frame()

CWM_F_Pest_resultModel$Parameter = rownames(CWM_F_Pest_resultModel)
CWM_F_Pest_resultModel = CWM_F_Pest_resultModel %>% left_join(importance_df)
print(CWM_F_Pest_resultModel)

CWM_F_Pest_final <- lmer(CWM_F_Pest ~ X2_F_hgt + X2_F_LA + X2_F_flowT + X2_F_seeds + 
                           (1|yr) + 
                           (1|Origin2:SR2:pot) + 
                           (1|res_sp_list) + 
                           (1|add_sp_list), 
                         data = pot_cwm_corr_CWM_F_Pest, control = lmerControl(optimizer = "bobyqa"))
vif(CWM_F_Pest_final)
MuMIn::r.squaredGLMM(CWM_F_Pest_final)


################################################################################
# plot
plot_data <- read.xlsx("Figure_4c_model_select.xlsx", sheet = "delta<4", colNames = T)


library(tidyverse)

# 定义标准预测变量
standard_predictors <- c("X2_F_LA", "X2_F_hgt", "X2_F_LDMC", "X2_F_germT", "X2_F_flowT", "X2_F_seeds")

# 按响应变量和预测变量重塑数据
df_standard <- plot_data %>%
  # 筛选出有 Estimate 的行（即非截距行，或根据需要调整）
  filter(Predictors %in% standard_predictors | Predictors == "(Intercept)") %>%
  # 为每个响应变量创建完整的预测变量组合
  complete(Responses, Predictors = c("(Intercept)", standard_predictors)) %>%
  # 按原始顺序排列
  arrange(Responses, match(Predictors, c("(Intercept)", standard_predictors))) %>%
  # 可选：重新组织列的顺序
  select(Responses, Predictors, Estimate, Std..Error, z.value, `Pr(>|z|)`, everything())


# re-name
df_standard$Responses[df_standard$Responses == "SQRTX2mass"] <- "Introduced assemblage biomass"
df_standard$Responses[df_standard$Responses == "SQRTRmass"] <- "Receipt native community biomass"
df_standard$Responses[df_standard$Responses == "SQRTMass"] <- "Whole community biomass"
df_standard$Responses[df_standard$Responses == "CWM_F_Disease"] <- "CWM of leaf infection infection"
df_standard$Responses[df_standard$Responses == "CWM_F_Pest"] <- "CWM of leaf defoliation"
df_standard$Responses[df_standard$Responses == "E_AP"] <- "Acid phosphatase activity"
df_standard$Responses[df_standard$Responses == "E_BG"] <- "β-1,4-glucosidase activity"
df_standard$Responses[df_standard$Responses == "E_NAG"] <- "β-1,4-N-acetylglucosaminidase activity"
df_standard$Responses[df_standard$Responses == "soilN"] <- "Soil TN content"
df_standard$Responses[df_standard$Responses == "Linsitu"] <- "Litter decomposition rate in field"

# re-name
unique(df_standard$Predictors)
df_standard$Predictors[df_standard$Predictors == "X2_F_hgt"] <- "Height"
df_standard$Predictors[df_standard$Predictors == "X2_F_LA"] <- "LA"
df_standard$Predictors[df_standard$Predictors == "X2_F_LDMC"] <- "LDMC"
df_standard$Predictors[df_standard$Predictors == "X2_F_germT"] <- "Germ data"
df_standard$Predictors[df_standard$Predictors == "X2_F_flowT"] <- "Flower data"
df_standard$Predictors[df_standard$Predictors == "X2_F_seeds"] <- "# of seeds"

df_standard <- subset(df_standard, Predictors != "(Intercept)")
df_standard$Responses = factor(df_standard$Responses, levels = (c("Soil TN content", "Litter decomposition rate in field",
                                                                  "Acid phosphatase activity", "β-1,4-N-acetylglucosaminidase activity", "β-1,4-glucosidase activity", 
                                                                  "CWM of leaf defoliation", "CWM of leaf infection infection",
                                                                  "Whole community biomass", "Receipt native community biomass", 
                                                                  "Introduced assemblage biomass")))

df_standard$Predictors = factor(df_standard$Predictors, levels = (c("Height", "LA", "LDMC", "Slope", 
                                                                       "Germ data", "Flower data", "# of seeds")))

# adding signal
df_standard$Significance <- case_when(
  df_standard$`Pr(>|z|)` < 0.001 ~ "***",
  df_standard$`Pr(>|z|)` < 0.01  ~ "**",
  df_standard$`Pr(>|z|)` < 0.05  ~ "*",
  #df_standard$`Pr(>|z|)` < 0.1   ~ "†",  
  TRUE ~ "")


df_standard <- df_standard %>%
  mutate(
    Std_Coefficient_removed = as.numeric(ifelse(`Pr(>|z|)` <= 0.1, Estimate, 0)),
    label_text = ifelse(`Pr(>|z|)` <= 0.90,
                        sprintf("%.2f\n%s", Estimate, Significance),NA))


ggplot(df_standard, aes(x = Predictors, y = Responses, fill = Std_Coefficient_removed)) +
  geom_tile(color = "#F8F8F8", linewidth = 0.5, na.rm = TRUE) + 
  geom_text(aes(label = label_text), 
            color = "black", size = 4, na.rm = TRUE) + 
  scale_fill_gradientn(limits = c(-0.6, 0.65), na.value = "#F8F8F8", name = NULL, 
                       colours = (RColorBrewer::brewer.pal(11,"RdBu"))) +
  theme_bw() + 
  #scale_y_discrete(expand = c(0, 0), labels = x_label_mapping) +
  scale_x_discrete(expand = c(0, 0), labels = y_label_mapping) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        panel.border = element_blank(),
        strip.background = element_rect(color="white", fill="white", size=0.5, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"),
        legend.position = "right") + # bottom
  labs(x = NULL, y = NULL,  fill = "Coefficient") +
  guides(fill = guide_colorbar( title.position = "top", title.hjust = 0.5,barwidth = 1, barheight = 14))

# 11.54 X 4.42

