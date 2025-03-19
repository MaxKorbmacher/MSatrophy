# Play with multivariate Gaussian process regression (in GAMLSS)
# Max Korbmacher, 17 Feb 2025
# last change: 14 March 2025
#
# I short disclaimer on the harmonization strategy:
# Cross-sectional HC data were harmonized separately from longitudinal MS data to avoid data leakage.
#
#
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
#
# PREP ####
#remotes::install_github("pvarvia/mgpr", ref = "main")
#library(mgpr)
#install.packages("zeallot")
# library(recipes)
# library(workflows)
# library(glmnet)
# library(zeallot)
library(ggplot2)
library(gridExtra)
library(gamlss)
library(gamlss.ggplots)
library(dplyr)
library(neuroCombat)
library(longCombat)
library(tidyverse) 
library(mgcv)
library(lme4)
library(lmerTest)
library(lmtest)
library(pscl)
library(MKmisc)
library(survey)
library(MatchIt)
library(effects)
# get data
dt = read.csv("/Users/max/Documents/Local/Data/Lifespan/long_thalamus.csv")
msBGO = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
demo_OSL = read.csv("/Users/max/Documents/Local/Data/Oslo/Oslo_demographics.csv")
oasis = read.delim("/Users/max/Documents/Local/Data/Lifespan/oasis_subcortical.txt")
oasis_demo = read.csv("/Users/max/Documents/Local/Data/Lifespan/oasis_demo.csv")
# prep dat
oasis$eid = factor(gsub("/","",oasis$Measure.volume))
names(oasis_demo)[names(oasis_demo) == "ID"] = "eid"
names(oasis_demo)[names(oasis_demo) == "M.F"] = "sex"
oasis_demo$sex = ifelse(oasis_demo$sex == "F", 0, 1)
names(oasis_demo)[names(oasis_demo) == "Age"] = "age"
oasis = merge(oasis,oasis_demo, by= "eid")
dt$session = ifelse(dt$session == "ADNI Screening",1,dt$session)
dt$session = ifelse(dt$session == "tp2",1,dt$session)
dt$session = ifelse(dt$session == "m00",1,dt$session)
dt$session = as.numeric(dt$session)
dt$scanner = as.numeric(factor(dt$scanner))
dt$sex = ifelse(dt$sex == 0,"Female",dt$sex)
dt$sex = ifelse(dt$sex == 1,"Male",dt$sex)
dt$sex = ifelse(dt$sex == "Female", 0,dt$sex)
dt$sex = ifelse(dt$sex == "F", 0,dt$sex)
dt$sex = ifelse(dt$sex != 0, 1,dt$sex)
hc = dt%>%filter(diagnosis == "HC")
hc = hc%>% filter(session < 2)
hc=na.omit(hc)
# 
# NOTE! What is treated here as the scanner in the Bergen data (i.e., the OFAMS clinical trial) 
# data is actually the hospital, not the scanner type/model.
msOSL = dt%>%filter(diagnosis == "MS") %>% select(eid, scanner ,session,age,sex,Left.Thalamus,Right.Thalamus,EstimatedTotalIntraCranialVol)
msBGO = msBGO %>% select(eid,session,age,sex,Left.Thalamus,Right.Thalamus,EstimatedTotalIntraCranialVol)
msBGO$scanner = ifelse(msBGO$eid > 100 & msBGO$eid < 200, "a", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 200 & msBGO$eid < 300, "b", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 300 & msBGO$eid < 400, "c", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 400 & msBGO$eid < 500, "d", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 500 & msBGO$eid < 600, "e", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 600 & msBGO$eid < 700, "f", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 700 & msBGO$eid < 800, "g", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 800 & msBGO$eid < 900, "h", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 900 & msBGO$eid < 1000, "i", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1000 & msBGO$eid < 1100, "j", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1100 & msBGO$eid < 1200, "k", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1200 & msBGO$eid < 1300, "l", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1300 & msBGO$eid < 1400, "m", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1400 & msBGO$eid < 1500, "n", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1500 & msBGO$eid < 1600, "o", msBGO$eid)
msBGO$scanner = ifelse(msBGO$eid > 1600, "p", msBGO$eid)
ms = rbind(msOSL, msBGO)
ms$thalamus = ms$Left.Thalamus + ms$Right.Thalamus
ms$sex = factor(ms$sex)
ms$eid = factor(ms$eid)
#
# check age groups
# hc$AgeGroup <- cut(hc$age, 
#                    breaks = c(-Inf
#                               ,5 ,10 ,15,20,25,30,35,40,45,50,55,60 ,65,70,75,80,85
#                               , Inf), 
#                    
#                    labels = c("0-4 years"
#                               ,"5-9 years","10-14 years","15-19 years","20-24 years"
#                               ,"25-29 years","30-34 years","35-39 years","40-44 years"
#                               ,"45-49 years","50-54 years","55-59 years","60-64 years"
#                               ,"65-69 years","70-74 years","75-79 years","80-84 years"
#                               ,"85+ years"),
#                    right = FALSE)
# table(hc$AgeGroup)
#
# harmonise data
hc$eid = factor(hc$eid)
hc$sex = factor(hc$sex)
hc$data = factor(hc$data)
covars = hc %>% dplyr::select(eid,sex,scanner,age)
hc = neuroCombat(t(hc%>%dplyr::select(EstimatedTotalIntraCranialVol, Left.Thalamus, Right.Thalamus)),batch=hc$scanner,mod=model.matrix(~covars$age+covars$sex), mean.only = T)
hc = data.frame(t(hc$dat.combat))
hc = cbind(covars,hc)
# viz thal x age
#plot(hc$age,hc$thalamus)

# full thalamic vol
hc$thalamus = hc$Right.Thalamus+hc$Left.Thalamus
#
# split by sex
M = hc%>%filter(sex=="1")
Fem = hc%>%filter(sex=="0")
#
# FIT MODELS ####
# fit the model using the Lambda Mu and Sigma (LMS) method
#fit = gamlss(thalamus ~ cs(age)+sex+scanner, sigma.formula = ~ cs(age), data = hc, family = BCTo)
#fitM = gamlss(thalamus ~ cs(age), sigma.formula = ~ cs(age), data = M, family = NO) # NO family is normal
#fitF = gamlss(thalamus ~ cs(age), sigma.formula = ~ cs(age), data = Fem)
#fit=lms(thalamus, age, data = hc)
#
#
#
#
#
# Models were fit and from here just loaded
# fitM=lms(thalamus, age, data = M)
# fitF=lms(thalamus, age, data = Fem)
# save models
# saveRDS(fitM, "/Users/max/Documents/Local/MS/GMresults/male_model.rds")
# saveRDS(fitF, "/Users/max/Documents/Local/MS/GMresults/female_model.rds")
# load models
fitM = readRDS("/Users/max/Documents/Local/MS/GMresults/male_model.rds")
fitF = readRDS("/Users/max/Documents/Local/MS/GMresults/female_model.rds")
# run diagnostics / sanity checks
# fittedPlot(fitF, x=Fem$age)
# plot(fitM)
# wp(fitM,xlim.all = 4,xvar = age,show.given = T)
# dtop(fitM, xvar = M$age, n.inter = 9)
# qstats = Q.stats(fitM, xvar = hc$age, n.inter = 9)
# print(qstats, digits = 3)
#
plot(M$age,fitted(fitM))
plot(Fem$age,fitted(fitF))
#
#plot(hc$age, qNO(pBCTo(hc$thalamus,mu,sigma,nu,tau)))
#
# plot centiles
#centiles.split(fitM)
#centiles(fitM,xvar=age, cent=c(90), points=T)
#centiles(fit,plot = T,legend = F)
centiles(fitM,plot = T,legend = F)
centiles(fitF,plot = T,legend = F)
# centile fans
# centiles.fan(fitM,M$age,cent=c(1,2.5,10,25,50,75,90,97.5,99), 
#              ylab="Thalamus volume", xlab="Age",median = T)
#
# NOTE THAT ALL THE MODELS USE THE Box-Cox t distribution
# centiles.fan(fit,hc$age,cent=c(0.5,2.5,50,97.5,99.5), 
#              ylab="Thalamus volume", xlab="Age",
#              median = T,ylim = c(0,30000))
pdf(file='/Users/max/Documents/Local/MS/GMresults/CentilePlots.pdf', width = 10, height = 5)  # Set up PDF output
par(mfrow = c(1, 2))  # Arrange plots in 1 row, 2 columns
# Plot for Males
centiles.fan(fitM, M$age, cent = c(0.5, 2.5, 50, 97.5, 99.5), 
             ylab = "Males' Thalamus Volume", xlab = "Age",
             median = TRUE, ylim = c(0, 30000), 
             main = "Centile Curves for Males")
# Plot for Females
centiles.fan(fitF, Fem$age, cent = c(0.5, 2.5, 50, 97.5, 99.5), 
             ylab = "Females' Thalamus Volume", xlab = "Age",
             median = TRUE, ylim = c(0, 30000), 
             main = "Centile Curves for Females")
dev.off()  # Close PDF device
#
# estimate Z-values
zest = function(model, data){
  # estimate the fitted paramters for BCTo
  mu = predict(model, newdata = data, type = "response", what = "mu")
  sigma = predict(model, newdata = data, type = "response", what = "sigma")
  nu = predict(model, newdata = data, type = "response", what = "nu")
  tau = predict(model, newdata = data, type = "response", what = "tau")
  # estimate Z-scores
  qNO(pBCT(data$thalamus,mu,sigma,nu,tau))
}
Fem$Z = zest(fitF,Fem)
M$Z = zest(fitM,M)

# for predictions in MS data first remove NA in thalamus
ms = ms[!is.na(ms$thalamus),]
ms = ms[!is.na(ms$sex),]
ms = ms[!is.na(ms$age),]
#
# cross-sectional comparisons of MS data with OASIS data
ms_cross = ms %>% filter(session == 1) %>% select(eid,age,sex,thalamus,EstimatedTotalIntraCranialVol,scanner)

covars = ms_cross %>% dplyr::select(eid,sex,age)
ms_cross = neuroCombat(t(ms_cross%>%dplyr::select(EstimatedTotalIntraCranialVol, thalamus)),batch=ms_cross$scanner,mod=model.matrix(~covars$age+covars$sex), mean.only = T)
ms_cross = data.frame(t(ms_cross$dat.combat))
ms_cross = cbind(covars,ms_cross)
ms_cross$diagnosis = "MS"
oasis$thalamus = oasis$Right.Thalamus.Proper+oasis$Left.Thalamus.Proper
oasis$EstimatedTotalIntraCranialVol = oasis$eTIV
oasis_cross = oasis %>% select(eid,age,sex,thalamus,EstimatedTotalIntraCranialVol)
oasis_cross$diagnosis = "HC"
oasis_cross$scanner = 112233
oasis_cross = oasis_cross[!grepl("_MR2",oasis_cross$eid),] # exclude follow-ups
names(oasis_cross)
cross = rbind(ms_cross, oasis_cross)
# match
test = matchit(as.factor(diagnosis) ~ age + sex, data = cross,
        method = "nearest", distance = "glm")
cross = match.data(test)
cross = cross %>% select(-c(distance, weights, subclass))
#



hist((cross %>% filter(diagnosis == "HC"))$thalamus)
hist((cross %>% filter(diagnosis == "MS"))$thalamus)

# Fc = cross %>% filter(sex == 0)
# Mc = cross %>% filter(sex == 1)
# Mc$Z = zest(fitM, Mc)
# Fc$Z = zest(fitF, Fc)
# cross = rbind(Fc,Mc)


cross %>% group_by(diagnosis)%>%summarise(C = cor(thalamus,age),Mean = mean(thalamus))
cross %>% group_by(diagnosis)%>%summarise(Correlation = cor(Z,age),Mean = mean(Z))
#
#
# check differences
m_c = lm(thalamus~sex+age*diagnosis, data = cross)
m_c = lm(Z~sex+age*diagnosis, data = cross)
summary(m_c)
plot(Effect(c("age","diagnosis"),m_c))
#
#
#
# 
# get longitudinal data from HCs
dat = na.omit(dt)
dat$thalamus = dat$Left.Thalamus + dat$Right.Thalamus
dat = dat %>% filter(diagnosis == "HC")
table(dat$session)

dat %>% filter(data == "Rockland")

dat = dat %>% select(-data)
ms$diagnosis = "MS"
dat = rbind(dat, ms)
#
# harmonise data once again (fresh start, not double harmonized!)
# this is done to account for the longitudinal data!
# first, exclude cases for which data are available for a single subject only at the specific site
dat =  filter(dat, !scanner %in% c(12,45,46,1201))
# longitudinal Combat
covars = dat %>% dplyr::select(sex,age, diagnosis)
dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
           features = c("EstimatedTotalIntraCranialVol", "thalamus"),
           formula = "age * sex + diagnosis", ranef = "(1|eid)", data = dat)
dat = dat$data_combat
dat = cbind(covars,dat)
#
# We can also add the Z-scores to the data frames and check whether they add to predictions
Males = dat %>% dplyr::filter(sex == 1)
Females = dat %>% filter(sex == 0)
Males$Z = zest(fitM, Males)
Females$Z = zest(fitF, Females)
dat = rbind(Males, Females)


dat %>% filter(diagnosis=="")

#
#
# CLASSIFY ####
#
# Now we can check what is contributing most in predicting diagnosis
m <- glmer(factor(diagnosis) ~ Z + thalamus.combat + EstimatedTotalIntraCranialVol.combat + 
             age + factor(sex) + (1 | eid), data = dat, 
           family = binomial)#, control = glmerControl(optimizer = "bobyqa"),nAGQ = 10)
m2 <- glmer(factor(diagnosis) ~ Z + 
             age + factor(sex) + (1 | eid), data = dat, 
           family = binomial)#, control = glmerControl(optimizer = "bobyqa"),nAGQ = 10)
m3 <- glmer(factor(diagnosis) ~ 
              age + factor(sex) + (1 | eid), data = dat, 
            family = binomial)
dat$pred = ifelse(predict(m,type = "response") >0.5,"MS","HC")
dat$pred2 = ifelse(predict(m2,type = "response") >0.5,"MS","HC")
dat$pred3 = ifelse(predict(m3,type = "response") >0.5,"MS","HC")
table(dat$pred, dat$diagnosis)
table(dat$pred2, dat$diagnosis)
table(dat$pred3, dat$diagnosis)
# conclusion: the resulting perfect predictions are potentially due to the number of time points.
# We try cross-sectional predictions:
cross = dat %>% filter(session == 1)
cross$diagnosis = factor(cross$diagnosis)

m4 = glm(diagnosis ~ Z + thalamus.combat + EstimatedTotalIntraCranialVol.combat + 
             age + factor(sex), data = cross, 
           family = binomial)
m5 = glm(diagnosis ~ thalamus.combat + EstimatedTotalIntraCranialVol.combat + 
           age + factor(sex), data = cross, 
         family = binomial)
m6 = glm(diagnosis ~ thalamus.combat + EstimatedTotalIntraCranialVol.combat + 
           age + factor(sex), data = cross, 
         family = binomial)
m7 = glm(diagnosis ~ Z + EstimatedTotalIntraCranialVol.combat + 
           age + factor(sex), data = cross, 
         family = binomial)
m8 = glm(diagnosis ~ Z + thalamus.combat + 
           age + factor(sex), data = cross, 
         family = binomial)
anova(m4,m5,m6,m7,m8, test = "Chisq") # first model is best
lrtest(m4,m5,m6,m7,m8) # first model is best
summary(m4)
pR2(m4)
pR2(m6) # but the difference to the other models is small
cross$pred = factor(ifelse(fitted(m4) > 0.5,"MS","HC"))
HLgof.test(fit = cross$pred, obs = cross$diagnosis)
regTermTest(m4, "Z")
regTermTest(m4, "thalamus.combat")
regTermTest(m4, "EstimatedTotalIntraCranialVol.combat")
regTermTest(m4, "age")
regTermTest(m4, "factor(sex)")
#
table(cross$pred, cross$diagnosis)
#
# LMER ####
m7 <- lmer(Z ~ diagnosis * age +  thalamus.combat + EstimatedTotalIntraCranialVol.combat + 
             factor(sex) + (1 | eid), data = dat)
summary(m7)
effectsize::standardize_parameters(m7)

m7 <- lmer(thalamus.combat ~ diagnosis + EstimatedTotalIntraCranialVol.combat + 
             age + factor(sex) + (1 | eid), data = dat)
summary(m7)

#
# MORE DESCRIPTIVES ####
#
# get centiles
FemCent = (centiles.pred(fitF, xname="age", xvalues= Fem$age, cent= c(5,25,50,75,95)))
MCent = (centiles.pred(fitF, xname="age", xvalues= M$age, cent= c(5,25,50,75,95)))
#
# you can get regular quantiles as well (but not useful here)
# get_quantiles_vals = function(data, vector_of_quantiles){
# # create a set of functions to estimate different percentiles
# p <- vector_of_quantiles
# p_names <- map_chr(p, ~paste0(.x*100, "%"))
# p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
#   set_names(nm = p_names)
# # apply the functions to data
# quantiles = data%>%mutate(agefr=floor(age))%>%group_by(agefr)%>%summarise_at(vars(thalamus), funs(!!!p_funs))
# return(quantiles)
# }
# #get_quantiles_vals(hc,c(.1, 0.25, 0.5, 0.75, 0.9, 0.99))
# Fquant = get_quantiles_vals(Fem,c(0.005, .025, 0.5, 0.975, 0.995))
# Mquant = get_quantiles_vals(M,c(0.005, .025, 0.5, 0.975, 0.995))
#
#
# Plot MS patient values
p1 = fitted_centiles(fitM,cent = c(0.5,2.5,50,97.5,99.5), points = F,line.col=c("grey","blue","black","blue","grey"))
pdat = p1$data
ggp = ggplot(pdat, aes(x=x)) +
  geom_line(aes(y=c.0.5)) +
  geom_line(aes(y=c.2.5)) +
  geom_line(aes(y=c.50)) +
  geom_line(aes(y=c.97.5)) +
  geom_line(aes(y=c.99.5)) +
  geom_ribbon(aes(x = x,ymin = c.0.5,ymax = c.2.5),fill = "#1b98e0", alpha=0.4)+
  geom_ribbon(aes(x = x,ymin = c.2.5,ymax = c.50),fill = "blue", alpha=0.4)+
  geom_ribbon(aes(x = x,ymin = c.50,ymax = c.97.5),fill = "blue", alpha=0.4)+
  geom_ribbon(aes(x = x,ymin = c.97.5,ymax = c.99.5),fill = "#1b98e0", alpha=0.4)+
  theme_bw() +ylab("Thalamus volume") + xlab("Age") + ylim(0,30000)+
  ggtitle("Males' Thalamic Volumes in MS")
# msgam = gam(thalamus~s(age,k=4, bs="cr")+sex+EstimatedTotalIntraCranialVol+s(eid, bs = "re"), data = ms)
# summary(msgam)
# plot(msgam)
# Mlmer = lmer(thalamus~age+factor(sex)+EstimatedTotalIntraCranialVol+(1|eid), data = ms)
# summary(Mlmer)


# Mgam = gam(thalamus~s(age,k=4, bs="cr")+EstimatedTotalIntraCranialVol+s(eid, bs = "re"), data = M)
# saveRDS(Mgam, "/Users/max/Documents/Local/MS/GMresults/male_gam.rds")
# Fgam = gam(thalamus~s(age,k=4, bs="cr")+EstimatedTotalIntraCranialVol+s(eid, bs = "re"), data = Fem)
# saveRDS(Fgam, "/Users/max/Documents/Local/MS/GMresults/female_gam.rds")

ggp + geom_point(data = ms%>%filter(sex=="1"), aes(x=age,y=thalamus),color = "black", alpha = 0.5) +
  geom_line(data = ms%>%filter(sex=="1"), aes(x=age,y=thalamus,group = eid),color = "black",alpha = 0.5) +
  geom_line(data = Mgam, aes(x=age,y=preds,group = eid),color = "black",alpha = 0.5)
#
#
# For future modelling:
# 1) load model
#my_model <- readRDS("model.rds").
# 2) make function to obtain Z-scores
zest = function(model, data){
  # estimate the fitted paramters for BCTo
  mu = predict(model, newdata = data, type = "response", what = "mu")
  sigma = predict(model, newdata = data, type = "response", what = "sigma")
  nu = predict(model, newdata = data, type = "response", what = "nu")
  tau = predict(model, newdata = data, type = "response", what = "tau")
  # estimate Z-scores
  qNO(pBCT(data$thalamus,mu,sigma,nu,tau))
}
# 3) for females, use model for females, and vice versa
# e.g., zest(Model_Females,Female_Persons_Data)
