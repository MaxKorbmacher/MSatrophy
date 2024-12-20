# Analyse GM atrophy among OFAMS clinical trial MS patients
# Max Korbmacher, 20.12.2024
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
savepath = "/Users/max/Documents/Local/MS/results/" # define res path
# 0. Prep ####
# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,lme4,lmerTest,effects,effectsize,interactions,gamm4,
               ggseg,ggtext,ggpubr,MuMIn)
# data
df = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
# join
# 1. General trend of GM vol degeneration ####
df$session = factor(df$session) # make session a factor
df$geno = factor(df$geno) # make genotype a factor
df$sex = factor(df$sex) # make sex a factor
levels(df$sex) = c("Male","Female")
df$TIV = df$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df$TotalVol = df$TotalGrayVol/1000000
# from mm^3 to dm^3 or liters
vol.m = lmer(TotalVol~session+TIV+sex+age+(1|eid),df)
summary(vol.m) 
plot(Effect("session",vol.m))
Effect("session",vol.m)
vol.m = lmer(TotalVol~TIV+sex+age+(1|eid),df)
summary(vol.m) 
plot(Effect("age",vol.m))
#plot(Effect("geno",vol.m))
#plot(Effect(c("geno","sex"),vol.m)) # seem like no interaction of sex and genotype

vol.m0 = lmer(TotalVol~session+TIV+age*sex+(1|eid),df)
interact_plot(vol.m0,pred =  age, modx = sex) # no interaction of genotype and age

# introduce a non-linear trend of age (yet similar to linear trend...)
vol.gam = gamm4(TotalVol~TIV+sex+s(age),random=~(1|eid),data=df)
plot(vol.gam$gam,pages=1)
#plot(vol.gam$mer,pages=1) # underlying mixed model
#anova(vol.gam$gam)
#
###### IMPORTANT: WE SEE AN OVERALL QUITE LINEAR DEGENERATION.
###### THAT MEANS THAT ONE CAN ALSO USE THE DATA UP TO 2 YEARS TO MODEL DEGENERATIVE PROCESSES

# 2. General association of GM vol and EDSS ####

edss.m = lmer(edss~TotalVol+sex+TIV+sex+age+(1|eid),df)
summary(edss.m) 
plot(Effect("TotalVol",edss.m),xlab = "Total brain volume",ylab="EDSS",main="Adjusted association of EDS and total volume")
effectsize::standardize_parameters(edss.m) # clear association between edss & vol
r.squaredGLMM(edss.m)

edss.m1 = lmer(edss~TotalVol+sex+TIV+age+geno+(1|eid),df)
summary(edss.m1) 
effectsize::standardize_parameters(edss.m1) # clear association between edss & vol
r.squaredGLMM(edss.m1)



edss.m2 = lmer(edss~TotalVol*sex+TIV+age+(1|eid),df)
summary(edss.m2) 
effectsize::standardize_parameters(edss.m2) # clear association between edss & vol
interact_plot(edss.m2,pred =  TotalVol, modx = sex) # interaction of sex

edss.m3 = lmer(edss~TotalVol*geno+sex+TIV+age+(1|eid),df)
summary(edss.m3) 
interact_plot(edss.m3,pred =  TotalVol, modx = geno) # interaction of genotype

# edss.m3 = lmer(edss~TotalVol*geno+geno*sex+TIV+age+(1|eid),df)
# summary(edss.m3) 
# interact_plot(edss.m3,pred =  sex, modx = geno) 
# effectsize::standardize_parameters(edss.m3)

#
r.squaredGLMM(edss.m)
r.squaredGLMM(edss.m1)
r.squaredGLMM(edss.m2)
r.squaredGLMM(edss.m3)

anova(edss.m1,edss.m3) # the two best performing models are not different from each other
# we pick hence the simpler model, edss.m1
# contains sex,age,geno,vol, no interaction effects
# 3. Regional age associations ####
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
df %>% select(starts_with("rh") & ends_with("volume")) %>% names())

reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+geno+(1|eid)"))
  mod = lmer(f,df)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = length(res$p)*res$p
plot_df = res
plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  labs(title="Regional volume loss") + 
  theme_void()
p1 = p1+labs(fill="Annual loss in mm<sup>3</sup>") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )

# 4. Regional EDSS associations ####
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+geno+(1|eid)"))
  mod = lmer(f,df)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = length(res$p)*res$p
plot_df = res
plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
p2 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  labs(title="Regional association between EDSS and volume") + 
  theme_void()
p2 = p2+labs(fill="mm<sup>3</sup> per one EDSS point") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
ggarrange(p1,p2,ncol=2)
