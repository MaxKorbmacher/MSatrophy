# Analyse GM atrophy among people with multiple sclerosis in the OFAMS (Norwegian multicenter) clinical trial and local data from the Oslo University Hospital
# Max Korbmacher, 20.12.2024
# last change: 14 March 2025
#
# Note: the data were NOT harmonized across datasets, but across scanners within datasets when attempting to replicate the findings from one dataset to another.
#       This was done to provide unbiased replications between Oslo and Bergen data.
#       However, the datasets were harmonized together when analysing them together!
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
savepath = "/Users/max/Documents/Local/MS/GMresults/" # define results/save/oputput path
#
################################### #
# Contents:
################################### #
# 0. Prep
# 1. General trend of GM vol degeneration
# 2. General association of GM vol and EDSS
# 3. Regional age associations
  # 3.1 Unstandardized
    # 3.1.1 Cortical
    # 3.1.1 SubCortical
  # 3.2 Standardized
    # 3.2.1 Cortical
    # 3.2.1 SubCortical
# 4. Regional EDSS associations
  # 4.1 Unstandardized
    # 4.1.1 Cortical
    # 4.1.1 SubCortical
  # 4.2 Standardized
    # 4.2.1 Cortical
    # 4.2.1 SubCortical
# 5. Create outputs
  # 5.1 Merge plots
  # 5.2 Save tables
################################### #
#
# 0. Prep ####
# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,lme4,lmerTest,effects,effectsize,interactions,gamm4,
               ggseg,ggtext,ggpubr,MuMIn,dplyr,ggplot2,standardize)
# data
lifespan = read.csv("/Users/max/Documents/Local/Data/Lifespan/cortical_subcortical.csv")
df = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
msOSL = read.csv("/Users/max/Documents/Local/Data/Oslo/MS_Einar.csv")
subOSL = read.delim("/Users/max/Documents/Local/Data/Lifespan/MS_subcortical.txt")
demoOSL = read.csv("/Users/max/Documents/Local/Data/Oslo/Oslo_demographics.csv")
# standardize and join datasets
names(subOSL)[names(subOSL) == "Measure.volume"] = "eid"
subOSL$eid = gsub("/","",subOSL$eid)
msOSL = merge(msOSL, subOSL, by = c("eid"))
names(demoOSL)[names(demoOSL) == "Individual_TimePoint"] = "eid"
msOSL = merge(msOSL, demoOSL %>% select(EDSS, eid), by = c("eid"))
names(msOSL)[names(msOSL) == "EDSS"] = "edss"
df$scanner = ifelse(df$eid > 100 & df$eid < 200, "a", df$eid)
df$scanner = ifelse(df$eid > 200 & df$eid < 300, "b", df$eid)
df$scanner = ifelse(df$eid > 300 & df$eid < 400, "c", df$eid)
df$scanner = ifelse(df$eid > 400 & df$eid < 500, "d", df$eid)
df$scanner = ifelse(df$eid > 500 & df$eid < 600, "e", df$eid)
df$scanner = ifelse(df$eid > 600 & df$eid < 700, "f", df$eid)
df$scanner = ifelse(df$eid > 700 & df$eid < 800, "g", df$eid)
df$scanner = ifelse(df$eid > 800 & df$eid < 900, "h", df$eid)
df$scanner = ifelse(df$eid > 900 & df$eid < 1000, "i", df$eid)
df$scanner = ifelse(df$eid > 1000 & df$eid < 1100, "j", df$eid)
df$scanner = ifelse(df$eid > 1100 & df$eid < 1200, "k", df$eid)
df$scanner = ifelse(df$eid > 1200 & df$eid < 1300, "l", df$eid)
df$scanner = ifelse(df$eid > 1300 & df$eid < 1400, "m", df$eid)
df$scanner = ifelse(df$eid > 1400 & df$eid < 1500, "n", df$eid)
df$scanner = ifelse(df$eid > 1500 & df$eid < 1600, "o", df$eid)
df$scanner = ifelse(df$eid > 1600, "p", df$eid)
df$data = "OFAMS"
df$sex = factor(df$sex) # make sex a factor
levels(df$sex) = c("M","F")
names(df)[names(df) == "Left.Thalamus"] = "Left.Thalamus.Proper"
names(df)[names(df) == "Right.Thalamus"] = "Right.Thalamus.Proper"
df$session = as.numeric(factor(df$session))
msOSL = msOSL[,names(msOSL)[names(msOSL) %in% names(df)]]
msOSL$session = ifelse(grepl("_01",msOSL$eid) == T, 1,0) + ifelse(grepl("_02",msOSL$eid) == T, 2,0) + 
  ifelse(grepl("_03",msOSL$eid) == T, 3,0)+ ifelse(grepl("_04",msOSL$eid) == T, 4,0) +
  ifelse(grepl("_05",msOSL$eid) == T, 5,0)+ ifelse(grepl("_06",msOSL$eid) == T, 6,0) +
  ifelse(grepl("_07",msOSL$eid) == T, 7,0)+ ifelse(grepl("_08",msOSL$eid) == T, 8,0) +
  ifelse(grepl("_09",msOSL$eid) == T, 9,0)
msOSL$eid = substr(msOSL$eid,1,7)
df = df[,names(df)[names(df) %in% names(msOSL)]]

# longitudinal Combat
LC = function(dat){covars = dat %>% dplyr::select(sex,age,edss,data)
features = c(dat %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
             dat %>% select(starts_with("rh") & ends_with("volume")) %>% names(),
             names(dat)[grepl(c("halamus|allidum|mygdala|
                                campus|utamen|audate|
                                CC|Cerebellum.Cortex"),names(dat))],
             dat %>% select(TotalGrayVol, EstimatedTotalIntraCranialVol,
                            ends_with("Hippocampus")) %>% names())
dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
                 features = features,
                 formula = "age * sex", ranef = "(1|eid)", data = dat)
dat = dat$data_combat
dat = cbind(covars,dat)
colnames(dat) = gsub('.combat','',colnames(dat))
return(dat)
}
df =  filter(df, !scanner %in% c(501, 702, 802, 803, 1409, 1410))
df0=df
df = rbind(LC(msOSL),LC(na.omit(df)))
#
# For analysis of both datasets together, the harmonisation must be made together
df2 = LC(rbind(msOSL,na.omit(df0)))
#
#
#
# note: from here onwards, findings are being examined in OFAMS data and then replicated in independent Oslo data
#
#
# 1. General trend of GM vol degeneration ####
df$session = factor(df$session) # make session a factor
#df$geno = factor(df$geno) # make genotype a factor
df$sex = factor(df$sex) # make sex a factor
df$TIV = df$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df$TotalVol = df$TotalGrayVol/1000000
df2$TIV = df2$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df2$TotalVol = df2$TotalGrayVol/1000000
# from mm^3 to dm^3 or liters
# OFAMS
vol.m = lmer(TotalVol~age*sex+TIV+(1|eid),df%>%filter(data == "OFAMS"))
#Effect("session",vol.m)
plot(Effect("age",vol.m))
summary(vol.m)$coefficients
effectsize::standardize_parameters(vol.m)
interact_plot(vol.m,pred =  age, modx = sex) # interaction of sex

# Oslo
vol.m = lmer(TotalVol~age*sex+TIV+(1|eid),df%>%filter(data == "MS"))
#summary(vol.m) 
plot(Effect("age",vol.m))
summary(vol.m)$coefficients
effectsize::standardize_parameters(vol.m)
#plot(Effect("geno",vol.m))
#plot(Effect(c("geno","sex"),vol.m)) # seem like no interaction of sex and genotype
#vol.m0 = lmer(TotalVol~session+TIV+age*sex+(1|eid),df)
interact_plot(vol.m,pred =  age, modx = sex) # no interaction of genotype and age

# introduce a non-linear trend of age (yet similar to linear trend...)
vol.gam = gamm4(TotalVol~TIV+sex+s(age),random=~(1|eid),data=df%>%filter(data == "MS"))
plot(vol.gam$gam,pages=1)
vol.gam = gamm4(TotalVol~TIV+sex+s(age),random=~(1|eid),data=df%>%filter(data == "OFAMS"))
plot(vol.gam$gam,pages=1)

#plot(vol.gam$mer,pages=1) # underlying mixed model
#anova(vol.gam$gam)
#
###### IMPORTANT: WE SEE AN OVERALL QUITE LINEAR DEGENERATION.
###### THAT MEANS THAT ONE CAN ALSO USE THE DATA UP TO 2 YEARS TO MODEL DEGENERATIVE PROCESSES IN CASE WE LIMIT THE RESEARCH TO THE OMICS SAMPLING PERIOD

# 2. General association of GM vol and EDSS ####
edss.m = lmer(edss~TotalVol+sex+TIV+age+(1|eid),df%>%filter(data == "OFAMS"))
summary(edss.m) 
plot(Effect("TotalVol",edss.m),xlab = "Total brain volume",ylab="EDSS",main="Adjusted association of EDS and total volume")
effectsize::standardize_parameters(edss.m) # clear association between edss & vol
r.squaredGLMM(edss.m)
edss.m = lmer(edss~TotalVol+sex+TIV+age+(1|eid),df%>%filter(data == "MS"))
summary(edss.m) 
plot(Effect("TotalVol",edss.m),xlab = "Total brain volume",ylab="EDSS",main="Adjusted association of EDS and total volume")
effectsize::standardize_parameters(edss.m) # clear association between edss & vol
r.squaredGLMM(edss.m)
#
# one can also control for the interaction of brain volume and sex
edss.m = lmer(edss~TotalVol*sex+TIV+age+(1|eid),df%>%filter(data == "OFAMS"))
summary(edss.m) 
plot(Effect("TotalVol",edss.m),xlab = "Total brain volume",ylab="EDSS",main="Adjusted association of EDS and total volume")
interact_plot(edss.m,pred =  TotalVol, modx = sex) # interaction of sex
effectsize::standardize_parameters(edss.m) # clear association between edss & vol
r.squaredGLMM(edss.m)
edss.m = lmer(edss~TotalVol*sex+TIV+age+(1|eid),df%>%filter(data == "MS"))
summary(edss.m) 
plot(Effect("TotalVol",edss.m),xlab = "Total brain volume",ylab="EDSS",main="Adjusted association of EDS and total volume")
effectsize::standardize_parameters(edss.m) # clear association between edss & vol
r.squaredGLMM(edss.m)
interact_plot(edss.m,pred =  TotalVol, modx = sex) # interaction of sex
#
#
#
# edss.m3 = lmer(edss~TotalVol*geno+sex+TIV+age+(1|eid),df)
# summary(edss.m3) 
# interact_plot(edss.m3,pred =  TotalVol, modx = geno) # interaction of genotype

# edss.m3 = lmer(edss~TotalVol*geno+geno*sex+TIV+age+(1|eid),df)
# summary(edss.m3) 
# interact_plot(edss.m3,pred =  sex, modx = geno) 
# effectsize::standardize_parameters(edss.m3)
#
#
# r.squaredGLMM(edss.m)
# r.squaredGLMM(edss.m1)
# r.squaredGLMM(edss.m2)
# r.squaredGLMM(edss.m3)
# 
# anova(edss.m1,edss.m3) # the two best performing models are not different from each other
# we pick hence the simpler model, edss.m1
# contains sex,age,vol, no interaction effects
# 3. Regional age associations ####
  # 3.1 Unstandardized ####
ROIs = names(df)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(df))]
change = function(data_frame){
  # 3.1.1 Cortical #### #
ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
n_cort_ROIs = length(ROIs)
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,data_frame)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
plot_df = res
plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
plot_df01 = plot_df
plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional volume loss") + 
  theme_void()
p1 = p1+labs(fill="mm<sup>3</sup>") + #mm<sup>3</sup>/year
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
# 3.1.2 Sub-Cortical #### #
# subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
ROIs = ROIs[!grepl(c("CC"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,data_frame)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res2 = list_rbind(reslist)
res2$p.corr = (length(res2$p)+(length(c(data_frame %>% dplyr::select(starts_with("lh") & ends_with("volume")) %>% names(),
                                      data_frame %>% dplyr::select(starts_with("rh") & ends_with("volume")) %>% names())
)))*res2$p # Bonferroni correction considering all ROIs
#cerebellum = res2[grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
res=rbind(res,res2)
res2 = res2[!grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
plot_df1 = res2
plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
#brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
coronal_brain_aseg = as_tibble(aseg) %>%
  filter(side == "coronal", !grepl("\\d", label))
#
# Note: the painful steps above leading to this merger are crucial.
# This will be useful for later use of ggseg (aseg). Enjoy!
plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
plot_df02 = plot_df1
#
plot_df1$age.beta = ifelse(plot_df1$p.corr > .05, NA,plot_df1$age.beta)
#
p1.1 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = age.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional volume loss") + 
  theme_void()
p1.1 = p1.1+labs(fill="mm<sup>3</sup>") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
p = ggarrange(p1,p1.1,nrow=1,widths=c(2,.5))
return(list(res,p))
}
#
# Apply the function which was established above
#
#age.plot = ggarrange(p1,p1.1,ncol=2,widths=c(2,.5),labels=c("a","b"))
age.plot = ggarrange(change(df%>%filter(data == "MS"))[[2]],change(df%>%filter(data == "OFAMS"))[[2]],nrow=2,labels=c("a","b"))
age.plot = annotate_figure(age.plot, top = text_grob("Annual regional volume loss",face = "bold", size = 17))
ggsave(paste(savepath,"age_plot_unstandardized.pdf",sep=""),age.plot, width = 14, height = 6)
#safe table of the results
unstd = rbind(change(df%>%filter(data == "MS"))[[1]],change(df%>%filter(data == "OFAMS"))[[1]])
unstd$data = c(replicate(nrow(unstd)/2,"Oslo"),replicate(nrow(unstd)/2,"OFAMS"))
write.csv(file = paste(savepath,"age_beta_unstandardized.csv",sep=""),unstd)
#
#
# create output table with all values in it 
#age.tab = rbind(plot_df01 %>% select(region,hemi,age.beta,SE,std.beta,CI_low,CI_high, t,p,p.corr),plot_df02 %>% select(region,hemi,age.beta,SE,std.beta,CI_low,CI_high,t,p,p.corr))
#
  # 3.2 Standardized ####
plist = the_reslist = list()
data_list = list(df%>%filter(data == "MS"),df%>%filter(data == "OFAMS"),df2)
for (i in 1:length(data_list)){
  data_frame = data_list[[i]]
  # 3.1.1 Cortical #### #
  ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
  n_cort_ROIs = length(ROIs)
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res = list_rbind(reslist)
  res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
  plot_df = res
  plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
  plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
  plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
  plot_df01 = plot_df
  plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
  p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.7,0.1)) +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  # 3.1.2 Sub-Cortical #### #
  # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
  ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
  ROIs = ROIs[!grepl(c("CC"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res2 = list_rbind(reslist)
  res2$p.corr = (length(res2$p)+(length(c(data_frame %>% dplyr::select(starts_with("lh") & ends_with("volume")) %>% names(),
                                          data_frame %>% dplyr::select(starts_with("rh") & ends_with("volume")) %>% names())
  )))*res2$p # Bonferroni correction considering all ROIs
  #cerebellum = res2[grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  res=rbind(data.frame(res),data.frame(res2))
  res2 = res2[!grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  plot_df1 = res2
  plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
  plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
  #brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
  coronal_brain_aseg = as_tibble(aseg) %>%
    filter(side == "coronal", !grepl("\\d", label))
  #
  # Note: the painful steps above leading to this merger are crucial.
  # This will be useful for later use of ggseg (aseg). Enjoy!
  plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
  plot_df02 = plot_df1
  #
  plot_df1$age.beta = ifelse(plot_df1$p.corr > .05, NA,plot_df1$age.beta)
  #
  p1.1 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.7,0.1)) +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1.1 = p1.1+labs(fill="Std.Coeff.") +
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  p = ggarrange(p1,p1.1,nrow=1,widths=c(2,.5))
  the_reslist[[i]] = res
  plist[[i]] = p
}
#
# Before plotting it all together, we add a plot based on the data of the healthy controls which was imported
#
# We take a subset of the lifespan data including only
#       1 time point
#       HCs 
#       and age and sex matches (at MS study onset),
#       which are harmonized.
# 
# 1) filter
Nm = (lifespan %>% filter(session != 2 & session != 3 & diagnosis == "HC"))
Nm$TIV = Nm$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df_list = list(Nm,df)
df_list[[1]]$session = 1
df_list[[2]]$diagnosis = "MS"
commons = Reduce(function(x, y){intersect(x, names(y))}, df_list, init = names(df_list[[1]]))
Nm = rbind(df_list[[1]],df_list[[2]][commons])
Nm = Nm %>% filter(session == 1)
Nm = na.omit(Nm)
# 2) match data frames for sex and age at MS study onset
test = matchit(as.factor(diagnosis) ~ age + sex, data = Nm,
               method = "full", distance = "glm")
Nm = match.data(test)
Nm = Nm %>% select(-c(distance, weights, subclass))
table(Nm$diagnosis)
#
factor(Nm$diagnosis)

# Following structure was applied:
# ROIs.region.    age.beta       CI_high       CI_low            p       p.corr  data
min(df$age)
max(df$age)

min(na.omit(lifespan$age))
max(na.omit(lifespan$age))

hist(df$age)
hist(lifespan$age)

#

names(lifespan)


lmlist = list(lifespan)
for (i in 1:length(lm_list)){
  data_frame = lm_list[[i]]
  # 3.1.1 Cortical #### #
  ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
  n_cort_ROIs = length(ROIs)
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+age"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res = list_rbind(reslist)
  res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
  plot_df = res
  plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
  plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
  plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
  plot_df01 = plot_df
  plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
  p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.7,0.1)) +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  # 3.1.2 Sub-Cortical #### #
  # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
  ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
  ROIs = ROIs[!grepl(c("CC"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res2 = list_rbind(reslist)
  res2$p.corr = (length(res2$p)+(length(c(data_frame %>% dplyr::select(starts_with("lh") & ends_with("volume")) %>% names(),
                                          data_frame %>% dplyr::select(starts_with("rh") & ends_with("volume")) %>% names())
  )))*res2$p # Bonferroni correction considering all ROIs
  #cerebellum = res2[grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  res=rbind(data.frame(res),data.frame(res2))
  res2 = res2[!grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  plot_df1 = res2
  plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
  plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
  #brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
  coronal_brain_aseg = as_tibble(aseg) %>%
    filter(side == "coronal", !grepl("\\d", label))
  #
  # Note: the painful steps above leading to this merger are crucial.
  # This will be useful for later use of ggseg (aseg). Enjoy!
  plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
  plot_df02 = plot_df1
  #
  plot_df1$age.beta = ifelse(plot_df1$p.corr > .05, NA,plot_df1$age.beta)
  #
  p1.1 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.7,0.1)) +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1.1 = p1.1+labs(fill="Std.Coeff.") +
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  p = ggarrange(p1,p1.1,nrow=1,widths=c(2,.5))
  the_reslist[[i]] = res
  plist[[i]] = p
}












# Plot
age.plot = ggarrange(plotlist = plist,ncol = 1,labels=c("a","b","c"))
age.plot = annotate_figure(age.plot, top = text_grob("Ageing effects on regional brain volumes",face = "bold", size = 17))
ggsave(paste(savepath,"age_plot_standardized.pdf",sep=""),age.plot, width = 14, height = 6)
#safe table of the results
unstd = list_rbind(the_reslist)
unstd$data = c(replicate(nrow(unstd)/length(data_list),"Oslo"),
               replicate(nrow(unstd)/length(data_list),"OFAMS"),
               replicate(nrow(unstd)/length(data_list),"both"))
write.csv(file = "age_beta_standardized.csv",unstd)
#
#
#
#
# For comparison: OFAMS data without the 10-year follow-up ####
# this QC steps allows for a closer mimicking of the follow-up times between datasets
plist2 = the_reslist2 = list()
data_list2 = list(df%>%filter(data == "MS"),
                  df%>%filter(data == "OFAMS")%>%filter(session!=16),
                  df2%>%filter(session!=16))
for (i in 1:length(data_list2)){
  data_frame = data_list2[[i]]
  # 3.1.1 Cortical #### #
  ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
  n_cort_ROIs = length(ROIs)
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res = list_rbind(reslist)
  res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
  plot_df = res
  plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
  plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
  plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
  plot_df01 = plot_df
  plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
  p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  # 3.1.2 Sub-Cortical #### #
  # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
  ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
  ROIs = ROIs[!grepl(c("CC"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res2 = list_rbind(reslist)
  res2$p.corr = (length(res2$p)+(length(c(data_frame %>% dplyr::select(starts_with("lh") & ends_with("volume")) %>% names(),
                                          data_frame %>% dplyr::select(starts_with("rh") & ends_with("volume")) %>% names())
  )))*res2$p # Bonferroni correction considering all ROIs
  #cerebellum = res2[grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  res=rbind(data.frame(res),data.frame(res2))
  res2 = res2[!grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
  plot_df1 = res2
  plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
  plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
  #brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
  coronal_brain_aseg = as_tibble(aseg) %>%
    filter(side == "coronal", !grepl("\\d", label))
  #
  # Note: the painful steps above leading to this merger are crucial.
  # This will be useful for later use of ggseg (aseg). Enjoy!
  plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
  plot_df02 = plot_df1
  #
  plot_df1$age.beta = ifelse(plot_df1$p.corr > .05, NA,plot_df1$age.beta)
  #
  p1.1 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = age.beta),color="black")+
    #scale_fill_viridis_c(option = "cividis", direction = -1)+
    scale_fill_gradient2(low = "blue",mid = "white",high="red") +
    #labs(title="Regional volume loss") + 
    theme_void()
  p1.1 = p1.1+labs(fill="Std.Coeff.") +
    theme(
      plot.title = element_markdown(),
      legend.title = element_markdown()
    )
  p = ggarrange(p1,p1.1,nrow=1,widths=c(2,.5))
  the_reslist2[[i]] = res
  plist2[[i]] = p
}
age.p = ggarrange(plotlist = plist,ncol = 1,labels=c("a","b","c"))
age.p = annotate_figure(age.p, top = text_grob("Annual regional volume loss",face = "bold", size = 17))
ggsave(paste(savepath,"age_plot_SHORT_standardized.pdf",sep=""),age.p, width = 14, height = 6)
#safe table of the results
unstd_SHORT = list_rbind(the_reslist)
unstd_SHORT$data = c(replicate(nrow(unstd_SHORT)/length(data_list),"Oslo"),
               replicate(nrow(unstd_SHORT)/length(data_list),"OFAMS"),
               replicate(nrow(unstd_SHORT)/length(data_list),"both"))
write.csv(file = paste(savepath,"age_beta_SHORT_standardized.csv",sep=""),unstd_SHORT)
#
#
#
#
#
#
#
#
# 4. Regional EDSS associations ####
  # 4.1 UnStandardized ####
    # 4.1.1 Cortical ####
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  edss.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],edss.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
plot_df = res
plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
plot_df01 = plot_df
#
# Note: UN-corrected p-values used here!!
#
plot_df$edss.beta = ifelse(plot_df$p > .05, NA,plot_df$edss.beta)
#
p2 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = edss.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p2 = p2+labs(fill="mm<sup>3</sup>/EDSS") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
    # 4.1.2 Sub-Cortical ####
ROIs = names(df)[2:46]
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  edss.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],edss.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = length(res$p)*res$p # p.adjust(res$p, method="BH")
plot_df1 = res
plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
#brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
coronal_brain_aseg = as_tibble(aseg) %>%
  filter(side == "coronal", !grepl("\\d", label))
plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
plot_df02 = plot_df1
#
# !!!!!!!
#
plot_df1$edss.beta = ifelse(plot_df1$p > .05, NA,plot_df1$edss.beta) # we use only uncorrected p-vals < 0.05 !!!!!
#
#
#
p2.2 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = edss.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p2.2 = p2.2+labs(fill="mm<sup>3</sup>/EDSS") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
  # 4.2 Standardized ####
    # 4.2.1 Cortical ####
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low)
}
res = list_rbind(reslist)
plot_df01 = merge(plot_df01,res,by="ROIs.region.")
plot_df = plot_df01
plot_df$std.beta = ifelse(plot_df$p > .05, NA,plot_df$std.beta)
p03 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p03 = p03+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
    # 4.2.2 Sub-Cortical ####
ROIs = names(df)[2:46]
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low)
}
res = list_rbind(reslist)
plot_df02 = merge(plot_df02,res,by="ROIs.region.")
plot_df = plot_df02
plot_df$std.beta = ifelse(plot_df$p > .05, NA,plot_df$std.beta)
p04 = ggplot(plot_df) + geom_brain(atlas = aseg,side = "coronal",aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p04 = p04+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
# make a table containing outputs of cortical and subcortical volume associations with 
#edss.tab = rbind(plot_df01 %>% select(region,hemi,edss.beta,SE,t,p,p.corr),plot_df02 %>% select(region,hemi,edss.beta,SE,t,p,p.corr))
edss.tab = rbind(plot_df01 %>% select(region,hemi,edss.beta,SE,std.beta,CI_high,CI_low,t,p,p.corr),plot_df02 %>% select(region,hemi,edss.beta,SE,std.beta,CI_high,CI_low,t,p,p.corr))
#
# 5. Create outputs ####
#
  # 5.1 Merge plots ####
#
# We see some strong associations between age and local volumes
#ggarrange(p1,p2,p1.1,p2.2,ncol=2,nrow=2,widths=c(1,1), heights = c(1.5,0.7))
#ggarrange(p1,p1.1,p2,p2.2,ncol=2,nrow=2,widths=c(1,.5), heights = c(1,1))
#
#
#
# Start with unstandardized plots
#
edss.plot = ggarrange(p2,p2.2,ncol=2,widths=c(2,.5),labels=c("c","d"))
edss.plot = annotate_figure(edss.plot, top = text_grob("Associations of EDSS and regional brain volume",face = "bold", size = 17))
plot1 = ggarrange(age.plot, edss.plot, nrow=2)
ggsave(paste(savepath,"reg_plot.pdf",sep=""),plot1, width = 14, height = 6)
#
#
# Now, the standardized plots
edss.plot = ggarrange(p03,p04,ncol=2,widths=c(2,.5),labels=c("c","d"))
edss.plot = annotate_figure(edss.plot, top = text_grob("Associations of EDSS and regional brain volume",face = "bold", size = 17))
plot1 = ggarrange(age.plot, edss.plot, nrow=2)
ggsave(paste(savepath,"reg_plot_standardized.pdf",sep=""),plot1, width = 14, height = 6)
#
  # 5.2 Save tables ####
write.csv(age.tab,file = paste(savepath,"age_tab.csv",sep=""))
write.csv(edss.tab,file = paste(savepath,"edss_tab.csv",sep=""))
#
#
#
#
# Done.
# test = read.csv(paste(savepath,"age_tab.csv",sep=""))
# test = test %>% filter(p<.05)
# 
# test = read.csv(paste(savepath,"edss_tab.csv",sep=""))
# test = test %>% filter(p<.05)
