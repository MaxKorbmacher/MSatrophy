# Analyse GM atrophy among people with multiple sclerosis in the OFAMS (Norwegian multicenter) clinical trial and local data from the Oslo University Hospital
# Max Korbmacher, 20.12.2024
# last change: 22 April 2025
#
# Note: the data were NOT harmonized across datasets, but across scanners within datasets when attempting to replicate the findings from one dataset to another.
#       This was done to provide unbiased replications between Oslo and OFAMS data.
#       However, the datasets were harmonized together for scanner differences when analysing them together!
#       The harmonisation was done longitudinally for longitudinal data and cross-sectional for cross sectional data.
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
# install.packages("devtools")
# library(devtools)
# devtools::install_github("jcbeer/longCombat")
# install_github("jfortin1/neuroCombatData")
# install_github("jfortin1/neuroCombat_Rpackage")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,lme4,lmerTest,effects,effectsize,interactions,gamm4,
               ggseg,ggtext,ggpubr,MuMIn,dplyr,ggplot2,standardize,longCombat,
               neuroCombat,haven,reshape2,simr,data.table)
# data
lifespan = read.csv("/Users/max/Documents/Local/Data/Lifespan/cortical_subcortical.csv")
df = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
OSLl  = read.delim("/Users/max/Documents/Local/Data/Oslo/lh.volume.UKBB.txt")
OSLr  = read.delim("/Users/max/Documents/Local/Data/Oslo/rh.volume.UKBB.txt")
OSLsub  = read.delim("/Users/max/Documents/Local/Data/Oslo/subcorticalstats.UKBB.txt")
demoOSL = read.csv("/Users/max/Documents/Local/Data/Oslo/Oslo_demographics.csv")
edss = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/edss.sav")
edss = reshape2::melt(edss,id.vars = c("Patnr"))
names(edss) = c("eid","session","edss")
edss$session = edss$session
edss$session = factor(edss$session)
levels(edss$session) = c(1,7,13,19,25)

# standardize and join datasets
colnames(OSLsub)[1] = colnames(OSLr)[1] = colnames(OSLl)[1] = "eid"
OSLl = OSLl %>% filter(grepl("long",eid))
OSLl$eid = substr(OSLl$eid,1,10)
OSLr = OSLr %>% filter(grepl("long",eid))
OSLr$eid = substr(OSLr$eid,1,10)
OSLsub = OSLsub %>% filter(grepl("long",eid))
OSLsub$eid = substr(OSLsub$eid,1,10)
msOSL = merge(OSLr,OSLl %>% select(-eTIV,-BrainSegVolNotVent),by="eid")
msOSL = merge(msOSL,OSLsub%>% select(-BrainSegVolNotVent),by="eid")
df$session = df$session+1
#names(subOSL)[names(subOSL) == "Measure.volume"] = "eid"
#subOSL$eid = gsub("/","",subOSL$eid)
#msOSL = merge(msOSL, subOSL, by = c("eid"))
names(demoOSL)[names(demoOSL) == "Individual_TimePoint"] = "eid" 
#
#
#
#
#
# NOTE!
# We filter for RRMS
demoOSL = demoOSL %>% filter(Subtype_MS == "RRMS")
#
#
#
#
#
#
msOSL = merge(msOSL, demoOSL %>% select(eid,EDSS,Scanner, Gender, Age), by = c("eid"))
names(msOSL)[names(msOSL) == "EDSS"] = "edss"
names(msOSL)[names(msOSL) == "Scanner"] = "scanner"
names(msOSL)[names(msOSL) == "Gender"] = "sex"
names(msOSL)[names(msOSL) == "Age"] = "age"
df$scanner = ifelse(df$eid > 100 & df$eid < 200, "a", df$eid)
df$scanner = ifelse(df$eid > 200 & df$eid < 300, "b", df$scanner)
df$scanner = ifelse(df$eid > 300 & df$eid < 400, "c", df$scanner)
df$scanner = ifelse(df$eid > 400 & df$eid < 500, "d", df$scanner)
df$scanner = ifelse(df$eid > 500 & df$eid < 600, "e", df$scanner)
df$scanner = ifelse(df$eid > 600 & df$eid < 700, "f", df$scanner)
df$scanner = ifelse(df$eid > 700 & df$eid < 800, "g", df$scanner)
df$scanner = ifelse(df$eid > 800 & df$eid < 900, "h", df$scanner)
df$scanner = ifelse(df$eid > 900 & df$eid < 1000, "i", df$scanner)
df$scanner = ifelse(df$eid > 1000 & df$eid < 1100, "j", df$scanner)
df$scanner = ifelse(df$eid > 1100 & df$eid < 1200, "k", df$scanner)
df$scanner = ifelse(df$eid > 1200 & df$eid < 1300, "l", df$scanner)
df$scanner = ifelse(df$eid > 1300 & df$eid < 1400, "m", df$scanner)
df$scanner = ifelse(df$eid > 1400 & df$eid < 1500, "n", df$scanner)
df$scanner = ifelse(df$eid > 1500 & df$eid < 1600, "o", df$scanner)
df$scanner = ifelse(df$eid > 1600, "p", df$scanner)
df$data = "OFAMS"
df$sex = factor(df$sex) # make sex a factor
levels(df$sex) = c("M","F")
# names(df)[names(df) == "Left.Thalamus"] = "Left.Thalamus.Proper"
# names(df)[names(df) == "Right.Thalamus"] = "Right.Thalamus.Proper"
#df$session = as.numeric(factor(df$session))
msOSL = msOSL[,names(msOSL)[names(msOSL) %in% names(df)]]
msOSL$session = ifelse(grepl("_01",msOSL$eid) == T, 1,0) + ifelse(grepl("_02",msOSL$eid) == T, 2,0) + 
  ifelse(grepl("_03",msOSL$eid) == T, 3,0)+ ifelse(grepl("_04",msOSL$eid) == T, 4,0) +
  ifelse(grepl("_05",msOSL$eid) == T, 5,0)+ ifelse(grepl("_06",msOSL$eid) == T, 6,0) +
  ifelse(grepl("_07",msOSL$eid) == T, 7,0)+ ifelse(grepl("_08",msOSL$eid) == T, 8,0) +
  ifelse(grepl("_09",msOSL$eid) == T, 9,0)
msOSL$eid = substr(msOSL$eid,1,7)
msOSL$data = "MS"
msBGO = df
#df = df[,names(df)[names(df) %in% names(msOSL)]]
# longitudinal Combat
LC = function(dat){
  features = c(dat %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
               dat %>% select(starts_with("rh") & ends_with("volume")) %>% names(),
               names(dat)[grepl(c("halamus|allidum|mygdala|
                                campus|utamen|audate|
                                CC|Cerebellum.Cortex"),names(dat))],
               dat %>% select(TotalGrayVol, EstimatedTotalIntraCranialVol,
                              ends_with("Hippocampus")) %>% names())
  dat = dat %>% select(all_of(features),eid,data,sex,session,scanner, age, edss)%>%na.omit
  covars = dat %>% dplyr::select(sex,edss,data,age)
  dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
                   features = features,
                   formula = "age + sex", ranef = "(1|eid)", data = dat)
  dat = dat$data_combat
  dat = cbind(covars,dat)
  colnames(dat) = gsub('.combat','',colnames(dat))
  return(dat)
}
#
#
#
#
#df =  filter(df, !scanner %in% c(501, 702, 802, 803, 1409, 1410,103, 104, 207, 214, 216, 303, 601, 701, 703, 708, 709, 801, 806, 808, 809, 813, 901, 1202, 1412, 1413, 1414))
msOSL =  filter(msOSL, !scanner %in% c("NMR1"))
msOSL$scanner = as.numeric(factor(msOSL$scanner))
msOSL$sex = ifelse(msOSL$sex == "male","M",msOSL$sex)
msOSL$sex = ifelse(msOSL$sex == "female","F",msOSL$sex)
df$scanner = as.numeric(factor(df$scanner))
df = merge(df,edss,by=c("eid","session"),all = T)
df$edss = coalesce(df$edss.x,df$edss.y)
df = df %>% select(-X,-X.1)
X = as.data.table(df)
keys = df %>% select(eid,session,sex,data,geno_HLA1501_1,scanner) %>% names
X1 = X %>% filter(session == 1)
X2 = X %>% filter(session == 24)
X3 = X %>% filter(session == 145)
X1 = X1[,lapply(.SD,base::mean),keys]
X2 = X2[,lapply(.SD,base::mean),keys]
X3 = X3[,lapply(.SD,base::mean),keys]
df = rbind(df %>% filter(session != 1 & session != 24 & session != 145),X1)
df = rbind(df,X2)
df = rbind(df,X3)
df0=df
df = rbind(LC(msOSL),LC((df)))
OFAMS_N_SCANS = df%>%filter(data=="OFAMS")%>%nrow
#
# For analysis of both datasets together, the harmonisation must be made together
msOSL = msOSL[,names(msOSL)[names(msOSL) %in% names(df0)]]
df0 = df0[,names(df0)[names(df0) %in% names(msOSL)]]
df2 = LC(rbind(msOSL,(df0)))
#
#
# note: from here onwards, findings are being examined in OFAMS data and then replicated in independent OUH data
#
#
# 1. General trend of GM vol degeneration ####
#df$session = factor(df$session) # make session a factor
#df$geno = factor(df$geno) # make genotype a factor
df$sex = factor(df$sex) # make sex a factor
df$TIV = df$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df$TotalVol = df$TotalGrayVol/1000000
df2$TIV = df2$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
df2$TotalVol = df2$TotalGrayVol/1000000
  # 3.2 Standardized ####
# add UKB longitudinal data to check for random effects
UKBlong = lifespan[lifespan$eid %in% (lifespan %>% filter(data == "UKB" & session == 2))$eid,]
UKBlong = UKBlong %>% filter(age < max(df$age)) # limits maximum age to that of the MS cohort
UKBlong$TIV = UKBlong$EstimatedTotalIntraCranialVol
UKBlong$sex = ifelse(ifelse(UKBlong$sex == "Female",0,UKBlong$sex) == "Male",1,ifelse(UKBlong$sex == "Female",0,UKBlong$sex))
plist = the_reslist = list()
UKBlong$TotalGrayVol = UKBlong$EstimatedTotalIntraCranialVol # placeholder variable
UKBlong$edss = 0 # placeholder variable
UKBlong = LC(UKBlong)
df2$TIV = df2$EstimatedTotalIntraCranialVol
UKBlong$TIV = UKBlong$TotalGrayVol
data_list = list(df%>%filter(data == "MS"),df%>%filter(data == "OFAMS"),df2, UKBlong)


for (i in 1:length(data_list)){
  data_frame = data_list[[i]]
  # 3.1.1 Cortical #### #
  ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
  #n_cort_ROIs = length(ROIs)
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
  res$p.corr = (80+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
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
    theme_void() + theme(legend.position="none")
  # p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
  #   theme(
  #     plot.title = element_markdown(),
  #     legend.title = element_markdown()
  #   )
  # 3.1.2 Sub-Cortical #### #
  # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
  ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
  ROIs = ROIs[!grepl(c("CC"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
  for(region in 1:(length(ROIs))){
    f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
    mod = lmer(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res2 = list_rbind(reslist)
  res2$p.corr = 80*res2$p # Bonferroni correction considering all ROIs
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
  plot_df1$age.beta = ifelse(plot_df1$p > .05, NA,plot_df1$age.beta)
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
#       and match age span,
#       which are harmonized.
# 
# 1) filter
Nm = (lifespan %>% filter(session != 2 & session != 3 & diagnosis == "HC"))
(lifespan %>% filter(diagnosis == "HC")) %>% nrow()
Nm$TIV = Nm$EstimatedTotalIntraCranialVol/1000000 # transform into liters / dm3
# df_list = list(Nm,df)
# df_list[[1]]$session = 1
# df_list[[2]]$diagnosis = "MS"
# commons = Reduce(function(x, y){intersect(x, names(y))}, df_list, init = names(df_list[[1]]))
# Nm = rbind(df_list[[1]],df_list[[2]][commons])
# Nm = Nm %>% filter(session == 1)
Nm = na.omit(Nm)
# 2) harmonize
covars = Nm %>% dplyr::select(eid,sex,scanner,age,data)
covars$sex = ifelse(covars$sex == "F" | covars$sex == "Female",0,1)
datasets = covars$data
covars$data = as.numeric(factor(Nm$data))
Nm = neuroCombat(t(Nm%>%dplyr::select(EstimatedTotalIntraCranialVol,starts_with("Left"),starts_with("Right"), starts_with("lh"),starts_with("rh"))),batch=as.numeric(factor(Nm$scanner)),mod=model.matrix(~covars$age+covars$sex), mean.only = T)
Nm = data.frame(t(Nm$dat.combat))
Nm = cbind(covars,Nm)
#
# 3) match age span
Nm$datasets = datasets
Nm = Nm %>% filter(age > min(df$age) & age < max(df$age))
# plot((Nm %>% filter(age > min(df$age) & age < max(df$age)))$age,(Nm %>% filter(age > min(df$age) & age < max(df$age)))$Right.Thalamus.Proper)
# cor((Nm %>% filter(age > min(df$age) & age < max(df$age)))$age,(Nm %>% filter(age > min(df$age) & age < max(df$age)))$Right.Thalamus.Proper)
# plot((df%>%filter(session==1))$age,(df%>%filter(session==1))$Right.Thalamus.Proper)
#
#
# # If wanting to match:
# test = matchit(as.factor(diagnosis) ~ age + sex, data = Nm,
#                method = "full", distance = "glm")
# summary(test)
# 
# Nm = match.data(test)
# Nm = Nm %>% select(-c(distance, weights, subclass))
# table(Nm$diagnosis)
# hist((Nm %>% filter(diagnosis == "HC"))$age)
# #
# factor(Nm$diagnosis)
#
# Plot
res_list = p_list = list()
lm_list = list(Nm)
for (i in 1:length(lm_list)){
  data_frame = lm_list[[i]]
  # 3.1.1 Cortical #### #
  ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
  #n_cort_ROIs = length(ROIs)
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+age"))
    mod = lm(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,4]
    reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
  }
  res = list_rbind(reslist)
  res$p.corr = (80+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
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
    theme_void() + theme(legend.position="none")
  # p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
  #   theme(
  #     plot.title = element_markdown(),
  #     legend.title = element_markdown()
  #   )
  # 3.1.2 Sub-Cortical #### #
  # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
  ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
  ROIs = ROIs[!grepl(c("CC"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
  for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+age"))
    mod = lm(f,data_frame)
    age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,4]
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
  res_list[[i]] = res
  p_list[[i]] = p
}
#
# Add the cross-sectional results to the longitudinal results
plist[[5]] = p_list[[1]]
the_reslist[[5]] = res_list[[1]]
#
# Plot all together
age.plot = ggarrange(plotlist = plist[1:4],ncol = 1,labels=c("OUH MS data","OFAMS MS data","Combined MS data", "UK Biobank HC data"),hjust = c(0,0,0))

# age.plot1 = ggarrange(plotlist = plist[1:3],ncol = 1,labels=c("OUH data","OFAMS data","Combined data"),hjust = c(0,0,0))
# age.plot1 = annotate_figure(age.plot1, top = text_grob("Multiple sclerosis",face = "bold", size = 17))
# 
# age.plot2 = ggarrange(plotlist = plist[4:5],ncol = 1,labels=c("UK Biobank","Lifespan sample (cross-sectional)"),hjust = c(0,0,0))
# age.plot2 = annotate_figure(age.plot2, top = text_grob("Healthy controls",face = "bold", size = 17))
# 
# age.plot = ggarrange(age.plot1,age.plot2,ncol=1,heights = c(3,2.05))

ggsave(paste(savepath,"age_plot_standardized.pdf",sep=""),age.plot, width = 10, height = 8)
#safe table of the results
unstd = list_rbind(the_reslist)
unstd$data = c(replicate(nrow(unstd)/(length(data_list)+1),"OUH"),
               replicate(nrow(unstd)/(length(data_list)+1),"OFAMS"),
               replicate(nrow(unstd)/(length(data_list)+1),"both"),
               replicate(nrow(unstd)/(length(data_list)+1),"UKB"),
               replicate(nrow(unstd)/(length(data_list)+1),"cross_sectional"))
unstd = unstd %>% filter(data != "cross_sectional")
unstd$p.corr = ifelse(unstd$p.corr > 1,1,unstd$p.corr)
unstd_copy = unstd
unstd_copy[2:6] = round(unstd_copy[2:6],2)
write.csv(file = paste(savepath,"age_beta_standardized.csv",sep=""),unstd_copy)
#
#
#
#
# Key info to report ####
# report the average ageing of participants in the longitudinal studies
paste("OFAMS participants (N=85, scans=",OFAMS_N_SCANS,") were aged 38.9±8.3 (range: 19-58) years at baseline and 49.6±8.6 at the final follow-up, with 65.9% being females at the first visit, and ",
      round((table(df %>% filter(data=="OFAMS") %>% filter(session=="145")%>%select(sex))[1]/df %>% filter(data=="OFAMS") %>% filter(session=="145")%>%nrow)*100,1),"% females at the final visit.")


df %>% filter(data=="OFAMS") %>% filter(session=="145")%>%select(sex)

paste("The Oslo sample participants (N=", nrow(msOSL%>%filter(session==min(as.numeric(msOSL$session)))),
      ", scans=",nrow(msOSL),") were, on average, aged ", 
      msOSL %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(mean(Age),1)),
      "±",
      msOSL %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(sd(Age),1)),
      " (range: ",
      msOSL %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(min(Age),1)),
      "-",
      msOSL %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(max(Age),1)),
      ") years at baseline and ",
      msOSL %>% group_by(eid) %>% summarize(Age = max(age)) %>% reframe(M = round(mean(Age),0)),
      "±",
      msOSL %>% group_by(eid) %>% summarize(Age = max(age)) %>% reframe(M = round(sd(Age),0)),
      " at the last follow-up, with ",
      round((table((msOSL %>% filter(session == 1))$sex)[1]/msOSL %>% filter(session == 1) %>% nrow)*100,1),
      "% being females at the first visit, and ",
      round((table((msOSL %>% group_by(eid) %>% filter(session == max(session)) %>%ungroup)$sex))[1]/(msOSL %>% group_by(eid) %>% filter(session == max(session)) %>%ungroup%>%nrow)*100,1),
      "% being females at the last visit",sep="")
paste("The healthy control ageing sample participants (N=", nrow(UKBlong%>%filter(session==min(as.numeric(UKBlong$session)))),
      ", scans=",nrow(UKBlong),") were, on average, aged ", 
      UKBlong %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(mean(Age),1)),
      "±",
      UKBlong %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(sd(Age),1)),
      " (range: ",
      UKBlong %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(min(Age),1)),
      "-",
      UKBlong %>% group_by(eid) %>% summarize(Age = min(age)) %>% reframe(M = round(max(Age),1)),
      ") years at baseline and ",
      UKBlong %>% group_by(eid) %>% summarize(Age = max(age)) %>% reframe(M = round(mean(Age),0)),
      "±",
      UKBlong %>% group_by(eid) %>% summarize(Age = max(age)) %>% reframe(M = round(sd(Age),0)),
      " at the single follow-up, with ",
      round((table((UKBlong %>% filter(session == 1))$sex)[1]/UKBlong %>% filter(session == 1) %>% nrow)*100,1),
      "% being females.",sep="")
lifespan$sex = factor(lifespan$sex)
levels(lifespan$sex) = c(0,0,1,0,0,1,1)
paste("Finally, the cross-sectional lifespan sample (N=", nrow(Nm),
      ") were, on average, aged ", 
      Nm%>% summarize(Age = round(mean(na.omit(age)),2)),
      "±",
      Nm%>% summarize(Age = round(sd(na.omit(age)),2)),
      " (range: ",
      Nm%>% summarize(Age = min(round(na.omit(age),2))),
      "-",
      Nm%>% summarize(Age = max(round(na.omit(age),2))),
      ") years,with ",
      round((table((Nm)$sex)[1]/Nm%>% nrow)*100,1),
      "% being females.",sep="")
#
# check all the HC at a single time point
table(Nm$datasets)
nrow(Nm)
Nm %>% group_by(datasets) %>% summarize(M = mean(age), SD = sd(age), Min = min(age), Max = max(age))

# some demographic background, here the age distributions of pwMS vs UKB:
# hist(df$age,xlab = "MS Age", main = "")
# hist(UKBlong$age, xlab = "UKB Age", main = "")
#
#
# msOSL %>% group_by(eid) %>% summarize(span = max(age)-min(age)) %>% summarize(ISI = mean(span))
# (msOSL %>% group_by(eid) %>% summarize(span = max(age)-min(age)) %>% reframe(range = range(span)))[1,]
# (msOSL %>% group_by(eid) %>% summarize(span = max(age)-min(age)) %>% reframe(range = range(span)))[2,]
#
# Assemble key information for the "age_plot_standardized.pdf" figure legend
paste("Figure X. Volume changes in relapsing-remitting multiple sclerosis and age-associations in healty controls. Panel a) presents volume changes in the Oslo sample (baseline N = ", 
      nrow(msOSL%>%filter(session==min(as.numeric(msOSL$session)))), 
      "). Panel b) presents volume changes in the OFAMS sample (baseline N = 85). Panel c) presents volume changes in both MS samples combined. Panel d) presents volume changes in the UK Biobank sample of healthy controls (N = ",
      UKBlong %>% filter(session==1) %>%nrow(),") Panel e) presents volume-age associations in a cross-section sample of HC (N = ",
      nrow(Nm)[1],"). The age range across time points for the longitudinal sample of pwMS was ",
      round(range(df$age)[1],2)," to ",round(range(df$age)[2],2)," years, for the longitudinal  healthy control sample ",
      round(range(UKBlong$age)[1],2)," to ",round(range(UKBlong$age)[2],2),", and for the cross-sectional sample ", round(range(Nm$age)[1],2)," to ", round(range(Nm$age)[2],2), " years.", sep="")
#
# We also highlight the strongest effects for the results section
unstd %>% group_by(data) %>%filter(age.beta == min(age.beta))
unstd %>% group_by(data) %>%filter(ROIs.region. == "Left.Thalamus" | ROIs.region. == "Right.Thalamus")
#unstd %>% group_by(data) %>%filter(ROIs.region. == "Left.Putamen" | ROIs.region. == "Right.Putamen")
unstd %>% group_by(data) %>%filter(ROIs.region. == "lh_superiorfrontal_volume" | ROIs.region. == "rh_superiorfrontal_volume")
unstd %>% group_by(data) %>%filter(ROIs.region. == "lh_rostralmiddlefrontal_volume" | ROIs.region. == "rh_rostralmiddlefrontal_volume")
unstd %>% group_by(data) %>%filter(ROIs.region. == "lh_parsorbitalis_volume" | ROIs.region. == "rh_parsorbitalis_volume")
unstd %>% group_by(data) %>%filter(ROIs.region. == "lh_parstriangularis_volume" | ROIs.region. == "rh_parstriangularis_volume")
#
# Comparison: global effects
mod = lmer(TotalVol~TIV+sex+age+(1|eid),df%>%filter(data=="OFAMS"))
effectsize::standardize_parameters(mod)
summary(mod)
mod = lmer(age~TIV+sex+TotalVol+(1|eid),df%>%filter(data=="MS"))
effectsize::standardize_parameters(mod)
summary(mod)$coefficients[4,5]


# MEDIAN EFFECTS for in-text presentation
unstd%>%group_by(data)%>%summarise(Md=median(age.beta),MAD=mad(age.beta))
#
#
# we get a good idea of the strongest effects in MS: putamen, thalamus, superior frontal cortex
# HC age-associations highlight frontal and temporal areas, as well as cerebellum, thalamus and putamen
#
# Check also the ROIs defined by significant effects which replicate across MS datasets
signi = (unstd %>% filter(data=="OFAMS" & p.corr < 0.05))$ROIs.region.[
  (unstd %>% filter(data=="OFAMS" & p.corr < 0.05))$ROIs.region. %in%
  (unstd %>% filter(data=="OUH" & p.corr < 0.05))$ROIs.region.]
#
# show which ROIs present in both hemisphers (here indicated by 2)
table(gsub("Right.","",gsub("Left.","",gsub("lh_","",gsub("rh_","",signi)))))

# Present results in-text
unstd[unstd$ROIs.region. %in% signi,] %>% filter(data == "OUH" & p.corr < 0.05) %>% arrange(age.beta)
unstd[unstd$ROIs.region. %in% signi,] %>% filter(data == "OFAMS" & p.corr < 0.05) %>% arrange(age.beta)
#
#
# Make a function for the estimation of the standard error (SE) and extraction of the effect size (ES)
getit = function(the_data_frame,variable){
  fun=formula(paste(variable,"~age+sex+EstimatedTotalIntraCranialVol+(1|eid)"))
  the_data_frame[variable] = scale(the_data_frame[variable])
  the_data_frame$EstimatedTotalIntraCranialVol = scale(the_data_frame$EstimatedTotalIntraCranialVol)
  the_data_frame$age = scale(the_data_frame$age)
  model = lmer(fun,the_data_frame)
  ES = summary(model)$coefficients[2] # coeff
  SE = summary(model)$coefficients[2,2] # se
  return(c(ES,SE))
}
# Make another function for the estimation of differences in model coefficients
# by Clogg, C. C., Petkova, E., & Haritou, A. (1995). Statistical methods for comparing regression coefficients between models. American Journal of Sociology, 100(5), 1261-1293. 
UKBlong$Left.Thalamus = UKBlong$Left.Thalamus.Proper
UKBlong$Right.Thalamus = UKBlong$Right.Thalamus.Proper
Zdiff = function(ES1,ES2,SE1,SE2){
  Z = (ES1-ES2)/(sqrt( (SE1^2) + (SE2^2) ))
  return(Z)
}
ES_MS = ES_UKB = Zscores = c()
for (i in 1:length(signi)){
  Zscores[i]=Zdiff(getit(df2,signi[i])[1],getit(UKBlong,signi[i])[1],getit(df2,signi[i])[2], getit(UKBlong,signi[i])[2])
  ES_UKB[i] = getit(UKBlong,signi[i])[1]
  ES_MS[i] = getit(df2,signi[i])[1]
}
Zres = data.frame(ROI=signi,Z=Zscores,p=pnorm(Zscores, mean = 0, sd = 1, lower.tail = TRUE), ES_UKB, ES_MS)
for (i in 1:length(signi)){
  Zscores[i]=Zdiff(getit(data_list[[1]],signi[i])[1],getit(UKBlong,signi[i])[1],getit(data_list[[1]],signi[i])[2], getit(UKBlong,signi[i])[2])
  ES_UKB[i] = getit(UKBlong,signi[i])[1]
  ES_MS[i] = getit(data_list[[1]],signi[i])[1]
}
Zres_Oslo = data.frame(ROI=signi,Z=Zscores,p=pnorm(Zscores, mean = 0, sd = 1, lower.tail = TRUE), ES_UKB, ES_MS)
for (i in 1:length(signi)){
  Zscores[i]=Zdiff(getit(data_list[[2]],signi[i])[1],getit(UKBlong,signi[i])[1],getit(data_list[[2]],signi[i])[2], getit(UKBlong,signi[i])[2])
  ES_UKB[i] = getit(UKBlong,signi[i])[1]
  ES_MS[i] = getit(data_list[[2]],signi[i])[1]
}
Zres_Bergen = data.frame(ROI=signi,Z=Zscores,p=pnorm(Zscores, mean = 0, sd = 1, lower.tail = TRUE), ES_UKB, ES_MS)

# show significant ones only
Zres %>% filter(p<.05)
Zres_Oslo %>% filter(p<.05)
Zres_Bergen %>% filter(p<.05)

fix_table = function(Ztab){
ntab = Ztab #%>% filter(p<.05)
# consider cortical stats only
ntab = ntab %>% filter(grepl("volume",ROI))
ntab = ntab[order(ntab$ROI),] #ROIntab = ntab[order(ntab$ROI),] # order the data frame
ntab$hemi = ifelse(grepl("lh_",ntab$ROI)==T,"left","right")
ntab$region = (brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",ntab$ROI)][1:34])[1:nrow(ntab)]
ntab2 = ntab%>%filter(hemi=="right")
ntab2$ROI = gsub("rh_","",ntab2$ROI)
ntab2$region = (brain_regions(dk)[gsub("lh_","",brain_labels(dk)) %in% gsub("_volume","",ntab2$ROI)][1:34])[1:nrow(ntab2)]
ntab$region = c((na.omit(ntab))$region, ntab2$region)
ntab$Z = ifelse(ntab$p>=0.05,0,ntab$Z)
return(ntab)
}
ZpOSL = ggplot(fix_table(Zres_Oslo)) + geom_brain(atlas = dk, aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  theme_void()
ZpBGO = ggplot(fix_table(Zres_Bergen)) + geom_brain(atlas = dk, aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  theme_void()
Zp = ggplot(fix_table(Zres)) + geom_brain(atlas = dk, aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  theme_void()
Zmap = ggarrange(ZpOSL,ZpBGO,Zp, ncol=1,labels=c("OUH data","OFAMS data","Combined data"),
                 hjust = c(0,0,0),common.legend = T, legend = "right")
Zmap = annotate_figure(Zmap, top = text_grob("Volume changes in MS compared to 20 years older healthy controls",face = "bold", size = 17))
ggsave(paste(savepath,"Zmap.pdf",sep=""),Zmap, width = 12, height = 6)
#
fix_table = function(Ztab){
  ntab = Ztab #%>% filter(p<.05)
  # consider cortical stats only
  ntab = ntab %>% filter(grepl("volume",ROI))
  ntab = ntab[order(ntab$ROI),] #ROIntab = ntab[order(ntab$ROI),] # order the data frame
  ntab$hemi = ifelse(grepl("lh_",ntab$ROI)==T,"left","right")
  ntab$region = (brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",ntab$ROI)][1:34])[1:nrow(ntab)]
  ntab2 = ntab%>%filter(hemi=="right")
  ntab2$ROI = gsub("rh_","",ntab2$ROI)
  ntab2$region = (brain_regions(dk)[gsub("lh_","",brain_labels(dk)) %in% gsub("_volume","",ntab2$ROI)][1:34])[1:nrow(ntab2)]
  ntab$region = c((na.omit(ntab))$region, ntab2$region)
  #ntab$Z = ifelse(ntab$p>=0.05,0,ntab$Z)
  return(ntab)
}
Zfix = rbind(fix_table(Zres_Oslo), fix_table(Zres_Bergen), fix_table(Zres))
Zfix$data = c(replicate(nrow(fix_table(Zres_Oslo)), "Oslo"),replicate(nrow(fix_table(Zres_Bergen)), "OFAMS"),replicate(nrow(fix_table(Zres)), "combined"))

Zfix_copy = rbind((Zres_Oslo), (Zres_Bergen), (Zres))
Zfix_copy[2:5] = round(Zfix_copy[2:5],2)
Zfix_copy$data = c(replicate(nrow((Zres_Oslo)), "Oslo"),replicate(nrow((Zres_Bergen)), "OFAMS"),replicate(nrow((Zres)), "combined"))
write.csv(Zfix_copy,paste(savepath,"Ztable.csv",sep=""))
#
# Closer examine the effects
#Zfix[Zfix$ROI %in% signi,] %>% filter(data =="Oslo")
#Zfix[Zfix$ROI %in% signi,] %>% filter(data =="OFAMS")

# show only Bonferroni corrected
Zfix$p.cor = Zfix$p*(length(signi)+80)
#Zfix[Zfix$ROI %in% signi,] %>% filter(data =="Oslo" & p.cor<0.05)
#Zfix[Zfix$ROI %in% signi,] %>% filter(data =="OFAMS" & p.cor<0.05)
Zfix[Zfix$ROI %in% signi,] %>% filter(data =="combined" & p.cor<0.05)
Zfix[Zfix$ROI %in% signi,] %>% filter(data =="OFAMS" & p.cor<0.05)
#
# Check whether thalamus and superior frontal cortex degenerate faster in pwMS compared to HC
# !! 
# For simplicity, we add the volumes of our ROIs into a single volume score
# !!
#
# Z1 = Zdiff(getit(df2,"Thalamus")[1],getit(UKBlong,"Thalamus")[1],getit(df2,"Thalamus")[2], getit(UKBlong,"Thalamus")[2])
#
#
#
# !! 
# For simplicity, we add the volumes of our ROIs into a single volume score for the top-associations
# !!
#
#
# 
# df2$Thalamus = df2$Left.Thalamus.Proper+df2$Right.Thalamus.Proper
# UKBlong$Thalamus = UKBlong$Left.Thalamus.Proper+UKBlong$Right.Thalamus.Proper
# Z1 = Zdiff(getit(df2,"Thalamus")[1],getit(UKBlong,"Thalamus")[1],getit(df2,"Thalamus")[2], getit(UKBlong,"Thalamus")[2])
# #
# # superior frontal cortex
# df2$SFC = df2$lh_superiorfrontal_volume+df2$rh_superiorfrontal_volume
# UKBlong$SFC = UKBlong$lh_superiorfrontal_volume+UKBlong$rh_superiorfrontal_volume
# Z2 = Zdiff(getit(df2,"SFC")[1],getit(UKBlong,"SFC")[1],getit(df2,"SFC")[2], getit(UKBlong,"SFC")[2])
# #
# paste("MS vs HC Thalamic degeneration Z = ",round(Z1,3),", p = ", round(pnorm(Z1, mean = 0, sd = 1, lower.tail = TRUE),3),sep="")
# paste("MS vs HC Superior Frontal Cortex degeneration Z = ",round(Z2,3),", p = ", pnorm(Z2, mean = 0, sd = 1, lower.tail = TRUE),sep="")
# #
# #
# # Plot age associations in these regions
# plotfunc = function(the_data_frame, variable, ylim,data_decriptor){
#   fun=formula(paste(variable,"~age+sex+EstimatedTotalIntraCranialVol+(1|eid)"))
#   mod=lmer(fun,the_data_frame)
#   plot = plot(Effect("age",mod),xlab = "Age",ylab=variable,
#               main=paste("Adjusted association of Age and ", variable, " volume in ", data_decriptor,sep=""),
#               xlim = c(20,70),ylim = ylim)
#   return(plot)
# }
# key_ass = ggarrange(plotfunc(df2,"Thalamus",c(13000,17000),"MS"),plotfunc(UKBlong,"Thalamus",c(13000,17000),"HC"),
#           plotfunc(df2,"SFC",c(35000,50000),"HC"),plotfunc(UKBlong,"SFC",c(35000,50000),"MS"))
# ggsave(paste(savepath,"key_age_ROI_associations.pdf",sep=""),key_ass, width = 12, height = 10)
# 
#
#
# # For comparison: OFAMS data without the 10-year follow-up ####
# # this QC steps allows for a closer mimicking of the follow-up times between datasets
# plist2 = the_reslist2 = list()
# data_list2 = list(df%>%filter(data == "MS"),
#                   df%>%filter(data == "OFAMS")%>%filter(session!=16),
#                   df2%>%filter(session!=16))
# for (i in 1:length(data_list2)){
#   data_frame = data_list2[[i]]
#   # 3.1.1 Cortical #### #
#   ROIs = c(data_frame %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
#            data_frame %>% select(starts_with("rh") & ends_with("volume")) %>% names())
#   #n_cort_ROIs = length(ROIs)
#   reslist=list()
#   for(region in 1:length(ROIs)){
#     f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
#     mod = lmer(f,data_frame)
#     age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
#     CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
#     CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
#     p = summary(mod)$coefficients[4,5]
#     reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
#   }
#   res = list_rbind(reslist)
#   res$p.corr = (length(res$p)+14)*res$p # Bonferroni correction considering also the 14 subcortical areas
#   plot_df = res
#   plot_df = plot_df[order(plot_df$ROIs.region.),] # order the data frame
#   plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
#   plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
#   plot_df01 = plot_df
#   plot_df$age.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$age.beta)
#   p1 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = age.beta),color="black")+
#     #scale_fill_viridis_c(option = "cividis", direction = -1)+
#     scale_fill_gradient2(low = "blue",mid = "white",high="red") +
#     #labs(title="Regional volume loss") + 
#     theme_void()
#   p1 = p1+labs(fill="Std.Coeff.") + #mm<sup>3</sup>/year
#     theme(
#       plot.title = element_markdown(),
#       legend.title = element_markdown()
#     )
#   # 3.1.2 Sub-Cortical #### #
#   # subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
#   ROIs = names(data_frame)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),names(data_frame))]
#   ROIs = ROIs[!grepl(c("CC"),ROIs)]
#   ROIs = ROIs[order(ROIs)]
#   reslist=list()
#   for(region in 1:length(ROIs)){
#     f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
#     mod = lmer(f,data_frame)
#     age.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
#     CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
#     CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
#     p = summary(mod)$coefficients[4,5]
#     reslist[[region]] = data.frame(ROIs[region],age.beta,CI_high,CI_low,p)
#   }
#   res2 = list_rbind(reslist)
#   res2$p.corr = (length(res2$p)+(length(c(data_frame %>% dplyr::select(starts_with("lh") & ends_with("volume")) %>% names(),
#                                           data_frame %>% dplyr::select(starts_with("rh") & ends_with("volume")) %>% names())
#   )))*res2$p # Bonferroni correction considering all ROIs
#   #cerebellum = res2[grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
#   res=rbind(data.frame(res),data.frame(res2))
#   res2 = res2[!grepl(".Cerebellum.Cortex",res2$ROIs.region.),]
#   plot_df1 = res2
#   plot_df1 = plot_df1[order(plot_df1$ROIs.region.),] # order the data frame
#   plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
#   #brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
#   coronal_brain_aseg = as_tibble(aseg) %>%
#     filter(side == "coronal", !grepl("\\d", label))
#   #
#   # Note: the painful steps above leading to this merger are crucial.
#   # This will be useful for later use of ggseg (aseg). Enjoy!
#   plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
#   plot_df02 = plot_df1
#   #
#   plot_df1$age.beta = ifelse(plot_df1$p.corr > .05, NA,plot_df1$age.beta)
#   #
#   p1.1 = ggplot(plot_df1) + geom_brain(atlas = aseg, side = "coronal",aes(fill = age.beta),color="black")+
#     #scale_fill_viridis_c(option = "cividis", direction = -1)+
#     scale_fill_gradient2(low = "blue",mid = "white",high="red") +
#     #labs(title="Regional volume loss") + 
#     theme_void()
#   p1.1 = p1.1+labs(fill="Std.Coeff.") +
#     theme(
#       plot.title = element_markdown(),
#       legend.title = element_markdown()
#     )
#   p = ggarrange(p1,p1.1,nrow=1,widths=c(2,.5))
#   the_reslist2[[i]] = res
#   plist2[[i]] = p
# }
# age.p = ggarrange(plotlist = plist,ncol = 1,labels=c("a","b","c"))
# age.p = annotate_figure(age.p, top = text_grob("Annual regional volume loss",face = "bold", size = 17))
# ggsave(paste(savepath,"age_plot_SHORT_standardized.pdf",sep=""),age.p, width = 14, height = 6)
# #safe table of the results
# unstd_SHORT = list_rbind(the_reslist)
# unstd_SHORT$data = c(replicate(nrow(unstd_SHORT)/length(data_list),"Oslo"),
#                replicate(nrow(unstd_SHORT)/length(data_list),"OFAMS"),
#                replicate(nrow(unstd_SHORT)/length(data_list),"both"))
# write.csv(file = paste(savepath,"age_beta_SHORT_standardized.csv",sep=""),unstd_SHORT)
# #
# #
# #
# #
# #
# #
# #
# #
# 4. Regional EDSS associations ####
#   # 4.1 UnStandardized ####
  # 4.2 Standardized ####
    # 4.2.1 Cortical ####
# EDSS
edss.tab = edss.plot = list()
for (i in 1:3){
  df = data_list[[i]]
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res = list_rbind(reslist)
plot_df = res[order(res$ROIs.region.),] # order the data frame
plot_df$region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df$ROIs.region.)][1:34]
plot_df$hemi = ifelse(grepl("lh_",plot_df$ROIs.region.)==T,"left","right")
plot_df01 = plot_df
#plot_df$std.beta = ifelse(plot_df$p*46 > .05, NA,plot_df$std.beta) # Bonf corrected are displayed
p03 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.2,0)) +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p03 = p03+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
    # 4.2.2 Sub-Cortical ####
ROIs = names(df)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC"),ROIs)]
#ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+edss+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res = list_rbind(reslist)
plot_df1 = res[order(res$ROIs.region.),] # order the data frame
plot_df1$label = brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))]
#brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",plot_df1$ROIs.region.)][1:34]
coronal_brain_aseg = as_tibble(aseg) %>%
  filter(side == "coronal", !grepl("\\d", label))
plot_df1 = merge(plot_df1,coronal_brain_aseg,by="label") # merge the data with the atlas labels
plot_df02 = plot_df1
plot_df02$hemi = ifelse(grepl("lh_",plot_df02$ROIs.region.)==T,"left","right")
#plot_df02$std.beta = ifelse(plot_df02$p*43 > .05, NA,plot_df02$std.beta) # Bonferroni correction by nb regions
# plot it
p04 = ggplot(plot_df02) + geom_brain(atlas = aseg,side = "coronal",aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.2,0)) +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p04 = p04+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
edss.plot[[i]] = ggarrange(p03,p04,ncol=2,widths=c(2,.5))#,labels=c("c","d"))
plot_df02 = plot_df02 %>% select(names(plot_df01))
edss.tab[[i]] = rbind(plot_df01,plot_df1[,colnames(plot_df1) %in% colnames(plot_df01)])
}
edss.tab[[1]]$p.cor = edss.tab[[1]]$p*82
edss.tab[[2]]$p.cor = edss.tab[[2]]$p*82

r1 = edss.tab[[1]] %>% filter(p<0.05/82) %>% pull(ROIs.region.)
r2 = edss.tab[[2]] %>% filter(p<0.05/82) %>% pull(ROIs.region.)

# THESE ARE THE REGIONS WHICH SURVE MULTIPLE COMPARISON:
r1[r1 %in% r2]

# More lenient?
r1 = edss.tab[[1]] %>% filter(p<0.05) %>% pull(ROIs.region.)
r2 = edss.tab[[2]] %>% filter(p<0.05) %>% pull(ROIs.region.)

# THESE ARE THE REGIONS WHICH SURVE MULTIPLE COMPARISON:
r1[r1 %in% r2]

# These are the model parameters
edss.tab[[1]][edss.tab[[1]]$ROIs.region. %in% r1[r1 %in% r2],]
edss.tab[[2]][edss.tab[[1]]$ROIs.region. %in% r1[r1 %in% r2],]


# Oslo
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+edss+age+(1|eid),data_list[[1]])
c(a$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
# OFAMS
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+edss+age+(1|eid),data_list[[2]])
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
#
# 5. Create outputs ####
#
  # 5.1 Merge plots ####
#
#
# standardized plots
#edss.plot = annotate_figure(edss.plot, top = text_grob("Associations of EDSS and regional brain volume",face = "bold", size = 17))
#plot1 = ggarrange(age.plot, edss.plot, nrow=2)
#edss.plot0 = ggarrange(plotlist = edss.plot,ncol=1,labels=c("Oslo cohort","OFAMS cohort","Combined cohort"))
#ggsave(paste(savepath,"EDSS_plot_standardized.pdf",sep=""),edss.plot0, width = 13, height = 6)
#
# # Merge tables
# # show significant correlates of EDSS
# edss.tab[[2]][edss.tab[[2]]$ROIs.region. %in% signi,] %>% filter(p<.05/(length(signi)+80))
# edss.tab[[3]][edss.tab[[3]]$ROIs.region. %in% signi,] %>% filter(p<.05/(length(signi)+80))
# # more lenient alpha level
# edss.tab[[2]][edss.tab[[2]]$ROIs.region. %in% signi,] %>% filter(p<.005)
# edss.tab[[3]][edss.tab[[3]]$ROIs.region. %in% signi,] %>% filter(p<.005)
# 
# 
# edss.tab[[1]]$p.corr = edss.tab[[1]]$p*(length(signi)+80)
# edss.tab = edss.tab[[1]]
# edss.tab[edss.tab$ROIs.region. %in% signi,] %>% filter(p.corr<.05)#filter(p<.05/(length(signi)+83))
# edss.tab[edss.tab$ROIs.region. %in% signi,] %>% filter(p<.05)
# edss.tab %>% filter(p.corr<.05)#filter(p<.05/(length(signi)+83))
# edss.tab_copy = edss.tab
# edss.tab_copy[2:5] = round(edss.tab_copy[2:5],2)
# edss.tab_copy[8] = round(edss.tab_copy[8],2)
# 
#   # 5.2 Save tables ####
# # write.csv(age.tab,file = paste(savepath,"age_tab.csv",sep=""))
# write.csv(edss.tab_copy,file = paste(savepath,"edss_tab.csv",sep=""))
# #
#
# Some demographics (Age & EDSS)

# First visit OUH
demodesc = data_list[[1]] %>% 
  group_by(eid) %>% 
  summarise(session = min(session)) %>% 
  slice_max(session, n = length(unique(data_list[[1]]$eid)))
demodesc = merge(demodesc,data_list[[1]],by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% (data_list[[1]])$eid,]
paste("OUH first visit descriptives")
table(demodesc$sex)
paste("Age: ",round(mean(demodesc$age),1),"±",round(sd(demodesc$age),1),sep="")
paste("EDSS: ",round(mean(demodesc$edss),1),"±",round(sd(demodesc$edss),1),sep="")


# Final visit OUH
demodesc = data_list[[1]] %>% filter(session != 1) %>%
  group_by(eid) %>% 
  summarise(session = max(session)) %>% 
  slice_max(session, n = length(unique(data_list[[1]]$eid)))
demodesc = merge(demodesc,data_list[[1]],by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% (data_list[[1]] %>% filter(session == 1))$eid,]
paste("OUH Last visit descriptives")
table(demodesc$sex)
paste("Age: ",round(mean(demodesc$age),1),"±",round(sd(demodesc$age),1),sep="")
paste("EDSS: ",round(mean(demodesc$edss),1),"±",round(sd(demodesc$edss),1),sep="")

# First visit OFAMS
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
demodesc = data_list[[2]] %>% 
  group_by(eid) %>% 
  summarise(session = min(session)) %>% 
  slice_max(session, n = length(unique(data_list[[1]]$eid)))
demodesc = merge(demodesc,data_list[[2]],by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% (data_list[[2]])$eid,]
demodesc = demodesc[demodesc$eid %in% demo10[ifelse(demo10$Missing_reason == "Trakk seg",F,T),]$Patnr,]
paste("OFAMS first visit descriptives")
table(demodesc$sex)
paste("Age: ",round(mean(demodesc$age),1),"±",round(sd(demodesc$age),1),sep="")
paste("EDSS: ",round(mean(demodesc$edss),1),"±",round(sd(demodesc$edss),1),sep="")

# Last visit OFAMS
demodesc = data_list[[2]] %>% 
  group_by(eid) %>% 
  summarise(session = max(session)) %>% 
  slice_max(session, n = length(unique(data_list[[2]]$eid)))
demodesc = merge(demodesc,data_list[[2]],by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% demo10[ifelse(demo10$Missing_reason == " ",T,F),]$Patnr,]

paste("OFAMS Last visit descriptives")
table(demodesc$sex)
paste("Age: ",round(mean(demodesc$age),1),"±",round(sd(demodesc$age),1),sep="")
paste("EDSS: ",round(mean(demodesc$edss),1),"±",round(sd(demodesc$edss),1),sep="")

#
#
# 6. Overlap between faster ageing regions and EDSS-related regions ####
# overlap between significant ageing regions and EDSS
EDSS.signi = edss.tab %>% filter(p<.05/43) %>% pull(ROIs.region.)
EDSS.signi[EDSS.signi %in% signi]


# the following regions overlap between
edss_rel = edss.tab %>% filter(p<.05/43) %>% pull(ROIs.region.)
# let's look at the overlap between faster ageing than older folks regions
Z_vars = (Zres%>% filter(p<.05) %>% pull(ROI))
Z_vars[Z_vars %in% edss_rel]
# this can be plotted (if there is anything)
plot_df999=plot_df
plot_df999$std.beta = ifelse(plot_df$ROIs.region. == "lh_insula_volume",NA,plot_df$std.beta)
plot_df999$std.beta = ifelse(plot_df999$ROIs.region. == "lh_rostralmiddlefrontal_volume",NA,plot_df999$std.beta)
p033 = ggplot(plot_df999) + geom_brain(atlas = dk,aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.2,0)) +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p033 = p033+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
p033
# as well as the overlap of generally ageing regions and EDSS correlated regions
age_vars = unstd %>% filter(data=="both" & p.corr<0.05) %>% pull(ROIs.region.)
age_vars[age_vars %in% edss_rel]
# this can be plotted
plot_df999=plot_df
plot_df999$std.beta = ifelse(plot_df$ROIs.region. == "lh_insula_volume",NA,plot_df$std.beta)
plot_df999$std.beta = ifelse(plot_df999$ROIs.region. == "lh_entorhinal_volume",NA,plot_df999$std.beta)
p033 = ggplot(plot_df999) + geom_brain(atlas = dk,aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-0.2,0)) +
  #labs(title="Regional association between EDSS and volume") + 
  theme_void()
p033 = p033+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
p033
#
#
# Check fatigue and PASAT in the OFAMS data
# Fatigue ####
fati = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/fatigue.sas7bdat") # fatigue scores
fati_OSL = read.csv("/Users/max/Documents/Local/Data/Oslo/fatigue_Oslo.csv")
#
# prep eid and session
fati_OSL$session = substr(fati_OSL$eid,9,11)
fati_OSL$eid = substr(fati_OSL$eid,1,7)
fati_OSL = fati_OSL %>% select(eid,session,fatigue)
fati_OSL$session = as.numeric(fati_OSL$session)
fati_OSL = merge(msOSL,fati_OSL,by=c("eid","session"))
fati.cop = fati_OSL

# fix session factor levels
fati$session = factor(fati$VISIT)
fati$session = ifelse(fati$session == "Baseline",1,fati$session)
fati$session = ifelse(fati$session == 4,7,fati$session)
fati$session = ifelse(fati$session == 2,13,fati$session)
fati$session = ifelse(fati$session == 3,25,fati$session)
fati$eid = as.numeric(fati$patno)
fati$fatigue = fati %>% select(A,     B,     C,     D,     E,     F,     G,     H,     I) %>% rowMeans()
fati = fati %>% select(eid,session,fatigue)
# start with fatigue
df = merge(msBGO,fati,by = c("eid","session"))
df = df[df$eid %in% demo10[ifelse(demo10$Missing_reason == "Trakk seg",F,T),]$Patnr,]

paste("OFAMS first visit descriptives")
paste("FSS: ",round(mean(df%>% filter(session==25)%>%pull(fatigue)),1),"±",round(sd(df%>% filter(session==25)%>%pull(fatigue)),1),sep="")
paste("FSS: ",round(mean(df%>% filter(session==1)%>%pull(fatigue)),1),"±",round(sd(df%>% filter(session==1)%>%pull(fatigue)),1),sep="")
paste("FSS: ",round(mean(df%>% filter(session==1)%>%pull(fatigue)),1),"±",round(sd(df%>% filter(session==1)%>%pull(fatigue)),1),sep="")
df%>% filter(session==1) %>% nrow
df%>% filter(session==25) %>% nrow

# first visit OUH
demodesc = fati_OSL %>% 
  group_by(eid) %>% 
  summarise(session = min(session)) %>% 
  slice_max(session, n = length(unique(fati_OSL$eid)))
demodesc = merge(demodesc,fati_OSL,by = c("eid","session"))
paste("FSS: ",round(mean(na.omit(demodesc$fatigue)),1),"±",round(sd(na.omit(demodesc$fatigue)),1),sep="")

# last visit OUH
demodesc = fati_OSL %>% filter(session != 1) %>%
  group_by(eid) %>% 
  summarise(session = max(session)) %>% 
  slice_max(session, n = length(unique(fati_OSL$eid)))
demodesc = merge(demodesc,fati_OSL,by = c("eid","session"))
paste("FSS: ",round(mean(na.omit(demodesc$fatigue)),1),"±",round(sd(na.omit(demodesc$fatigue)),1),sep="")

# longitudinal Combat
LC = function(dat){
  features = c(dat %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
             dat %>% select(starts_with("rh") & ends_with("volume")) %>% names(),
             names(dat)[grepl(c("halamus|allidum|mygdala|
                                campus|utamen|audate|
                                CC|Cerebellum.Cortex"),names(dat))],
             dat %>% select(TotalGrayVol, EstimatedTotalIntraCranialVol,
                            ends_with("Hippocampus")) %>% names())
  dat = dat %>% select(all_of(features),eid,data,sex,session,scanner, age, fatigue)%>%na.omit
  covars = dat %>% dplyr::select(sex,fatigue,data,age)
  dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
                 features = features,
                 formula = "age + sex", ranef = "(1|eid)", data = dat)
  dat = dat$data_combat
  dat = cbind(covars,dat)
  colnames(dat) = gsub('.combat','',colnames(dat))
  return(dat)
}
df = LC(df)
fati_OSL = fati_OSL %>% filter(scanner!=2)
fati_OSL = LC(fati_OSL)

# start with BGO sample
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
           df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
    mod = lmer(f,df)
    std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
  }
res1 = list_rbind(reslist)
res1 %>% filter(p<0.05/(length(signi)+83))
res1[res1$ROIs.region. %in% signi,] %>% filter(p<0.05)
  ROIs = names(df)
  ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
  #ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
  ROIs = ROIs[order(ROIs)]
  reslist=list()
for(region in 1:length(ROIs)){
    f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
    mod = lmer(f,df)
    std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
    CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
    CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
    p = summary(mod)$coefficients[4,5]
    reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
  }
res = list_rbind(reslist)
BGO.fattab = rbind(res1,res)
BGO.fattab_copy = BGO.fattab
BGO.fattab_copy[2:5] = round(BGO.fattab_copy[2:5],2)

# now, the OSL sample
df=fati_OSL
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res1 = list_rbind(reslist)
res1 %>% filter(p<0.05/(length(signi)+83))
res1[res1$ROIs.region. %in% signi,] %>% filter(p<0.05)
ROIs = names(df)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
#ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
Osl.fat.tab = rbind(res1,res)
Osl.fat.tab_copy = Osl.fat.tab
Osl.fat.tab_copy[2:5] = round(Osl.fat.tab_copy[2:5],2)

Osl.fat.tab %>% filter(p<0.05/80)
Osl.fat.tab[Osl.fat.tab$ROIs.region. %in% signi,] %>% filter(p<0.05)
BGO.fattab %>% filter(p<0.05/(length(signi)+80))
BGO.fattab[BGO.fattab$ROIs.region. %in% signi,] %>% filter(p<0.05)

fat.tab_copy = rbind(BGO.fattab_copy,Osl.fat.tab_copy)
fat.tab_copy$Data = c(replicate(nrow(BGO.fattab_copy),"OFAMS"),replicate(nrow(BGO.fattab_copy),"Oslo"))
#
#
#
# There is no overlap between the effects of fatigue in brain volumes between samples.
# Hence, we analyse them together.
df = merge(msBGO,fati,by = c("eid","session"))
df$scanner = as.numeric(factor(df$scanner))
fati.cop = fati.cop %>% filter(scanner!=2)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
cort_ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
df$data = "OFAMS"
fati.cop$data = "Oslo"
fati_all=
  rbind(df%>%select(eid,session,scanner,age,sex,data,fatigue,TotalGrayVol, EstimatedTotalIntraCranialVol,
                    all_of(ROIs),all_of(cort_ROIs)),
        fati.cop%>%select(eid,session,scanner,age,sex,data,TotalGrayVol, EstimatedTotalIntraCranialVol,
                          fatigue,all_of(ROIs),all_of(cort_ROIs)))
fati_all = LC(fati_all)
# 
#
#
# fatigue in both samples together
df=fati_all
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res1 = list_rbind(reslist)
res1 %>% filter(p<0.05/(length(signi)+83))
res1[res1$ROIs.region. %in% signi,] %>% filter(p<0.05)
ROIs = names(df)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
#ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
all.fattab = rbind(res1,res)
all.fattab %>% filter(p<0.05/(80))
all.fattab[all.fattab$ROIs.region. %in% signi,] %>% filter(p<0.05)
all.fattab_copy = all.fattab
all.fattab_copy[2:5] = round(all.fattab_copy[2:5],2)
write.csv(all.fattab_copy,file = paste(savepath,"fatigue_tab.csv",sep=""))
#
#
# Check also global effects
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid),fati_all)
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid),fati_all%>%filter(data=="OFAMS"))
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+fatigue+age+(1|eid),fati_all%>%filter(data=="Oslo"))
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
#
#
# PASAT ####
# zscored PASAT
PASAT_OSL = read.csv("/Users/max/Documents/Local/Data/Oslo/Database_cognitive_phenotypes_MS_Oslo.csv")
#
# OSL
PASAT_OSL$eid = gsub("MS","MS_",PASAT_OSL$subject_id)
PASAT_OSL$session = PASAT_OSL$tpoint
PASAT_OSL$PASAT = PASAT_OSL$MACFIMS_PASAT3_zscore
PASAT_OSL = PASAT_OSL %>% select(eid,session,PASAT)
PASAT_OSL = merge(PASAT_OSL,msOSL,by=c("eid","session"))

demodesc = PASAT_OSL %>% 
  group_by(eid) %>% 
  summarise(session = min(session)) %>% 
  slice_max(session, n = length(unique(PASAT_OSL$eid)))
demodesc = merge(demodesc,PASAT_OSL,by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% (PASAT_OSL)$eid,]
paste("OUH first visit descriptives")
paste("PASAT: ",round(mean(na.omit(demodesc$PASAT)),1),"±",round(sd(na.omit(demodesc$PASAT)),1),sep="")
paste("N = ", length(demodesc$PASAT),sep="")

demodesc = PASAT_OSL %>% 
  group_by(eid) %>% 
  summarise(session = max(session)) %>% 
  slice_max(session, n = length(unique(PASAT_OSL$eid)))
demodesc = merge(demodesc,PASAT_OSL,by = c("eid","session"))
demodesc = demodesc[demodesc$eid %in% (PASAT_OSL)$eid,]
paste("OUH last visit descriptives")
paste("PASAT: ",round(mean(na.omit(demodesc$PASAT)),1),"±",round(sd(na.omit(demodesc$PASAT)),1),sep="")
paste("N = ", length(demodesc$PASAT),sep="")

# OFAMS
pasat1 = demo10%>%dplyr::select(Patnr,BL_PASATcorrect,PASAT_24M, PASAT_OFAMS10)
pasat1 = reshape2::melt(pasat1, id.vars = c("Patnr"))
names(pasat1) = c("eid","session","PASAT")
pasat1$session = ifelse(pasat1$session == "BL_PASATcorrect",1,0)+ifelse(pasat1$session == "PASAT_24M",25,0)+ifelse(pasat1$session == "PASAT_OFAMS10",145,0)
pasat1$session = pasat1$session+1
df = merge(msBGO,pasat1,by = c("eid","session"),all = T)
df = df[df$eid %in% demo10[ifelse(demo10$Missing_reason == "Trakk seg",F,T),]$Patnr,]
df%>%filter(session==2)%>%summarise(M=mean(na.omit(PASAT)), SD=sd(na.omit(PASAT)),N=length(unique(eid)))
df%>%filter(session==146)%>%summarise(M=mean(na.omit(PASAT)), SD=sd(na.omit(PASAT)),N=length(unique(eid)))
df%>%filter(session==2) %>% select(PASAT) %>% na.omit() %>%nrow
df%>%filter(session==146) %>% select(PASAT) %>% na.omit() %>%nrow

df$scanner = as.numeric(factor(df$scanner))
# longitudinal Combat
LC = function(dat){
  features = c(dat %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
               dat %>% select(starts_with("rh") & ends_with("volume")) %>% names(),
               names(dat)[grepl(c("halamus|allidum|mygdala|
                                campus|utamen|audate|
                                CC|Cerebellum.Cortex"),names(dat))],
               dat %>% select(TotalGrayVol, EstimatedTotalIntraCranialVol,
                              ends_with("Hippocampus")) %>% names())
  dat = dat %>% select(all_of(features),eid,data,sex,session,scanner, age, PASAT)%>%na.omit
  covars = dat %>% dplyr::select(sex,PASAT,data,age)
  dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
                   features = features,
                   formula = "age + sex", ranef = "(1|eid)", data = dat)
  dat = dat$data_combat
  dat = cbind(covars,dat)
  colnames(dat) = gsub('.combat','',colnames(dat))
  return(dat)
}
df = LC(df)
PASAT_OSL = LC(PASAT_OSL)

# first, estimate effects in OFAMS data
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res1 = list_rbind(reslist)
ROIs = names(df)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res = list_rbind(reslist)
PASAT.tab = rbind(res1,res)
PASAT.tab_copy = PASAT.tab
PASAT.tab_copy[2:5] = round(PASAT.tab_copy[2:5],2)
#
#
# Second, doPASAT_OSL the same in OSL data
ROIs = c(PASAT_OSL %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         PASAT_OSL %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,PASAT_OSL)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res1 = list_rbind(reslist)
ROIs = names(PASAT_OSL)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,PASAT_OSL)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res = list_rbind(reslist)
PASAT.tab_OSL = rbind(res1,res)
PASAT.tab_OSL_copy = PASAT.tab_OSL
PASAT.tab_OSL_copy[2:5] = round(PASAT.tab_OSL_copy[2:5],2)


PASAT.tab %>% filter(p<0.05/(length(signi)+83)) #OFAMS
PASAT.tab_OSL %>% filter(p<0.05/(length(signi)+83)) #Oslo

PASAT.tab_OSL[PASAT.tab_OSL$ROIs.region. %in% signi,] %>% filter(p<0.005)
PASAT.tab[PASAT.tab$ROIs.region. %in% signi,] %>% filter(ROIs.region. == "rh_lateralorbitofrontal_volume")

PASAT.tab[PASAT.tab$ROIs.region. %in% signi,] %>% filter(p<0.005)
PASAT.tab_OSL[PASAT.tab_OSL$ROIs.region. %in% signi,] %>% filter(ROIs.region. == "Right.Thalamus")



PASAT.tab_copy = rbind(PASAT.tab_copy,PASAT.tab_OSL_copy)

PASAT.tab_OSL %>% filter(p<0.005) 
PASAT.tab %>% filter(p<0.005) 

#
#
# Now, analyse both datasets together.
PASAT = rbind(df,PASAT_OSL)
ROIs = c(PASAT %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         PASAT %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist = list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,PASAT)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res1 = list_rbind(reslist)
ROIs = names(PASAT)
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid)"))
  mod = lmer(f,PASAT)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low,p)
}
res = list_rbind(reslist)
PASAT.all.tab = rbind(res1,res)
PASAT.all.tab %>% filter(p<0.05/(length(signi)+83))
PASAT.all.tab[PASAT.all.tab$ROIs.region. %in% signi,] %>% filter(p<0.05)

# Check also global effects
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid),PASAT)
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid),PASAT%>%filter(data=="OFAMS"))
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
mod = lmer(TotalGrayVol~EstimatedTotalIntraCranialVol+sex+PASAT+age+(1|eid),PASAT%>%filter(data=="MS"))
c(effectsize::standardize_parameters(mod)$Std_Coefficient[4],
  effectsize::standardize_parameters(mod)$CI_high[4],
  effectsize::standardize_parameters(mod)$CI_low[4],
  summary(mod)$coefficients[4,5])
#
#
#
# Power ####
################### for highlighted AGE associations:
## Superior frontal
model = lmer(lh_superiorfrontal_volume~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(rh_superiorfrontal_volume~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(lh_superiorfrontal_volume~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
model = lmer(rh_superiorfrontal_volume~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
## pars orbitalis 
model = lmer(lh_parsorbitalis_volume~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(rh_parsorbitalis_volume~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(lh_parsorbitalis_volume~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
model = lmer(rh_parsorbitalis_volume~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)

## thalamus
model = lmer(Left.Thalamus~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(Right.Thalamus~age+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(Left.Thalamus~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
model = lmer(Right.Thalamus~age+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)

################### for highlighted EDSS associations:
model = lmer(lh_parsorbitalis_volume~edss+EstimatedTotalIntraCranialVol+sex+age+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(rh_parsorbitalis_volume~edss+EstimatedTotalIntraCranialVol+sex+age+(1|eid),data_list[[1]])
powerSim(model,nsim=100)
model = lmer(lh_parsorbitalis_volume~edss+EstimatedTotalIntraCranialVol+sex+age+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
model = lmer(rh_parsorbitalis_volume~edss+EstimatedTotalIntraCranialVol+sex+age+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
#
model = lmer(lh_superiorfrontal_volume~edss+TIV+sex+(1|eid),data_list[[1]])
powerSim(model,nsim=1000)
model = lmer(lh_superiorfrontal_volume~edss+TIV+sex+(1|eid),data_list[[2]])
powerSim(model,nsim=100)
#
############# Fatigue
df = merge(msBGO,fati,by = c("eid","session"))
df$scanner = as.numeric(factor(df$scanner))
model = lmer(rh_caudalmiddlefrontal_volume~fatigue+EstimatedTotalIntraCranialVol+sex+(1|eid),df)
powerSim(model,nsim=1000)
model = lmer(lh_precuneus_volume~fatigue+EstimatedTotalIntraCranialVol+sex+(1|eid),df)
powerSim(model,nsim=1000)

model = lmer(rh_lateralorbitofrontal_volume~PASAT+EstimatedTotalIntraCranialVol+sex+(1|eid),PASAT%>%filter(data=="MS"))
powerSim(model,nsim=1000)
model = lmer(Left.Thalamus~PASAT+EstimatedTotalIntraCranialVol+sex+(1|eid),PASAT%>%filter(data=="OFAMS"))
powerSim(model,nsim=1000)

#
# Done.