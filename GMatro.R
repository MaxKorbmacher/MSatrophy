# Analyse GM atrophy among OFAMS clinical trial MS patients
# Max Korbmacher, 20.12.2024
# last change: 03 Jan 2025
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
vol.m = lmer(TotalVol~session+TIV+sex+(1|eid),df)
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
###### THAT MEANS THAT ONE CAN ALSO USE THE DATA UP TO 2 YEARS TO MODEL DEGENERATIVE PROCESSES IN CASE WE LIMIT THE RESEARCH TO THE OMICS SAMPLING PERIOD

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
# contains sex,age,vol, no interaction effects
# 3. Regional age associations ####
  # 3.1 Unstandardized ####
    # 3.1.1 Cortical ####
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
n_cort_ROIs = length(ROIs)
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,df)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = (length(res$p)+12)*res$p # Bonferroni correction considering also the 12 subcortical areas
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
p1 = p1+labs(fill="mm<sup>3</sup>/year") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
#
    # 3.1.2 Sub-Cortical ####
# subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
ROIs = names(df)[2:46]
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,df)
  age.beta = summary(mod)$coefficients[4]
  SE = summary(mod)$coefficients[4,2]
  t = summary(mod)$coefficients[4,4]
  p = summary(mod)$coefficients[4,5]
  reslist[[region]] = data.frame(ROIs[region],age.beta,SE,t,p)
}
res = list_rbind(reslist)
res$p.corr = (length(res$p)+n_cort_ROIs)*res$p # Bonferroni correction considering all ROIs
plot_df1 = res
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
p1.1 = p1.1+labs(fill="mm<sup>3</sup>/year") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
# create output table with all values in it 
#age.tab = rbind(plot_df01 %>% select(region,hemi,age.beta,SE,std.beta,CI_low,CI_high, t,p,p.corr),plot_df02 %>% select(region,hemi,age.beta,SE,std.beta,CI_low,CI_high,t,p,p.corr))
#
  # 3.2 Standardized ####
    # 3.1.1 Cortical ####
ROIs = c(df %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
         df %>% select(starts_with("rh") & ends_with("volume")) %>% names())
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low)
}
res = list_rbind(reslist)
plot_df01 = merge(plot_df01,res,by="ROIs.region.")
plot_df = plot_df01
plot_df$std.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$std.beta)
#
p01 = ggplot(plot_df) + geom_brain(atlas = dk,aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional volume loss") + 
  theme_void()
p01 = p01+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
#
    # 3.1.2 Sub-Cortical ####
# subcortical: thalamus, pallidum, amygdala, hippocampus, putamen, accumbens area, and caudate nucleus
ROIs = names(df)[2:46]
ROIs = ROIs[grepl(c("halamus|allidum|mygdala|campus|utamen|audate|CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[!grepl(c("CC|Cerebellum.Cortex"),ROIs)]
ROIs = ROIs[order(ROIs)]
reslist=list()
for(region in 1:length(ROIs)){
  f = formula(paste(ROIs[region],"~TIV+sex+age+(1|eid)"))
  mod = lmer(f,df)
  std.beta = effectsize::standardize_parameters(mod)$Std_Coefficient[4]
  CI_high = effectsize::standardize_parameters(mod)$CI_high[4]
  CI_low = effectsize::standardize_parameters(mod)$CI_low[4]
  reslist[[region]] = data.frame(ROIs[region],std.beta,CI_high,CI_low)
}
res = list_rbind(reslist)
plot_df = merge(plot_df02,res,by="ROIs.region.")
plot_df02 = plot_df
plot_df$std.beta = ifelse(plot_df$p.corr > .05, NA,plot_df$std.beta)
#
p02 = ggplot(plot_df) + geom_brain(atlas = aseg,side = "coronal",aes(fill = std.beta),color="black")+
  #scale_fill_viridis_c(option = "cividis", direction = -1)+
  scale_fill_gradient2(low = "blue",mid = "white",high="red") +
  #labs(title="Regional volume loss") + 
  theme_void()
p02 = p02+labs(fill="Std. Beta") +
  theme(
    plot.title = element_markdown(),
    legend.title = element_markdown()
  )
#
# create output table with all values in it 
age.tab = rbind(plot_df01 %>% select(region,hemi,age.beta,SE,std.beta,CI_high,CI_low,t,p,p.corr),plot_df02 %>% select(region,hemi,age.beta,SE,std.beta,CI_high,CI_low,t,p,p.corr))

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
res$p.corr = length(res$p)*res$p # p.adjust(res$p, method="BH")
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
age.plot = ggarrange(p1,p1.1,ncol=2,widths=c(2,.5),labels=c("a","b"))
age.plot = annotate_figure(age.plot, top = text_grob("Annual regional volume loss",face = "bold", size = 17))
edss.plot = ggarrange(p2,p2.2,ncol=2,widths=c(2,.5),labels=c("c","d"))
edss.plot = annotate_figure(edss.plot, top = text_grob("Associations of EDSS and regional brain volume",face = "bold", size = 17))
plot1 = ggarrange(age.plot, edss.plot, nrow=2)
ggsave(paste(savepath,"reg_plot.pdf",sep=""),plot1, width = 14, height = 6)
#
#
# Now, the standardized plots
age.plot = ggarrange(p01,p02,ncol=2,widths=c(2,.5),labels=c("a","b"))
age.plot = annotate_figure(age.plot, top = text_grob("Annual regional volume loss",face = "bold", size = 17))
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
