# Data wrangling chaos: put multimodal data together
# Max Korbmacher, 20.12.2024
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# load pkgs
pacman::p_load(haven,rlist, dplyr,tidyr,reshape2)
# load demographics
test=read_sas("/Users/max/Documents/Local/MS/data/screening.sas7bdat")
demg = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/MRI_A_E_D.sav")
demo = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/Demographics.sav") # Demographics for inclusion age
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
edss = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/edss.sav") # contains Expanded Disability Status Scale (EDSS) scores
# load brain stats
icv1 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats/aparc_stats_area_lh_20210816145055.csv")
icv2 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220214_update/manual_edits_batch_20210913/aparc_stats_area_lh_20220214103330.csv")
#icv3 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_wmedits_rerun_20211201/aparc_stats_area_lh_20220209170109.csv")
icv4 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_no3T_pre-m120/stats_cross_20220209/manual_edits_batch_20210913/aparc_stats_area_lh_20220209165901.csv")
icv5 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210127_stats/aparc_stats_area_lh_20210210092415.csv")
icv6 = read.csv("/Users/max/Documents/Local/freesurfer/freesurfer_7.1.1_recon-all_3T_m120/batch_20210304/stats/aparc_stats_area_lh_20210315160653.csv")
# read in already established cortical stats
cortical = read.csv("/Users/max/Documents/Local/MS/data/final_cort_dat.csv")
# load FS stats
filenames = list.files("/Users/max/Documents/Local/MS/tabular_subcort/", pattern = ".csv", full.names = T)
ldf = lapply(filenames, read.csv)
ldf = list.rbind(ldf)
# Lesion info
lesion_count = read.csv("/Users/max/Documents/Local/MS/data/lesion_count.csv")
lesion_vol = read.csv("/Users/max/Documents/Local/MS/data/lesion_vol.csv")

# read in serum, geno markers
blood = read.csv("/Users/max/Documents/Local/MS/data/Copy of Paired_serum-MRI.csv")
geno = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/hla_drb1_15_Juni_2011_alle pas.sav") # genotype 

###
# save the merged FreeSurfer table
write.csv(ldf,"/Users/max/Documents/Local/MS/data/FSsubc_data.csv",row.names = F)
#
# add age to the demographics
test$age = as.numeric(difftime(test$SC_DATEOFVISIT,test$DATEOFBIRTH)/365)
test0 = test %>% select(patno,Sex)
names(test0) = c("eid","sex")
#test = test %>% select(age,Sex)
#
#
#
#
########### AGE CALC
#tmp=ldf%>% filter(grepl("baseline",eid)) # select the ids we have available
#
tmp=ldf # create tmp df
session=ifelse(grepl("m24",tmp$eid),"24",tmp$eid)
session=ifelse(grepl("m120",session),"144",session)
session=ifelse(grepl("m12",session),"12",session)
session=ifelse(grepl("m11",session),"11",session)
session=ifelse(grepl("m10",session),"10",session)
session=ifelse(grepl("m9",session),"9",session)
session=ifelse(grepl("m8",session),"8",session)
session=ifelse(grepl("m7",session),"7",session)
session=ifelse(grepl("m6",session),"6",session)
session=ifelse(grepl("m5",session),"5",session)
session=ifelse(grepl("m4",session),"4",session)
session=ifelse(grepl("m3",session),"3",session)
session=ifelse(grepl("m2",session),"2",session)
session=ifelse(grepl("m1",session),"1",session)
tmp$session=as.numeric(ifelse(grepl("baseline",session),"0",session))
tmp$eid = as.numeric(gsub("\\_.*","",gsub("sub-","",tmp$eid))) # standardize ids
test$eid = test$patno
demg$eid=demg$Patnr # equal to test$eid for merging
#
# identify age for each visit using demg$Visitnr and demg$Visit_date
months = levels(factor(demg$Visitnr))
age_list = list()
for (i in 1:length(months)){
  #new = full_join(demg %>% filter(Visitnr==months[i]),test, by = c("eid","Visitnr"))
  new = demg %>% filter(Visitnr==months[i])
  new = merge(new,test,by="eid")
  new$age = as.numeric(difftime(new$Visit_date,new$DATEOFBIRTH)/365)
  age_list[[i]] = new %>% select(eid,Visitnr,age)
  #print(new$age)
}
age_list=list.rbind(age_list)
age_list$session = as.numeric(age_list$Visitnr)
new = full_join(tmp,age_list,by=c("eid","session"))
#
# quick and dirty solution to make sure all age fields are filled
## one might want to remove NA age for baseline (unfortunately, these guys are then lost, there is no way of recovering this data)
#new = new %>% filter(!(eid %in% (new %>% filter(session==0 & is.na(age)==T) %>% select(eid))$eid & session==0))
#new %>% group_by(session) %>% summarize(NA_sum = sum(is.na(age)), Rows = length(age))
# now, we add the number of months at each visit to the baseline age (in years)
# this is not the perfect solution, but a feasible solution.
for (i in 1:length(new$age)){
  if (is.na(new$age[i])==TRUE){
    if (new %>% filter(eid == eid[i] & session == 0) %>% nrow() > 0){
      a0 = new %>% filter(eid == eid[i] & session == 0) %>% select(age)
      new$age[i] = as.numeric(a0$age)+new[i,]$session/12
      }}}
new = new %>% select(-Visitnr)
#new = na.omit(new)
new = cbind(new,cortical %>% select(-eid,-session,-age)) # merge cortical and subcortical stats
new = new %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
new$session = factor(new$session)
#
# add demographics and biol info (edss, geno, NfL, CHL, etc)
#
### edss
demo10$M144_EDSSscore = as.numeric(unlist(ifelse(demo10$EDSS_score_10 == " ",NA,gsub(",",".",demo10$EDSS_score_10))))
edss10 = demo10 %>% select(Patnr,M144_EDSSscore)
edss = merge(edss,edss10,by="Patnr")
edss = melt(edss, id.vars="Patnr")
names(edss) = c("eid","session","edss")
edss$session = factor(edss$session)
levels(edss$session) = c(1,6,12,18,24,144)
new = full_join(new,edss,by=c("eid","session"))
new = new %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
#
## blood
blood = blood %>% select(Sub,Visit.nr,NfL..pg.ml.,CH3L.1..mg.ml..mean)
names(blood) = c("eid","session","NfL","CH3L1")
blood$session =factor(blood$session)
new = full_join(new,blood,by=c("eid","session"))
new = new %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
#
## risk gene carriership
names(geno) = c("eid","geno_HLA1501_1")
new = full_join(new,geno,by="eid")
#
## Lesions
lesion_vol = melt(lesion_vol,id.vars="eid")
names(lesion_vol) = c("eid","session","lesion_vol")
levels(lesion_vol$session) = c(0,1,2,3,4,5,6,7,8,9,12,24,120,15,18)
lesion_count = melt(lesion_count,id.vars="eid")
names(lesion_count) = c("eid","session","lesion_count")
levels(lesion_count$session) = c(0,1,2,3,4,5,6,7,8,9,12,24,120,15,18)
new = full_join(new,lesion_count,by=c("eid","session"))
new = full_join(new,lesion_vol,by=c("eid","session"))
new = new %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)
#
## sex
new = full_join(new,test0,by="eid")
#
# filter withdrawn and null (AGAIN :-) )
new = new %>% filter(!eid == 403) %>% filter(!eid == 807) %>% filter(!eid == 1106)%>% filter(!eid == 1408)

# save csv
write.csv(new,"/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
