# DEAP model 2017

library(DescTools)
library(tidyverse)
library(caret)
library(broom)
library(foba)
library(bartMachine)
library(LiblineaR)
library(elasticnet)
library(gamm4)
library(MuMIn)
library(corrplot)
library(reshape2)
library(ggplot2)
library("FactoMineR")
library("readxl")
library(rjson)
library(stargazer)
library(knitr)
library(R.matlab)
library(tableone)
library('parallel')
library('doParallel')
library('yaml')
library('lme4')
library(data.table)
options(scipen=999)

###SET THESE BEFORE RUNNING#####

ivs=readLines("/home/max/Documents/linear_mixed_model_abcd/finalfullnames_adhd.txt")
dvs=c('cbcl_scr_syn_attention_r')
covs=c("race.4level","sex","high.educ.bl","household.income.bl","age",'cbcl_scr_syn_internal_r','ehi_ss_score')
spec_covs=c('tfmri_mid_all_beta_mean.motion','tfmri_nback_all_beta_mean.motion','tfmri_sst_all_beta_mean.motion')
re_covs=c("rel_family_id","mri_info_device.serial.number")
all_meds<-sprintf("medinv_plus_rxnorm_med%s_p",c(01:15))
allvars_butivsmots=c('src_subject_id','eventname','fsqc_qc',dvs,re_covs,covs)
allvars=c(allvars_butivsmots,spec_covs,"smri_vol_subcort.aseg_intracranialvolume",'tfmri_sst_beh_perform.flag',all_meds,ivs)

#####SCRIPT BEGINS#### NO FURTHER CHANGES NEEDED########

data =  readRDS( paste0("/home/max/Documents/linear_mixed_model_abcd/nda2.0.1.Rds"))
backup_data=data
data=backup_data
data = data[c(allvars)]
data <-data[ which(data$eventname=='baseline_year_1_arm_1'),]

#####make names more manageable#####
colnames(data)<-sub("tfmri_mid_all_", "mid_", colnames(data))
colnames(data)<-sub("_beta_subcort.aseg", "_aseg", colnames(data))
colnames(data)<-sub("_beta_cort.destrieux", "_dest", colnames(data))
colnames(data)<-sub("tfmri_nback_all_", "nb_", colnames(data))
colnames(data)<-sub("_beta_subcort.aseg", "_aseg", colnames(data))
colnames(data)<-sub("_beta_cort.destrieux", "_dest", colnames(data))
colnames(data)<-sub("tfmri_sst_all_", "sst_", colnames(data))
colnames(data)<-sub("_beta_subcort.aseg", "_aseg", colnames(data))
colnames(data)<-sub("_beta_cort.destrieux", "_dest", colnames(data))
col=which( colnames(data)=="nb_2.back.vs.0.back_dest_g.and.s.frontomargin.lh")
ivs=colnames(data[c(col:length(data))])

#scale icv because its huge
data$smri_vol_subcort.aseg_intracranialvolume<-scale(data$smri_vol_subcort.aseg_intracranialvolume)

#add puberty
puberty=read.csv("/home/max/Documents/linear_mixed_model_abcd/puberty.csv")
colnames(puberty)[2]='puberty'
data=merge(puberty,data,by="src_subject_id",all.y=TRUE)
covs=c(covs,'puberty')

#make function to pick out if each medication variable contains a stimulant med for each subject
myfunction<-function(x){
  grepl("adderall|concerta|methylphenidate|ritalin|focalin|strattera|amphetamine|quillivant|guanfacine|evekeo|atomoxetine|lisdexamfetamine|dexedrine|dynavel|adzenys|metadate|kapvay|clonidine|intuniv|daytrana|methylin|dextrostat|zenzedi|tenex|catapres|aptensio|cotempla|daytrana|quillichew|bupropion|wellbutrin|norpramin|desipramine|impiprmine|tofranil|nortriptyline|aventyl|pamelor",x,ignore.case = TRUE)
}

#run the function then consolidate outputs

v1<-rowSums(sapply(data[all_meds],myfunction))
stimulant_med_yesno<-v1>0
stimulant_med_yesno<-as.factor(stimulant_med_yesno)
data<-cbind(data,stimulant_med_yesno)
covs<-c(covs,'stimulant_med_yesno')

#######split into specific modalities#########
mid_beta_mean.motion<-data$mid_beta_mean.motion
nb_beta_mean.motion<-data$nb_beta_mean.motion
sst_beta_mean.motion<-data$sst_beta_mean.motion
smri_vol_subcort.aseg_intracranialvolume<-data$smri_vol_subcort.aseg_intracranialvolume
allmot<-rowMeans(cbind(mid_beta_mean.motion,nb_beta_mean.motion,sst_beta_mean.motion),na.rm=TRUE)

#split into separate files for each task
run_names<-list('smri','mid','nb','sst')
for (l in run_names)
{
  assign(l,data[c('src_subject_id','eventname','fsqc_qc',dvs,re_covs,covs)])
  n=data[c(ivs)]
  if (l=='mid'){mid<-cbind(mid,mid_beta_mean.motion,(n[, grep("mid_", names(n)) ]))}
  if (l=='nb'){nb<-cbind(nb,nb_beta_mean.motion,(n[, grep("nb_", names(n)) ]))}
  if (l=='sst'){sst<-cbind(sst,data$'tfmri_sst_beh_perform.flag',sst_beta_mean.motion,(n[, grep("sst_", names(n)) ]))}
  if (l=='smri'){smri<-cbind(smri,allmot,smri_vol_subcort.aseg_intracranialvolume,(n[, grep("smri", names(n))]))}
}
spec_covs=c("smri_vol_subcort.aseg_intracranialvolume",'mid_beta_mean.motion','nb_beta_mean.motion','sst_beta_mean.motion')

#combine separate task files into a list
runs=list(smri,mid,nb,sst)


#remove fmri/mri missing data
mid <- mid[complete.cases(mid[c(ivs[grep("mid_", ivs) ])]), ]
nb <- nb[complete.cases(nb[c(ivs[grep("nb_", ivs) ])]), ]
sst <- sst[complete.cases(sst[c(ivs[grep("sst_", ivs) ])]), ]
smri <- smri[complete.cases(smri[c(ivs[grep("smri_", ivs) ])]), ]

#########SMRI QC###########
mid <- mid[mid$fsqc_qc == 'accept', ]
nb <- nb[nb$fsqc_qc == 'accept', ]
sst <- sst[sst$fsqc_qc == 'accept', ]
smri <- smri[smri$fsqc_qc == 'accept', ]

#########FMRI QC###########
QCed_final<-read_excel("/home/max/Documents/linear_mixed_model_abcd/ABCD_QC_voxel_vertex.xlsx")
QCed_final$nBack<-gsub("INV","NDAR_INV", QCed_final$nBack)
QCed_final$MID<-gsub("INV","NDAR_INV", QCed_final$MID)
QCed_final$SST<-gsub("INV","NDAR_INV", QCed_final$SST)
names(QCed_final)[1]<-"src_subject_id"
names(QCed_final)[2]<-"src_subject_id"
names(QCed_final)[3]<-"src_subject_id"
QCed_nb<-QCed_final[1]
QCed_mid<-QCed_final[2]
QCed_sst<-QCed_final[3]
QCed_nb<-QCed_nb[complete.cases(QCed_nb), ]
QCed_mid<-QCed_mid[complete.cases(QCed_mid), ]
QCed_sst<-QCed_sst[complete.cases(QCed_sst), ]

mid <- mid[(mid$src_subject_id %in% QCed_mid$src_subject_id),]
nb <- nb[(nb$src_subject_id %in% QCed_nb$src_subject_id),]
sst <- sst[(sst$src_subject_id %in% QCed_sst$src_subject_id),]

rm(QCed_final,QCed_mid,QCed_nb,QCed_sst)

#remove phillips
philips_exclude<-read.delim("/home/max/Documents/linear_mixed_model_abcd/philips.txt",header = FALSE)
names(philips_exclude)[1]<-"src_subject_id"
mid <- mid[!(mid$src_subject_id %in% philips_exclude$src_subject_id),]
nb <- nb[!(nb$src_subject_id %in% philips_exclude$src_subject_id),]
sst <- sst[!(sst$src_subject_id %in% philips_exclude$src_subject_id),]

mid <- mid[complete.cases(mid[c(dvs)]), ]
nb <- nb[complete.cases(nb[c(dvs)]), ]
sst <- sst[complete.cases(sst[c(dvs)]), ]
smri <- smri[complete.cases(smri[c(dvs)]), ]

mid <- mid[complete.cases(mid[c(covs)]), ]
nb <- nb[complete.cases(nb[c(covs)]), ]
sst <- sst[complete.cases(sst[c(covs)]), ]
smri <- smri[complete.cases(smri[c(covs,'allmot','smri_vol_subcort.aseg_intracranialvolume')]), ]

#remove bad SST subjects
sst_exclude<-read.csv("/home/max/ADHD/Overlap_Flagged_3.0_Release.csv",header = TRUE)
names(sst_exclude)[1]<-"src_subject_id"
sst_exclude <- sst_exclude[sst_exclude$eventname == 'baseline_year_1_arm_1', ]
sst_exclude2<-read.csv("/home/max/ADHD/Violators_3.0_Release.csv",header = TRUE)
names(sst_exclude2)[1]<-"src_subject_id"
sst_exclude2 <- sst_exclude2[sst_exclude2$eventname == 'baseline_year_1_arm_1', ]
sst <- sst[!(sst$src_subject_id %in% sst_exclude$src_subject_id),]
sst <- sst[!(sst$src_subject_id %in% sst_exclude2$src_subject_id),]
#names(sst)[17]<- 'tfmri_sst_beh_perform.flag'
names(sst)[16]<- 'tfmri_sst_beh_perform.flag'
sst <- sst[sst$tfmri_sst_beh_perform.flag == 'acceptable', ]
ssrt<-read.csv("/home/max/ADHD/SSRT_Release_3.0.csv",header = TRUE)
ssrt <- ssrt[ssrt$eventname == 'baseline_year_1_arm_1', ]
names(ssrt)[1]<-"src_subject_id"
sst$tfmri_sst_beh_perform.flag<-NULL

smri <- smri[complete.cases(smri), ]
mid <- mid[complete.cases(mid), ]
nb <- nb[complete.cases(nb), ]
sst <- sst[complete.cases(sst), ]

######winsorize mri###########
smri_win=smri
mid_win=mid
sst_win=sst
nb_win=nb

start_col=which(colnames(nb)=="nb_2.back.vs.0.back_dest_g.and.s.frontomargin.lh")
smri_win[,start_col:length(smri)]<-as.data.frame(sapply(smri[,start_col:length(smri)],Winsorize,probs = c(0.05, 0.95),type=7))
mid_win[,start_col:length(mid)]<-as.data.frame(sapply(mid[,start_col:length(mid)],Winsorize,probs = c(0.05, 0.95),type=7))
sst_win[,start_col:length(sst)]<-as.data.frame(sapply(sst[,start_col:length(sst)],Winsorize,probs = c(0.05, 0.95),type=7))
nb_win[,start_col:length(nb)]<-as.data.frame(sapply(nb[,start_col:length(nb)],Winsorize,probs = c(0.05, 0.95),type=7))

#user_data = data
runs=list(smri,mid,nb,sst)
runs_win=list(smri_win,mid_win,nb_win,sst_win)
covs_ls=list()
ivs_ls=list()

#####demographics######
range(data$age)
table(data$high.educ.bl)
table(data$household.income.bl)
range(data$cbcl_scr_syn_attention_r)
table(data$race.4level)
table(data$sex)
table(data$stimulant_med_yesno)

for (r in 1:4) 
{
user_data=runs_win[[r]]
#ivs_ls[[r]]=readLines( paste0("/home/max/Documents/linear_mixed_model_abcd/ivs_",run_names[[r]],".txt"))
col=which( colnames(user_data)=="stimulant_med_yesno")
col2=col+2
col3=col+3
if (r==1){ivs_ls[[r]]=names(user_data[col3:length(user_data)])}
if (r>1){ivs_ls[[r]]=names(user_data[col2:length(user_data)])}
if (r==1){covs_ls[[r]]=c(covs,'smri_vol_subcort.aseg_intracranialvolume','allmot')}
if (r>1){covs_ls[[r]]=c(covs,spec_covs[c(r)])}
}

#####build mixed model function for abcd####
mixed_model_hcp=function(x,y,covs,data){
stat_holder <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_holder[k] <- as.double()

form_cov_only <- formula(paste(y, "~", paste(covs, collapse="+")))
form <- formula(paste(y, "~", x, "+", paste(covs, collapse="+")))
model <- gamm4(form, data=data, random =~(1|mri_info_device.serial.number/rel_family_id) )
model2 <- gamm4(form_cov_only, data=data, random =~(1|mri_info_device.serial.number/rel_family_id))
r2_delta = round(as.numeric(r.squaredLR(model$mer,model2$mer)),5)
sg<-summary(model$gam)

for (statnum in 1:4){
stat_holder[1,statnum]<-sg$p.table[2,statnum]
}
stat_holder[1,5]<-r2_delta
return(stat_holder)
}      

mixed_model_hcp_nocov=function(x,y,data){
  stat_holder <- data.frame()
  stat_names <- c('B','SE','t','p','R2')
  for (k in stat_names) stat_holder[k] <- as.double()
  
  nullmod <- formula(paste(y, "~", 1 ))
  form_cov_only <- formula(paste(y, "~", x ))
  model <- gamm4(form_cov_only, data=data, random =~(1|mri_info_device.serial.number/rel_family_id) )
  nlmodel <- gamm4(nullmod, data=data, random =~(1|mri_info_device.serial.number/rel_family_id) )
  r2_delta = round(as.numeric(r.squaredLR(model$mer,nlmodel$mer)),5)
  sg<-summary(model$gam)
  
  for (statnum in 1:4){
    stat_holder[1,statnum]<-sg$p.table[2,statnum]
  }
  stat_holder[1,5]<-r2_delta
  return(stat_holder)
}    

###################make ordered factor############

#get only regions sig in en for shortening me analyses
ivs1=readLines("/home/max/Documents/linear_mixed_model_abcd/nb_nocov.txt")
ivs2=readLines("/home/max/Documents/linear_mixed_model_abcd/smri_nocov.txt")
ivs3=readLines("/home/max/Documents/linear_mixed_model_abcd/nb_cov.txt")
ivs4=readLines("/home/max/Documents/linear_mixed_model_abcd/nb_medcov.txt")
ivs5=readLines("/home/max/Documents/linear_mixed_model_abcd/nb_factor.txt")
ivs6=readLines("/home/max/Documents/linear_mixed_model_abcd/smri_factor.txt")

colnames(runs_win[[3]])[c(17:length(runs_win[[3]]))]<-gsub(".", "_", colnames(runs_win[[3]][c(17:length(runs_win[[3]]))]), fixed = TRUE)
colnames(runs_win[[1]])[c(17:length(runs_win[[1]]))]<-gsub(".", "_", colnames(runs_win[[1]][c(17:length(runs_win[[1]]))]), fixed = TRUE)
  

#########run mixed model##########
y <- dvs #'ksads_14_853_p'
x <- ivs3[2]
covs <- covs_ls[[3]]

stat_list_nb<-(lapply(x, mixed_model_hcp, covs=covs,y=y, data=runs_win[[3]]))

stat_list_nb_nocov<-(lapply(ivs1[1:46], mixed_model_hcp_nocov,y=y, data=runs_win[[3]]))
stat_list_smri<-(lapply(ivs2, mixed_model_hcp_nocov,y=y, data=runs_win[[1]]))
stat_list_nb_standcov<-(lapply(ivs3, mixed_model_hcp,y=y, covs=covs_ls[[3]][-c(9)], data=runs_win[[3]]))
stat_list_nb_medcov<-(lapply(ivs4, mixed_model_hcp,y=y, covs=covs_ls[[3]], data=runs_win[[3]]))

stat_list_nb_factor<-(lapply(ivs5, mixed_model_hcp_nocov,y=y, data=runs_win[[3]]))
stat_list_smri_factor<-(lapply(ivs6, mixed_model_hcp_nocov,y=y, data=runs_win[[1]]))

stat_matrix_nb <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_nb[k] <- as.double()
for (i in 1:length(stat_list_nb_nocov)){stat_matrix_nb[i,]=stat_list_nb_nocov[[i]]}
for (n in 1:length(ivs1)) row.names(stat_matrix_nb)[n] <- ivs1[n]

stat_matrix_smri <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_smri[k] <- as.double()
for (i in 1:length(stat_list_smri)){stat_matrix_smri[i,]=stat_list_smri[[i]]}
for (n in 1:length(ivs2)) row.names(stat_matrix_smri)[n] <- ivs2[n]

stat_matrix_nb_standcov <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_nb_standcov[k] <- as.double()
for (i in 1:length(stat_list_nb_standcov)){stat_matrix_nb_standcov[i,]=stat_list_nb_standcov[[i]]}
for (n in 1:length(ivs3)) row.names(stat_matrix_nb_standcov)[n] <- ivs3[n]

stat_matrix_nb_medcov <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_nb_medcov[k] <- as.double()
for (i in 1:length(stat_list_nb_medcov)){stat_matrix_nb_medcov[i,]=stat_list_nb_medcov[[i]]}
for (n in 1:length(ivs4)) row.names(stat_matrix_nb_medcov)[n] <- ivs4[n]

write.csv(stat_matrix_nb,"/home/max/Documents/ABCD_ADHD/stat_matrix_nb.csv")
write.csv(stat_matrix_nb_standcov,"/home/max/Documents/ABCD_ADHD/stat_matrix_nb_standcov.csv")
write.csv(stat_matrix_nb_medcov,"/home/max/Documents/ABCD_ADHD/stat_matrix_nb_medcov.csv")
write.csv(stat_matrix_smri,"/home/max/Documents/ABCD_ADHD/stat_matrix_smri.csv")

stat_matrix_nb_factor <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_nb_factor[k] <- as.double()
for (i in 1:length(stat_list_nb_factor)){stat_matrix_nb_factor[i,]=stat_list_nb_factor[[i]]}
for (n in 1:length(ivs5)) row.names(stat_matrix_nb_factor)[n] <- ivs5[n]

stat_matrix_smri_factor <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_smri_factor[k] <- as.double()
for (i in 1:length(stat_list_smri_factor)){stat_matrix_smri_factor[i,]=stat_list_smri_factor[[i]]}
for (n in 1:length(ivs6)) row.names(stat_matrix_smri_factor)[n] <- ivs6[n]

write.csv(stat_matrix_nb_factor,"/home/max/Documents/ABCD_ADHD/stat_matrix_nb_factor.csv")
write.csv(stat_matrix_smri_factor,"/home/max/Documents/ABCD_ADHD/stat_matrix_smri_factor.csv")

#create dset for matlab
for (r in 1:4){
  names(runs_win[[r]])[7]<-'mri_info_device_serial_number'
  names(runs_win[[r]])[1]<-'id_redcap'
  runs_win[[r]]=runs_win[[r]][c('id_redcap','cbcl_scr_syn_attention_r',ivs_ls[[r]])]
}

write.csv(runs_win[[1]],'/home/max/Documents/ABCD_ADHD/smri_win_RR.csv',row.names = FALSE)
write.csv(runs_win[[2]],'/home/max/Documents/ABCD_ADHD/mid_win_RR.csv',row.names = FALSE)
write.csv(runs_win[[3]],'/home/max/Documents/ABCD_ADHD/nb_win_RR.csv',row.names = FALSE)
write.csv(runs_win[[4]],'/home/max/Documents/ABCD_ADHD/sst_win_RR.csv',row.names = FALSE)

test=read.csv('/home/max/Documents/ABCD_ADHD/TESTofChanges_fixed_nb_noresid.csv')

#residualize dset thing for matlab, without meds
for (r in 1:4){
  form_cov_only <- formula(paste('cbcl_scr_syn_attention_r', "~",paste(covs_ls[[r]][c(1:8,10:length(covs_ls[[r]]))], collapse="+")))
  model=gamm4(form_cov_only, data=runs_win[[r]], random =~(1|mri_info_device.serial.number/rel_family_id) )
  runs_win[[r]][c('residualized_cbcl_scr_syn_attention_r')] <- as.numeric(model$gam$residuals)
  names(runs_win[[r]])[7]<-'mri_info_device_serial_number'
  names(runs_win[[r]])[1]<-'id_redcap'
  runs_win[[r]]=runs_win[[r]][c('id_redcap','residualized_cbcl_scr_syn_attention_r',ivs_ls[[r]])]
}

write.csv(runs_win[[1]],'/home/max/Documents/ABCD_ADHD/smri_win_RR_standardcovs.csv',row.names = FALSE)
write.csv(runs_win[[2]],'/home/max/Documents/ABCD_ADHD/mid_win_RR_standardcovs.csv',row.names = FALSE)
write.csv(runs_win[[3]],'/home/max/Documents/ABCD_ADHD/nb_win_RR_standardcovs.csv',row.names = FALSE)
write.csv(runs_win[[4]],'/home/max/Documents/ABCD_ADHD/sst_win_RR_standardcovs.csv',row.names = FALSE)

#residualize dset thing for matlab, with meds
for (r in 1:4){
  form_cov_only <- formula(paste('cbcl_scr_syn_attention_r', "~",paste(covs_ls[[r]], collapse="+")))
  model=gamm4(form_cov_only, data=runs_win[[r]], random =~(1|mri_info_device.serial.number/rel_family_id) )
  runs_win[[r]][c('residualized_cbcl_scr_syn_attention_r')] <- as.numeric(model$gam$residuals)
  names(runs_win[[r]])[7]<-'mri_info_device_serial_number'
  names(runs_win[[r]])[1]<-'id_redcap'
  runs_win[[r]]=runs_win[[r]][c('id_redcap','residualized_cbcl_scr_syn_attention_r',ivs_ls[[r]])]
}

write.csv(runs_win[[1]],'/home/max/Documents/ABCD_ADHD/smri_win_RR_medcovary.csv',row.names = FALSE)
write.csv(runs_win[[2]],'/home/max/Documents/ABCD_ADHD/mid_win_RR_medcovary.csv',row.names = FALSE)
write.csv(runs_win[[3]],'/home/max/Documents/ABCD_ADHD/nb_win_RR_medcovary.csv',row.names = FALSE)
write.csv(runs_win[[4]],'/home/max/Documents/ABCD_ADHD/sst_win_RR_medcovary.csv',row.names = FALSE)

library(gtools)
for (r in 1:4){
runs_win[[r]]$cbcl_scr_syn_attention_r<-quantcut(runs_win[[r]]$cbcl_scr_syn_attention_r, q=3)
runs_win[[r]]$cbcl_scr_syn_attention_r<-as.numeric(runs_win[[r]]$cbcl_scr_syn_attention_r)
}

#create dset for matlab with factor scores
for (r in 1:4){
  names(runs_win[[r]])[7]<-'mri_info_device_serial_number'
  names(runs_win[[r]])[1]<-'id_redcap'
  runs_win[[r]]=runs_win[[r]][c('id_redcap','cbcl_scr_syn_attention_r',ivs_ls[[r]])]
}

write.csv(runs_win[[1]],'/home/max/Documents/ABCD_ADHD/smri_win_RR_factor.csv',row.names = FALSE)
write.csv(runs_win[[2]],'/home/max/Documents/ABCD_ADHD/mid_win_RR_factor.csv',row.names = FALSE)
write.csv(runs_win[[3]],'/home/max/Documents/ABCD_ADHD/nb_win_RR_factor.csv',row.names = FALSE)
write.csv(runs_win[[4]],'/home/max/Documents/ABCD_ADHD/sst_win_RR_factor.csv',row.names = FALSE)

#create dset for matlab ksads
for (r in 1:4){
  names(runs_win[[r]])[7]<-'mri_info_device_serial_number'
  names(runs_win[[r]])[1]<-'id_redcap'
  runs_win[[r]]=runs_win[[r]][c('id_redcap','ksads_14_853_p',ivs_ls[[r]])]
}

write.csv(runs_win[[1]],'/home/max/Documents/ABCD_ADHD/smri_win_RR_ksads.csv',row.names = FALSE)
write.csv(runs_win[[2]],'/home/max/Documents/ABCD_ADHD/mid_win_RR_ksads.csv',row.names = FALSE)
write.csv(runs_win[[3]],'/home/max/Documents/ABCD_ADHD/nb_win_RR_ksads.csv',row.names = FALSE)
write.csv(runs_win[[4]],'/home/max/Documents/ABCD_ADHD/sst_win_RR_ksads.csv',row.names = FALSE)


#----------------Enet ANALYSES------------#

#############Set up Train/Test Split#########

#form_cov_only <- formula(paste('cbcl_scr_syn_attention_r', "~","+",paste(covs_ls[[r]], collapse="+")))
#model=gamm4(form_cov_only, data=runs[[r]], random =~(1|mri_info_device.serial.number/rel_family_id) )
#runs[[r]][c('residualized_cbcl_scr_syn_attention_r')] <- as.numeric(residuals(model))

### split into test and training sets (80/20)
trainIndex=list()
dataTrain=list()
dataTest=list()
for (r in 2:length(runs)){

set.seed(123)
trainIndex[[r]] <- createDataPartition(runs[[r]]$cbcl_scr_syn_attention_r, p = .8, 
                                  list = FALSE, 
                                  times = 1)

dataTrain[[r]] <- runs[[r]][ trainIndex[[r]],]
dataTest[[r]]  <- runs[[r]][-trainIndex[[r]],]

#residualize covariates
form_cov_only <- formula(paste('cbcl_scr_syn_attention_r', "~",paste(covs_ls[[r]], collapse="+")))
model1=gamm4(form_cov_only, data=dataTrain[[r]], random =~(1|mri_info_device.serial.number/rel_family_id) )
dataTrain[[r]][c('residualized_cbcl_scr_syn_attention_r')] <- as.numeric(model1$gam$residuals)
model2=gamm4(form_cov_only, data=dataTest[[r]], random =~(1|mri_info_device.serial.number/rel_family_id) )
dataTest[[r]][c('residualized_cbcl_scr_syn_attention_r')] <- as.numeric(model2$gam$residuals)
}

#########data processing##################

#winsorize for elastic net
for (r in 2:length(runs)){
  dataTrain[[r]][,ivs_ls[[r]]]<-as.data.frame(sapply(dataTrain[[r]][,ivs_ls[[r]]],Winsorize,probs = c(0.05, 0.95),type=7))
  dataTest[[r]][,ivs_ls[[r]]]<-as.data.frame(sapply(dataTest[[r]][,ivs_ls[[r]]],Winsorize,probs = c(0.05, 0.95),type=7))
}

trainTransformed=list()
testTransformed=list()

for (r in 2:length(runs)){
preProcValues <- preProcess(dataTrain[[r]][c(ivs_ls[[r]],dvs)], method = c("center", "scale"))
trainTransformed[[r]] <- predict(preProcValues, dataTrain[[r]])
testTransformed[[r]] <- predict(preProcValues, dataTest[[r]])
}

#######set up internal 5 fold cv#########
for (r in 2:length(runs)){
trainTransformed[[r]]$cv=createFolds(trainTransformed[[r]]$cbcl_scr_syn_attention_r, k = 5, list = FALSE, returnTrain = FALSE)
}

dataTrain_in=list(list(),list(),list(),list(),list())
dataVal=list(list(),list(),list(),list(),list())

for (r in 2:length(runs)){
for (i in 1:5){
  dataTrain_in[[r]][[i]] <- trainTransformed[[r]][ trainTransformed[[r]]$cv!=i,]
  dataVal[[r]][[i]]  <- trainTransformed[[r]][ trainTransformed[[r]]$cv==i,]
}
}

##########ML using caret############
enFit_nb=list()
enFit_all=list(list(),list(),list(),list(),list())
pred=list(list(),list(),list(),list(),list())
val_stats=list(list(),list(),list(),list(),list())
fitControl <- trainControl(
  method = "repeatedcv", #cv
  number = 20,
  repeats = 1)

#grid <- expand.grid(size=c(5,10,20,50), k=c(1,2,3,4,5))
myGrid <- expand.grid(
  fraction = seq(0.05, 1, length = 20),
  lambda = seq(0.0001, 1, length = 20)
)

set.seed(111)

for (i in 1:5){
  enFit_nb[[i]] <- train(cbcl_scr_syn_attention_r~ ., 
                      data = subset(dataTrain_in[[3]][[i]],select = c(ivs_ls[[3]],dvs)),
                      method = "enet",
                      trControl = fitControl,
                      tuneGrid = myGrid)
}

max_r2=function(r2_vec)
  {
  #gets maximum r-square from a caret internal cv
  #eg call: nb_maxr2=max_r2(enFit_nb[[1]]$results$Rsquared)
  mr2=max(r2_vec)
  if (exists('mr2_best')==FALSE){mr2_best=0}
  if (mr2>mr2_best){mr2_best=mr2}
  return(mr2_best)
}

nb_maxr2=list()
for (i in 1:5){nb_maxr2[i]=max_r2(enFit_nb[[i]]$results$Rsquared)}

#sapply(enFit_nb[[1:5]]$results$Rsquared , max_r2)

pred <- list()
for (i in 1:5){
  pred[[i]] <- predict(enFit_nb[[i]], dataVal[[3]][[i]])
  val_stats[[i]]=postResample(pred = pred[[i]], obs = dataVal[[3]][[i]]$cbcl_scr_syn_attention_r)
}

#calculate/organize validation scores
val_scores<-data.frame(postResample(pred = pred[[1]], obs = dataVal[[3]][[1]]$cbcl_scr_syn_attention_r))
for (i in 2:5) {
  val_scores[i]<-postResample(pred = pred[[i]], obs = dataVal[[3]][[i]]$cbcl_scr_syn_attention_r)
}
val_scores=transpose(val_scores)
names(val_scores)=c("RMSE","Rsquared","MAE")
row.names(val_scores)=c("Fold 1","Fold 2","Fold 3","Fold 4","Fold 5")

#what does this do and why was it here?
#obs<-dataVal[[3]][[4]]$cbcl_scr_syn_attention_r

par(mar=c(1,1,1,1))
plot(enFit_nb[[4]], main="Model Accuracies")
varimp_mars <- varImp(enFit_nb[[4]])
plot(varimp_mars, main="Variable Importance", top = 20)


####################DEMOGRPAHICS######################
library(ggplot2)
#standard cbcl
barfill <- "#4271AE"
barlines <- "#1F3552"
ggplot(data=runs_win[[1]], aes(runs_win[[1]]$cbcl_scr_syn_attention_r)) + 
  geom_histogram(colour = barlines, fill = barfill,binwidth=1) +
  scale_x_continuous(name = "ADHD Symptomatology")+
  scale_y_continuous(name = "Participants") +
  theme(axis.text.x = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45),
        axis.text.y = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45))+
  theme(axis.title.y = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 90))+
  theme(axis.title.x = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 00))

#factors
barfill <- "#4271AE"
barlines <- "#1F3552"
ggplot(data=runs_win[[1]], aes(runs_win[[1]]$cbcl_scr_syn_attention_r)) + 
  geom_histogram(colour = barlines, fill = barfill,binwidth=1) +
  scale_x_continuous(name = "ADHD Symptomatology", breaks=1:3,
                     labels=c("Low", "Medium","High"))+
  scale_y_continuous(name = "Participants") +
  theme(axis.text.x = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45),
        axis.text.y = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45))+
  theme(axis.title.y = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 90))+
  theme(axis.title.x = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 00))

#ksads bar graph
barfill <- "#4271AE"
barlines <- "#1F3552"
ggplot(data=runs_win[[1]], aes(runs_win[[1]]$ksads_14_853_p)) + 
  geom_bar(colour = barlines, fill = barfill) +
  scale_x_continuous(name = "ADHD Symptomatology", breaks=0:1,
                     labels=c("Control", "ADHD"))+
  scale_y_continuous(name = "Participants") +
  theme(axis.text.x = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45),
        axis.text.y = element_text(family='URWBookman',face="bold", color=barlines, 
                                   size=rel(1.5), angle=45))+
  theme(axis.title.y = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 90))+
  theme(axis.title.x = element_text(family='URWBookman',face="bold",
                                    size = rel(1.3), color="black", angle = 00))

###########demographics###############
for (r in 1:4) runs[[r]]=merge(runs[[r]],backup_data[c('src_subject_id','hisp')])
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
demographics<-data.frame(Age=rep(NA,8),Sex1=rep(NA,8),Hisp=rep(NA,8),Race1=rep(NA,8),Rac2=rep(NA,8),Race3=rep(NA,8),Race4=rep(NA,8),Income1=rep(NA,8),Income2=rep(NA,8),Income3=rep(NA,8),Education1=rep(NA,8),Education2=rep(NA,8),Education3=rep(NA,8),Education4=rep(NA,8),Education5=rep(NA,8))

for (r in 1:4) 
{
  data=runs[[r]]
  demographics[r,1]<-mean(data$age)
  demographics[(r+4),1]<-sd(data$age)
  demographics[r,2]<-table(data['sex'])[1]
  demographics[r,3]<-table(data['hisp'])[2]
  for (i in 1:4) {demographics[r,(3+i)]<-table(data['race.4level'])[i]}
  for (i in 1:3) {demographics[r,(7+i)]<-table(data['household.income.bl'])[i]}
  for (i in 1:5) {demographics[r,(10+i)]<-table(data['high.educ.bl'])[i]}
  for (i in 1:2) {demographics[(r+4),2]<-table(data['sex'])[1]/length(data$sex)}
  for (i in 1:2) {demographics[(r+4),3]<-table(data['hisp'])[2]/length(data$hisp)}
  for (i in 1:4) {demographics[(r+4),(3+i)]<-table(data['race.4level'])[i]/length(data$race.4level)}
  for (i in 1:3) {demographics[(r+4),(7+i)]<-table(data['household.income.bl'])[i]/length(data$household.income.bl)}
  for (i in 1:5) {demographics[(r+4),(10+i)]<-table(data['high.educ.bl'])[i]/length(data$high.educ.bl)}
}

#for (coll in 2:14){
#for (i in 5:8) {demographics[i,coll]<-percent(demographics[i,coll])}}
is.num <- sapply(demographics, is.numeric)
demographics[is.num] <- lapply(demographics[is.num], round, 2)

library(data.table)
demographics_t <- transpose(demographics)
row.names(demographics_t)<-names(demographics)
colnames(demographics_t)<-c('smri','mid','nb','sst','%smri','%mid','%nb','%sst')

write.csv(demographics_t,"/home/max/Dropbox/ABCD/ADHD/TP/RR/demographics_RR.csv")

################cov only model for paper#############
mixed_model_hcp_nocov=function(x,y,data){
  stat_holder <- data.frame()
  stat_names <- c('B','SE','t','p','R2')
  for (k in stat_names) stat_holder[k] <- as.double()
  
  nullmod <- formula(paste(y, "~", 1 ))
  form_cov_only <- formula(paste(y, "~", x ))
  model <- gamm4(form_cov_only, data=data, random =~(1|mri_info_device.serial.number/rel_family_id) )
  nlmodel <- gamm4(nullmod, data=data, random =~(1|mri_info_device.serial.number/rel_family_id) )
  r2_delta = round(as.numeric(r.squaredLR(model$mer,nlmodel$mer)),5)
  sg<-summary(model$gam)
  
  for (statnum in 1:4){
    for (rw in 2:(length(sg$p.table)/4)){
    stat_holder[(rw-1),statnum]<-sg$p.table[(rw),statnum]
    }
  }
  for (rw in 2:(length(sg$p.table)/4)){
  stat_holder[(rw-1),5]<-r2_delta
  }
  return(list(stat_holder,sg))
}  


y <- dvs #'ksads_14_853_p'
x <- covs_ls[[r]][c(10)]
#covs <- #covs_ls[[3]][2:10]
  
slco <-(lapply(x, mixed_model_hcp_nocov, y, data=runs_win[[r]]))

stat_matrix_covs <- data.frame()
stat_names <- c('B','SE','t','p','R2')
for (k in stat_names) stat_matrix_covs[k] <- as.double()
c=1
for (i in 1:length(x)){
  for (rw in 1:length(slco[[i]][[1]][,1])){
  stat_matrix_covs[c,]=slco[[i]][[1]][rw,]
  c=c+1
  }
}
c=1
for (i in 1:length(x)){
  for (rw in 1:length(slco[[i]][[1]][,1])){
  row.names(stat_matrix_covs)[c] <- rownames(slco[[i]][[2]]$p.table)[[(rw+1)]]
  c=c+1
  }
}

write.csv(stat_matrix_covs,"/home/max/Dropbox/ABCD/ADHD/TP/RR/cov_lme_RR_sst.csv")

###################BEHAVIORAL CORRELATES ANALYSIS#################
######set up behavior analyses

######set up behavior analyses##################
behav_data<-read.csv("/home/max/ADHD/full_sampel_mri_qc_for_bader_03162020.csv")
nb_win<-merge(nb,behav_data[c('src_subject_id','tfmri_sst_all_beh_total_meanrt','dprime_2back')],by='src_subject_id')
sst_win<-merge(sst,behav_data[c('src_subject_id','tfmri_sst_all_beh_total_meanrt','dprime_2back')],by='src_subject_id')

beh<-c('dprime_2back','tfmri_sst_all_beh_total_meanrt')

#########run behav mixed model##########
y <- dvs #'ksads_14_853_p'
x1 <- beh[1]
x2 <- beh[2]

#stat_list_nb<-(lapply(x, mixed_model_hcp, covs=covs,y=y, data=runs_win[[3]]))

mixed_model_hcp_nocov(x1,y,nb_win)
mixed_model_hcp_nocov(x2,y,sst_win)

#########test random effects################
y <- dvs #'ksads_14_853_p'
x <- 'nb_2.back.vs.0.back_dest_s.intrapariet.and.p.trans.lh'
covs <- covs_ls[[3]]

form <- formula(paste(y, "~", x, "+", paste(covs, collapse="+")))
model <- gamm4(form, data=runs_win[[3]], random =~(1|mri_info_device.serial.number/rel_family_id) )
model2 <- gamm4(form, data=runs_win[[3]], random =~(1|mri_info_device.serial.number) )
r2_delta = round(as.numeric(r.squaredLR(model$mer,model2$mer)),5)

summary(model$mer)
summary(model2$mer)
