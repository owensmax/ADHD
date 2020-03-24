# DEAP model 2017

###SET THESE BEFORE RUNNING#####

ivs=readLines("/home/max/Documents/linear_mixed_model_abcd/finalfullnames_adhd.txt")
dvs= 'ksads_14_853_p' #'cbcl_scr_syn_attention_r'
covs=c("race.4level","sex","high.educ.bl","household.income.bl","age",'cbcl_scr_syn_internal_r','ehi_ss_score')
spec_covs=c("smri_vol_subcort.aseg_intracranialvolume",'tfmri_mid_all_beta_mean.motion','tfmri_nback_all_beta_mean.motion','tfmri_sst_all_beta_mean.motion')
re_covs=c("rel_family_id","mri_info_device.serial.number")
allvars_butivsmots=c('src_subject_id','eventname','fsqc_qc','iqc_dmri_good_ser',dvs,re_covs,covs)
allvars=c(allvars_butivsmots,spec_covs,ivs)

phil_exclude=FALSE
pubs=TRUE
dti=FALSE
output_name="/home/max/Documents/linear_mixed_model_abcd/LME_output_diagnosis.txt"

#####SCRIPT BEGINS#### NO FURTHER CHANGES NEEDED########

library(gamm4)
library(rjson)
library(stargazer)
library(knitr)
library(MuMIn)
library(R.matlab)
library(tableone)
library('parallel')
library('doParallel')
library('yaml')
library('lme4')
library('rjson')
library('parallel')

#data =  readRDS( paste0("/home/max/Documents/linear_mixed_model_abcd/nda2.0.1.Rds"))
#backup_data=data
data=backup_data

data = data[c(allvars)]
data <- data[complete.cases(data[c(covs)]), ]
data <- data[complete.cases(data[c(dvs)]), ]
data <-data[ which(data$eventname=='baseline_year_1_arm_1'),]
data <- data[data$fsqc_qc == 'accept', ]

if(dti==TRUE)
{
  data <- data[data$iqc_dmri_good_ser > 0, ]
}

if (pubs==TRUE)
{
puberty=read.csv("/home/max/Documents/linear_mixed_model_abcd/puberty.csv")
colnames(puberty)[2]='puberty'
data=merge(puberty,data,by="src_subject_id")
covs=c(covs,'puberty')
}

#data <- data[complete.cases(data), ]

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

#####data checks######
range(data$age)
table(data$high.educ.bl)
table(data$household.income.bl)
range(data$cbcl_scr_syn_attention_r)
table(data$race.4level)
table(data$sex)
table(data$puberty)

#######split into specific modalities#########
#runs=list(smri,mid,nb,sst)
run_names<-list('smri','mid','nb','sst')

mid_beta_mean.motion<-data$mid_beta_mean.motion
nb_beta_mean.motion<-data$nb_beta_mean.motion
sst_beta_mean.motion<-data$sst_beta_mean.motion
smri_vol_subcort.aseg_intracranialvolume<-data$smri_vol_subcort.aseg_intracranialvolume
allmot<-rowMeans(cbind(mid_beta_mean.motion,nb_beta_mean.motion,sst_beta_mean.motion),na.rm=TRUE)
allmot<-allmot[complete.cases(allmot)]

for (l in run_names)
{
  assign(l,data[c('src_subject_id','eventname','fsqc_qc','iqc_dmri_good_ser',dvs,re_covs,covs)])
  n=data[c(ivs)]
  if (l=='mid'){mid<-cbind(mid,mid_beta_mean.motion,(n[, grep("mid_", names(n)) ]))}
  if (l=='nb'){nb<-cbind(nb,nb_beta_mean.motion,(n[, grep("nb_", names(n)) ]))}
  if (l=='sst'){sst<-cbind(sst,sst_beta_mean.motion,(n[, grep("sst_", names(n)) ]))}
  if (l=='smri'){smri<-cbind(smri,smri_vol_subcort.aseg_intracranialvolume,(n[, grep("smri", names(n))]))}
}
spec_covs=c("smri_vol_subcort.aseg_intracranialvolume",'mid_beta_mean.motion','nb_beta_mean.motion','sst_beta_mean.motion')

#remove phillips
philips_exclude<-read.delim("/home/max/Documents/linear_mixed_model_abcd/philips.txt",header = FALSE)
names(philips_exclude)[1]<-"src_subject_id"
mid <- mid[!(mid$src_subject_id %in% philips_exclude$src_subject_id),]
nb <- nb[!(nb$src_subject_id %in% philips_exclude$src_subject_id),]
sst <- sst[!(sst$src_subject_id %in% philips_exclude$src_subject_id),]
  
#remove fmri/mri missing data
mid <- mid[complete.cases(mid), ]
nb <- nb[complete.cases(nb), ]
sst <- sst[complete.cases(sst), ]
smri <- smri[complete.cases(smri), ]

#########qc removals###########
library("readxl")
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

######winsorize mri###########

#winsorize
library('DescTools')
col=which( colnames(data)=="nb_2.back.vs.0.back_dest_g.and.s.frontomargin.lh")
for (i in col:length(smri))
{
  smri[c(i)]<-Winsorize(smri[,i], probs = c(0.05, 0.95),type=7)
}

for (i in col:length(mid))
{
  mid[c(i)]<-Winsorize(mid[,i], probs = c(0.05, 0.95),type=7)
}

for (i in col:length(sst))
{
  sst[c(i)]<-Winsorize(sst[,i], probs = c(0.05, 0.95),type=7)
}

for (i in col:length(nb))
{
  nb[c(i)]<-Winsorize(nb[,i], probs = c(0.05, 0.95),type=7)
}

#user_data = data
runs=list(smri,mid,nb,sst)
#counter=2

#####added stuff for automation######
percent <- function(x, digits = 0, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}
demographics<-data.frame(Age=rep(NA,9),Sex=rep(NA,9),Race=rep(NA,9),Income=rep(NA,9),Education=rep(NA,9),Sex_percent=rep(NA,9),race_percent=rep(NA,9),Income_percent=rep(NA,9),Education_percent=rep(NA,9))
demographics[1,1]<-mean(data$age)
demographics[2,1]<-sd(data$age)
for (i in 1:2) {demographics[i,2]<-table(data['sex'])[i]}
for (i in 1:4) {demographics[i,3]<-table(data['race.4level'])[i]}
for (i in 1:3) {demographics[i,4]<-table(data['household.income.bl'])[i]}
for (i in 1:5) {demographics[i,5]<-table(data['high.educ.bl'])[i]}
for (i in 1:2) {demographics[i,6]<-table(data['sex'])[i]/length(data$sex)}
for (i in 1:4) {demographics[i,7]<-table(data['race.4level'])[i]/length(data$race.4level)}
for (i in 1:3) {demographics[i,8]<-table(data['household.income.bl'])[i]/length(data$household.income.bl)}
for (i in 1:5) {demographics[i,9]<-table(data['high.educ.bl'])[i]/length(data$high.educ.bl)}

for (i in 6:9) {demographics[,i]<-percent(demographics[,i])}
is.num <- sapply(demographics, is.numeric)
demographics[is.num] <- lapply(demographics[is.num], round, 0)

library(data.table)
demographics_t <- transpose(demographics)
row.names(demographics_t)<-names(demographics)

write.csv(demographics_t,"/home/max/Dropbox/ABCD/ADHD/demographics.csv")


sg<-summary(model$gam)
for (statnum in 1:4){
  stat_list[[splitnum]][[statnum]][simnum,varnum]<-sg$p.table[2,statnum]
  p=p+1
}
stat_list[[splitnum]][[5]][simnum,varnum]<-r2_delta

########end added stuff#############


for (r in 1:4) 
{
  user_data=runs[[r]]
  #ivs=readLines( paste0("/home/max/Documents/linear_mixed_model_abcd/ivs_",run_names[[r]],".txt"))
  col=which( colnames(user_data)=="puberty")
  col2=col+2
  ivs=names(user_data[col2:(col2+1)])
  if (r>1){covs=covs[-length(covs)]}
  covs=c(covs,spec_covs[c(r)])
  #counter+1
  
  for (d in ivs[1:2])
  {
    for (i in dvs)
    {
      data<-user_data
    
#input your covariates here in cov.fixed and randome effect variables in rand.var
inputs = list(dep.var.=list(d),ind.var.=list(i),usercov.=list(),cov.fixed=list(covs),smo.var=list(""),log.var=list(""),int.var=list(""),sq.var=list(""),gr.var=list(),ws.var=list(""),fl.var=list(),rand.var=list(re_covs))

# could filter with a json file like this in the fl.var section: "/home/max/Documents//linear_mixed_models/filterSets_ABCD_711b.json"

# User defined model code

##########################
##  user customization  ##
##########################
TEST = FALSE
# Change the coding from factor levels to numbers
#levels(data$physical_activity1_y) = c("", 0:7);
#data$physical_activity1_y = as.numeric(as.character(data$physical_activity1_y))

# Change the reference category of a factor level variable
# data$sex <- factor(data$sex, levels=c("M","F"))

#print(toJSON(inputs))


extract.variables = function(a){
  rt = c()
  for (l in 1:length(a) ){ 
    if(length(a[[l]]) > 0){
      for(e in 1:length(a[[l]])){
        if(unlist(a[[l]][[e]]) != "")
          rt = c(rt, unlist(a[[l]][[e]]))
      } 
    }
  }
  rt_inster = c()
  for( item in 1:length(rt)){
    if(!is.character(rt[item])){
      next
    }
    if(length(unlist(strsplit(rt[item],"[*]"))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"[*]")))
    }
    else if(length(unlist(strsplit(rt[item],"[+]"))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"[+]")))
    }
    
    else if(length(unlist(strsplit(rt[item],"^2", fixed=TRUE))) > 1){
      rt_inster = c(rt_inster, unlist(strsplit(rt[item],"^2", fixed=TRUE)))
    }
    else if( length(regmatches(rt[item], gregexpr("(?<=\\().*?(?=\\))", rt[item], perl=T))[[1]] ) > 0 ){
      
      rt_inster = c(rt_inster, regmatches(rt[item], gregexpr("(?<=\\().*?(?=\\))", rt[item], perl=T))[[1]])
    }else{
      rt_inster = c(rt_inster, rt[item])
    }
  }
  rt = rt_inster
  return(rt);
}

if(length(unlist(inputs[['dep.var.']])) == 0 ){
  stop("Dependent variable is empty.")
}



varList.initial = extract.variables(inputs)
#if smooth by variables, need to split first
if(sum(grepl("by =",varList.initial))>=1){
  s.by.split = varList.initial[grepl("by =",varList.initial)]
  split.vars = unlist(strsplit(s.by.split, ", by = "))
  varList.initial = c(varList.initial,split.vars)
  varList.initial = varList.initial[!duplicated(varList.initial)]
}

vars.in.data = varList.initial[varList.initial %in% names(data)]
vars.keep = c("src_subject_id",'eventname',"rel_family_id","mri_info_device.serial.number",vars.in.data)
vars.keep = vars.keep[!duplicated(vars.keep)]

data = data[,vars.keep]

##################
##  functions   ##
##################

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

sep.vars = function(x){
  x = gsub(" ", "", x, fixed = TRUE)
  x = unlist(strsplit(x,"+",fixed=T))
  return(x)
}

censor =  function(x, fraction=.005){
  if(length(fraction) != 1 || fraction < 0 || fraction > 1){
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction/2, 1-fraction/2), na.rm = T)
  x[ x < lim[1] ] <- NA
  x[ x > lim[2] ] <- NA
  x
}

#########################
##  data  extraction   ##
#########################

trigger.warning = F

#exit script silently if @inputs is empyt
empty = T;
for (key in names(inputs)){
  #print(inputs[[key]]);
  #print(length(inputs[[key]]))
  if(length(inputs[[key]]) != 0 && inputs[[key]] != ""){
    empty = F;
  }
}

if(empty){
  options(warn=-1)
  opt <- options(show.error.messages=FALSE) 
  on.exit(options(opt)) 
  stop() 
}

#By using Rserve gamm4 is already loaded

### censor/windsorize first
wsVar = unlist(inputs[['ws.var']])
wsVar = sep.vars(wsVar)
if(length(wsVar)>0){
  for(ii in 1:length(wsVar)){
    if(class(data[,wsVar[ii]]) == "numeric"){
      data[,wsVar[ii]] = censor(data[,wsVar[ii]])
    }
  }
}

dependendVar   = unlist(inputs[['dep.var.']]);
dependendVar.name = NULL

if(length(dependendVar) == 0 ){
  stop("Dependent variable is empty.")
}


##if y is log-transformed...
if(substring(dependendVar,1,4) == "log("){
  dependendVar.name = substring(dependendVar,5,nchar(dependendVar)-1)
  if(sum(data[[dependendVar.name]] <= 0, na.rm=T) > 0){
    data[[dependendVar.name]][data[[dependendVar.name]] <= 0] = NA
    trigger.warning = T
    warning.logging.0 = paste0("1 or more log transformed variable contains values <=0. All <=0 values replaced with NA.")
  }
  data$Y.log = log(data[[dependendVar.name]])
  new.name = paste0("log.",dependendVar.name)
  names(data)[names(data) == "Y.log"] = new.name
  dependendVar = new.name
}

if( dependendVar %in% names(data)){
  if(is.factor(data[[dependendVar]])){
    if(nlevels(data[[dependendVar]]) > 2){
      stop("Categorical variables with more than 2 levels are not supported as dependent variables. \nConsider converting your categorical variable into a continuous variable.")
    }
  }
} else {
  stop(paste("Dependent variable <", dependendVar ,">does not exist in the database"));
}


if( dependendVar %in% names(data)){
  if(is.factor(data[[dependendVar]])){
    if(nlevels(data[[dependendVar]]) > 2){
      stop("Categorical variables with more than 2 levels are not supported as dependent variables. \nConsider converting your categorical variable into a continuous variable.")
    }
  }
} else {
  stop(paste("Dependent variable <", dependendVar ,">does not exist in the database"));
}


independendVar = unlist(inputs[['ind.var.']]);
usercovVar     = paste(unlist(inputs[['usercov.']]),  sep='+')

smoothVar.all = unlist(inputs[['smo.var']])
smoothVar.all = sep.vars(smoothVar.all)

logVar = unlist(inputs[['log.var']])
logVar = sep.vars(logVar)
#check if 0's in logged vars. if so, add 0.0001
strip.log = substring(logVar,5,nchar(logVar)-1)
if(length(strip.log)>0){
  if(sum(data[,strip.log]<=0 , na.rm=T) > 0){
    trigger.warning = T
    warning.logging.0 = paste0("1 or more log transformed variable contains values <=0. All <=0 values replaced with NA.")
    for(ii in 1:length(strip.log)){
      log.var_i = strip.log[ii]
      if(sum(data[,log.var_i]<=0 , na.rm=T) > 0){
        #data[[log.var_i]][data[[log.var_i]] == 0] = 0.0001
        data[[log.var_i]][data[[log.var_i]] <= 0] = NA
        
      }
    }
  }
}

interactionVar = unlist(inputs[['int.var']])

sqVar = unlist(inputs[['sq.var']])
sqVar = sep.vars(sqVar)
sqVar = substring(sqVar,1,nchar(sqVar)-2)

sqVar_SQUARED = NULL
if(length(sqVar)>0){
  for(ii in 1:length(sqVar)){
    sqVar_SQUARED[ii] = paste0(sqVar[ii],"_SQUARED")
    data[,sqVar_SQUARED[ii]] = data[,sqVar[ii]]^2
  }
}

groupVar = unlist(inputs[['gr.var']])
if(length(groupVar) > 0){
  if(is.character(groupVar) & nchar(groupVar) ==0){
    groupVar = NULL 
  } 
}

subsetVar   = unlist(inputs[['fl.var']]);
if(length(subsetVar) > 0){
  if(is.character(subsetVar) & nchar(subsetVar) ==0){
    subsetVar = NULL 
  } 
}
#Interacting grouping variable with independent variable
if(length(groupVar)>0){
  #may need to change it instead of searching for the string, it strips string first and looks for exact match
  is.smooth = grepl(independendVar, smoothVar.all)
  is.log =    independendVar %in% substring(logVar,5,nchar(logVar)-1)
  is.square = independendVar %in% sqVar
  ## If independent var is smooth
  if( sum(is.smooth) > 0 ){
    #replace s(independendVar) with s(independendVar, by = groupVar)
    smoothVar.all[is.smooth] = paste0("s(", independendVar,",by=",groupVar, ")")
  } else if(is.square){
    #else if squared, add var*groupvar and var_SQUARED*groupvar
    interactionVar = c(interactionVar, paste0(independendVar,"*",groupVar), paste0(independendVar,"_SQUARED*",groupVar) )
  } else if(is.log){    
    # else if log, interact with log(var)
    log.independent = logVar[independendVar == substring(logVar,5,nchar(logVar)-1)]
    interactionVar = c(interactionVar, paste0(log.independent,"*",groupVar) )
  } else {
    # else make normal interaction
    interactionVar = c(interactionVar, paste0(independendVar,"*",groupVar) )
  }
}  #may need to do another if with ^2 independent variables
smoothVarInt.ind = grepl(",by=", smoothVar.all)
smoothVarInt = smoothVar.all[smoothVarInt.ind]
smoothVar =    smoothVar.all[!smoothVarInt.ind]

if(0 %in% nchar(sqVar)) sqVar = character()

print(sqVar)
print(paste("length sqVar",length(sqVar)))
#nestVar = c("Site", "FamilyID")
#usercovVar =  usercovVar[!(usercovVar %in% nestVar)]

#TODO: seperate Site and Familiy to another catagory of random effect

#if(include.random.site){
#  inputs[['cov.fixed']][[which(unlist(inputs[['cov.fixed']]) == "abcd_site")]] = NULL
#}
#inputs[['cov.fixed']][[which(unlist(inputs[['cov.fixed']]) == "rel_family_id")]] = NULL

covfixedVar    = paste(unlist(inputs[['cov.fixed']]), sep='+')

smoothVarInt.stripped.term1 = ""
smoothVarInt.stripped.term2 = ""


if(length(smoothVarInt)>0){
  smoothVarInt.stripped.term1 = unlist(lapply( strsplit(smoothVarInt,","), function(x)x[[1]]))
  smoothVarInt.stripped.term1 = substring(smoothVarInt.stripped.term1,3,nchar(smoothVarInt.stripped.term1))
  smoothVarInt.stripped.term2 = unlist(lapply( strsplit(smoothVarInt,"by="), function(x)x[[2]]))
  smoothVarInt.stripped.term2 = substring(smoothVarInt.stripped.term2,1,nchar(smoothVarInt.stripped.term2)-1)
}

### if usercovVar's have been transformed, then they are stored in usercovVar as well as the 
### transformed vars (smoothVar, logVar, etc.), in which they will need to be removed from
### usercovVar before putting into formula
print(paste("before remove",usercovVar))
print(c( substring(smoothVar,3,nchar(smoothVar)-1),
         smoothVarInt.stripped.term1,
         smoothVarInt.stripped.term2,
         substring(logVar,5,nchar(logVar)-1),
         sqVar ))
cov.ind.remove = usercovVar %in% c( substring(smoothVar,3,nchar(smoothVar)-1),
                                    smoothVarInt.stripped.term1,
                                    smoothVarInt.stripped.term2,
                                    substring(logVar,5,nchar(logVar)-1),
                                    sqVar )
if(sum(cov.ind.remove)>0 ){
  usercovVar = usercovVar[!cov.ind.remove]
}
print(paste("after remove",usercovVar))

#remove covfixedVar when it is a smooth, smooth interaction (first variable), or log -- should usually only be age of covfixedVar
#ok to keep it in as interaction, squared, or smooth interaction term
covfixedVar.ind.remove = covfixedVar %in% c( substring(smoothVar,3,nchar(smoothVar)-1),
                                             smoothVarInt.stripped.term1,
                                             substring(logVar,5,nchar(logVar)-1))
if(sum(covfixedVar.ind.remove)>0 ){
  covfixedVar = covfixedVar[!covfixedVar.ind.remove]
}


form_arr = c(independendVar, usercovVar,covfixedVar, smoothVar.all, logVar, interactionVar, sqVar, sqVar_SQUARED);
#form_arr = c(independendVar, usercovVar,covfixedVar, smoothVar, logVar, interactionVar, sqVar, sqVar_SQUARED);
#form_arr = c(independendVar, usercovVar,covfixedVar, smoothVar, logVar, interactionVar, sqVar, sqVar_SQUARED, groupVar);

### similarly for the independent variable...
#if independent variable is a smooth variable, log variable, or squared variable remove independendendVar from form_arr
if(independendVar %in% c(  substring(smoothVar,3,nchar(smoothVar)-1),
                           smoothVarInt.stripped.term1,
                           substring(logVar,5,nchar(logVar)-1),
                           sqVar ) ){
  form_arr = c(usercovVar,covfixedVar, smoothVar.all, logVar, interactionVar, sqVar, sqVar_SQUARED);
  #form_arr = form_arr[form_arr != independendVar]
}
form_arr = form_arr[form_arr!=""]
#########################
##  remove duplicates  ##
#########################
#take out duplicate variables
form_arr = form_arr[!duplicated(form_arr)]
formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
#get variables involve in the formula
varList = all.vars(as.formula(formulastr));
varList.independent = varList[-1]

#trigger.warning = F
#take out duplicated variables according to their values
if(length(varList.independent)>1){
  var.combin = combn(varList.independent,2)
  identical.values = c()
  for(ii in 1:ncol(var.combin)){
    vars_i = var.combin[,ii]
    #if two variables perfectly correlated, store both variables
    if(cor(as.numeric(data[,vars_i[1]]),as.numeric(data[,vars_i[2]]), use = "complete.obs") %in% c(1,-1)){
      # identical.values[ii] = list(vars_i)
      identical.values = c(identical.values,list(vars_i))
    }
  }
  vars.duplicate.remove = c()
  if(length(identical.values) > 0){
    for(ii in 1:length(identical.values)){
      vars_both_i = identical.values[[ii]]
      #if independendVar is in vars_both, remove other variable, otherwise doesn't matter which is removed
      if(independendVar %in% vars_both_i){
        vars.duplicate.remove[ii] = vars_both_i[vars_both_i != independendVar]
      } else{
        vars.duplicate.remove[ii] = vars_both_i[2]
      }
    }
    trigger.warning = T
    warning.duplicates = paste0("perfectly correlated variables detected, removed variable '", vars.duplicate.remove, "'")
    form_arr = form_arr[!form_arr %in% vars.duplicate.remove]
    formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
    varList = all.vars(as.formula(formulastr));
    varList.independent = varList[-1]
  }
}

print(varList)
print(formulastr)

#########################
##  data manipulation  ##
#########################
#data = data[data$eventname == "baseline_year_1_arm_1",]
print(summary(data[[independendVar]]))
#if independent variable has 5 or less unique values change it to character/factor variable
categorical.independent = FALSE
#if( length(table(data[[independendVar]])) < 6 ){
#  data[[independendVar]] = as.character(data[[independendVar]])
#  categorical.independent = TRUE
#} else{
#  data[[independendVar]] = as.numeric(as.character(data[[independendVar]]))
#}

if(class(data[[independendVar]]) != "numeric"){
  categorical.independent = TRUE
}

#user defined covariates
#for(ucov in unlist(inputs[['usercov.']])){
#  data[[ucov]] = as.numeric(as.character(data[[ucov]]))
#}

print(summary(data[[independendVar]]))
#determine if logistic regression or not
categorical.dependent = FALSE
data[[dependendVar]][data[[dependendVar]] == ""] = NA
if( length(table(data[[dependendVar]])) == 2 ){
  data[[dependendVar]] = as.factor(data[[dependendVar]])
  categorical.dependent = TRUE
} else{
  data[[dependendVar]] = as.numeric(as.character(data[[dependendVar]]))
}
#data[[dependendVar]] = as.numeric(as.character(data[[dependendVar]]))

#if("demo_prnt_marital_v2" %in% colnames(data)){
# Type of household
# 1 = Married, 6 = Living with a partner
# 2 = Widowed, 3 = Divorced, 4 = Separated, 5 = Never married
#  marital_v2.old = data$demo_prnt_marital_v2
#  data$demo_prnt_marital_v2 = NA
#  data$demo_prnt_marital_v2[marital_v2.old %in% c(1,6)] = 1
#  data$demo_prnt_marital_v2[marital_v2.old %in% 2:5] = 0
#}

# if("gender" %in% colnames(data)){
#   data = data[data$gender %in% c("M","F"),]
# }
#if("race.ethnicity" %in% colnames(data)){
#  data$race.ethnicity = as.factor(data$race.ethnicity)
#}
#income
##NDA using demo_prtnr_income_v2
#if("demo_comb_income_v2" %in% colnames(data)){
#  data$hhinc = NA
#  data$hhinc[data$demo_comb_income_v2 %in% 1:6]  = "[<50K]"
#  data$hhinc[data$demo_comb_income_v2 %in% 7:8]  = "[>=50K&<100K]"
#  data$hhinc[data$demo_comb_income_v2 %in% 9:10] = "[>=100K]"
#  data$demo_comb_income_v2 = data$hhinc
#}
if("household.income.bl" %in% names(data)){
  data$household.income.bl = as.character(data$household.income.bl)
  data$household.income.bl[data$household.income.bl == "[>=50K & <100K]"] = "[>=50K& <100K]"
  data$household.income.bl = as.factor(data$household.income.bl)
}



#data = data[c("src_subject_id","rel_family_id","abcd_site",varList)]
##################
##  subset data ##
##################
if(length(subsetVar)>0){
  json_data = rjson::fromJSON(file = subsetVar);
  subset = data.frame(src_subject_id=unlist(lapply(json_data[[1]]$set,function(d){ d[1] })), eventname=unlist(lapply(json_data[[1]]$set,function(d){d[2]})));
  data = merge(subset, data, all.x = T, all.y = F)
}


##################################################
##  identify variables of different timepoints  ##
##################################################
####before complete.cases, need to figure out which variables are in which events, and potentially transform


# var.timepoints.with.data = matrix(NA, nrow = length(varList), ncol =  length(unique(data$eventname))+1)
# var.timepoints.with.data[,1] = varList
var.timepoints.with.data = list()
for(ii in 1:length(varList)){
  var_i = varList[ii]
  valid.data.per.timepoint = aggregate(data[var_i], list(data$eventname), function(x)sum(!is.na(x)))
  timepoints.with.data = as.character(valid.data.per.timepoint$Group.1[valid.data.per.timepoint[var_i] != 0])
  # var.timepoints.with.data[ii,(1:length(timepoints.with.data))+1] = timepoints.with.data
  var.timepoints.with.data[[ii]] = timepoints.with.data
}
n.timepoints.per.var = unlist(lapply(var.timepoints.with.data,length))

dep.timepoints = unlist(var.timepoints.with.data[dependendVar == varList])
indep.timepoints = unlist(var.timepoints.with.data[independendVar == varList])

#if DEPENDENT is follow-up only (no baseline visits), and INDEPENDENT is baseline only, shift DEPENDENT rows to baseline
if(sum(dep.timepoints != "baseline_year_1_arm_1") >= 1 & sum(dep.timepoints == "baseline_year_1_arm_1") == 0 & 
   sum(indep.timepoints != "baseline_year_1_arm_1") == 0 & sum(indep.timepoints == "baseline_year_1_arm_1") == 1){
  
  data.dep = data[data$eventname == dep.timepoints , c("src_subject_id","eventname",dependendVar)]
  data.dep$eventname[data.dep$eventname == dep.timepoints] = indep.timepoints
  data[[dependendVar]] = NULL
  data = merge(data, data.dep, by = c("src_subject_id","eventname"), all = T)
}

###################################################


data = data[complete.cases(data),]

print(dim(data))

###################################################
##  remove columns from model if 1 unique value  ##
###################################################

#less.than.2.levels = sapply(data , function(x) length(unique(x)) ) < 2
less.than.2.levels = sapply(data[,varList] , function(x) length(unique(x)) ) < 2

#trigger.warning = F
if(sum(less.than.2.levels) >0 ){
  trigger.warning = T
  #vars.to.remove = names(data)[less.than.2.levels]
  vars.to.remove = names(less.than.2.levels)[less.than.2.levels]
  
  
  #stop script if dependent variable has less than 2 unique values
  if(dependendVar %in% vars.to.remove){
    stop(paste0("'", dependendVar,"'" ," variable has <2 unique values"))
  }
  if(independendVar %in% vars.to.remove){
    stop(paste0("'", independendVar,"'" ," variable has <2 unique values"))
  }
  
  #remove vars.to.remove from form_arr, formulastr, varList, and varList.independent
  form_arr = form_arr[!form_arr %in% vars.to.remove]
  formulastr = paste(dependendVar," ~ ",paste(form_arr,collapse='+'))
  
  varList = all.vars(as.formula(formulastr));
  varList.independent = varList[-1]
  
  warning.1.unique.value = paste0("'", vars.to.remove,"'" ," variable has <2 unique values - removed from model")
}

########################################################################
##  remove independendVar (& groupVar) from formula to get delta R^2  ##
########################################################################

run.effect.size = T

if(run.effect.size){
  
  
  logVar.stripped    = substring(logVar,5,nchar(logVar)-1)
  smoothVar.stripped = substring(smoothVar,3,nchar(smoothVar)-1)
  
  if(length(interactionVar)>0){
    ##INTERACTION
    interaction.stripped.list = strsplit(interactionVar,"*",fixed=T)
  } else{
    interaction.stripped.list = character()
  }
  
  
  if(independendVar %in% smoothVar.stripped){
    ##SMOTH
    smoothVar.remove = smoothVar[smoothVar.stripped == independendVar]
    form_arr2 = form_arr[!(form_arr %in% smoothVar.remove)]
  } else if(independendVar %in% smoothVarInt.stripped.term1){
    ##SMOTH INTERACTION
    smoothVarInt.remove = smoothVarInt[smoothVarInt.stripped.term1 == independendVar]
    form_arr2 = form_arr[!(form_arr %in% smoothVarInt.remove)]
  } else if(independendVar %in% logVar.stripped){
    ##LOG
    logVar.remove = logVar[logVar.stripped == independendVar]
    form_arr2 = form_arr[!(form_arr %in% logVar.remove)]
    if(length(groupVar)>0){
      logVarinteraction.remove = c(paste0(logVar.remove,"*",groupVar) , groupVar)
      form_arr2 = form_arr2[!(form_arr2 %in% logVarinteraction.remove)]
    }
  } else if(independendVar %in% sqVar){
    ##SQUARED
    sqVar.remove = c(independendVar, paste0(independendVar,"_SQUARED"))
    form_arr2 = form_arr[!(form_arr %in% sqVar.remove)]
    if(length(groupVar)>0){
      sqVarinteraction.remove = c(paste0(sqVar.remove,"*",groupVar), groupVar)
      form_arr2 = form_arr2[!(form_arr2 %in% sqVarinteraction.remove)]
    }
  } else if(independendVar %in% unlist(interaction.stripped.list)){
    ##INTERACTION
    independ.location = c()
    for(ii in 1:length(interaction.stripped.list)){
      int.stripped.list_i = interaction.stripped.list[[ii]]
      independ.location[ii] = sum(independendVar == int.stripped.list_i)
    }
    interaction.remove = interactionVar[independ.location == 1]
    interactionVars.remove = unlist(strsplit(interaction.remove,"*",fixed=T))
    all.interactionVars.remove = c(interaction.remove,interactionVars.remove)
    form_arr2 = form_arr[!(form_arr %in% all.interactionVars.remove)]
  } else {
    #NO INDEPENDENT VAR TRANSFORMATION
    form_arr2 = form_arr[!(form_arr %in% independendVar)]
  }
  
  if(length(form_arr2) > 0){
    formulastr2 = paste(dependendVar," ~ ",paste(form_arr2,collapse='+'))
  } else{
    formulastr2 = paste(dependendVar," ~ 1")
  }
  
  
}

########################################
##  make sure log doesn’t create NAs  ##
########################################

# if(length(logVar.stripped)>0){
#   for(ii in 1:length(logVar.stripped)){
#     min.log_i = min(data[[logVar.stripped[ii]]])
#     if(min.log_i < 0)  
#       data[[logVar.stripped[ii]]] = data[[logVar.stripped[ii]]] - min.log_i + 1
#   }
# }



print("Before calling the model")



######################
##  random effects  ##
######################
include.random.scanner = "mri_info_device.serial.number" %in% unlist(inputs$rand.var)
include.random.site    = "abcd_site" %in% unlist(inputs$rand.var)
include.random.subject = "src_subject_id" %in% unlist(inputs$rand.var)

#IF less than 2 unique "eventname" & "include.random.subject" toggled on, 
#   need to toggle off random subject effect and produce warning (or error)
if(length(unique(data$eventname)) <2 & include.random.subject){
  include.random.subject = FALSE
  trigger.warning = T
  warning.1.unique.eventname = paste0("Your model data contains a single data point per subject (event name “baseline_year_1_arm_1”) only. To prevent problems with model convergence suggest to remove SUBJECT from the list of random effects. It is not required with this particular cross-sectional analysis.")
}
#IF 2+ unique "eventname" & "include.random.subject" toggled off, 
#   need to toggle on random subject effect and produce warning (or error)
if(length(unique(data$eventname)) > 1 & !include.random.subject){
  # include.random.subject = TRUE
  trigger.warning = T
  warning.2.unique.eventname = paste0("Your model data contains multiple data points per subject (event name “baseline_year_1_arm_1”). Suggest to include SUBJECT from the list of random effects.")
}



formula.random.str = c()
if(include.random.site & !include.random.subject){
  formula.random.str = "(1|abcd_site/rel_family_id)"
} else if(include.random.scanner & !include.random.subject){
  formula.random.str = "(1|mri_info_device.serial.number/rel_family_id)"
} else{
  if(include.random.site) formula.random.str[1] = "(1|abcd_site)"
  if(include.random.scanner) formula.random.str[1] = "(1|mri_info_device.serial.number)"
  formula.random.str[2] = "(1|rel_family_id)"
  if(include.random.subject) formula.random.str[3] = "(1|src_subject_id)"
}
formula.random.str = formula.random.str[!is.na(formula.random.str)]
formula.random.str = paste0("~",paste(formula.random.str, collapse = "+"))
formula.random = formula.random.str

#if(include.random.site == T & include.random.subject == T){
#  formula.random = "~(1|abcd_site)+(1|rel_family_id)+(1|src_subject_id)"
#} else if(include.random.site == F & include.random.subject == T){
#  formula.random = "(1|rel_family_id)+(1|src_subject_id)"
#} else if(include.random.site == T & include.random.subject == F){
#  formula.random = "~(1|abcd_site/rel_family_id)"
#} else if(include.random.site == F & include.random.subject == F){
#  formula.random = "~(1|rel_family_id)"
#}

#################
##  run model  ##
#################

if(categorical.dependent){ #logistic regression
  #model  = gamm4(as.formula(formulastr) , random = as.formula(formula.random), family = binomial , data = data)
  #model2 = gamm4(as.formula(formulastr2) , random = as.formula(formula.random), family = binomial , data = data)
  #intercept model
  #model3 = gamm4(as.formula(paste0(dependendVar, " ~1")) , random = as.formula(formula.random), family = binomial, data = data)
  
  formula_list = c(formulastr, formulastr2, paste0(dependendVar, " ~1"));
  cl <- makeCluster(3)
  registerDoParallel(3)
  model_list = foreach(ii= 1:length(formula_list)) %dopar% {
    gamm4(as.formula(formula_list[ii]), random = as.formula(formula.random), family = binomial, data = data)
  }
  print(model_list)
  model = model_list[[1]]
  model2 = model_list[[2]]
  model3 = model_list[[3]]
  stopCluster(cl)
  
  
} else{ #linear regression
  formula_list = c(formulastr, formulastr2, paste0(dependendVar, " ~1"));
  cl <- makeCluster(3)
  registerDoParallel(3)
  model_list = foreach(ii= 1:length(formula_list),.packages='gamm4') %dopar% {
    gamm4(as.formula(formula_list[ii]), random = as.formula(formula.random), data = data)
  }
  print(model_list)
  model = model_list[[1]]
  model2 = model_list[[2]]
  model3 = model_list[[3]]
  stopCluster(cl)
}



print(summary(model))
print(summary(model2$gam))

########################
##  other statistics  ##
########################
#N subjects, aic, bic, & r^2
n = summary(model$gam)$n
aic = round(AIC(model$mer),2)
bic = round(BIC(model$mer),2)
#r2   = round(summary(model$gam)$r.sq,5)
#r2_2 = round(summary(model2$gam)$r.sq,5)
r2       = round(as.numeric(r.squaredLR(model$mer,model3$mer)),5)
r2_delta = round(as.numeric(r.squaredLR(model$mer,model2$mer)),5)
aic2 = round(AIC(model2$mer),2)
bic2 = round(BIC(model2$mer),2)


#compute R^2 and delta R^2
#model.r2 = r2beta(model$mer, method = 'nsj')
#model2.r2 = r2beta(model2$mer, method = 'nsj')
#r2dt = r2dt(model.r2,model2.r2)

#r2 = round(model.r2$Rsq[model.r2$Effect ==  "Model"],5)
#r2_delta = round(r2dt$d,5)

#change name in sink to change where output is printed
sink(output_name,append=TRUE)
print(formula.random.str)
print(capture.output(summary(model$gam)))
print(paste0("r2 = ",r2))
print(paste0("r2_delta = ",r2_delta))
sink()

  }
}
}
