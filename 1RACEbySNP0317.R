rm(list=ls())
library(geepack)
library(dplyr)
setwd("/group/research/mvp001/Xuan/IL6R")

# snps
load("/group/research/mvp001/Xuan/SNPRelease4/results-snp.RData")
colnames(res)[1]='id'
colnames(res)[-1]=paste0('x',1:127)

# race
phe.type = 'B'
gen.type = 'B'

# load(paste0('/scratch/scratch10/mvp001/Nogues_Isabella/IL6R_rs2228145/DEMO/demo_',phe.type,gen.type,'.RData'))
# demo=DEMO[,c('mvp001_id','race.white', 'race.black')]
# colnames(demo)=c('id','white','black')
if (phe.type == "A") { 
  code.file.list <- list.files(
    "/data/data1/mvp001/clinical_data/20201125", "PheWasCode", full.names = TRUE) 
  count.file <- code.file.list[3]
}else if (phe.type == "B"){
  code.file.list <- list.files(
    "/data/data1/mvp001/clinical_data/20210721", "PheWasCode", full.names = TRUE) 
  count.file <- code.file.list[6]
} 
data = fread(count.file,data.table=FALSE)

demo.var=c("mvp001_id","race")
demo = as.data.table(data[,demo.var])
demo[, race := ifelse(race == "N", NA_character_, race)]
demo = na.omit(demo)
demo = demo[, list(mvp001_id, white = ifelse(race == "W", 1, 0), black = ifelse(race == "B", 1, 0))]
demo = as.data.frame(demo)
colnames(demo)[1]='id'

# demo.var=c("mvp001_id","race")
# demo = as.data.table(data[,demo.var])
# demo[, race := ifelse(race == "N", NA_character_, race)]
# demo[,hispanic := ifelse(hispanic == 9, NA_integer_,hispanic)]
# demo$hispanic = as.numeric(demo$hispanic == 1)
# demo = na.omit(demo)
# demo = demo[, list(mvp001_id,
#                    white = ifelse(race == "W" & hispanic == 0, 1, 0),
#                    black = ifelse(race == "B" & hispanic == 0, 1, 0))]
# demo = as.data.frame(demo)
# colnames(demo)[1]='id'

dat=left_join(demo,res,by='id')
dat = na.omit(dat)


## logistic regression
system.time({
formula = as.formula(paste0("white ~", paste( colnames(dat)[-(1:3)],collapse='+') ))
fit= geeglm(formula, family ="binomial", id=1:dim(dat)[1], corstr="independence", data=dat)
pred=predict(fit,dat,type='response')
})

system.time({
formula = as.formula(paste0("black ~", paste( colnames(dat)[-(1:3)],collapse='+') ))
fit.black= geeglm(formula, family ="binomial", id=1:dim(dat)[1], corstr="independence", data=dat)
pred.black=predict(fit.black,dat,type='response')
})

# choose cutoff
cutoff=seq(0.99,0,-0.01)
ppvs=sapply(1:length(cutoff),function(k){ppv.fun(dat$white,as.numeric(pred>cutoff[k]))})
out=data.frame(cutoff,ppvs)
sel=min(out$cutoff[out$ppvs>=0.95])
prediction.white=as.numeric(pred>sel)

accs=sapply(1:length(cutoff),function(k){acc.fun(dat$white,as.numeric(pred>cutoff[k]))})
out.white=data.frame(cutoff,t(accs))
colnames(out.white)[-1]=c('tpr','fpr','ppv','npv')

ppvs=sapply(1:length(cutoff),function(k){ppv.fun(dat$black,as.numeric(pred.black>cutoff[k]))})
out=data.frame(cutoff,ppvs)
sel=min(out$cutoff[out$ppvs>=0.95])
prediction.black=as.numeric(pred.black>sel)

accs=sapply(1:length(cutoff),function(k){acc.fun(dat$black,as.numeric(pred.black>cutoff[k]))})
out.black=data.frame(cutoff,t(accs))
colnames(out.black)[-1]=c('tpr','fpr','ppv','npv')

out=data.frame(dat$id,prediction.white,prediction.black)
save(out,file=paste0('predRACEbySNP_',phe.type,gen.type,'.rda') )

out=list('out.white'=out.white,'out.black'=out.black)
save(out,file=paste0('predRACEbySNP_',phe.type,gen.type,'ACC.rda') )




ppv.fun=function(label,prediction){
  tp=sum(label*prediction)
  fp=sum((1-label)*prediction)
  fn=sum(label*(1-prediction))
  tn=sum((1-label)*(1-prediction))
  ppv=tp/(tp+fp)
  # npv=tn/(tn+fn)
  # tpr=tp/(tp+fn)
  # fpr=fp/(fp+tn)
  # out=c(tpr,fpr,ppv,npv)
}

acc.fun=function(label,prediction){
  tp=sum(label*prediction)
  fp=sum((1-label)*prediction)
  fn=sum(label*(1-prediction))
  tn=sum((1-label)*(1-prediction))
  ppv=tp/(tp+fp)
  npv=tn/(tn+fn)
  tpr=tp/(tp+fn)
  fpr=fp/(fp+tn)
  out=c(tpr,fpr,ppv,npv)
}

