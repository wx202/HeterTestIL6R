rm(list=ls())
timestart<-Sys.time()
setwd('/Users/xuanwang/Desktop/Projects/IL6R_Analyses/Xuan')
library(readr)
library("dplyr")
library('tidyr')
source("~/Desktop/Projects/IL6R_Analyses/Xuan/test_heterogeneity.R")

{
dat1=read_csv(paste0('download0616/ICD_results.BB_phewas_result_RaceWhite.csv'))
dat2=read_csv(paste0('download0616/ICD_results.BB_phewas_result_RaceBlack.csv'))
mvp.icd=rbind(dat1,dat2)
mvp.icd=mvp.icd[,c("phenotype","white","Estimate","Std. Error","Pr(>|z|)","mean.y","sample.size")]
mvp.white=mvp.icd[mvp.icd$white==1,]
mvp.white$white=NULL
colnames(mvp.white)=c("phe.codes","Estimate_White","Std.err_White",
                      "p.value_White", "mean.y_White", "sample.size_White")
mvp.black=mvp.icd[mvp.icd$white==0,]
mvp.black$white=NULL
colnames(mvp.black)=c("phe.codes","Estimate_Black","Std.err_Black",
                      "p.value_Black", "mean.y_Black","sample.size_Black")
mvp.icd=left_join(mvp.black,mvp.white,by='phe.codes')

## lab 
dat1=read_csv(paste0('download0616/lab_results.BB_phewas_result_RaceWhite.csv'))
dat2=read_csv(paste0('download0616/lab_results.BB_phewas_result_RaceBlack.csv'))
mvp.lab=rbind(dat1,dat2)
mvp.lab=mvp.lab[,c("phenotype","white","Estimate","Std. Error","Pr(>|t|)","mean.y","sample.size")]
mvp.white=mvp.lab[mvp.lab$white==1,]
mvp.white$white=NULL
colnames(mvp.white)=c("phe.codes","Estimate_White","Std.err_White",
                      "p.value_White", "mean.y_White", "sample.size_White")
mvp.black=mvp.lab[mvp.lab$white==0,]
mvp.black$white=NULL
colnames(mvp.black)=c("phe.codes","Estimate_Black","Std.err_Black",
                      "p.value_Black", "mean.y_Black","sample.size_Black")
mvp.lab=left_join(mvp.black,mvp.white,by='phe.codes')

## infection.code
dat1=read_csv(paste0('download0616/infection.code_results.BB_phewas_result_RaceWhite.csv'))
dat2=read_csv(paste0('download0616/infection.code_results.BB_phewas_result_RaceBlack.csv'))
mvp.infection=rbind(dat1,dat2)
mvp.infection=mvp.infection[,c("phenotype" ,"white" , "Estimate","Std. Error","Pr(>|z|)",
                   "mean.y","sample.size")]
mvp.white=mvp.infection[mvp.infection$white==1,]
mvp.white$white=NULL
colnames(mvp.white)=c("phe.codes","Estimate_White","Std.err_White",
                      "p.value_White", "mean.y_White", "sample.size_White")
mvp.black=mvp.infection[mvp.infection$white==0,]
mvp.black$white=NULL
colnames(mvp.black)=c("phe.codes","Estimate_Black","Std.err_Black",
                      "p.value_Black", "mean.y_Black","sample.size_Black")
mvp.infection=left_join(mvp.black,mvp.white,by='phe.codes')

}

aa=gsub('phe_','',mvp.icd$phe.codes)
aa=gsub('_','.',aa)
aa=as.numeric(aa)
ind=which(!aa%in%floor(aa[aa%%1!=0]) )
mvp.icd=mvp.icd[ind,]

mvp.lab=unique(mvp.lab)
aa=sort(c(mvp.lab$phe.codes[grepl('median', mvp.lab$phe.codes)]))
mvp.lab=mvp.lab[mvp.lab$phe.codes%in%aa,]

{
#### heter test
outcome='lab'
if (outcome=='icd'){dat=mvp.icd}
if (outcome=='lab'){dat=mvp.lab}
if (outcome=='infection'){dat=mvp.infection}
if (outcome=='lab'){
  dat=na.omit(dat)
  ind=which(dat$Std.err_Black>0 & dat$Std.err_White>0)
  dat=dat[ind,]
  # dat=dat[dat$Estimate_Black>-1 & dat$Estimate_Black<1 & dat$Estimate_White>-1 & dat$Estimate_White<1,]
}else{
  dat=na.omit(dat)
  ind=which(dat$mean.y_Black>5e-3 & dat$mean.y_White>5e-3 )
  dat=dat[ind,]
  # dat=dat[dat$Std.err_Black>0 & dat$Std.err_White>0 &
  #           dat$Estimate_Black>-1 & dat$Estimate_Black<1 & dat$Estimate_White>-1 & dat$Estimate_White<1,]
}

beta_hat <- rbind(dat$Estimate_Black, dat$Estimate_White)
sigma_hat <- rbind(dat$Std.err_Black, dat$Std.err_White)
pvl_res <- test_heterogeneity_beta(beta_hat, sigma_hat, level=0.1)
dat$phe.codes[pvl_res$het_effect]

# dat$select=0
# dat$select[pvl_res$het_effect]=1

tau <- 0.005
res_test <- adapt_test(pvl_res$p_het, pvl_res$p_mean, tau = tau, level = 0.1)
dat$phe.codes[res_test$select]

# tau <- 0.005
# tau_mean <- 0.005
# # which(pvl_res$p_mean < tau_mean)
# res_test <- adapt_test(pvl_res$p_het, ifelse(pvl_res$p_mean < tau_mean, 0.2, 0.8),
#                        tau = tau, level = 0.05)

res <- data.frame(dat,res_test$p.mean, res_test$p.heter)
res$select=0
res$select[res_test$select]=1
res=res[order(-res$select),] 
res.new=res

res.all=data.frame(unique(sort(c(res.new$phe.codes))))
names(res.all)='phe.codes'
ver=20.1; version='B'
mvpicd=read_csv(paste0('Isabella/ICD_White&Black-MVP_Phe-v',ver,'-Gen-v',ver,'_Summary.csv'))
map=mvpicd[,c("phenotype.names","phe.codes")]
res.all=data.frame('phenotype.names'=map$phenotype.names[match(res.all$phe.codes,map$phe.codes)],res.all)
# res.all=left_join(res.all,res.old,by='phe.codes')
res.all=left_join(res.all,res.new,by='phe.codes')
res.all=res.all[order(-res.all$select),]
# res.all=res.all[order(-res.all$select.y),]

if (outcome %in% c('icd','infection')){
res.all$ES_Black=exp(res.all$Estimate_Black)
res.all$ESlow_Black=exp(res.all$Estimate_Black-1.96*res.all$Std.err_Black)
res.all$ESupp_Black=exp(res.all$Estimate_Black+1.96*res.all$Std.err_Black)
res.all$ES_White=exp(res.all$Estimate_White)
res.all$ESlow_White=exp(res.all$Estimate_White-1.96*res.all$Std.err_White)
res.all$ESupp_White=exp(res.all$Estimate_White+1.96*res.all$Std.err_White)
}
if (outcome=='lab'){
  res.all$ES_Black= res.all$Estimate_Black
  res.all$ESlow_Black= (res.all$Estimate_Black-1.96*res.all$Std.err_Black)
  res.all$ESupp_Black= (res.all$Estimate_Black+1.96*res.all$Std.err_Black)
  res.all$ES_White= res.all$Estimate_White
  res.all$ESlow_White= (res.all$Estimate_White-1.96*res.all$Std.err_White)
  res.all$ESupp_White= (res.all$Estimate_White+1.96*res.all$Std.err_White)
}

write.csv(res.all,file=paste0("out/",ver,'_',outcome,"_test_",Sys.Date(),".csv")  )
}

# aa <- read_csv("out/icd_test_2022-04-22.csv")
# # aa <- read_excel("out/lab_test_2022-04-13.xls")
# # aa <- read_excel("out/infection_test_2022-04-13.xls")
# aa=aa[aa$phe.codes=='phe_288_2',]
# m1=aa$Estimate_Black.x;m2=aa$Estimate_White.x;n1=aa$sample.size_Black.x;n2=aa$sample.size_White.x
# s1=aa$Std.err_Black.x;s2=aa$Std.err_White.x
# t.test2(m1,m2,s1,s2,n1,n2)
# m1=aa$Estimate_Black.y;m2=aa$Estimate_White.y;n1=aa$sample.size_Black.y;n2=aa$sample.size_White.y
# s1=aa$Std.err_Black.y;s2=aa$Std.err_White.y
# t.test2(m1,m2,s1,s2,n1,n2)
# 
# 
# t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
# {
#   if( equal.variance==FALSE ) 
#   {
#     se <- sqrt( (s1^2/n1) + (s2^2/n2) )
#     # welch-satterthwaite df
#     df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
#   } else
#   {
#     # pooled standard deviation, scaled by the sample sizes
#     se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
#     df <- n1+n2-2
#   }      
#   t <- (m1-m2-m0)/se 
#   dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
#   names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
#   return(dat) 
# }
# 
# ##MGB
#   # mvp$phe.codes=gsub('phe_','',mvp$phe.codes)
#   # mvp$phe.codes=gsub('_','.',mvp$phe.codes)
#   # mvp$phe.codes=as.numeric(mvp$phe.codes)
#   # mvp.black=mvp[,c("phe.codes","Estimate_Black","p.value_Black","sample.size_Black")]
#   # colnames(mvp.black)=c('code','est','p','size')
#   # mvp.black$race=0
#   # mvp.white=mvp[,c("phe.codes","Estimate_White","p.value_White","sample.size_White")]
#   # colnames(mvp.white)=c('code','est','p','size')
#   # mvp.white$race=1
#   
#   # mgb.black=read_csv('/Users/xuanwang/Desktop/Projects/IL6R_Analyses/Xuan/results_MGB_0313/result_IL6R_phewas_Black_pca_no.csv')
#   # mgb.black=mgb.black[,c("phecode","beta","p.value","sample.size")]
#   # colnames(mgb.black)=c('code','est','p','size')
#   # mgb.black$race=0
#   # 
#   # mvp.black=mvp.black[match(mgb.black$code,mvp.black$code),]
#  

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 
