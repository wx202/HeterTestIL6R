## run code ##
## /data/data1/mvp001/clinical_data/20180425
library(bit64)
library(data.table)
library(geepack)
library(sandwich)
library(Matrix)
library(dplyr)

skew <- function(x)
{
  z = (x - mean(x)) / max(sd(x), 1e-10)
  mean(z ** 3)
}


## parameters
args = commandArgs(trailingOnly=TRUE)
print(args)
# quit(save='no')

#Lab, ICD, MAP, or infection.code
outcome = args[1]
#AA, BA, BB
#phe.type = args[2]
#gen.type = args[3]
#Usually defined as 2
#threshold = as.integer(args[3])
#1- 50
part = as.integer(args[2])
#5, no
#PC = args[5]
#afr, asn, eur, his
#race = args[6]

dat_num = as.integer(args[3]) #(code: 1-3; MAP: 1-2, lab: 1-2)

white = as.integer(args[4])


code.file.list.A <- list.files(
  "/data/data1/mvp001/clinical_data/20201125", 
  "PheWasCode", full.names = TRUE)

code.file.list.B <- list.files(
  "/data/data1/mvp001/clinical_data/20210721", 
  "PheWasCodesV2", full.names = TRUE)


phe288.file.list.A <- list.files(
  "/data/data1/mvp001/clinical_data/20220207", 
  "192", full.names = TRUE)

phe288.file.list.B <- list.files(
  "/data/data1/mvp001/clinical_data/20220207", 
  "201", full.names = TRUE)


lab.file.list <- list.files(
  "/data/data1/mvp001/clinical_data/20211101", 
  "PheWasLab", full.names = TRUE)

# To be changed once MAP results are available
MAP.file.list <- list.files(
  "/data/data1/mvp001/clinical_data/20210712", 
  "MAPMRM", full.names = TRUE)

infection.code.file.list <- list.files(
  "/data/data1/mvp001/clinical_data/20210713", 
  "SeriousInfections", full.names = TRUE)


## load outcome data
for (phe.type in c("A","B")){
	for(gen.type in c("A","B")){

		if(outcome=='ICD' | outcome =='caseControl'){
			if (phe.type == "A"){
				code.dat = fread(code.file.list.A[dat_num],data.table = F)
			} else {
				code.dat = fread(code.file.list.B[dat_num],data.table = F)
			}
			code = code.dat[,grep('phe',colnames(code.dat))]
			mvp001.ID =code.dat$mvp001_id
			code$mvp001_id=NULL
			fam = "binomial"
		}else if(outcome=='phe288'){
			if (phe.type == "A"){
				phe288.dat = fread(phe288.file.list.A[dat_num],data.table = F)
			} else {
				phe288.dat = fread(phe288.file.list.B[dat_num],data.table = F)
			}
			phe288 = phe288.dat[,grep('phe',colnames(phe288.dat))]
			mvp001.ID =phe288.dat$mvp001_id
			phe288$mvp001_id=NULL
			fam = "binomial"
		}else if(outcome=='MAP'){
			MAP.dat = fread(MAP.file.list[dat_num],data.table = F)
			MAP = MAP.dat[,grep('phe',colnames(MAP.dat))]
			mvp001.ID = MAP.dat$mvp001_id
			MAP$mvp001_id=NULL
			fam = "quasibinomial"
		}else if(outcome =='lab'){
			lab = fread(lab.file.list[dat_num],data.table = F)
			dt.nm = grep('_dt',colnames(lab))
			if (length(dt.nm) > 0) lab = lab[,-dt.nm]
			mvp001.ID =lab$mvp001_id
			lab$mvp001_id=NULL
			fam = "gaussian"
		}else if(outcome=='infection.code'){
			if (phe.type == "A"){
				code.dat = fread(code.file.list.A[dat_num],data.table = F)
			} else {
				code.dat = fread(code.file.list.B[dat_num],data.table = F)
			}
			infection.code.dat = fread(infection.code.file.list[dat_num],data.table = F)
			infection.code = matrix(0,nrow(code.dat),2)
			infection.code[,1] = as.numeric(code.dat$mvp001_id %in% infection.code.dat$mvp001_id)
			#IDs = intersect(code.dat$mvp001_id, infection.code.dat$mvp001_id)
			infection.code[which(code.dat$mvp001_id %in% infection.code.dat$mvp001_id),2] = infection.code.dat[na.omit(match(code.dat$mvp001_id,infection.code.dat$mvp001_id)),2]
			#infection.code[which(code.dat$mvp001_id %in% infection.code.dat$mvp001_id),2] = 
#infection.code.dat[intersect(code.dat$mvp001_id, infection.code.dat$mvp001_id),2]
		
			fam = c("binomial","poisson")
			mvp001.ID = code.dat$mvp001_id
		}else{
		#Generate error if neither code nor lab have been provided
		stop('outcome must be code, MAP, infection.code, or lab!')
		}

		## run the algorithm

		if(outcome=='ICD'| outcome=='caseControl'){
			phenotype.list = colnames(code)
		}else if(outcome=='MAP'){
			phenotype.list = colnames(MAP)
		}else if(outcome=='infection.code'){
			phenotype.list = c("serious_infections_binary","serious_infections_count")
		}else if(outcome == "lab"){
			phenotype.list = colnames(lab)
		} else {
			phenotype.list = colnames(phe288)
		}


		part.seq.n = round(length(phenotype.list)/50)
		if(outcome=='ICD'| outcome =='caseControl' | outcome =='MAP'){
		  if(part<50){
		    part.seq = 1:part.seq.n + part.seq.n*(part-1)
		  }else{
		    part.seq = (1+part.seq.n*(part-1)):length(phenotype.list)
		  }
		}else{
		    part.seq = seq_along(phenotype.list)
		}


		#Summary of regression analyses


		#for (race in c("afr","asn","eur","his")){
		#for (race in c("afr","asn","eur","his","all")){
			race = "all"
			print(race)
			res = NULL
			snp.name = 'IL6R'
			PC = ifelse(race == "all","no",5)
			#adjust for either race or PCs
			#load demographic + SNP data from Scratch space/folder
			if(tolower(PC)=='no'){ # adjust by race, instead of PCs
			  covariates = c(snp.name,'visitmons','totalcodes','age','sex1')
			  #covariates = c(snp.name,'visitmons','totalcodes','age','sex1','race.white','race.black')
			  output.file.white <- sprintf(
			    "/group/research/mvp001/Xuan/IL6R/output/%s_results.%s%s_phewas_result_RaceWhite_File%d_part%d.csv",outcome,phe.type,gen.type,dat_num, part)
			  output.file.black <- sprintf(
			    "/group/research/mvp001/Xuan/IL6R/output/%s_results.%s%s_phewas_result_RaceBlack_File%d_part%d.csv",outcome,phe.type,gen.type,dat_num, part)
			  load(paste0('/scratch/scratch10/mvp001/Nogues_Isabella/IL6R_rs2228145/DEMO/demo_',phe.type,gen.type,'.RData'))
			  DEMO$race.white=NULL; DEMO$race.black=NULL
			  load("/group/research/mvp001/Xuan/IL6R/predRACEbySNP_BB.rda")
			  colnames(out)=c('mvp001_id','race.white','race.black')
			  DEMO=left_join(DEMO,out,by='mvp001_id')
			  DEMO=na.omit(DEMO)
			  white.patients = DEMO[,"race.white"] == 1 & DEMO[,"race.black"] == 0
                          black.patients = DEMO[,"race.white"] == 0 & DEMO[,"race.black"] == 1
			  demo.white = DEMO[white.patients,]; demo.black = DEMO[black.patients,]
			}else{
			  PC = as.integer(PC)
			  #covariates = c(snp.name,'visitmons','totalcodes','age','sex1', paste('PC',1:PC,sep=''))
			  covariates = c(snp.name,'visitmons','totalcodes','age','sex1')
			  output.file <- sprintf(
			    "/scratch/scratch10/mvp001/Nogues_Isabella/IL6R_rs2228145/%s/results.PC.%s%s/%s/phewas_result_pc%d_File%d_part%d.csv", outcome,phe.type,gen.type, toupper(race),PC, dat_num, part)
			  load(paste0('/scratch/scratch10/mvp001/Nogues_Isabella/IL6R_rs2228145/DEMO/demo_',race,'_',phe.type,gen.type,'.RData'))
			}

			formula = as.formula(paste0(" y ~", paste(covariates,collapse='+'),sep='' ))
			demo.white = demo.white[,c("mvp001_id", covariates)]
			demo.black = demo.black[,c("mvp001_id", covariates)]
			if (white == 1) { 
				output.file= output.file.white; demo = demo.white
			  } else {
				output.file= output.file.black; demo = demo.black
			  }

			for (i in part.seq) {
			  # read a specific phenotype
			  if (outcome == "ICD") {
			    outcome.dat = data.frame(mvp001_id = mvp001.ID, y = ifelse(code[,i] >= 2, 1, 0))
 			  }else if (outcome == "phe288") {
			    outcome.dat = data.frame(mvp001_id = mvp001.ID, y = ifelse(phe288[,i] >= 2, 1, 0))			 
			  }else if (outcome == "infection.code") {
			    outcome.dat = data.frame(mvp001_id = mvp001.ID, y = as.matrix(infection.code)[,i]) 
			  } else if (outcome == 'MAP'){
			    outcome.dat = data.frame(mvp001_id = mvp001.ID,y=MAP[,i])
			  } else if (outcome == "lab") {
			    #count.name = gsub('_mean',"_count",phenotype.list[i])
			    #count.name = gsub('_min',"_count",count.name)
			    #count.name = gsub('_max',"_count",count.name)
			    #count.tmp = count[,count.name]
			    #count.tmp = lab[,count.name]
			    #outcome.dat = data.frame(mvp001_id=mvp001.ID[count.tmp>0],y=lab[count.tmp > 0, i])
			    outcome.dat = data.frame(mvp001_id=mvp001.ID,y=lab[, i])
			  } else if (type == 'caseControl'){
			    outcome.dat = data.frame(mvp001.ID, y=code[,i])
			    outcome.dat = outcome[outcome$y!=1,] # filter out the ones
			    outcome.dat$y = ifelse(outcome$y >= 2, 1, 0)
			  }

			  
			  # inner join phenotype outcome with covariates
			  ID = intersect(outcome.dat$mvp001_id, demo$mvp001_id)
			  if(length(ID) > 0){
				  data = cbind(outcome.dat[match(ID, outcome.dat$mvp001_id),], demo[match(ID,demo$mvp001_id),])
				  data = na.omit(data)
				  mean.y = mean(data$y)
				  sd.y =  sd(data$y)
				  skew.y = skew(data$y)
				  num.levels.y = length(unique(data$y))

				    if (num.levels.y <= 2) {
				      data$y = ifelse(data$y > mean(data$y), 1, 0)   
				    } else if (skew.y > (-Inf) && all(data$y > 0.0)) {
				      data$y = log(data$y)
				    }
				 
				    if (outcome == "infection.code") { FAMILY = fam[i]} else {FAMILY = fam}
				    print(FAMILY)
				    model <- geeglm(formula, family = FAMILY, data = data,
				      corstr = "independence", id = 1:dim(data)[1])

			    # extract coefficients estimate and se from fitted model
				    if(!is.na(model$coefficients[snp.name])){
				      tmp = summary(model)$coefficients[snp.name,]
				      colnames(tmp)[dim(tmp)[2]] = "p.value"
				      res = rbind(res,c(phenotype = phenotype.list[i], 'white'=white,tmp,
				         sample.size = nrow(data),
					      mean.y = mean.y, sd.y = sd.y, skew.y = skew.y,
					      num.levels.y = num.levels.y))
			  	      }
			    	}
			    # cat(phenotype.list[i]," ",!is.na(model$coefficients[snp.name])," ",length(ID) > 0, '\n')
			   }
			
			#print(warnings())
			# write final table
			write.table(res,  file = output.file, sep = ",", quote = TRUE, col.names = TRUE, row.names = FALSE)

			}
	#}
}






quit(save='no')






