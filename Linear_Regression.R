options(repos="https://repo.bppt.go.id/cran/")

library(R6)
PD <- R6Class(public = list(
  dire=NULL,Migration=NULL,monthoutsample=NULL,
  res = NULL, Cutoff=NULL, month=NULL, inSample=NULL,outSample=NULL,
  TODR = NULL,MEF = NULL, Cluster = NULL, LMTable = NULL, code=NULL,TPCA=NULL,Reg=NULL,
  Tbacktest=NULL,  TStepwise=NULL,scaledTS=NULL,bestmodels=NULL,MulMatrix=NULL,wfReg=NULL,IAS39=NULL,
  initialize=function(dir="C:/Users/NR925AL/Documents/Mitsui/December 2019/PD",monthoutsample=12){
    options(java.parameters="- Xmx4000m")
    require(xlsx);require(lubridate); require(data.table)
    require(tidyverse); require(stringr); require(foreach)
    require(PCAmixdata);require(ClustOfVar)
    require(RcppEigen); require(caret); require(leaps); require(BB)
    require(setRNG);require(car);  require(lmtest)
    require(gridExtra);require(nortest)
    setwd(dir)
    self$monthoutsample<-monthoutsample
    self$dire=dir
    invisible(self)
  },
  ODR=function(x=j,filename="datas",dir=self$dire){
    
    
    ETFinish<-data.frame(fread(paste0(dir,"/Input/DataETFinish.csv")))
    ETFinish$CON_DATE<-as.Date(ETFinish$CON_DATE,"%Y-%m-%d")
    ETFinish<-ETFinish %>% filter(CON_DATE <="2019-12-31")
    ETFinish<-ETFinish[,c(1,24,43)]
    ETFinish[ETFinish$Tgl_Ditarik!="",'UPD_STAT']<-"R"
    CrossNo<-data.frame(fread(paste0(dir,"/Input/DataCrossNo_Dec2019.csv")))
    ETFinish[is.na(match(ETFinish$Contract,CrossNo$Contract))==FALSE,"UPD_STAT"]<-"R&R"
    
    Nominativ<-data.frame(fread(paste0(dir,"/Input/",filename,x,".csv")))
    Nominativ<-Nominativ[,c(1,3,7,27)]
    Nominativ[,1]<-as.Date(Nominativ[,1])
    Nominativ<-Nominativ[order(Nominativ$CONTRACT,Nominativ$ASOF),];rownames(Nominativ)<-seq_len(nrow(Nominativ))
    
    Nominativ[Nominativ$OneObligor==0,'mia']<-2
    Nominativ[Nominativ$OneObligor>0&Nominativ$OneObligor<=30,'mia']<-3
    Nominativ[Nominativ$OneObligor>30&Nominativ$OneObligor<=60,'mia']<-4
    Nominativ[Nominativ$OneObligor>60&Nominativ$OneObligor<=90,'mia']<-5
    Nominativ[Nominativ$OneObligor>90,'mia']<-6
    
    k<-length(Nominativ)
    for (i in 1:12){
      Nominativ[,k+i]<-lead(Nominativ$mia,i)
      a<-lead(Nominativ$CONTRACT,i)
      a[is.na(a)]<-0
      Nominativ[a!=Nominativ$CONTRACT,k+i]<-0
    }
    Nominativ$mia12<-Nominativ[,k+12]
    Nominativ[apply(Nominativ[,k+0:12],1,max)==6,'mia12']<-6
    Nominativ<-Nominativ[,-c(k+1:12)]
    Nominativ$UPD_STAT<-ETFinish[match(Nominativ[,"CONTRACT"],ETFinish$Contract),"UPD_STAT"]
    Nominativ[is.na(Nominativ$UPD_STAT),"UPD_STAT"]<-0
    Nominativ[Nominativ$mia12==0&(Nominativ$UPD_STAT=="F"|Nominativ$UPD_STAT=="T"),'mia12']<-1
    Nominativ[Nominativ$mia12==0&(Nominativ$UPD_STAT=="W"|Nominativ$UPD_STAT=="R&R"|Nominativ$UPD_STAT=="R"),'mia12']<-6
    Nominativ[Nominativ$mia12==0,"mia12"]<-1
    
    Migration<-Nominativ %>% group_by(ASOF,mia,mia12) %>% summarise(Outstanding=sum(PRINCIPALCONT))
    Migration<-data.frame(spread(Migration,key = mia12,value =Outstanding))
    Migration[is.na(Migration)]<-0
    Migration$Outstanding<-rowSums(Migration[,3:8])
    Migration<-Migration[,c(1,2,9,3:8)]
    
    ODR<-Migration[Migration$mia!=6,]
    ODR<-head(data.frame(ODR %>% group_by(ASOF) %>% summarize(Outstanding=sum(Outstanding),Default=sum(X6))),-12)
    
    ODR$ODR<-ODR$Default/ODR$Outstanding
    ODR$MAODR<-rowMeans(sapply(0:11,function(x)lag(ODR$ODR,x)))
    
    ODR$Logit<-log(ODR$MAODR/(1-ODR$MAODR))
    ODR<-ODR[,c(1,6)]
    
    ODR<-ODR[is.na(ODR[,2])==FALSE,]
    print(colnames(ODR)[2])
    colnames(ODR)<-c("ASOF","ODR")
    
    self$Migration<-Migration
    self$TODR<-ODR
    
    invisible(self)
  },
  MEFSelection=function(x=1,ODR=self$TODR,dir=self$dire){
    #Read MEF and MEF expected relationship
    
    MEF<-data.frame(fread(paste0(dir,"/Input/PredictedMEF.csv"),header =  TRUE))
    MEF$Month<-as.Date(MEF$Month,"%m/%d/%Y")
    MEFE<-fread(paste0(dir,"/Input/MEFER.csv"))
    #change MEF into numeric  
    MEF[,-1]<-apply(MEF[,-1],2,function(x)as.numeric(x))
    
    #remove mef with NA, and align MEF starting date with ODR
    MEF<-MEF[MEF$Month>=range(ODR$ASOF)[1],]
    MEF<-MEF[,sapply(MEF,function(x)anyNA(x)==FALSE)]
    MEF<-MEF[,sapply(MEF,function(x){any(is.infinite(x))==FALSE})]
    
    
    colnames(MEF)[nchar(colnames(MEF))<3]<-paste(colnames(MEF)[nchar(colnames(MEF))<3],'_',sep='')
    cond<-apply(MEF,2,function(x)(sum(is.na(x),na.rm = TRUE)/length(x)<0.01))
    MEF<-MEF[,cond]
    
    
    self$inSample<-tail(ODR$ASOF,1) %m-% months(self$monthoutsample)

    #Combine ODR and MEF to have correlation
    InSampleODR<-ODR[(ODR$ASOF<=self$inSample),]
    InSampleODR<-InSampleODR[order(InSampleODR$ASOF),]
    colnames(InSampleODR)[1]<-"Month"
    rownames(InSampleODR)<-as.character(c(1:length(InSampleODR$Month)))
    ODRMEF<-plyr::join(InSampleODR,MEF,by="Month")
    
    #Combine ODR and MEF to have correlation
    Correlation<-cor(ODRMEF[,2:length(ODRMEF)],method=c("pearson"))
    Correlation<-data.frame(Correlation)
    MEFCor<-Correlation[-c(1),]
    MEFid<-row.names(MEFCor)

    #Select k highest correlation MEF with odr per mefbase
    selection<-function(k=1,x){
      listdf<-list()
      listdt<-list()
      for(i in 1 :nrow(MEFE)){
        
        if(unlist(MEFE[i,4])=="Positive"){
          a<-MEFCor[substr(MEFid,1,nchar(MEFE[i,2]))==unlist(MEFE[i,2]),1:4]
          a<-a[a[,k]>=0.00001,]
          if(nrow(a)==0){
            listdf[[i]]<-listdt[[i]]<-NULL
          }else{
            listdt[[i]]<-row.names(a[order(a[,k],decreasing=TRUE
            )[1:min(x,nrow(a))],])
            listdf[[i]]<-row.names(a)
          }
          
        }  else if(unlist(MEFE[i,4])=="Negative"){
          a<-MEFCor[substr(MEFid,1,nchar(MEFE[i,2]))==unlist(MEFE[i,2]),1:4]
          a<-a[a[,k]<=(-0.00001),]
          if(nrow(a)==0){
            listdf[[i]]<-listdt[[i]]<-NULL
          }else{
            listdt[[i]]<-row.names(a[order(a[,k],decreasing=FALSE
            )[1:min(x,nrow(a))],])
            listdf[[i]]<-row.names(a)
          }
        }  else{
          a<-MEFCor[substr(MEFid,1,nchar(MEFE[i,2]))==unlist(MEFE[i,2]),1:4]
          a<-abs(a[abs(a[,1])>=0.0001,])
          if(nrow(a)==0){
            listdf[[i]]<-listdt[[i]]<-NULL
          }else{
            listdt[[i]]<-row.names(a[order(a[,k],decreasing=TRUE
            )[1:min(x,nrow(a))],])
            listdf[[i]]<-row.names(a)
          }
        }
      }
      listdf<-list(listdf,listdt)
      return(listdf)
    }
    MEFc<-unlist(selection(1,x)[[2]])
    MEFc<-MEFc[MEFc!="NA"]
    Correl<-data.frame("Correlation"=Correlation[,1])
    row.names(Correl)<-colnames(Correlation)
    
    #compile data to be put in environment
    listdf<-list(MEFc,MEF,Correl,data.frame(MEFE))
    
    self$MEF<-listdf
    invisible(self)
  },
  Clustering=function(data=self$MEF){
    
    #cluster MEF until all second eigen values is below 2
    
    #get MEF data
    SelectedMEF<-data[[2]][,match(data[[1]],colnames(data[[2]]))]
    a<-SelectedMEF
    
    #Perform correlation matrix
    cluster<-list()
    cormat<-cor(a)
    i=1
    
    #count eigen values
    max<-max(eigen(cormat)$values[-1])
    listdf<-list()
    
    #cluster and cut MEF into i number of cluster    
    clusvar<-hclustvar(a)
    cluster[[i]]<-cutreevar(clusvar,i,matsim=FALSE)$cluster
    listReturn<-list(cluster[[i]],unlist(listdf),cormat)
    
    #re-do cluster until all second eigen value is below 2
    while (max>1) {
      i=i+1
      cluster[[i]]<-cutreevar(clusvar,i,matsim=FALSE)$cluster
      listdf<-list()
      for(l in 1:i){
        if(length(data.frame(a[,which(cluster[[i]]==l)]))==1){0
        }else{listdf[[l]]<-eigen(cor(a[,which(cluster[[i]]==l)]))$values[-1]}
      }
      max<-max(unlist(listdf))
      print(max)
      listReturn<-list(data.frame(do.call("cbind",cluster)),listdf,cormat)
    }
  
      
    print("Success")
    self$Cluster<-listReturn
    invisible(self)
  },
  Stepwise=function(ODRate=self$TODR,MEF=self$MEF,maxcoefs=min(10,length(self$Cluster[[2]]))){
    #perform stepwise regression
    
    
    #combine ODR and all MEF
    SelectedODRate<-ODRate[order(ODRate$ASOF),]
    SelectedODRate<-SelectedODRate[,c(1,2)]
    colnames(SelectedODRate)[1]<-"Month"
    SelectedODRate<-SelectedODRate[SelectedODRate$Month<=self$inSample,]
    SelectedMEF<-MEF[[2]][,match(MEF[[1]],colnames(MEF[[2]]))]
    cond<-apply(SelectedMEF,2,function(x){anyNA(x)==FALSE})
    SelectedMEF$Month=NULL
    SelectedMEF$Month<-MEF[[2]]$Month
    SelectedMEF<-SelectedMEF[SelectedMEF$Month<=self$inSample,cond]
    ODRMEF<-plyr::join(SelectedMEF,SelectedODRate,by="Month")
    ODRMEF<-ODRMEF[,c(length(ODRMEF),1:(length(ODRMEF)-2))]
    colnames(ODRMEF)[1]<-"ODR"
    
    
    #Prepare data
    MEFCode<-data.frame(self$MEF[[4]])
    old.seed<-setRNG(list(kind="Mersenne-Twister",normal.kind="Inversion",seed=1234))
    listModel<-list()
    listResults<-list()
    listbestfactor<-list()
    
    #perform 20 loops of stepwise regression
    for(k in 1:20){
      x=1
      maxcoef<-maxcoefs
      print(k)
      train.control<-trainControl(method='cv',number=10)
      step.model<-train(ODR~.,data=ODRMEF,
                        method="leapSeq",
                        tuneGrid=data.frame(nvmax=1:maxcoef),
                        trControl=train.control
      )
      a<-names(coef(step.model$finalModel,step.model$bestTune[1,1]))
      
      #check double appearance of coefficient in choosen model
      for(i in seq_len(nrow(MEFCode))){
        if(length(a[substr(a,1,nchar(MEFCode[i,2]))==MEFCode[i,2]])>1){
          result=FALSE
          break
        }else{result=TRUE}
      }
      #if there are more than 1 model withe the same base MEF, redo stepwise with maxcoef-1
      while(result==FALSE){
        x=x+1  
        maxcoef=maxcoef-1
        print(paste(k,x,"maxcoef :",maxcoef))
        train.control<-trainControl(method='cv',number=10)
        step.model<-train(ODR~.,data=ODRMEF,
                          method="leapSeq",
                          tuneGrid=data.frame(nvmax=1:maxcoef),
                          trControl=train.control
        )
        a<-names(coef(step.model$finalModel,step.model$bestTune[1,1]))
        for(i in seq_len(nrow(MEFCode))){
          if(length(a[substr(a,1,nchar(MEFCode[i,2]))==MEFCode[i,2]])>1){
            result=FALSE
            break
          }else{result=TRUE}
        } 
      }
      
      
      
      #Get results
      listResults[[k]]<-step.model$results
      listModel[[k]]<-(coef(step.model$finalModel,step.model$bestTune[1,1]))
      print((coef(step.model$finalModel,step.model$bestTune[1,1])))
      listbestfactor[[k]]<-step.model$bestTune
      
    }
    
    #take out model coefficients
    comparison<-do.call('rbind',listbestfactor)
    modelnum<-unique((comparison))
    Comparison<-do.call('rbind',listResults)
    DTComparison<-data.table(Comparison)
    modelfactor<-sapply(listModel,function(x)length(x))
    models<-data.table(do.call('rbind',lapply(listModel,function(x)names(x))))
    model<-list()
    for(i in seq_len(length(modelnum[,1]))){
      model[[i]]<-unique(models[modelfactor==(modelnum[i,1]+1),][,2:(modelnum[i,1]+1)])
    }
    
    result<-list(model,listModel,listResults)
    self$TStepwise<-result
  },
  MAPE=function(a){
    #function to calculate insample MAPE
    return(mean(abs((a$fitted.values-a$model[,1])/a$model[,1])))
  },
  selectODR=function(ODRate=self$TODR,MEF=self$MEF){
    #function to get combined ODR and MEF
    
    SelectedODRate<-ODRate[order(ODRate$ASOF),]
    colnames(SelectedODRate)[1]<-"Month"
    SelectedMEF<-MEF[[2]]
    SelectedMEF$Month<-MEF[[2]]$Month
    SelectedMEF<-SelectedMEF[SelectedMEF$Month<=self$inSample,]
    ODRMEF<-plyr::join(SelectedMEF,SelectedODRate,by="Month")
    ODRMEF<-ODRMEF[,c(1,length(ODRMEF),2:(length(ODRMEF)-1))]
    colnames(ODRMEF)[2]="ODR"
    
    return(ODRMEF)
  },
  RegMod=function(ODRate=self$TODR,MEF=self$MEF,Clusteringv=self$Cluster,Stepwisev=self$TStepwise){
    
    #regression, test and order according to MAPE and test results
    
    #initiate numbering and lists
    m=1
    x=1
    lists<-list()
    listmodel<-list()
    
    
    #cut ODR into insample ODR
    
    ODRate<-ODRate[order(ODRate$ASOF),]
    row.names(ODRate)<-seq_len(length(ODRate$ASOF))
    selectedODR<-self$selectODR()
    
    #pick MEF for outsample prediction
    predictionODR<-data.table(MEF[[2]])[Month>self$inSample]
    
    #re-regress for stepwise models
    for(l in seq_len(length(Stepwisev[[1]]))){
      #take out coefficients
      stepwisemef<-Stepwisev[[1]][[l]]
      #regress
      RegMod<-lm(ODR~.,data=selectedODR[selectedODR$Month<=self$inSample,c('ODR',unlist(stepwisemef))])
      #predict
      prediction<-data.table(logPredict=predict.lm(RegMod,MEF[[2]]),
                             ASOF=MEF[[2]]$Month)[ASOF>self$inSample]
      
      prediction<-prediction[c(1:self$monthoutsample,nrow(prediction)),]
      prediction$predODR<-1/(1+exp(-prediction$logPredict))
      
      listmodel[[m]]<-RegMod
      
      #build data table from regression
      lists[[m]]<-data.frame(coefficient=names(RegMod$coefficients),
                             Beta=RegMod$coefficients,
                             ODR1Year=prediction[self$monthoutsample+1,predODR],
                             modelnum=m)
      prediction<-prediction[1:self$monthoutsample]
      prediction$logtest<-ODRate[match(prediction$ASOF,ODRate$ASOF),'ODR']
      prediction$test<-1/(1+exp(-prediction$logtest))

      
      #Add statistical test columns
      lists[[m]]$vif<-any(if(length(RegMod$coefficients)<=5){0}else{vif(RegMod)}>5)==FALSE
      lists[[m]]$pvalue<-any(summary(RegMod)$coef[,4]>0.05)==FALSE
      lists[[m]]$rsquare<-summary(RegMod)$adj.r.squared>=0.5
      lists[[m]]$ResNorm<-shapiro.test(RegMod$residuals)$p.value>=0.05
      lists[[m]]$htest<-bptest(RegMod)$p.value>=0.05 
      a<-MEF[[4]][match(substr(names(RegMod$coefficients),1,4),MEF[[4]][,2]),4]
      c<-lists[[m]]$Beta[!is.na(a)]
      b<-rep("TRUE",length(lists[[m]][!is.na(a),2]))
      a<-a[!is.na(a)]
      b[a=="Positive"]<-c[a=="Positive"]>0
      b[a=="Negative"]<-c[a=="Negative"]<0
      
      lists[[m]]$signTest<-sum(b==FALSE)==0
      lists[[m]]$lastscore<-sum(as.numeric(lists[[m]][1,c('vif','pvalue','rsquare','ResNorm','htest','signTest')]))/6 
      lists[[m]]$lastlogMAPE<-mean(abs((prediction$logPredict-prediction$logtest)/prediction$logtest))
      lists[[m]]$lastMAPE<-mean(abs((prediction$predODR-prediction$test)/prediction$test))
      lists[[m]]$vifscore<-if(length(RegMod$coefficients)<=2){0}else{c(0,vif(RegMod))}
      lists[[m]]$pvaluescore<-summary(RegMod)$coef[,4]
      lists[[m]]$rsquarescore<-summary(RegMod)$adj.r.squared
      lists[[m]]$ResNormscore<-shapiro.test(RegMod$residuals)$p.value
      lists[[m]]$htestscore<-bptest(RegMod)$p.value 
      lists[[m]]$expectedRelation<-c(NA,a)
      m=m+1
      x=x+1
    }
    MEFs<-Clusteringv[[1]]
    MEFs$MEF<-row.names(MEFs) 
    MEFComb<-list()
    
    #create combinations of all MEF
    #create clusters to be combined
    for(l in 2:(length(MEFs)-1)){
      listdf<-list()
      for(i in seq_len(l)){
        listdf[[i]]<-as.factor(MEFs[MEFs[,l]==i,'MEF'])
      }
      MEFComb[[l-1]]<-listdf
      
      
    }
    
    #Combine MEF
    if (length(MEFComb)>=3){
      MEFComb<-list(data.table(MEFs$MEF),CJ(MEFComb[[1]][[1]],MEFComb[[1]][[2]]),  
                    CJ(MEFComb[[2]][[1]],MEFComb[[2]][[2]],MEFComb[[2]][[3]]),
                    CJ(MEFComb[[3]][[1]],MEFComb[[3]][[2]],MEFComb[[3]][[3]],MEFComb[[3]][[4]]))
    }else if (length(MEFComb)==2){
      MEFComb<-list(data.table(MEFs$MEF),CJ(MEFComb[[1]][[1]],MEFComb[[1]][[2]]),  
                    CJ(MEFComb[[2]][[1]],MEFComb[[2]][[2]],MEFComb[[2]][[3]]))  
    }else if (length(MEFComb)==1){
      MEFComb<-list(data.table(MEFs$MEF),CJ(MEFComb[[1]][[1]],MEFComb[[1]][[2]]))
    }
    #Discard rows of MEF with combination containing same base MEF
    
    
    for(i in seq_len(length(MEFComb))){
      
      MEFComb[[i]]<-MEFComb[[i]][apply(MEFComb[[i]],1,function(x){any(table(substr(x,1,4))>1)})==FALSE]
      print(nrow(MEFComb[[i]]))
      if(nrow(MEFComb[[i]])==0)MEFComb[[i]]=NULL
    }
    
    
    #Regression 
    for(l in 1:length(MEFComb)){
      
      Comb<-MEFComb[[l]]
      if(l ==1){row<-nrow(Comb)}else{row<-row+nrow(Comb)}
      print(row)
      Comb[,names(Comb):=lapply(.SD,as.character)]
      start<-Sys.time()
      for(i in 1:nrow(MEFComb[[l]])){
        RegMod<-lm(ODR~.,data=selectedODR[selectedODR$Month<=self$inSample,c('ODR',unlist(Comb[i,]),recursive=TRUE)])
        prediction<-data.table(logPredict=predict.lm(RegMod,MEF[[2]]),
                               ASOF=MEF[[2]]$Month)[ASOF>self$inSample]
        
        prediction<-prediction[c(1:self$monthoutsample,nrow(prediction)),]
        prediction$predODR<-1/(1+exp(-prediction$logPredict))
        
        
        listmodel[[m]]<-RegMod
        #build data table from regression
        lists[[m]]<-data.frame(coefficient=names(RegMod$coefficients),
                               Beta=RegMod$coefficients,
                               ODR1Year=prediction[self$monthoutsample+1,predODR],
                               modelnum=m)
        prediction<-prediction[1:self$monthoutsample]
        prediction$logtest<-ODRate[match(prediction$ASOF,ODRate$ASOF),'ODR']
        prediction$test<-1/(1+exp(-prediction$logtest))
        
        
        #Add statistical test columns
        lists[[m]]$vif<-any(if(length(RegMod$coefficients)<=2){0}else{vif(RegMod)}>5)==FALSE
        lists[[m]]$pvalue<-any(summary(RegMod)$coef[,4]>0.05)==FALSE
        lists[[m]]$rsquare<-summary(RegMod)$adj.r.squared>=0.5
        lists[[m]]$ResNorm<-shapiro.test(RegMod$residuals)$p.value>=0.05
        lists[[m]]$htest<-bptest(RegMod)$p.value>=0.05 
        a<-MEF[[4]][match(substr(names(RegMod$coefficients),1,4),MEF[[4]][,2]),4]
        
        c<-lists[[m]]$Beta[!is.na(a)]
        b<-rep("TRUE",length(lists[[m]][!is.na(a),2]))
        a<-a[!is.na(a)]
        b[a=="Positive"]<-c[a=="Positive"]>0
        b[a=="Negative"]<-c[a=="Negative"]<0
        
        lists[[m]]$signTest<-sum(b==FALSE)==0
        lists[[m]]$lastscore<-sum(as.numeric(lists[[m]][1,c('vif','pvalue','rsquare','ResNorm','htest','signTest')]))/6
        lists[[m]]$lastlogMAPE<-mean(abs((prediction$logPredict-prediction$logtest)/prediction$logtest))
        lists[[m]]$lastMAPE<-mean(abs((prediction$predODR-prediction$test)/prediction$test))
        lists[[m]]$vifscore<-if(length(RegMod$coefficients)<=2){0}else{c(0,vif(RegMod))}
        lists[[m]]$pvaluescore<-summary(RegMod)$coef[,4]
        lists[[m]]$rsquarescore<-summary(RegMod)$adj.r.squared
        lists[[m]]$ResNormscore<-shapiro.test(RegMod$residuals)$p.value
        lists[[m]]$htestscore<-bptest(RegMod)$p.value 
        lists[[m]]$expectedRelation<-c(NA,a)
        x=x+1
        m=m+1
        if(x%%100==0){print(paste(m,Sys.time()));gc();object.size(lists)}
      }
      
      
      end<-Sys.time()
      print(end-start);gc()
    }
    
    #combine regression result into a table
    LmTable<-rbindlist(lists)
    LmTable<-LmTable[order(-LmTable$lastscore,LmTable$lastlogMAPE)]
    
    self$Reg<-listmodel
    self$LMTable<-LmTable
  },
  CV=function(lmtable=self$LMTable,regress=self$Reg){
    
    mape <- function(actual, predicted){
      actual=1/(1+exp(-actual));predicted=1/(1+exp(-predicted))
      mean(abs((actual - predicted)/actual))
    }
    mapeexpSummary <- function (data,
                                lev = NULL,
                                model = NULL) {
      out <- mape(data$obs, data$pred)  
      names(out) <- "MAPEEXP"
      out
    }
    
    
    modelnum<-unique(unlist(lmtable[lmtable$lastscore==1&lmtable$lastMAPE<0.1,'modelnum']))
    x<-1
    result<-list()
    for ( i in modelnum){
    passedmodel<-regress[[i]]
    train.control<-trainControl(method='cv',number=floor(nrow(passedmodel$model)/6),summaryFunction = mapeexpSummary)
    cvtest<-train(ODR~.,data=passedmodel$model,
                      method="lm",
                      metric = 'MAPEEXP',
                      trControl=train.control)
    result[[x]]<-cvtest$result
    if(x%%100==0){print(paste(x,Sys.time(),"Out of",length(modelnum)))}
    x=x+1
    }
  result<-do.call('rbind',result)
  result$modelnum=modelnum
  lmtable[,c('cvMAPE','sdMAPE')]<-result[match(lmtable$modelnum,result$modelnum),c('MAPEEXP','MAPEEXPSD')]
  self$LMTable<-lmtable
  gc(gc)},
  run=function(j=1,dir=self$dire){
    for(x in j){
      #wrapper function to run for each segment(25-36)
      print(x)
      self$ODR(x)
      self$MEFSelection(3)
      self$Clustering()
      self$Stepwise()
      self$RegMod()
      self$CV()
      fwrite(self$LMTable,paste(dir,"/Output/lmtable",x,".csv",sep=""))
      fwrite(self$Migration,paste(dir,"/Output/Migration",x,".csv",sep=""))
      gc(gc)
    }
  }
))
self<-PD$new()
self$run(c(9))
