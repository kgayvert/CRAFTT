#########################################
# note - this code can only be run when the command line GSEA is installed on the user's system,
#    	which can be obtained from: http://www.broadinstitute.org/gsea/downloads.jsp

# User-specified path to the folder containing the command line GSEA and user directory
library(igraph)
library(data.table)

g=readRDS(paste0("network.rds")

myTF="CMYC"

# To identify best drug for myTF inhibition, repeat this all other drugs in dataset 
# (this code loads precalculated)
setwd(usr.dir)
myTF.ind<-which(mydrug.TF.list==myTF)
myTF.FWER<-c(mydrug.FWER[myTF.ind],readRDS(paste0(myTF,".FWER.rds")))
myTF.NES<-c(mydrug.NES[myTF.ind],readRDS(paste0(myTF,".NES.rds")))
myTF.ES<-c(mydrug.ES[myTF.ind],readRDS(paste0(myTF,".ES.rds")))
myTF.p<-c(mydrug.p[myTF.ind],readRDS(paste0(myTF,".p.rds")))
names(myTF.FWER)[1]<-names(myTF.NES)[1]<-names(myTF.ES)[1]<-names(myTF.p)[1]<-"mydrug"
druglist<-names(myTF.FWER)

##### STEP 2: Network Analysis #####
myTF.PL<-rep(NA,length(druglist))
names(myTF.PL)<-druglist
TFind=which(V(g)$name=="myTF"&V(g)$type=="TF")
for(i in 1:length(myTF.PL)){
  drugname = druglist[i]
  if(i%%100==0){ cat(paste0("Calculating path lengths for ",drugname," (",i,"/",length(druglist),")\n"))	}
  drug<-toupper(gsub("[^[:alnum:]]","",drugname))
  dind=which(V(g)$name==drug)
  if(length(dind)>0){
    PL[i]<-shortest.paths(g,dind,TFind,mode="out")		
  }
}

##### STEP 3: Prioritization #####
myTF.results<-data.table(Drug=druglist,PL=myTF.PL,NES=myTF.NES,FWER=myTF.FWER,ES=myTF.ES,p=myTF.p)
myTF.predictions<-myTF.results[myTF.results$FWER<0.1]
myTF.inhibition.predictions<-myTF.predictions[myTF.predictions$NES<0]; myTF.inhibition.predictions<-myTF.inhibition.predictions[order(myTF.inhibition.predictions$PL)]
myTF.reactivation.predictions<-myTF.predictions[myTF.predictions$NES>0]; myTF.reactivation.predictions<-myTF.reactivation.predictions[order(myTF.reactivation.predictions$PL)]
