#########################################
# note - this code can only be run when the command line GSEA is installed on the user's system,
#  		which can be obtained from: http://www.broadinstitute.org/gsea/downloads.jsp

# User-specified path to the folder containing the command line GSEA and user directory
args<-commandArgs(TRUE)
drug.dir<-args[1]
drugname<-gsub(".*/|.txt","",test)

geneset.name<-"ENCODE_ChIPseqTargets_Promoters.gmt"  	
GSEA.path<-"GSEA/"

usr.dir=getwd()

##### Data Acquistion and Prep #####
# Read in drug-induced differential expression data
X <- read.table(drug.dir)
expression <- as.matrix(X[,2],ncol=1)
order.ind<-order(expression,decreasing=FALSE)
dataset<- as.matrix(expression[order.ind,],ncol=1)
rownames(dataset) <- as.vector(X[,1])[order.ind]
drug.gct.name<-paste0(drugname,".gct")

# Convert the drug-induced expression data into gct format (if does not exist)
exist.dir<-list.files(pattern=paste0("^",drug.gct.name,"$"))
if(length(exist.dir)==0){
	col1=c("#1.2",length(dataset),"Gene",rownames(dataset))
	col2=c("",2,"DESCRIPTION",rep("na",length(dataset)))
	col3=c("","","DE",log2(dataset))
	col4=c("","","Original",rep(1,length(dataset)))
	table=cbind(col1,col2,col3,col4)
	write.table(table,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE,file=drug.gct.name)			
}


##### STEP 1: GSEA #####
# Run cmd line GSEA
system(paste0("java -Xmx2000m -cp ",GSEA.path,"build/gsea2.jar xtools.gsea.Gsea -res ",usr.dir,"/",drug.gct.name," -cls ",usr.dir,"cmap.cls -gmx ",usr.dir,geneset.name," -collapse false -permute gene_set -scoring_scheme classic -rpt_label ",drugname," -out ",usr.dir,"results/ -set_min 5 -set_max 1500 "))

# Read back in results
setwd("results")
newdir<-list.files(pattern=paste0("*",drugname,".Gsea*"))
num<-gsub(paste0(drugname,".Gsea."),"",newdir)
system(paste0("perl filter.report.pl ",drugname," ",num))
gsea.results<-read.table(paste0("gsea_report_",drugname,".txt"),sep="\t",stringsAsFactors=FALSE);colnames(gsea.results)<-gsea.results[1,];gsea.results<-gsea.results[2:dim(gsea.results)[1],]
mydrug.TF.list<-gsea.results[,1]
mydrug.ES<-gsea.results[,2]
mydrug.NES<-gsea.results[,3]
mydrug.p<-gsea.results[,4]
mydrug.FWER<-gsea.results[,6]

