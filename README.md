# CRAFTT
#### Computational drug-Repositioning Approach For Targeting Transcription factors


- Last Updated: 1/12/2015
- Correspondence to:  Katie Gayvert, kmg257 [at] cornell [dot] edu

##################
### Requirements #
##################
- R - tested on version  **3.2.2** (2015-08-14) -- **"Fire Safety"**
- R dependecies: data.table, igraph
- Other dependencies: command line installation of GSEA, available from http://www.broadinstitute.org/gsea/downloads.jsp

###########################################
### To Run Prediction Step (command line) #
###########################################
- Note: to run from command line, the geneset and GSEA program must be located in the current directory or the script must be updated to set the correct directories

R --vanilla sample_prediction_step.R /path/to/drug/pertubation/file
