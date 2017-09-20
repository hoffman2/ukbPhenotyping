####################################################################################################
#Counts ICD10
####################################################################################################
#Read in 3 character set(This is the list of codes from column "CODE_ONLY" 3D)
regen3Char <- fread("/GWD/appbase/projects/RD-TSci-PhewasUKB/josh/regeneronConcordance/3charCodes.ls",header=F)
#Read in 4 character set (This is the list of codes from column "CODE_ONLY" 4D)
regen4Char <- fread("/GWD/appbase/projects/RD-TSci-PhewasUKB/josh/regeneronConcordance/4charCodes.ls",header=F)
#Read in bulk primary and secondary HES codes
hesDataPrim <- fread("/GWD/appbase/projects/RD-TSci-PhewasUKB/josh/phenotypes/HESIN_26041_all.minus4226129.tsv",header=T,sep="\t",na.strings="")
hesDataSec <- fread("/GWD/appbase/projects/RD-TSci-PhewasUKB/josh/phenotypes/HESIN_SECONDARY_DIAG10_26041_all.minus508567.tsv",header=T,sep="\t",na.strings="")
#Convert diagnosis field name to be consistent across primary and secondary
names(hesDataPrim)[5]="diag"
names(hesDataSec)[4]="diag"
#Add a flag for whether prim/sec diagnosis
hesDataPrim[,type:="primary"]
hesDataSec[,type:="secondary"]
####################################################################################################
#Tabulate 3 digit codes: Stratified by primary/secondary counts
hesDataCombined<-rbind(hesDataPrim[,c(1,5,21)],hesDataSec[,c(1,4,6)])
hesDataCombined[,diag:=paste0("p_",substr(hesDataCombined$diag,1,3))] #Changed to 3 characters
#Remove duplicated HES codes based on eid,code, and whether primary or secondary
hesDataCombined <- hesDataCombined[!duplicated(hesDataCombined,by=c("eid","diag","type"))]
#Tabulate case counts based on the list of ICD10 codes read in above
{
combinedCount=NULL
finalCount=NULL
myCount=NULL
for(phenoCode in regen3Char$V1) {
	myCount=hesDataCombined[diag %like% phenoCode,][,.N,by=.(diag,type)]
	combinedCount=rbind(combinedCount,myCount)}
	finalCount=cbind(combinedCount[type=="primary",],combinedCount[type=="secondary",])
	print(finalCount)}
####################################################################################################
#3 digit codes: Combined primary/secondary counts
hesDataCombined<-rbind(hesDataPrim[,c(1,5,21)],hesDataSec[,c(1,4,6)])
hesDataCombined[,diag:=paste0("p_",substr(hesDataCombined$diag,1,3))] #Changed to 3 characters
#Remove duplicated HES codes based on eid,code, and whether primary or secondary
hesDataCombined <- hesDataCombined[!duplicated(hesDataCombined,by=c("eid","diag"))]
{
combinedCount=NULL
finalCount=NULL
myCount=NULL
for(phenoCode in regen3Char$V1) {
	myCount=hesDataCombined[diag %like% phenoCode,][,.N,by=.(diag)]
	combinedCount=rbind(combinedCount,myCount)}
	print(combinedCount)}
####################################################################################################
#Tabulate 4 digit codes: Stratified by primary/secondary counts
hesDataCombined<-rbind(hesDataPrim[,c(1,5,21)],hesDataSec[,c(1,4,6)])
hesDataCombined[,diag:=paste0("p_",substr(hesDataCombined$diag,1,4))] #Changed to 4 characters
#For generating startified primary secondary counts
hesDataCombined <- hesDataCombined[!duplicated(hesDataCombined,by=c("eid","diag","type"))]
#Tabulate case counts based on the list of ICD10 codes read in above
{
combinedCount=NULL
finalCount=NULL
myCount=NULL
for(phenoCode in regen4Char$V1) {
	myCount=hesDataCombined[diag %like% phenoCode,][,.N,by=.(diag,type)]
	combinedCount=rbind(combinedCount,myCount)}
	finalCount=cbind(combinedCount[type=="primary",],combinedCount[type=="secondary",])
	print(finalCount)}
####################################################################################################
#Tabulate 4 digit codes: Combined primary/secondary counts
hesDataCombined<-rbind(hesDataPrim[,c(1,5,21)],hesDataSec[,c(1,4,6)])
hesDataCombined[,diag:=paste0("p_",substr(hesDataCombined$diag,1,4))] #Changed to 4 characters
hesDataCombined <- hesDataCombined[!duplicated(hesDataCombined,by=c("eid","diag"))]
{
combinedCount=NULL
finalCount=NULL
myCount=NULL
for(phenoCode in regen4Char$V1) {
	myCount=hesDataCombined[diag %like% phenoCode,][,.N,by=.(diag)]
	combinedCount=rbind(combinedCount,myCount)}
	print(combinedCount)}
####################################################################################################
