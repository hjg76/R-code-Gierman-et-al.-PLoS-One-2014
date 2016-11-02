#removes medulla from kidney data
setwd("A:/R_stuff/test9_All_eQTL_106xVAR_Max8NA/R code and input files")

n=19 #number of part files (19)

for (i in 17:n){

f1 <- read.delim(file=paste("part.",i,".2log.txt",sep=""), header=FALSE) #store table in f1

s1 <- subset(f1, select=1:56) #store first 56 columns into s1

write.table(s1, file=paste("part.",i,".2log_rdw_crtx_max8NA.txt",sep=""), row.names=FALSE, col.names=FALSE)
}
