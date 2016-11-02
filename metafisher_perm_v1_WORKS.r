## script calculates combined fisher p for permutations
library(survcomp)
START <- date() # for outputfile naming
stime <- gsub(":",";",(substr(START,12,19))) #starttime in outputfile
cutoff=0.05

for (i in 0:19){ #parts
  
  hlist<-read.csv(paste("A:/R_stuff/u27_metafisher/overview/Overview_H25_p",i,".csv", sep=""))
  hlistsub<-subset(hlist, select =2)
  colnames(hlistsub) <- ("RefSeq|SNP|U133")
  rlist<-read.csv(paste("A:/R_stuff/u27_metafisher/overview/Overview_R56_p",i,".csv", sep=""))
  rlistsub<-subset(rlist, select =2)
  colnames(rlistsub) <- ("RefSeq|SNP|U133")
  fisherp<-data.frame()

  for (k in 1:20) { #20 perm files of 50P each

    setwd("A:/R_stuff/test22_BioX2_H25_output")
    hfiles<-list.files(pattern=paste("permatrixND_part",i,"_50P",sep=""))
    h<-read.csv(hfiles[k])
    hsub<-subset(h, select = 10:59)
    #colnames(hsub) <- c(1:50)
    hcombi<-cbind(hlistsub,hsub) #nb hcombi is dataframe
    setwd("A:/R_stuff/u26_BioX2_R56_output")
    rfiles<-list.files(pattern=paste("R56_permatrixND_part",i,"_50P",sep=""))
    r<-read.csv(rfiles[k])
    rsub<-subset(r, select = 10:59)
    #colnames(rsub) <- c("pgeno_R56","bgeno_R56","bage_R56","page_R56" ,"brace_R56","prace_R56","bsex_R56","psex_R56")
    rcombi<-cbind(rlistsub,rsub) #nb hcombi is dataframe
         
    if (nrow(hlistsub) != nrow(h) || nrow(rlistsub) != nrow(r)){
      stop("WARNING: UNEQUAL ROW nrs!")
      } else { #print("rows are equal")
        m1<-merge(hcombi,rcombi,by.x="RefSeq|SNP|U133",by.y="RefSeq|SNP|U133") #joins both tables on column with probeset identifiers
        x<-nrow(m1)
        
        for (j in 1:x){ #rows
          for (l in 2:51){ # 50P cols
            n<-(l+50)
              if(is.na(m1[j,l])==F && is.na(m1[j,n])==F){
                if ((m1[j,l])<cutoff && (m1[j,n])<cutoff){
                                        
                  p<-c(m1[j,l],m1[j,n])
                  fcp<-combine.test(p=p, method="fisher")                          
                  fisherp<-rbind(fisherp,fcp)
              }
            }
          }
        }
      }
      print (paste("NOW PARSING Perm50P file nr. ",k,sep=""))
      print(date())
    }
              
  print(paste("NOW PRINTING: fisher_perm_p",i,".csv",sep=""))
  print(date())
  setwd("A:/R_stuff/u30_metafisher_perm")
  write.csv(fisherp, file=paste("fisher_perm_p",i,".csv",sep=""))
}

print(warnings())
print(START) #START
print(date()) #END
  