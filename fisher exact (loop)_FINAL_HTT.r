# Simulates different scenarios for finding functional variants in 0-n SC and in X/54 HapMaps (X=0)
# Fisher exact (or chi square) for 2x2 table
# TABLE LAYOUT
#
#       n(SC)| N(ctrls)
#     	col1 | col2
# row1 | x1  | x3	|	Negatives
# row2 | x2  | x4	|	Positives (ie functional variant)
#
# Table can be enlarged to rown (x in col2 of row1 will be x(n+1) etc
# rown | xn  | etc  |   Condition n

n=9 # total of population n (max nr of SCs OR alleles)
cutoff = 0.05 # could be 0.05
tests=1211 #nr of tests for bonferroni corr
bonf_cutoff =(cutoff/tests) # adjusted P value cutoff
bf<-prettyNum(bonf_cutoff,digits=2)
cat("Bonf cutoff for P =",cutoff,": ",bonf_cutoff) #bonf cutoff
plist=NULL

for (i in 0:n){

x1=(n-i)	# negatives n (SC)
x2=i		# positives n (SC)
x4= 		#positives N (default 0) (HM54)
x3=54-x4   	#positives N (default 0) (HM54)

table <- matrix(c(x1, x2, x3, x4), nr=2, #nr= nr of datarows
dimnames=list(c("negatives", "positives"), c("SCs", "Ctrls"))) #headers

#NB only use Fisher when one of the nrs in contingency table is small (ie one cell < 5) chisq.test(table) # dont use chi square when one of the nrs in table is low
ft<-c(fisher.test(table)) #vectorizes fisher test output
fisher.test(table)
#ft_pval_adj <- (ft$p.value*bonf) #bonf correction using p value from fisher
plist <- c(plist, ft$p.value)
print(table)
cat("\n")
#cat("positives: ",i) #positives
cat("P-val=",format(ft$p.value, digits=2), fill=FALSE) #Unadj. Pval

if (ft$p.value < bonf_cutoff) {
	cat(" --> Bonf significant!", fill=FALSE, "\n")
	} else {
	cat("(not bonf signif.)", "\n")
	}
}
cat("Bonf cutoff for P =",cutoff,": ",format(bonf_cutoff, digits=2)) #bonf cutoff
#print(format(plist, digits=2), quote=FALSE)
plot(0:n, log10(plist), type="b", col="red", xlab="Nr. of positive SC alleles", ylab="log10 of P value", 
     main=c("Nr. of positives in HM54:",x4,"Bonf=",bf))

#,segments(0,log10(bonf_cutoff),n,log10(bonf_cutoff), col=par("fg"), lty=2)) #segments is not working yet