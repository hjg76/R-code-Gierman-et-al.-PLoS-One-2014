s1 <- read.delim(file="part.0_test10.snp.txt")
x <- subset(s1, select=1:56)
for( i in 1:
write.table(s1, file="s3.txt", row.names=FALSE, col.names=FALSE)