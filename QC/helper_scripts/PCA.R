#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

pdf("PCA.pdf")
par(mfrow=c(2,1),mar=c(4,4,2,2))
val <- read.table(args[1])
plot(c(seq(1,length(val$V1),by=1)),val$V1/sum(val$V1)*100,xlab="PC",ylab="Percent Variance Explained")
vec <- read.table(args[2])
plot(vec$V3,vec$V4,cex=0.5,xlab="PC 1",ylab = "PC 2")
dev.off()
