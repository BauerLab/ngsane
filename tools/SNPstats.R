library(VariantAnnotation)

args<-commandArgs(TRUE)

cat("********** read in file")
vcf <- readVcf(args[1], args[3])
calls <- geno(vcf)$GT

cat("********** calculate the distance")
e = apply(calls,2,function(x)colSums(x!=calls))
write.csv(e,args[4])

cat("********** plot")
pdf(file=args[2], width=7,height=7)
heatmap(e)

# plots heatmap but without weired color wharp, so that now the diagonal
# should indeed be zero
heatmap(e,Colv=NA,Rowv=NA, scale="none")

# doublecheck with a simple matrix plot
library(lattice)
levelplot(e)

dev.off()

