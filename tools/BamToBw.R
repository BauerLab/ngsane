

args <- commandArgs(trailingOnly = TRUE)

file_arg<-args[1]
sample_arg<-args[2]
stranded_arg<-args[3]
firststrand_arg<-args[4]
path<-args[5]

library(rtracklayer)
library(GenomicFeatures)


RNAbamTobw <- function(file, name, stranded=TRUE, firstStrand=TRUE,path) 
{
    

                rs <- readGappedAlignmentPairs(file)
                rs.count <- length(rs)/1e6
                rs.cov <- unlist(grglist(rs))
                if (stranded) {
                                rs.cov <- lapply(split(rs.cov, strand(rs.cov))[c("+", "-")], coverage)
                                if (!firstStrand) names(rs.cov) <- rev(names(rs.cov))
                                export(rs.cov[["+"]]/rs.count, BigWigFile(paste(path,name, "_+.bw", sep="")))
                                export(rs.cov[["-"]]/rs.count, BigWigFile(paste(path,name, "_-.bw", sep="")))
                                export(rs.cov[["-"]]/-rs.count, BigWigFile(paste(path,name, "_--.bw", sep="")))
                              } 
                                else export(coverage(rs.cov)/rs.count, BigWigFile(paste(path,name, ".bw", sep="")))


}


RNAbamTobw(file_arg,sample_arg,stranded_arg,firststrand_arg,path)

