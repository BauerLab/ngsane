

args <- commandArgs(trailingOnly = TRUE)

file_arg<-args[1] 
sample_arg<-args[2]
stranded_arg<-args[3]
path_arg<-args[4]
firststrand_arg<-args[5]
isPaired_arg<-args[6]

library(rtracklayer)
library(GenomicFeatures)

RNAbamTobw <- function(file, name, stranded=TRUE, firstStrand=TRUE, paired=TRUE,path) 
{
    require(rtracklayer)
    if (paired) {
				cat("[NOTE] Read gapped aligment paired\n")
                rs <- readGappedAlignmentPairs(file)
				        if (summary(rs)[1]<1) stop("[ERROR] Alignment file empty (is it really paired?)")
				        rs.count <- length(rs)/1e6
                rs.cov <- unlist(grglist(rs))
                if (stranded) {
                  					cat("[NOTE] Stranded\n")
                            rs.cov <- lapply(split(rs.cov, strand(rs.cov))[c("+", "-")], coverage)
                                if (firstStrand == TRUE) names(rs.cov) <- rev(names(rs.cov))
                                export(rs.cov[["+"]]/rs.count, BigWigFile(paste(path,"/",name, "_+.bw", sep="")))
                                export(rs.cov[["-"]]/rs.count, BigWigFile(paste(path,"/",name, "_-.bw", sep="")))
                                export(rs.cov[["-"]]/-rs.count, BigWigFile(paste(path,"/",name, "_--.bw", sep="")))
                              }
                else {
							cat("[NOTE] Unstranded\n")
							export(coverage(rs.cov)/rs.count, BigWigFile(paste(path,"/",name, ".bw", sep="")))
                }
	}    
    else {
				cat("[NOTE] Read gapped aligment unpaired\n")
                rs <- readGappedAlignments(file)
        				if (summary(rs)[1]<1)  stop("[ERROR] Alignment file empty")
                rs.count <- length(rs)/1e6
                rs.cov <- unlist(grglist(rs))
                if (stranded) {
                                cat("[NOTE] Stranded\n")
                                rs.cov <- lapply(split(rs.cov, strand(rs.cov))[c("+", "-")], coverage)
                                if (firstStrand == TRUE) names(rs.cov) <- rev(names(rs.cov))
                                export(rs.cov[["+"]]/rs.count, BigWigFile(paste(path,"/",name, "_+.bw", sep="")))
                                export(rs.cov[["-"]]/rs.count, BigWigFile(paste(path,"/",name, "_-.bw", sep="")))
                                export(rs.cov[["-"]]/-rs.count, BigWigFile(paste(path,"/",name, "_--.bw", sep="")))
                              } 
                                else {
                                  cat("[NOTE] Unstranded\n")
                                  export(coverage(rs.cov)/rs.count, BigWigFile(paste(path,"/",name, ".bw", sep="")))
                                }
                }
}

RNAbamTobw(file=file_arg,name=sample_arg,stranded=stranded_arg,firstStrand=firststrand_arg,paired=isPaired_arg, path=path_arg)

