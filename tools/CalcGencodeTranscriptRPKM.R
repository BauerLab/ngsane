
args <- commandArgs(trailingOnly = TRUE)
#annotation (gencode) file eg gencode.v14.annotation.gtf
annoF<-args[1]
#counts per transcript file from ht seq count eg gencode.v14.annotation.transcript
countsF<-args[2]
#sample name eg RNA3kChr16
sampleName<-args[3]

library(GenomicFeatures)
library(edgeR)

#as for per gene RPKM calc but replace with transcripts.

#read in annotation file
gencode.annotation.gtf<-read.table(file=annoF,sep="\t",stringsAsFactors=F,comment.char="#")
#name columns
colnames(gencode.annotation.gtf)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
#subset exons
gencode.annotation.exon.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="exon"),]

#find transcript_id for each exon
transcript_id<-sapply(gencode.annotation.exon.gtf$attribute,function(x) gsub(" transcript_id ","",strsplit(x,";")[[1]])[2])
#create granges of all exons and name by transcript_id
exons<-GRanges(seqnames=gencode.annotation.exon.gtf$seqname,ranges=IRanges(start=gencode.annotation.exon.gtf$start,end=gencode.annotation.exon.gtf$end,names=transcript_id),strand=gencode.annotation.exon.gtf$strand)
#find length of transcripts
transcript.length.gencode<-sum(width(split(exons, names(ranges(exons)))))
#read in ht-seq-counts
counts<-read.table(countsF,stringsAsFactors=F)
#ensure transcripts lengths by transcripts_ids are same order as count table (should be identical)
ordered.gene.length<-transcript.length.gencode[match(counts[,1],names(transcript.length.gencode))]
#calculate rpkms using edgeR
RPKM<-rpkm(matrix(counts[,2]), gene.length=transcript.length.gencode)
#write table of rpkms
write.csv(cbind("ENST"=counts[,1],"RPKM"=RPKM[,1]),file=paste(sampleName,".RPKM.csv",sep=""),quote=FALSE,row.names=FALSE)

sink(type = "message")
sessionInfo()



