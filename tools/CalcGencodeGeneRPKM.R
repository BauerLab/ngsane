
args <- commandArgs(trailingOnly = TRUE)
#annotation (gencode) file eg gencode.v14.annotation.gtf
annoF<-args[1]
#counts per gene file from ht seq count eg gencode.v14.annotation.gene
countsF<-args[2]
#sample name eg RNA3kChr16
sampleName<-args[3]

library(GenomicFeatures)
library(edgeR)

#read in annotation file
gencode.annotation.gtf<-read.table(file=annoF,sep="\t",stringsAsFactors=F,comment.char="#")
#name columns
colnames(gencode.annotation.gtf)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
#subset exons
gencode.annotation.exon.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="exon"),]
#find gene_id for each exon
gene_id<-sapply(gencode.annotation.exon.gtf$attribute,function(x) gsub("gene_id ","",strsplit(x,";")[[1]])[1])
#create granges of all exons and name by gene_id
exons<-GRanges(seqnames=gencode.annotation.exon.gtf$seqname,ranges=IRanges(start=gencode.annotation.exon.gtf$start,end=gencode.annotation.exon.gtf$end,names=gene_id),strand=gencode.annotation.exon.gtf$strand)
#extract collapsed exon lengths per gene_id
gene.length.gencode<-sapply(split(exons, names(ranges(exons))), function(x) sum(width(reduce(x))))
#read in ht-seq-counts
counts<-read.table(countsF,stringsAsFactors=F)
#ensure gene lengths by gene_ids are same order as count table (should be identical)
ordered.gene.length<-gene.length.gencode[match(counts[,1],names(gene.length.gencode))]
#calculate rpkms using edgeR
RPKM<-rpkm(matrix(counts[,2]), gene.length=ordered.gene.length)
#write table of rpkms
write.csv(cbind("ENSG"=counts[,1],"RPKM"=RPKM[,1]),file=paste(sampleName,".RPKM.csv",sep=""),quote=FALSE,row.names=FALSE)

#### calculate RPKM per median & average gene length as well

#caution overwriting some of the above values in workspace.

#find transcript_id for each exon
transcript_id<-sapply(gencode.annotation.exon.gtf$attribute,function(x) gsub(" transcript_id ","",strsplit(x,";")[[1]])[2])
#create granges of all exons and name by gene_id
exons<-GRanges(seqnames=gencode.annotation.exon.gtf$seqname,ranges=IRanges(start=gencode.annotation.exon.gtf$start,end=gencode.annotation.exon.gtf$end,names=transcript_id),strand=gencode.annotation.exon.gtf$strand)
## "length of a gene was calculated as the median length of all transcripts relating to that gene" Oshlack, Biology Direct 2009
#find length of transcripts
transcript.length.gencode<-sum(width(split(exons, names(ranges(exons)))))
# find median transcript length per gene_id
# find gene_id for transcript_id
gencode.annotation.transcript.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="transcript"),]
gene_id<-sapply(gencode.annotation.transcript.gtf$attribute,function(x) gsub("gene_id ","",strsplit(x,";")[[1]])[1])
names(gene_id)<-sapply(gencode.annotation.transcript.gtf$attribute,function(x) gsub(" transcript_id ","",strsplit(x,";")[[1]])[2])
#rename list of transcript lengths with gene_id
names(transcript.length.gencode)<-gene_id[match(names(transcript.length.gencode),names(gene_id))]
#calculate median transcript length per gene_id
xx<-tapply(transcript.length.gencode,names(transcript.length.gencode),function(x) median(x))
#tidy up
gene.length.gencode<-as.vector(xx)
names(gene.length.gencode)<-names(xx)
#overwrite above values for median table
counts<-read.table(countsF,stringsAsFactors=F)
ordered.gene.length<-gene.length.gencode[match(counts[,1],names(gene.length.gencode))]
RPKM<-rpkm(matrix(counts[,2]), gene.length=ordered.gene.length)
write.csv(cbind("ENSG"=counts[,1],"RPKM"=RPKM[,1]),file=paste(sampleName,".RPKM.median.csv",sep=""),quote=FALSE,row.names=FALSE)

#calculate average transcript length per gene_id
xx<-tapply(transcript.length.gencode,names(transcript.length.gencode),function(x) mean(x))
#tidy up
gene.length.gencode<-as.vector(xx)
names(gene.length.gencode)<-names(xx)
#overwrite above values for median table
counts<-read.table(countsF,stringsAsFactors=F)
ordered.gene.length<-gene.length.gencode[match(counts[,1],names(gene.length.gencode))]
RPKM<-rpkm(matrix(counts[,2]), gene.length=ordered.gene.length)
write.csv(cbind("ENSG"=counts[,1],"RPKM"=RPKM[,1]),file=paste(sampleName,".RPKM.average.csv",sep=""),quote=FALSE,row.names=FALSE)

sink(type = "message")
sessionInfo()
