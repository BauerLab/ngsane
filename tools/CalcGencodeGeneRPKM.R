
args <- commandArgs(trailingOnly = TRUE)
#annotation file name
annoF<-args[1]
#counts table per transcript from ht seq count
countsF<-args[2]
#sample name
sampleName<-args[3]
#annotation version identifier
annoversion<-args[4]
require(GenomicFeatures)
#read in gencode annotation
gencode.annotation.gtf<-read.table(file=annoF,sep="\t",stringsAsFactors=F,comment.char="#")
#give meaningful column names
colnames(gencode.annotation.gtf)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
#select all exons in annotation file to a dataframe
gencode.annotation.exon.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="exon"),]
 #find the name of the gene id  each exon is associated with
gene_id<-sapply(gencode.annotation.exon.gtf$attribute,function(x) gsub("gene_id ","",strsplit(x,";")[[1]])[1])
exons<-GRanges(seqnames=gencode.annotation.exon.gtf$seqname,ranges=IRanges(start=gencode.annotation.exon.gtf$start,end=gencode.annotation.exon.gtf$end,names=gene_id),strand=gencode.annotation.exon.gtf$strand)
#pull out all genes
gencode.annotation.gene.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="gene"),]
gencode.annotation.genes<-sapply(gencode.annotation.gene.gtf$attribute,function(x) gsub("gene_id ","",strsplit(x,";")[[1]])[1])
#find the sum of the lengths of all exons associated with each gene
##exon_length_per_gene_id<-sapply(unique(exons@ranges@NAMES),function(x) sum(width(reduce(exons[which(exons@ranges@NAMES==x)]))))
#speed up without lookup
exon_length_per_gene_id<-sapply(split(exons, names(ranges(exons))), function(x) sum(width(reduce(x))))
#reorder to match up with table
exon_length_per_gene_id<-exon_length_per_gene_id[match(gencode.annotation.genes,names(exon_length_per_gene_id))]
#add per gene exon length sums (and per kilobase)
gencode.annotation.gene.gtf.len<-cbind(gencode.annotation.gene.gtf,"exon_length_per_gene_id"=exon_length_per_gene_id,"exon_length_per_gene_id_per1000"=exon_length_per_gene_id/1000)
#read in per transcript counts
counts<-read.table(countsF,stringsAsFactors=F)
#calculate library size per million
lib_size_per_million<-sum(counts$V2)/1000000
# get all gene ids
ids<-sapply(gencode.annotation.gene.gtf.len$attribute,function(x) gsub("gene_id ","",strsplit(x,";")[[1]])[1])
#combine tables by order of gene ids
tab2<-cbind(counts,gencode.annotation.gene.gtf.len[match(counts$V1,ids),])
#calculate RPKMS for each gene
RPKMs<-(tab2[,2]/as.vector(tab2[,"exon_length_per_gene_id_per1000"]))/lib_size_per_million
#make a table that has RPKMs and Gene Name
RPKM<-cbind("RPKM"=RPKMs,"Gene"=tab2$V1)
#write to file with sample name and annotation version in filename
write.csv(RPKM,file=paste(sampleName,"RPKM",annoversion,".csv",sep=""),quote=FALSE,row.names=FALSE)

