
args <- commandArgs(trailingOnly = TRUE)
#annotation file name
annoF<-args[1]
#counts table per transcript from ht seq count
countsF<-args[2]
#sample name
sampleName<-args[3]
#annotation version identifier
annoversion<-args[4]

#read in gencode annotation
gencode.annotation.gtf<-read.table(file=annoF,sep="\t",stringsAsFactors=F,comment.char="#")
#give meaningful column names
colnames(gencode.annotation.gtf)<-c("seqname","source","feature","start","end","score","strand","frame","attribute")
#select all exons in annotation file to a dataframe
gencode.annotation.exon.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="exon"),]
#find the size of each exon
diffexon<-gencode.annotation.exon.gtf$end-gencode.annotation.exon.gtf$start
#find the name of the transcript each exon is associated with
names(diffexon)<-sapply(gencode.annotation.exon.gtf$attribute,function(x) gsub("transcript_id ","",strsplit(x,";")[[1]])[2])
#find the sum of the lengths of all exons associated with each transcript
exon_length_per_gene_id<-tapply(diffexon, names(diffexon), sum)
#pull out all transcripts
gencode.annotation.gene.gtf<-gencode.annotation.gtf[which(gencode.annotation.gtf$feature=="transcript"),]
#add per transcript exon length sums (and per kilobase)
gencode.annotation.gene.gtf.len<-cbind(gencode.annotation.gene.gtf,"exon_length_per_gene_id"=exon_length_per_gene_id,"exon_length_per_gene_id_per1000"=exon_length_per_gene_id/1000)
#read in per transcript counts
counts<-read.table(countsF,stringsAsFactors=F)
#calculate library size per million
lib_size_per_million<-sum(counts$V2)/1000000
# get all transcript ids
ids<-sapply(gencode.annotation.gene.gtf.len$attribute,function(x) gsub(" transcript_id ","",strsplit(x,";")[[1]])[2])
#combine tables by order of transcript ids
tab2<-cbind(counts,gencode.annotation.gene.gtf.len[match(counts$V1,ids),])
#calculate RPKMS for each transcript
RPKMs<-(tab2[,2]/as.vector(tab2[,"exon_length_per_gene_id_per1000"]))/lib_size_per_million
#make a table that has RPKMs and Transcript Name
RPKM<-cbind("RPKM"=RPKMs,"Transcript"=tab2$V1)
#write to file with sample name and annotation version in filename
write.csv(RPKM,file=paste(sampleName,"RPKMtranscript",annoversion,".csv",sep=""),quote=FALSE,row.names=FALSE)

