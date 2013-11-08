######################################################
#
# Visualising Cuffdiff output using cummeRbund: some
# functions and their usage for doing so.
#
######################################################
library("cummeRbund")

args<-commandArgs(TRUE)

# PolyATail filtered (read: high quality for protein coding genes, but only subset of individuals)
setwd('/datastore/cmis/bau04c/Documents/datahome/catomocolo/RNAseqPolyA/differential/cuffdiffStratObeseVsLean')
# Whole RNA data (read: captuerd everything that is expressed, e.g. ncRNA, for all individuals)
# setwd('/datastore/cmis/bau04c/Documents/datahome/catomocolo/RNAseq/differential/cuffdiffStratObeseVsLean')

# From commandline
# setwd(args[1])

cat('********** read in the data,  by changing the working directory') 
cuff <- readCufflinks()
cuff

pdf(file=args[2], width=10,height=7)

##########
# Diagnostic plots
#########

# dispersion plot
dispersionPlot(genes(cuff))

# squared coefficient of variation 
fpkmSCVPlot(genes(cuff))
fpkmSCVPlot(isoforms(cuff))

# produce pairwise scatterplots of FPKM for each sample
csScatterMatrix(genes(cuff),replicates=T)
csScatterMatrix(genes(cuff))

# produce density plots of the FPKMs for each sample
csDensity(genes(cuff),replicates=T)
csDensity(genes(cuff))

# boxplots of data for each condition
csBoxplot(genes(cuff),replicates=T)
csBoxplot(genes(cuff))

# cluster dendrogram
csDendro(genes(cuff),replicates=T)
csDendro(genes(cuff))

# PCA-plot
PCAplot(genes(cuff),"PC1","PC2")

###########
# Significant genes
###########

# produces a list of significant genes for at least 1 comparison
sigGeneIds <- getSig(cuff,alpha=0.05,level="genes")
length(sigGeneIds)

# a visual matrix of the number of significant genes for each 
# comparison of conditions
sigMatrix(cuff,alpha=0.05,level="genes",orderByDist=F)

# volcano
csVolcanoMatrix(genes(cuff))

############
# My genes
# @Rob These are the reserach questions:
#
# INFLAMATORY
# Are inflamatory genes more expressed in patient 14, 40, 57,50,(59), who have human 
# material in their digesta e.g. from bleeding or gut wall cell shedding. 
# Note some of the individuals are only in whole RNA.
#
# ADIPOSE
# Are adipose samples really from adipose tissue, i.e. do adipose samples 
# labelled as (<ID>ar, eg. 14ar) really express the adipose marker genes
# This is interesting because Paul found a large number of microbial material 
# in the adipose samples: which means a) Paul's analysis is wrong, b) adipose samples
# are acutally something else (e.g. gut), c) microbial material invades as deep as 
# adipose tissue 
#
# CALCIUM 
# These are the genes Desma is interrested in
#
#
############

#Myges
myIDsInf <- read.delim("/datastore/cmis/bau04c/Documents/datahome/catomocolo/doc/genesOfInterestInflammation/HGNCtoRefSeq.txt")
myGenesInf <- getGenes(cuff,t(myIDsInf[2]))
myIDsAdip <- read.delim("/datastore/cmis/bau04c/Documents/datahome/catomocolo/doc/genesOfInterestAdipokine/HGNCtoRefSeq.txt")
myGenesAdip <- getGenes(cuff,t(myIDsAdip[2]))
myIDsCalc <- read.delim("/datastore/cmis/bau04c/Documents/datahome/catomocolo/doc/genesOfInterest/HGNCtoRefSeq.txt")
myGenesCalc <- getGenes(cuff,t(myIdsCalc[2]))

# semi useless plotes because there are so many genes
expressionBarplot(myGenesInf)
expressionBarplot(myGenesAdip)
expressionBarplot(myGenesCalc)

#heatmap(fpkm(myGenesInf))

#matrix with filename -> replicates
#x<-(cbind(repFpkm(myGenesInf)[1],repFpkm(myGenesInf)[4] ,repFpkm(myGenesInf)[8]))
#qplot(rep_name,gene_id, data=x, fill=log(fpkm), geom="raster")
#ggplot(x, aes(factor(rep_name), factor(gene_id), fill = fpkm) + geom_raster())


# @Rob
# make matrix for genes of interest to plot as a heatmap
# here I was hoping that the 'replicates' which are actually the different individuals
# and tissue would nicely cluster and stand out somehow
gn <- unique(repFpkm(myGenesInf)[1])
sn <- unique(repFpkm(myGenesInf)[4])
fn <- replicates(cuff)[1] # I think this is always in the same order as sn
matInf<-matrix(t(repFpkm(myGenesInf)[8]),dim(gn)[1],dim(sn)[1],byrow=TRUE,dimnames=list(t(gn),t(fn)))
heatmap(matInf)#, Colv=NA)




#############
# Cluster genes
# @Rob not too sure this is actually usefull, clustering genes by
# their "expression" profile, pfff...
#############

#cluster genes of similar expression 
ic<-csCluster(genes(cuff),k=16)
csClusterPlot(ic)

csClusterPlot(csCluster(myGenesAdip,k=4))
csClusterPlot(csCluster(myGenesInf,k=16))
csClusterPlot(csCluster(myGenesCalc,k=16))


############# 
# My gene
# @Rob But this is pretty neat: finding genes that have a similar
# profile than a gene of interest (here CRNDE)
#############
mygene=getGene(cuff,"XLOC_016380") # CRNDE locus 1
mygene2=getGene(cuff,"XLOC_016381") # CRNDE Locus 2 

expressionPlot(isoforms(mygene),replicates=T)
expressionBarplot(isoforms(mygene),replicates=T)

#NR_034106     NR_034105

mySimilar<-findSimilar(cuff,"XLOC_016380",n=20)
expressionPlot(mySimilar,logMode=T,showErrorbars=F)


sink(type = "message")
sessionInfo()

