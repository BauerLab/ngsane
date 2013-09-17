library(Repitools)
library(ggplot2)
library(splines)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
chip<-args[1]
input<-args[2]
outputpath<-args[3]

ChIPQC <- function(rsChIP, rsINPUT, chrSizes, ipname, cname, outputpath, windowSize=1000, dataPoints=1000)
{	

	gn.windows <- genomeBlocks(chrSizes, chrs=seq(length(chrSizes)), windowSize)
	ChIPcounts <- countOverlaps(gn.windows, rsChIP)
	INPUTcounts <- countOverlaps(gn.windows, rsINPUT)
	
	gn.ordered = data.frame(cbind(ChIPcounts, INPUTcounts)[order(ChIPcounts),])
	
	gn.ordered$ChIPcounts = cumsum(gn.ordered$ChIPcounts)
	gn.ordered$INPUTcounts = cumsum(gn.ordered$INPUTcounts)
	
	
	gn.ordered$ChIPcounts = gn.ordered$ChIPcounts/gn.ordered$ChIPcounts[length(gn.windows)]
	gn.ordered$INPUTcounts = gn.ordered$INPUTcounts/gn.ordered$INPUTcounts[length(gn.windows)]	

	spaced <- c(round(seq(1,length(gn.windows), by=(length(gn.windows)-1)/dataPoints)))
	df1 <- data.frame("bin"=c(1: length(gn.ordered$ChIPcounts[spaced])),"value"=gn.ordered$ChIPcounts[spaced], "type"= ipname)
	df2 <- data.frame("bin"=c(1: length(gn.ordered$INPUTcounts[spaced])),"value"=gn.ordered$INPUTcounts[spaced], "type"= cname)
	df <- data.frame(rbind(df1,df2))
	maxdist <- which.max(abs(gn.ordered$ChIPcounts[spaced]-gn.ordered$INPUTcounts[spaced]))
	
	p1 <- ggplot(df, aes(x=bin, y=value, group=type)) 
	p1 <- p1 + geom_vline(xintercept = maxdist, colour="grey")
	p1 <- p1 + geom_line(aes(color=type))
	p1 <- p1 + xlab("Percentage of bins") + ylab("Percentage of reads") 
	p1 <- p1 + theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.text=element_text(size=8), legend.title=element_blank(), legend.background = element_rect(fill = "white", colour = NA))	
	p1 <- p1 + scale_x_continuous(breaks=seq(0,dataPoints, by=(dataPoints)/4), labels=c(0,25,50,75,100))
	
	# print plots
	pdf(paste(outputpath,"/", ipname, ".pdf", sep=""), width=4, height=4)
    p1
	dev.off()

    png(paste(outputpath,"/", ipname, ".png", sep=""), width=200, height=200, units = "px")
    p1
    dev.off()
}

rsChIP <- BAM2GRanges(chip)
rsINPUT <- BAM2GRanges(input)
#only count read starts, avoids reads being across multiple windows
rsChIP <- resize(rsChIP, 1, fix="start")
rsINPUT <- resize(rsINPUT, 1, fix="start")

# get the genome sizes form the bam header
bh <- scanBamHeader(chip)
chrSizes <- (bh[[names(bh)]][["targets"]])
    
ChIPQC(rsChIP,rsINPUT, chrSizes , basename(file_path_sans_ext(chip)), basename(file_path_sans_ext(input)), outputpath)
