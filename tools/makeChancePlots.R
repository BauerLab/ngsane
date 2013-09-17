library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(splines)
library(gridExtra)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
chip<-args[1]
input<-args[2]
outputpath<-args[3]

ChIPQC <- function(rsChip, rsInput, ipname, cname, outputpath, windowSize=1000, dataPoints=1000)
{	
	hg19.windows <- genomeBlocks(Hsapiens, chrs=1:24, windowSize)
	ChIPcounts <- countOverlaps(hg19.windows, rsChIP)
	INPUTcounts <- countOverlaps(hg19.windows, rsINPUT)
	
	hg19.ordered = data.frame(cbind(ChIPcounts, INPUTcounts)[order(ChIPcounts),])
	
	hg19.ordered$ChIPcounts = cumsum(hg19.ordered$ChIPcounts)
	hg19.ordered$INPUTcounts = cumsum(hg19.ordered$INPUTcounts)
	
	
	hg19.ordered$ChIPcounts = hg19.ordered$ChIPcounts/hg19.ordered$ChIPcounts[length(hg19.windows)]
	hg19.ordered$INPUTcounts = hg19.ordered$INPUTcounts/hg19.ordered$INPUTcounts[length(hg19.windows)]	

	spaced <- c(round(seq(1,length(hg19.windows), by=(length(hg19.windows)-1)/dataPoints)))
	df1 <- data.frame("bin"=c(1: length(hg19.ordered$ChIPcounts[spaced])),"value"=hg19.ordered$ChIPcounts[spaced], "type"= ipname)
	df2 <- data.frame("bin"=c(1: length(hg19.ordered$INPUTcounts[spaced])),"value"=hg19.ordered$INPUTcounts[spaced], "type"= cname)
	df <- data.frame(rbind(df1,df2))
	maxdist <- which.max(abs(hg19.ordered$ChIPcounts[spaced]-hg19.ordered$INPUTcounts[spaced]))
	
	p1 <- ggplot(df, aes(x=bin, y=value, group=type)) 
	p1 <- p1 + geom_vline(xintercept = maxdist, colour="grey")
	p1 <- p1 + geom_line(aes(color=type))
	p1 <- p1 + xlab("Percentage of bins") + ylab("Percentage of reads") 
	p1 <- p1 + theme(legend.position = c(0, 1), legend.justification = c(0, 1), legend.text=element_text(size=8), legend.title=element_blank(), legend.background = element_rect(fill = "white", colour = NA))	
	p1 <- p1 + scale_x_continuous(breaks=seq(0,dataPoints, by=(dataPoints)/4), labels=c(0,25,50,75,100))
	
	
	# get derivative
	x1 = seq(1, dataPoints-2)
	x2 = seq(3, dataPoints)
	slope1 = 1/2*(hg19.ordered$ChIPcounts[spaced[x2]]-hg19.ordered$ChIPcounts[spaced[x1]])
	slope2 = 1/2*(hg19.ordered$INPUTcounts[spaced[x2]]-hg19.ordered$INPUTcounts[spaced[x1]])
	
	df3 <- data.frame(bin=c(1: length(slope1)), value=slope1, type="ChIP")
	df4 <- data.frame(bin=c(1: length(slope1)), value=slope2, type="INPUT")
	df5 <- data.frame(rbind(df3,df4))
	
	p2 <- ggplot(df5, aes(x=bin, y=value, group=type)) + xlab("Percentage of bins") + ylab("dy/dx")
	#p2 <- p2 + geom_smooth(method = "lm",formula = y~bs(x, degree = 5),se = TRUE, alpha=0.5)
	p2 <- p2 + geom_vline(xintercept = maxdist, colour="grey")
	p2 <- p2 + geom_line(aes(color=type))
	p2 <- p2 + theme(legend.position="none")
	p2 <- p2 + scale_x_continuous(breaks=seq(0,dataPoints, by=(dataPoints)/4), labels=c(0,25,50,75,100))
	
	# print plots
	pdf(paste(outputpath,"/", ipname, ".pdf", sep=""), width=8, height=4)
	grid.arrange(p1, p2 , ncol=2)
	dev.off()

        png(paste(outputpath,"/", ipname, ".png", sep=""), width=800, height=400, units = "px")
        grid.arrange(p1, p2 , ncol=2)
        dev.off()
}

rsChIP <- BAM2GRanges(chip)
rsINPUT <- BAM2GRanges(input)
#only count read starts, avoids reads being across multiple windows
rsChIP <- resize(rsChIP, 1, fix="start")
rsINPUT <- resize(rsINPUT, 1, fix="start")

ChIPQC(rsChip,rsInput, basename(file_path_sans_ext(chip)), basename(file_path_sans_ext(input)), outputpath)

