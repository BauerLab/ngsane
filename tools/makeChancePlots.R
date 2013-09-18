library(Repitools)
library(ggplot2)
library(splines)
library(tools)
library(scales)

args <- commandArgs(trailingOnly = TRUE)
chip<-args[1]
input<-args[2]
chipname <-args[3]
inputname<-args[4]
outputpath<-args[5]

windowSize=1000
dataPoints=1000

rsChIP <- BAM2GRanges(chip)
rsINPUT <- BAM2GRanges(input)
#only count read starts, avoids reads being across multiple windows
rsChIP <- resize(rsChIP, 1, fix="start")
rsINPUT <- resize(rsINPUT, 1, fix="start")

# get the genome sizes form the bam header
bh <- scanBamHeader(chip)
chrSizes <- (bh[[names(bh)]][["targets"]])
    
gn.windows <- genomeBlocks(chrSizes, chrs=seq(length(chrSizes)), windowSize)
ChIPcounts <- countOverlaps(gn.windows, rsChIP)
INPUTcounts <- countOverlaps(gn.windows, rsINPUT)

gn.ordered = data.frame(cbind(ChIPcounts, INPUTcounts)[order(ChIPcounts),])

gn.ordered$ChIPcounts = cumsum(gn.ordered$ChIPcounts)
gn.ordered$INPUTcounts = cumsum(gn.ordered$INPUTcounts)


gn.ordered$ChIPcounts = gn.ordered$ChIPcounts/gn.ordered$ChIPcounts[length(gn.windows)]
gn.ordered$INPUTcounts = gn.ordered$INPUTcounts/gn.ordered$INPUTcounts[length(gn.windows)]	

spaced <- c(round(seq(1,length(gn.windows), by=(length(gn.windows)-1)/dataPoints)))
df1 <- data.frame("bin"=c(1: length(gn.ordered$ChIPcounts[spaced])),"value"=gn.ordered$ChIPcounts[spaced], "type"=chipname)
df2 <- data.frame("bin"=c(1: length(gn.ordered$INPUTcounts[spaced])),"value"=gn.ordered$INPUTcounts[spaced], "type"=inputname)
df <- data.frame(rbind(df1,df2))
maxdist <- which.max(abs(gn.ordered$ChIPcounts[spaced]-gn.ordered$INPUTcounts[spaced]))

p1 <- ggplot(df, aes(x=bin, y=value, group=type)) +
    geom_vline(xintercept = maxdist, colour="grey") +
    geom_line(aes(color=type)) +
    xlab("Percentage of bins") + ylab("Percentage of reads")  +
    theme(legend.position = "top", legend.text=element_text(size=8), legend.text=element_blank(),  legend.title=element_blank(), legend.background = element_rect(fill = "white", colour = NA)) +
    guides(colour = guide_legend(nrow = 2)) +
    scale_x_continuous(breaks=seq(0,dataPoints, by=(dataPoints)/4), labels=c(0,25,50,75,100))

png(paste(outputpath,"/", chipname, ".png", sep=""), width=250, height=300, units="px")
p1
dev.off()

# print plots
pdf(paste(outputpath,"/", chipname, ".pdf", sep=""), width=4, height=4.5)
p1
dev.off()
