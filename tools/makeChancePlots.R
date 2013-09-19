library(Repitools)
library(ggplot2)
library(scales)  

# function to plot enrichment 
plotfxn<-function(df, dataPoints, maxdist)
{
  p <- ggplot(df, aes(x=bin, y=value, group=type)) +
    geom_vline(xintercept = maxdist, colour="grey") +
    geom_line(aes(color=type)) +
    xlab("Percentage of bins") + ylab("Percentage of reads")  +
    theme(legend.position = c(0, 1), legend.justification = c(0, 1),  legend.text=element_text(size=8), legend.text=element_blank(),  legend.title=element_blank(), legend.background = element_rect(fill=alpha('blue', 0.35), colour = NA)) +
    scale_color_manual(values=c("#ce423c", "#4981b4")) +
    guides(colour = guide_legend(nrow = 2)) +
    scale_x_continuous(breaks=seq(0,dataPoints, by=(dataPoints)/4), labels=c(0,25,50,75,100)) 

  return (p)
}

args <- commandArgs(trailingOnly = TRUE)
chip<-args[1]
input<-args[2]
chipname <-args[3]
inputname<-args[4]
outputpath<-args[5]

windowSize=1000
dataPoints=1000

rs.ChIP <- BAM2GRanges(chip)
rs.INPUT <- BAM2GRanges(input)
#only count read starts, avoids reads being across multiple windows
rs.ChIP <- resize(rs.ChIP, 1, fix="start")
rs.INPUT <- resize(rs.INPUT, 1, fix="start")

# get the genome sizes form the bam header
rs.bh <- scanBamHeader(chip)
rs.chrSizes <- (rs.bh[[names(rs.bh)]][["targets"]])
    
gn.windows <- genomeBlocks(rs.chrSizes, chrs=seq(length(rs.chrSizes)), windowSize)
ChIPcounts <- countOverlaps(gn.windows, rs.ChIP)
INPUTcounts <- countOverlaps(gn.windows, rs.INPUT)

gn.ordered = data.frame(cbind(ChIPcounts, INPUTcounts)[order(ChIPcounts),])

gn.ordered$ChIPcounts = cumsum(gn.ordered$ChIPcounts)
gn.ordered$INPUTcounts = cumsum(gn.ordered$INPUTcounts)

gn.ordered$ChIPcounts = gn.ordered$ChIPcounts/gn.ordered$ChIPcounts[length(gn.windows)]
gn.ordered$INPUTcounts = gn.ordered$INPUTcounts/gn.ordered$INPUTcounts[length(gn.windows)]	

rs.spaced <- c(round(seq(1,length(gn.windows), by=(length(gn.windows)-1)/dataPoints)))
rs.df1 <- data.frame("bin"=c(1: length(gn.ordered$ChIPcounts[rs.spaced])),"value"=gn.ordered$ChIPcounts[rs.spaced], "type"=chipname)
rs.df2 <- data.frame("bin"=c(1: length(gn.ordered$INPUTcounts[rs.spaced])),"value"=gn.ordered$INPUTcounts[rs.spaced], "type"=inputname)
df <- data.frame(rbind(rs.df1,rs.df2))
maxdist <- which.max(abs(gn.ordered$ChIPcounts[rs.spaced]-gn.ordered$INPUTcounts[rs.spaced]))

p1 <- plotfxn(df, dataPoints, maxdist)

png(paste(outputpath,"/", chipname, ".png", sep=""), width=250, height=250, units="px")
p1
dev.off()

# print plots
pdf(paste(outputpath,"/", chipname, ".pdf", sep=""), width=4, height=4)
p1
dev.off()

# remove all variables not needed anymore
rm(list=setdiff(ls(), c("outputpath","chipname", "inputname", "df", "plotfxn", "dataPoints", "maxdist")))
# save session
save.image(file = paste(outputpath,"/", chipname, ".Rdata", sep=""))


