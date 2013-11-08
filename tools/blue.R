library("ggplot2")

args <- commandArgs(trailingOnly = TRUE)

pdf(file=args[2], width=15,height=7)

x<-read.delim(args[1], row.names = NULL ,header=T, quote="\"")
#head(x)

ggplot(x, aes(x=copy, y=perc)) +
  geom_line(aes(color=sample)) +
  labs(y = "frequency (%)", x = "copies of reads", title = "Read histogram")
dev.off()

sink(type = "message")
sessionInfo()
