library(ggplot2)
library(reshape)

args <- commandArgs(trailingOnly = TRUE)
pdf(file=args[2], width=9,height=7)

file=args[1]
x<-read.table(file, row.names = NULL ,header=T, quote="\"")
head(x)
distribution <- melt(x[,c(2:length(x))], id.vars="sample")
colnames(distribution)=c("sample","feature","value")

ggplot(distribution, aes(x = sample, y=value)) + 
  geom_bar(stat="identity", aes(fill = feature), position = "fill") + 
  scale_y_continuous("fraction") + 
  labs(title=args[3]) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  coord_flip()
#  facet_grid(. ~ type , space = "free", scales = "free_x")

ggplot(distribution, aes(x = sample, y=value)) + 
  geom_bar(stat="identity", aes(fill = feature)) +
  labs(title=args[3]) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  coord_flip()  
#  facet_grid(. ~ type , space = "free", scales = "free_x")


sink(type = "message")
sessionInfo()
