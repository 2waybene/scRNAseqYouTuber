#!/usr/bin/env Rscript
## Credit goes to https://www.biostars.org/p/355901/
library(ggplot2)

dt <- read.csv("Atg9b_expression.csv", header = TRUE)

dt
# set up `d` per your code
# add `sample`, `treatment`, and `expression` columns per this answer



table(dt$Treatment)
table(dt$Type)
p <- ggplot(dt, aes(group=Type)) + 
 # geom_bar(aes(x=Type, y=RPKM, fill=Treatment), width=0.75, position="dodge2", stat="identity") + 
  geom_bar(aes(x=Type, y=Count, fill=Treatment), width=0.75, position="dodge2", stat="identity") + 
 # facet_grid(~Treatment) +
  scale_colour_brewer(palette="Set1") + 
  scale_fill_brewer(palette="Set1") +
  #ggtitle('Shaw RNAseq Atg9b RPKM') +
  ggtitle('Shaw RNAseq Atg9b count') +
  xlab('Samples') +
  ylab('Expression') +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p)



