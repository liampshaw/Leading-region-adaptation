# Palindrome density
library(ggplot2)

dir.create("../figures/Supp/", showWarnings = FALSE, recursive = TRUE)


plasmid_table = read.csv('../data/TableS1_plasmid_metadata.csv', header=T,
                         row.names = 1)
plasmid_table = plasmid_table[plasmid_table$Leading.region.analysis==TRUE,]

d = read.csv('../data/palindrome_density.tsv', header=T, sep='\t')
d$PTU = plasmid_table[d$contig,"PTU"]

d$start.factor.kb = as.factor(d$start/1000)
pdf('../figures/Supp/figS4-palindrome-density-einverted.pdf', width=8, height=6)
ggplot(d[which(d$start<50000),], 
       aes(start.factor.kb, count))+
  geom_boxplot()+facet_wrap(~PTU, scales="free_y")+
  xlab("Position (kb)")+
  ylab("Number of palindromes per 5kb")+
  theme_bw()
dev.off()
