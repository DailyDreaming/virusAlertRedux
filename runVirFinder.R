#!/usr/bin/env Rscript
#usage: Rscript --vanilla /home/quokka/PycharmProjects/biotool/runVirFinder.R /home/quokka/PycharmProjects/data/SRR6172653.fastq 0.95 /home/quokka/PycharmProjects/tools/R
# /home/quokka/PycharmProjects/data/SRR6172653.fastq
# predResult <- VF.pred('/home/quokka/PycharmProjects/data/SRR6172653.fastq')
args = commandArgs(trailingOnly=TRUE)

#args[1] should contain path points to the fasta file
#args[2] is the significance cutoff

# .libPaths('/home/quokka/PycharmProjects/tools/R')
.libPaths(args[3])

# library("tidyverse")
library(VirFinder)

# fasta.inputfile ="data/SRR5383888.fa"
# threshold=0.95

fasta.inputfile <- args[1]
threshold <- args[2]

pred.result <- VF.pred(fasta.inputfile) %>%
  arrange(score) %>%
  mutate(fname = row_number())

write.table(pred.result, "virfinder_hits.csv", quote = F, sep=",", col.names=F)

ggplot(pred.result, aes(x=fname, y=score))+
  geom_point()+
  ylim(c(0,1))+
  geom_hline(yintercept = as.numeric(threshold), color="red")+
  theme_bw()+
  xlab("read index")+
  ylab("")

ggsave("virfinder_hits.png")

