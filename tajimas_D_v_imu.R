rm(list = ls())
library(tidyverse)
library(ape)
library(pegas)
library(PopGenome)
library(xlsx)
library(ggplot2)

lrf=c()
all_D=c()
gene2_D=c()
gene3_D=c()
gene4_D=c()
intergene_D=c()
all_p=c()
gene2_p=c()
gene3_p=c()
gene4_p=c()
intergene_p=c()

for (i in 0:19){
  rf=5*i
  print(rf)
  all_genome = read.dna(paste("H:/sim/genome_0_all_imu_",format(rf, digits=1, nsmall=0),".fas",sep=''), format = "fasta")
  gene2 = read.dna(paste("H:/sim/genome_0_gene2_imu_",format(rf, digits=1, nsmall=0),".fas",sep=''), format = "fasta")
  gene3 = read.dna(paste("H:/sim/genome_0_gene3_imu_",format(rf, digits=1, nsmall=0),".fas",sep=''), format = "fasta")
  gene4 = read.dna(paste("H:/sim/genome_0_gene4_imu_",format(rf, digits=1, nsmall=0),".fas",sep=''), format = "fasta")
  intergene = read.dna(paste("H:/sim/genome_0_intergene_imu_",format(rf, digits=1, nsmall=0),".fas",sep=''), format = "fasta")
  result_all=tajima.test(all_genome)
  result_gene2=tajima.test(gene2)
  result_gene3=tajima.test(gene3)
  result_gene4=tajima.test(gene4)
  result_intergene=tajima.test(intergene)
  
  all_D=c(all_D,result_all$D)
  gene2_D=c(gene2_D,result_gene2$D)
  gene3_D=c(gene3_D,result_gene3$D)
  gene4_D=c(gene4_D,result_gene4$D)
  intergene_D=c(intergene_D,result_intergene$D)
  
  all_p=c(all_p,result_all$Pval.normal)
  gene2_p=c(gene2_p,result_gene2$Pval.normal)
  gene3_p=c(gene3_p,result_gene3$Pval.normal)
  gene4_p=c(gene4_p,result_gene4$Pval.normal)
  intergene_p=c(intergene_p,result_intergene$Pval.normal)
  lrf=c(lrf,rf)
}
df_all=data.frame(rf=lrf,D=all_D, pnormal=all_p)
df_gene2=data.frame(rf=lrf,D=gene2_D, pnormal=gene2_p)
df_gene3=data.frame(rf=lrf,D=gene3_D, pnormal=gene3_p)
df_gene4=data.frame(rf=lrf,D=gene4_D, pnormal=gene4_p)
df_intergene=data.frame(rf=lrf,D=intergene_D, pnormal=intergene_p)

write.xlsx(df_all, file='H:/sim/tajimas_D_vs_imu.xlsx',
           sheetName = "all_D", append = FALSE)
write.xlsx(df_gene2, file='H:/sim/tajimas_D_vs_imu.xlsx',
           sheetName = "gene2", append = TRUE)
write.xlsx(df_gene3, file='H:/sim/tajimas_D_vs_imu.xlsx',
           sheetName = "gene3", append = TRUE)
write.xlsx(df_gene4, file='H:/sim/tajimas_D_vs_imu.xlsx',
           sheetName = "gene4", append = TRUE)
write.xlsx(df_intergene, file='H:/sim/tajimas_D_vs_imu.xlsx',
           sheetName = "intergene", append = TRUE)

ggplot(df_all, aes(x=rf, y=all_D))+
  geom_bar(stat='identity', fill="#999999")+
  xlab("RF") + ylab("Tajima's D: whole genome")
ggsave("H:/sim/tD_all_vs_imu.png")

ggplot(df_gene2, aes(x=rf, y=gene2_D))+
  geom_bar(stat='identity', fill="#E69F00")+
  xlab("RF") + ylab("Tajima's D: gene2")
ggsave("H:/sim/tD_gene2_vs_imu.png")

ggplot(df_gene3, aes(x=rf, y=gene3_D))+
  geom_bar(stat='identity', fill="#56B4E9")+
  xlab("RF") + ylab("Tajima's D: gene3")
ggsave("H:/sim/tD_gene3_vs_imu.png")

ggplot(df_gene4, aes(x=rf, y=gene4_D))+
  geom_bar(stat='identity', fill="#009E73")+
  xlab("RF") + ylab("Tajima's D: gene4")
ggsave("H:/sim/tD_gene4_vs_imu.png")

ggplot(df_intergene, aes(x=rf, y=intergene_D))+
  geom_bar(stat='identity', fill="#CC79A7")+
  xlab("RF") + ylab("Tajima's D: intergene")
ggsave("H:/sim/tD_intergene_vs_imu.png")
