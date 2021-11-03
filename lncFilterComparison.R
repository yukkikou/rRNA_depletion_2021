#identity lncRNA different in tA and rD
#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA/lncFilter") 
library(tidyverse)
library(ggpubr)


#TPM matrix
total_matrix <- read.csv("transcript_tpm_matrix.csv")
#lncRNA id
total_lnc_inter <- read.delim2("final_intersect_lncRNA.id",header = F)
lnc_matrix <- filter(total_matrix,transcript_id %in% total_lnc_inter$V1)
z_score_matrix_all <- t(scale(t(lnc_matrix[,-1])))
pheatmap(z_score_matrix_all)

#######################################################################
###trd=1
lnc_matrix_tmp1 <- filter(lnc_matrix,(C1>1&C2>1&EC3>1&RY3>1&EY3>1&EY4>1&M1>1&M2>1&M3>1&Y1>1&Y2>1&Y3>1))
dim(lnc_matrix_tmp1)
#heatmap
#both phase
z_score_matrix <- t(scale(t(lnc_matrix_tmp1[,-1])))
pheatmap(z_score_matrix)
#mycelium
lnc_matrix_tmp1_m <- lnc_matrix_tmp1[,c(1:4,7:9)]
z_score_matrix_m <- t(scale(t(lnc_matrix_tmp1_m[,-1])))
pheatmap(z_score_matrix_m)
#yeast
lnc_matrix_tmp1_y <- lnc_matrix_tmp1[,-c(2:4,7:9)]
z_score_matrix_y <- t(scale(t(lnc_matrix_tmp1_y[,-1])))
pheatmap(z_score_matrix_y)

#######################################################################
###trd=1 in rD
lnc_matrix_tmp1_rD <-  filter(lnc_matrix,((C1>1&C2>1&EC3>1)|(RY3>1&EY3>1&EY4>1)))
#lnc_matrix_tmp1_rD_loose <-  filter(lnc_matrix,((C1>1|C2>1|EC3>1)|(RY3>1|EY3>1|EY4>1)))
dim(lnc_matrix_tmp1_rD) #229
#dim(lnc_matrix_tmp1_rD_loose) #326

table(rownames(lnc_enrich) %in% lnc_matrix_tmp1_rD$transcript_id) #6
#table(rownames(lnc_enrich) %in% lnc_matrix_tmp1_rD_loose$transcript_id) #22 (all)


#######################################################################
ixd_matrix <- lnc_matrix[,-1]
rownames(ixd_matrix) <- lnc_matrix[,1]
ixd_matrix_0 <- ixd_matrix
ixd_matrix_0[ixd_matrix_0<1]<-0
z_score_matrix_filter <- t(scale(t(ixd_matrix_0)))
z_score_matrix_filter <- z_score_matrix_filter[which(rowSums(z_score_matrix_filter==0)==0),]
#dim(ixd_matrix)
#dim(ixd_matrix[which(rowSums(ixd_matrix)!=0),])
#dim(z_score_matrix_filter)
#pheatmap(z_score_matrix_filter,show_rownames = F)


lnc_targetoff_all <- filter(ixd_matrix_0,(C1==0&C2==0&EC3==0&RY3==0&EY3==0&EY4==0)&(M1>1|M2>1|M3>1|Y1>1|Y2>1|Y3>1))
dim(lnc_targetoff_all) #58

lnc_enrich <- filter(ixd_matrix_0,(C1>1|C2>1|EC3>1|RY3>1|EY3>1|EY4>1)&(M1==0&M2==0&M3==0&Y1==0&Y2==0&Y3==0))
dim(lnc_enrich) #22

#phase-related
#yeast enrich
yeast_enrich <- filter(ixd_matrix,(RY3>1&EY3>1&EY4>1)&(Y1<1&Y2<1&Y3<1)) #12 
#yeast targetoff
yeast_targetoff <- filter(ixd_matrix,(RY3<1&EY3<1&EY4<1)&(Y1>1&Y2>1&Y3>1)) #14
#mycelium enrich
mycelium_enrich <- filter(ixd_matrix,(C1>1&C2>1&EC3>1)&(M1<1&M2<1&M3<1)) #12
#mycelium targetoff
mycelium_targetoff <- filter(ixd_matrix,(C1<1&C2<1&EC3<1)&(M1>1&M2>1&M3>1)) #6

table(rownames(yeast_enrich) %in% rownames(lnc_enrich))
table(rownames(mycelium_enrich) %in% rownames(lnc_enrich))

######################################################################
#output
write.csv(lnc_targetoff_all,"lnc_targetoff.csv")
write_lines(rownames(lnc_targetoff_all),"lnc_targetoff_id.txt")
write.csv(lnc_enrich,"lnc_noval.csv")
write_lines(rownames(lnc_enrich),"lnc_noval_id.txt")
write_lines(lnc_matrix_tmp1_rD$transcript_id,"lnc_tmp1_id.txt")


