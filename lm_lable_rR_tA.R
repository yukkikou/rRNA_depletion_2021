#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpp)
library(ggpmisc)
library(corrplot)
library(patchwork)

#part 4E
co_tpm <- read.csv("comparison/s3_gene_tpm_matrix.csv")
#co_tpm <- read.csv("comparison/rD_tA/gene_tpm_matrix.csv")
head(co_tpm)
rD_gene_tpm_log10 <- geneExplog10(co_tpm,c(5:6,14))
head(rD_gene_tpm_log10)
hk_gene_tpm_log10 <- geneExplog10(co_tpm,c(10:13))
head(hk_gene_tpm_log10)
tA_gene_tpm_log10 <- geneExplog10(co_tpm,c(15:17))
head(tA_gene_tpm_log10)


#lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")
co_rD_tA <- coLinearLable(rD_gene_tpm_log10,tA_gene_tpm_log10,c(4,9),1,3,"E  ")
co_hk_tA <- coLinearLable(hk_gene_tpm_log10,tA_gene_tpm_log10,c(5,10),2,3,"")
co_rD_tA|co_hk_tA
ggsave(file="figure/final/r_tA_lm.pdf",width=210*0.9,height=95,units="mm")


#############################################################
#validation
rA_ab <- filter(rA_tpm_log10, abs(mean-mean.1)>1)
hA_ab <- filter(hA_tpm_log10, abs(mean-mean.1)>1)
ra <- lmpolt(rA_ab,1,3)
ha <- lmpolt(hA_ab,2,3)
ra|co_rD_tA
ha|co_tA_hk
ra|ha
exp(tA_gene_tpm_log10,"MSTRG.7192")
exp(rD_gene_tpm_log10,"MSTRG.2219")
##SUPPLEMENT
write.csv(rA_ab,"figure/final/supp/rA_abnormal_point_over_1.csv")
write.csv(hA_ab,"figure/final/supp/hA_abnormal_point_over_1.csv")


#############################################################
#function

coLinearLable <- function(rD,tA,range,xnlab,ynlab,tag){
  co_tpm_log10 <- cbind(rD,tA)
  co_tpm_log10 <- co_tpm_log10[,range]
  co_tpm_log10$trend <- sapply((co_tpm_log10$mean-co_tpm_log10$mean.1), 
                               function(x){if(x>1) 'up' else if(x<(-1)) 'down' else 'none'})
  #labeled point
  up <- filter(co_tpm_log10, mean-mean.1>2)
  up$gene <- rownames(up)
  down <- filter(co_tpm_log10, mean-mean.1<(-2))
  down$gene <- rownames(down)
  if (nrow(down)>3){
    down$minus <- down$mean-down$mean.1
    down <- down[order(down$minus),]
    down$gene[4:nrow(down)] <- ""    
  }
  lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")
  my.formula <- y ~ x
  lm.scattor <- lm(co_tpm_log10$mean.1 ~ co_tpm_log10$mean)
  #summary(lm.scattor)
  co <- ggplot(co_tpm_log10,aes(x=mean,y=mean.1,color=trend))+
    geom_point(size=1)+
    geom_smooth(inherit.aes = F, data=co_tpm_log10, aes(x=mean, y=mean.1), 
                method = "lm", se=FALSE, color="black", formula = my.formula)+
    stat_cor(inherit.aes = F, data=co_tpm_log10, aes(x=mean, y=mean.1),method="pearson",label.y = 4.9,size=3)+
    stat_regline_equation(inherit.aes = F, data=co_tpm_log10, aes(x=mean, y=mean.1),label.y = 4.6,size=3)+
    scale_color_manual(values = c("red","grey","blue"))+
    geom_point(inherit.aes = F, data=down, aes(x=mean, y=mean.1), 
               size = 1.2, color='darkred')+ 
    geom_point(inherit.aes = F, data=up, aes(x=mean, y=mean.1), 
               size = 1.2, color='darkblue')+ 
    ggrepel::geom_text_repel(inherit.aes = F, data = up, 
                             aes(x=mean, y=mean.1, label=gene), size=2.5)+
    ggrepel::geom_text_repel(inherit.aes = F, data = down, 
                             aes(x=mean, y=mean.1, label=gene), size=2.5)+
    #stat_poly_eq(formula = my.formula, 
    #             aes(label = paste(..rr.label.., ..p.value.label.., sep = "~`,`~")), 
    #             parse = TRUE) +
    geom_abline(coef = c(1, 1), col = "darkred",size=1)+
    theme_classic()+
    labs(y=paste(lab[ynlab],"log10(TPM+1)"),x=paste(lab[xnlab],"log10(TPM+1)"),tag=tag)+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
}

lmpolt <- function(co_tpm_log10,xnlab,ynlab){
  my.formula <- y ~ x
  lm.scattor <- lm(co_tpm_log10$mean.1 ~ co_tpm_log10$mean)
  summary(lm.scattor)
  lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")  
  co <- ggplot(co_tpm_log10,aes(x=mean,y=mean.1,label = rownames(co_tpm_log10)))+
    geom_point(size=1.3,colour="#E5C494")+
    theme_classic()+
    geom_text()+ylim(0,5)+xlim(0,5)+
    labs(y=paste(lab[ynlab],"log10(TPM+1)"),x=paste(lab[xnlab],"log10(TPM+1)"))+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  co
}

exp <- function(tpm_log10,gene){
  tpm_log10[rownames(tpm_log10)==gene,]
}
