#statistics of lncRNA reads ratio in libraries respectly.
#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA/lncRatio")
library(tidyverse)
library(Cairo)

tasample <- c("M1","M2","M3","Y1","Y2","Y3")
rasample <- c("C1","C2","EC3","EY3","EY4","RY3")

tA_ratio <- lncRatioCal("tA",tasample,"poly(A)")
rD_ratio <- lncRatioCal("rD",rasample,"Probe")

ratio <- rbind(tA_ratio,rD_ratio)
ratio$sample <- factor(ratio$sample,levels = c(tasample,rasample))
ratio$id <- factor(ratio$id,levels = c("poly(A)","Probe"))



ggplot(ratio,aes(x=id, y=ratio*100))+
  geom_boxplot(width=0.5) +
  geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white")+
  #geom_point(aes(color=phase))+
  geom_jitter(aes(fill=phase),width =0.2,shape = 21,size=2.5)+
  scale_fill_manual(values = c("#E69F00", "#0072B2"))+ 
  scale_color_manual(values=c("white","white","white"))+
  facet_wrap(~ group) +
  theme_classic() +
  labs(x=NULL, y="lncRNA reads ratio (%)",tag="A")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.88, 0.88),
        legend.text = element_text(size=7,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))

ggsave(file="../figure/final/lncRNA_ratio.pdf",width=210/2,height=297/3,units="mm")



###############################################################
#function
lncRatioCal <- function(directory,sample,id){
  trans_count_matrix <- read.csv(paste(directory,"\\","transcript_count_matrix.csv", sep=""))
  
  lnc_all <- read.delim2(paste(directory,"\\","final_all_lncRNA.id",sep=""),header = F)
  names(lnc_all) <- "lncID"
  lnc_inter <- read.delim2(paste(directory,"\\","final_intersect_lncRNA.id",sep=""),header = F)
  names(lnc_inter) <- "lncID"
  
  ##total counts of each library
  total_counts <- colSums(trans_count_matrix[,-1])
  
  ##lncRNA counts of each library
  #all
  lnc_all_counts_matrix <- filter(trans_count_matrix,transcript_id %in% lnc_all$lncID)
  lnc_all_counts <- colSums(lnc_all_counts_matrix[,-1])
  lnc_all_ratio <- as.data.frame(round(lnc_all_counts/total_counts,4))
  colnames(lnc_all_ratio) <- "ratio"
  lnc_all_ratio$group <- "all lncRNA"
  lnc_all_ratio$sample <- sample
  lnc_all_ratio$phase <- c(rep("mycelium",3),rep("yeast",3))
  rownames(lnc_all_ratio) <- NULL
  #intersect
  lnc_inter_counts_matrix <- filter(trans_count_matrix,transcript_id %in% lnc_inter$lncID)
  lnc_inter_counts <- colSums(lnc_inter_counts_matrix[,-1])
  lnc_inter_ratio <- as.data.frame(round(lnc_inter_counts/total_counts,4))
  colnames(lnc_inter_ratio) <- "ratio"
  lnc_inter_ratio$group <- "intercect lncRNA"
  lnc_inter_ratio$sample <- sample
  lnc_inter_ratio$phase <- c(rep("mycelium",3),rep("yeast",3))
  rownames(lnc_inter_ratio) <- NULL
  #bind
  lnc_ratio <- rbind(lnc_all_ratio,lnc_inter_ratio)
  lnc_ratio$id <- id
  lnc_ratio  
}
