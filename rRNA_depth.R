library(ggplot2)
library(patchwork)

rD_samples <- c("C1","C2","EC1","EC3","EY3","EY4","Y2","Y3")
hk_samples <- c("PY.1","PY.2","PY.3","PY.4")
#salist

rD_rRNA_28S <- preMatrix("rRNA/rD_depth_28S.txt",rD_samples)
hk_rRNA_28S <- preMatrix("rRNA/HK_depth_28S.txt",hk_samples)
rD_rRNA_18S <- preMatrix("rRNA/rD_depth_18S.txt",rD_samples)
hk_rRNA_18S <- preMatrix("rRNA/HK_depth_18S.txt",hk_samples)

s18 <- displot(rD_rRNA_18S,hk_rRNA_18S,"C ")
ggsave(file="figure/final/18S_depth.pdf",width=210/2,height=297/5*3,units="mm")
s28 <- displot(rD_rRNA_28S,hk_rRNA_28S,"D ")
ggsave(file="figure/final/28S_depth.pdf",width=210/2,height=297/5*3,units="mm")

s18|s28
ggsave(file="figure/final/18S_28S_depth.pdf",width=200,height=297/5*3,units="mm")

###############################################################
#function
preMatrix <- function(fileurl,samples){
  r <- read.delim2(fileurl,header = F,sep = "\t")
  colnames(r) <- c("chr","position",samples)
  r_ratio <- r[,-c(1,2)]
  r_colsum <- colSums(r_ratio)
  r_ratio <- r_ratio/r_colsum*10000
  r_ratio$position <- r$position
  r_ratio$chr <- r$chr
  r_ratio  
}

displot <- function(rD_rRNA_28S,hk_rRNA_28S,tag){
  y3 <- ggplot(rD_rRNA_28S,aes(x=position,y=Y3))+
    geom_bar(stat = "identity")+
    theme_classic()+ 
    ylim(0,100)+
    labs(x=NULL, y="RY3",tag=tag)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  ey3 <- ggplot(rD_rRNA_28S,aes(x=position,y=EY3))+
    geom_bar(stat = "identity")+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="EY3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  ey4 <- ggplot(rD_rRNA_28S,aes(x=position,y=EY4))+
    geom_bar(stat = "identity")+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="EY4")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  py3 <- ggplot(hk_rRNA_28S,aes(x=position,y=PY.1))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="PY.3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  py4 <- ggplot(hk_rRNA_28S,aes(x=position,y=PY.2))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="PY.4")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  py1 <- ggplot(hk_rRNA_28S,aes(x=position,y=PY.3))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="PY.1")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  py2 <- ggplot(hk_rRNA_28S,aes(x=position,y=PY.4))+
    geom_bar(stat = "identity")+
    theme_classic()+
    theme_classic()+
    ylim(0,100)+
    labs(x=NULL, y="PY.2")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
  
  y3/ey3/ey4/py1/py2/py3/py4
}

