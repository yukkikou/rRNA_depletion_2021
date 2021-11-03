#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA")

library(tidyverse)
library(doBy)
library(RColorBrewer)
library(ggpp)
library(ggpmisc)
library(corrplot)
library(patchwork)
#install.packages("plotrix")
library(plotrix)
library(ggpubr)
display.brewer.all()
brewer.pal(11,"RdBu") 
display.brewer.pal(8,"Set2")
display.brewer.pal(9,"YlOrRd")
brewer.pal(8,"Set2")
brewer.pal(9,"Reds")

###################################################################
##part 1
ratio <- read.csv("C:/Users/16103/Desktop/lncRNA/rRNA_ratio_R.csv")
rownames(ratio) <- ratio[,1]
ratio <- ratio[,-1]
r <- c("5S","5.8S","ITS1","ITS2","18S","28S","rRNA")
yeast_ratio <- data.frame("group"=c(rep("RY3",7),rep("EY3",7),rep("EY4",7)),
                          "rRNA"=rep(rownames(ratio),3),
                          "r3"=c(ratio$Y3,ratio$EY3,ratio$EY4))
yeast_ratio$group <- factor(yeast_ratio$group,level=c("RY3","EY3","EY4"))
yeast_ratio$rRNA <- factor(yeast_ratio$rRNA,levels = r)

ggplot(yeast_ratio,aes(x=rRNA,y=r3,fill=group))+
  geom_bar(stat = "identity",position=position_dodge(),width=0.6)+
  scale_fill_brewer(palette = "Set2")+
  geom_hline(yintercept = 1,colour = "grey11",linetype="dashed")+
  ylim(0,1.2)+
  labs(x=NULL, y="rRNA residual ratio (%)",tag="D")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.1, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))

ggsave(file="figure/final/yeast_rD_ratio.pdf",width=210/2,height=297/4.8,units="mm")


###################################################################
##part 2A
rr <- data.frame("group"=c(rep("Probe",3),rep("Ribo-Zero Kit",4)),
                 "r3"=c(0.82,0.76,0.79,1.31,1.50,11.05,4.37))
#Wilcoxon rank sum test
x <- c(0.82,0.76,0.79)
y <- c(1.31,1.50,11.05,4.37)
wilcox.test(x,y,paired = F, alternative = "less")
t.test(x,y,paired = F)

me <- summaryBy(r3~group,rr,FUN=mean)
sd <- tapply(rr$r3,rr$group,std.error)
tapply(rr$r3,rr$group,sd)/2
rr <- left_join(rr,me,by=c("group" = "group"))
sd <- data.frame("group"=c("Probe","Ribo-Zero Kit"),
                 "r3.sd"=sd)
rr <- left_join(rr,sd,by=c("group" = "group"))
rr$group <- factor(rr$group,levels = c("Probe","Ribo-Zero Kit"))

ggplot(rr,aes(x=group,ymin=r3.mean-r3.sd, ymax=r3.mean+r3.sd,fill=group)) +
  geom_bar(position=position_dodge(), aes(y=r3.mean), stat="identity",width=0.4)+
  geom_point(aes(y=r3),fill="black")+
  geom_errorbar(position=position_dodge(0.4), colour="black",size=0.7,width=.2)+
  scale_fill_brewer(palette = "Set2")+
  labs(y="rRNA residual ratio (%)",x=NULL,tag="A")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position="none")
ggsave(file="figure/final/yeast_rD_hk_r3.pdf",width=210/3,height=297/5,units="mm")


###################################################################
##part 2B
hk_ratio <- data.frame("group"=c(rep("PY-1",7),rep("PY-2",7),rep("PY-3",7),rep("PY-4",7)),
                          "rRNA"=rep(rownames(ratio),4),
                          "r3"=as.numeric(c(ratio$SRR5028789,ratio$SRR5028790,ratio$SRR6435878,ratio$SRR6435879)))
hk_ratio$rRNA <- factor(hk_ratio$rRNA,levels = r)

hk_ratio <- hk_ratio[,2:3]

me <- summaryBy(r3~rRNA,hk_ratio,FUN=mean)
sd <- tapply(hk_ratio$r3,hk_ratio$rRNA,std.error)
hk_ratio <- left_join(hk_ratio,me,by=c("rRNA" = "rRNA"))
sd <- data.frame("rRNA"=r,
                 "r3.sd"=sd)
hk_ratio <- left_join(hk_ratio,sd,by=c("rRNA" = "rRNA"))

ggplot(hk_ratio,aes(x=rRNA,ymin=r3.mean-r3.sd, ymax=r3.mean+r3.sd,fill=rRNA)) +
  geom_bar(position=position_dodge(), aes(y=r3.mean), stat="identity",width=0.6)+
  geom_point(aes(y=r3),fill="black")+
  geom_errorbar(position=position_dodge(0.4), colour="black",size=0.7,width=.2)+
  scale_fill_brewer(palette = "Set2")+
  theme_classic()+
  labs(x=NULL, y="rRNA residual ratio (%)",tag="B")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.1, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))
ggsave(file="figure/final/yeast_hk_ratio.pdf",width=210/3*2,height=297/4.8,units="mm")

###################################################################
#rRNA_depth.R for FigC&D


###################################################################
##part 2E
#CV
s3_gene_count <- read.csv("comparison/s3_gene_count_matrix.csv")
head(s3_gene_count)
#yeast
yeast_cv <- cal_CV(s3_gene_count,c(5:6,14),"Probe")
#mycelium
mycelium_cv <- cal_CV(s3_gene_count,c(2:4),"Probe")
#HK-yeast
kit_cv <- cal_CV(s3_gene_count,c(10:13),"Ribo-Zero Kit")
cv <- rbind(yeast_cv,kit_cv)
ggplot(cv,aes(x=group,ymin=cv.mean-cv.sd, ymax=cv.mean+cv.sd,fill=group)) +
  geom_bar(position=position_dodge(), aes(y=cv.mean), stat="identity",width=0.4)+
  geom_errorbar(position=position_dodge(0.4), colour="black",size=0.7,width=.2)+
  scale_fill_brewer(palette = "Accent")+
  theme_classic()+
  labs(y="Mean coefficient of variance",x=NULL,tag="E  ")+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position="none")
ggsave(file="figure/final/CV_rD_hk_r3.pdf",width=210/3,height=297/5,units="mm") 


###################################################################
##part 3B
#gene expression correlation  (lm)
#log10 (TPM + 1)
#merge togetger and requantity
co_tpm <- read.csv("comparison/s3_gene_tpm_matrix.csv")
#co_tpm <- read.csv("comparison/rD_tA/gene_tpm_matrix.csv")
head(co_tpm)
rD_gene_tpm_log10 <- geneExplog10(co_tpm,c(5:6,14))
head(rD_gene_tpm_log10)
hk_gene_tpm_log10 <- geneExplog10(co_tpm,c(10:13))
head(hk_gene_tpm_log10)

#lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")
co_rD_hk <- coLinear(rD_gene_tpm_log10,hk_gene_tpm_log10,c(4,10),1,2,"F")
ggsave(file="figure/final/rD_hk_lm.pdf",width=210/3,height=75,units="mm")

###################################################################
##part 3C (or delect)
#Q-Q plot

plot.new()

rD_hk_qq <- cbind(rD_gene_tpm_log10$mean,hk_gene_tpm_log10$mean)
rD_tA_qq <- cbind(rD_gene_tpm_log10$mean,tA_gene_tpm_log10$mean)
tA_hk_qq <- cbind(tA_gene_tpm_log10$mean,hk_gene_tpm_log10$mean)

pdf(file = "figure/final/s3_qq.pdf", width = 4, height = 4)
qqplot(rD_hk_qq[,1],rD_hk_qq[,2],col="blue",
       pch=1,cex=1.2,xlab = "Probe depletion", ylab = "Ribo-Zreo Kit",
       cex.axis=1.5,cex.lab=1.5)
lines(0:4.5, 0:4.5,lwd=3,col="black")
qqplot(rD_tA_qq[,1],rD_tA_qq[,2],col="blue",
       pch=1,cex=1.2,xlab = "Probe depletion", ylab = "poly(A) enrichment",
       cex.axis=1.5,cex.lab=1.5)
lines(0:4.5, 0:4.5,lwd=3,col="black")
qqplot(tA_hk_qq[,1],tA_hk_qq[,2],col="blue",
       pch=1,cex=1.2,xlab = "poly(A) enrichment", ylab = "Ribo-Zreo Kit",
       cex.axis=1.5,cex.lab=1.5)
lines(0:4.5, 0:4.5,lwd=3,col="black")
dev.off()


###################################################################
##part 3D
#corrplot of all samples
co_tpm_matrix <- co_tpm[,c(5:6,10:14)] 
rownames(co_tpm_matrix) <- co_tpm$gene_id
co_tpm_matrix <- log10(co_tpm_matrix+1)
co_cor <- cor(co_tpm_matrix)

pdf(file = "figure/final/rD_hk_correlation.pdf", width = 4, height = 3.5)
pheatmap(co_cor,display_numbers = TRUE,number_color = "black", fontsize_number = 6,
         color = colorRampPalette((brewer.pal(n = 9, name ="Blues")))(100),
         border_color = "white")
dev.off()

###################################################################
##part 4B
my_ratio <- data.frame("group"=c(rep("C1",7),rep("C2",7),rep("EC3",7)),
                          "rRNA"=rep(rownames(ratio),3),
                          "r3"=c(ratio$C1,ratio$C2,ratio$EC3))
my_ratio$group <- factor(my_ratio$group,level=c("C1","C2","EC3"))
my_ratio$rRNA <- factor(my_ratio$rRNA,levels = r)

ggplot(my_ratio,aes(x=rRNA,y=r3,fill=group))+
  geom_bar(stat = "identity",position=position_dodge(),width=0.6)+
  scale_fill_brewer(palette = "Set2")+
  ylab("rRNA residual ratio (%)")+
  geom_hline(yintercept = 0.4,colour = "grey11",linetype="dashed")+
  ylim(0,0.5)+
  labs(x=NULL, y="rRNA residual ratio (%)",tag="B")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.1, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))

ggsave(file="figure/final/mycelium_rD_ratio.pdf",width=210/2,height=297/4.8,units="mm")

###################################################################
##part 4C
rr <- data.frame("group"=c(rep("Probe-Mycelium",3),rep("Probe-Yeast",3),rep("Kit-Yeast",4)),
                 "r3"=c(0.36,0.39,0.23,0.82,0.76,0.79,1.31,1.50,11.05,4.37))
me <- tapply(rr$r3,rr$group,mean)
me <- data.frame("group"=names(me),"r3.mean"=me)
sd <- tapply(rr$r3,rr$group,std.error)
sd <- data.frame("group"=names(sd),"r3.sd"=sd)
rr <- left_join(rr,me,by=c("group" = "group"))
rr <- left_join(rr,sd,by=c("group" = "group"))
rr$group <- factor(rr$group,levels = rev(c("Probe-Mycelium","Probe-Yeast","Kit-Yeast")))

ggplot(rr,aes(x=group,ymin=r3.mean-r3.sd, ymax=r3.mean+r3.sd,fill=group)) +
  geom_bar(position=position_dodge(), aes(y=r3.mean), stat="identity",width=0.5)+
  geom_point(aes(y=r3),fill="black")+
  geom_errorbar(position=position_dodge(0.4), colour="black",size=0.7,width=.2)+
  #geom_hline(yintercept = 1,colour = "grey11",linetype="dashed")+
  scale_fill_brewer(palette = "Set2")+
  geom_text(aes(y = r3.mean, label = round(r3.mean,2),vjust = 0, hjust = -0.8), 
            show.legend = F,color="dimgrey")+
  labs(x=NULL, y="rRNA residual ratio (%)",tag="C")+
  theme_classic()+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position="none")
ggsave(file="figure/final/s3_rRNA_ratio.pdf",width=210/2,height=297/5,units="mm")

#Wilcoxon rank sum test
py <- c(0.82,0.76,0.79)
pm <- c(0.36,0.39,0.23)
ky <- c(1.31,1.50,11.05,4.37)

x <- list(pm,py,ky)

kruskal.test(x)

wilcox.test(pm,py,paired = F,alternative = "less")
wilcox.test(pm,ky,paired = F,alternative = "less")
wilcox.test(py,ky,paired = F,alternative = "less")

###################################################################
##part 4D: qualimap
###################################################################
##part 4E: lm_lable_rR_tA.R

###################################################################
##part 5A: lncRatio.R
#(##part 5B: info.R)
##part 5B: lnc.R (lncRNA_rlr.pdf)
##part 5C: lnc.R (lnc_m_y.pdf)
##part 5D: gel and qPCR

#############################################################
#function
cal_CV <- function(data,range,group){
  #for top 1000 gene
  rownames(data) <- data$gene_id
  data <- data[,range]
  data$mean <- apply(data,1,mean)
  data$sd <- apply(data,1,std.error)
  data$cv <- data$sd/data$mean
  data <- data[order(-data$mean),]
  data <- data[1:1000,]
  data_cv_mean <- mean(data$cv)
  data_cv_sd <- sd(data$cv)
  re <- data.frame("cv.mean"=data_cv_mean,
                   "cv.sd"=data_cv_sd,
                   "group"=group)
  re
}

geneExplog10 <- function(rD_gene_tpm,range){
  rownames(rD_gene_tpm) <- rD_gene_tpm$gene_id
  rD_gene_tpm_matrix <- rD_gene_tpm[,range]
  rD_gene_tpm_matrix_log10 <- log10(rD_gene_tpm_matrix+1)
  rD_gene_tpm_matrix_log10$mean <- apply(rD_gene_tpm_matrix_log10,1,mean)
  rD_gene_tpm_matrix_log10$sd <- apply(rD_gene_tpm_matrix_log10,1,std.error)
  rD_gene_tpm_matrix_log10
}

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

coLinear <- function(rD_gene_tpm_log10,hk_gene_tpm_log10,range,xnlab,ynlab,tag){
  co_tpm_log10 <- cbind(rD_gene_tpm_log10,hk_gene_tpm_log10)
  co_tpm_log10 <- co_tpm_log10[,range]
  #head(co_tpm_log10)
  #cor.test(co_tpm_log10$mean,co_tpm_log10$mean.1,method = 'pearson')
  lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")
  my.formula <- y ~ x
  lm.scattor <- lm(co_tpm_log10$mean.1 ~ co_tpm_log10$mean)
  summary(lm.scattor)
  co <- ggplot(co_tpm_log10,aes(x=mean,y=mean.1))+
    geom_point(size=1.3,colour="#E5C494")+
    geom_smooth(method = "lm", se=FALSE, color="black", formula = my.formula)+
    #stat_poly_eq(formula = my.formula, 
    #             aes(label = paste(..rr.label.., ..p.value.label.., sep = "~`,`~")), 
    #             parse = TRUE) +
    stat_cor(method="pearson",label.y = 4.9,size=3)+
    stat_regline_equation(label.y = 4.6,size=3)+
    geom_abline(coef = c(1, 1), col = "darkred",size=1)+
    theme_classic()+
    labs(y=paste(lab[ynlab],"log10(TPM+1)"),x=paste(lab[xnlab],"log10(TPM+1)"),tag=tag)+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=10),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
}

