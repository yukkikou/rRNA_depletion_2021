#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA/lncFilter/top10")
library(tidyverse)
library(doBy)
library(ggpubr)


qRe <- read_csv("hxy0727_qPCR.csv")
#View(qRe)

#caculate delta CT for detector
Mendo <- filter(qRe,Sample=="Mycelium"&Detector=="ACT")
Yendo <- filter(qRe,Sample=="Yeast"&Detector=="ACT")
#filter(qRe,Sample=="Mycelium"&Detector=="H2O")
#filter(qRe,Sample=="Yeast"&Detector=="H2O")
qRe$endo <- Mendo$AvgCt[1]
qRe$endo[qRe$Sample=="Yeast"] <- Yendo$AvgCt[1]
table(qRe$endo)
qRe$deltaCT <- qRe$Ct-qRe$endo

#find outlier
ggplot(qRe)+
  geom_bar(aes(x=Detector,y=deltaCT,fill=Sample),stat = "identity",position = "dodge")+
  theme_classic()

#remove outlier
qRe <- filter(qRe,Detector!="M6914")


#recaculate sd
deltaCTSD <- summaryBy(deltaCT~Detector+Sample,qRe,FUN=std.error)
deltaCTmean <- summaryBy(deltaCT~Detector+Sample,qRe,FUN=mean)

#aggregate(deltaCT~Detector, qRe, FUN=sd)
qRe <- left_join(qRe,deltaCTmean,by=c("Detector" = "Detector","Sample"="Sample"))
qRe <- left_join(qRe,deltaCTSD,by=c("Detector" = "Detector","Sample"="Sample"))

#View(qRe)
#unique(qRe$Detector)

#add error bar
ggplot(qRe,aes(x=Detector,ymin=deltaCT.mean-deltaCT.std.error, ymax=deltaCT.mean+deltaCT.std.error,fill=Sample)) +
  geom_bar(position=position_dodge(), aes(y=deltaCT.mean), stat="identity") +
  geom_errorbar(position=position_dodge(0.9), colour="black",size=0.5,width=.5)+
  scale_fill_manual(values = c('#999999','#E69F00'))+
  theme_classic()

#gene expression
qRe$Exp <- 2^(-qRe$deltaCT)
ExpSD <- summaryBy(Exp~Detector+Sample,qRe,FUN=std.error)
Expmean <- summaryBy(Exp~Detector+Sample,qRe,FUN=mean)
qRe <- left_join(qRe,Expmean,by=c("Detector" = "Detector","Sample"="Sample"))
qRe <- left_join(qRe,ExpSD,by=c("Detector" = "Detector","Sample"="Sample"))


#top 10 
top_id <- c("M2145","M2272","M2379","M3237","M4191","M4853","M5561","M6104","M6914","M761")
top_qRe <- filter(qRe,Detector %in% top_id)
top_qRe$Detector <- factor(top_qRe$Detector,top_id)
#top 5
top_5_id <- c("M4853","M2379","M6914","M4191","M2145")
top_qRe <- filter(qRe,Detector %in% top_5_id)
top_qRe$Detector <- factor(top_qRe$Detector,top_5_id)

#non_poly
non_poly_id <- c("M4853","M5561", "M5214","M55")
non_poly_qRe <- filter(qRe,Detector %in% non_poly_id)
non_poly_qRe$Detector <- factor(non_poly_qRe$Detector,non_poly_id)

ggplot(top_qRe,aes(x=Detector,ymin=Exp.mean-Exp.std.error, ymax=Exp.mean+Exp.std.error,fill=Sample)) +
  geom_bar(position=position_dodge(), aes(y=Exp.mean), stat="identity") +
  geom_errorbar(position=position_dodge(0.9), colour="black",size=0.5,width=.5)+
  scale_fill_manual(values = c('#999999','#E69F00'))+
  labs(x=NULL,tag="C")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.15, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))
ggsave(file="top_10_qPCR.pdf",width=210*0.5,height=297/3.5,units="mm")


tre <- sapply(unique(top_qRe$Detector),function(x){t.test(filter(top_qRe,Detector==x&Sample=="Mycelium")$Ct,
                                                y = filter(top_qRe,Detector==x&Sample=="Yeast")$Ct)},simplify=T)
tre[,9][3]
unique(top_qRe$Detector)


ggplot(non_poly_qRe,aes(x=Detector,ymin=Exp.mean-Exp.std.error, ymax=Exp.mean+Exp.std.error,fill=Sample)) +
  geom_bar(position=position_dodge(), aes(y=Exp.mean), stat="identity") +
  geom_errorbar(position=position_dodge(0.9), colour="black",size=0.5,width=.5)+
  scale_fill_manual(values = c('#999999','#E69F00'))+
  labs(x=NULL,tag="D")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.85, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))
ggsave(file="non_poly_qPCR.pdf",width=210/2,height=297/3.5,units="mm")

####################################################################################
##RNA-seq and qPCR
TPMData_rD_g <- TPMData_g[,c(1:6,10)]
head(TPMData_rD_g)
#non_polyA id: lnc_vali 

non_polyA <- geneTPMmean(TPMData_rD_g,lnc_vali,non_poly_qRe)
#top 10 id: top10gene
top5gene <- c("MSTRG.4853","MSTRG.2379","MSTRG.6914","MSTRG.4191","MSTRG.2145")
top_10 <- geneTPMmean(TPMData_rD_g,top10gene,top_qRe)
top_10_log <- geneTPMmean(TPMData_rD_g,top10gene,top_qRe)
top_5 <- geneTPMmean(TPMData_rD_g,top5gene,top_qRe)
dim(top_5)

ggplot(top_10_log,aes(x=exp,y=AvgCt))+
  geom_point(aes(color=Detector,shape = phase),size=2)+
  geom_smooth(method='lm',formula = y ~ x,color="black",fill = "wheat3")+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., ..p.value.label.., sep = "~`,`~")), 
               parse = TRUE) +
  labs(x="log2(TPM)",y="Ct",tag="D")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"))

ggsave(file="lm_R_qPCR.pdf",width=210/2,height=297/3.5,units="mm")

####################################################################################
#function
geneTPMmean <- function(TPMData_rD_g,id,qPCR){
  data <- filter(TPMData_rD_g,gene_id %in% id)
  data$Mexp <- apply(data[,2:4],1,mean)
  data$Yexp <- apply(data[,5:7],1,mean)
  exp <- data.frame("gene_id"=rep(data$gene_id,2),
                    "exp"=c(log2(data$Mexp),log2(data$Yexp)),
                    "phase"=c(rep("Mycelium",length(id)),rep("Yeast",length(id))))
  exp$Detector <- sapply(exp$gene_id, function(x){paste("M",strsplit(x,split = "\\.")[[1]][2],sep = "")})
  qtr <- qPCR[,c(1:2,5)]
  re <- left_join(exp,qtr,c("phase"="Sample","Detector"="Detector"))
  re <- unique(subset(re,select=c(gene_id,exp,phase,Detector,AvgCt)))
  na.omit(re)
}
