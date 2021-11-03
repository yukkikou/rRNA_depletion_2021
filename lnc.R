#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA") 
#do atfer DEseq2.R
#lnc transcript expression matrix

fi_rd <- read_lines("lncFilter/final_intersect_lncRNA.id")
#length(fi_rd) #514
fi_tmp1 <- read_lines("lncFilter/lnc_tmp1_id.txt")
#length(fi_tmp1) #229
lnc_noval <- read_lines("lncFilter/lnc_noval_id.txt")
#length(lnc_noval) #22
lnc_noval[lnc_noval %in% fi_tmp1]
#[1] "MSTRG.4830.1" "MSTRG.5214.1" "MSTRG.5214.3" "MSTRG.5561.4" "MSTRG.4853.1" "MSTRG.55.1"
non_pA_lnc <- filter(ixd_matrix_0,rownames(ixd_matrix) %in% lnc_noval[lnc_noval %in% fi_tmp1])
lnc_vali <- c("MSTRG.5214","MSTRG.5561","MSTRG.4853","MSTRG.55")

#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA/lncFilter") 
lnc_tmp1_res <- DElnc(dds_t,fi_tmp1,res_t,"lncRNA")
summary(lnc_tmp1_res) #115/229
#summary(lnc_res) #158/456
head(lnc_tmp1_res)

##########################
resOrdered <- lnc_tmp1_res[order(lnc_tmp1_res$padj),]
resSig <- subset(resOrdered,padj < 0.05)
resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
resSigUp <- subset(resSig,log2FoldChange > 1)
resSigDown <- subset(resSig,log2FoldChange < (-1))


table(lnc_noval %in% rownames(resOrdered)) #5
lnc_noval[lnc_noval %in% rownames(resSigUp)]
#[1] "MSTRG.5214.1" "MSTRG.5214.3" "MSTRG.5561.4" "MSTRG.4853.1" "MSTRG.55.1"

plot.new()
specGene(dds_g,"MSTRG.5214","#2b90d9","#f26d5b")
specGene(dds_g,"MSTRG.55","#2b90d9","#f26d5b")
specGene(dds_g,"MSTRG.5261","#2b90d9","#f26d5b")
specGene(dds_g,"MSTRG.4853","#2b90d9","#f26d5b")

t_52141 <- specGene(dds_t,"MSTRG.5214.1","#2b90d9","#f26d5b")
t_52143 <- specGene(dds_t,"MSTRG.5214.3","#2b90d9","#f26d5b")
t_55614 <- specGene(dds_t,"MSTRG.5561.4","#2b90d9","#f26d5b")
t_48531 <- specGene(dds_t,"MSTRG.4853.1","#2b90d9","#f26d5b")
t_551 <- specGene(dds_t,"MSTRG.55.1","#2b90d9","#f26d5b")

(t_52141+t_52143+t_55614+t_48531+t_551+plot_spacer())
ggsave(file="figure/final/lnc_m_y.pdf",width=210*0.9,height=297/3,units="mm")

#####################################
#top10
top10 <- rownames(resSigAll[1:10,])
top10gene <- c("MSTRG.4853","MSTRG.2379","MSTRG.6914","MSTRG.4191",
               "MSTRG.2145","MSTRG.5561","MSTRG.3237","MSTRG.6104",
               "MSTRG.761","MSTRG.2272")
write.csv(resSigAll[1:10,],"lncFilter/top10/DE-lncRNA-top10.csv")
M4853 <- specGene(dds_t,"MSTRG.4853.1","#2b90d9","#f26d5b")
M2379 <- specGene(dds_t,"MSTRG.2379.1","#2b90d9","#f26d5b")
M6419 <- specGene(dds_t,"MSTRG.6914.5","#2b90d9","#f26d5b")
M4191 <- specGene(dds_t,"MSTRG.4191.1","#2b90d9","#f26d5b")
M2145 <- specGene(dds_t,"MSTRG.2145.1","#2b90d9","#f26d5b")
M5561 <- specGene(dds_t,"MSTRG.5561.1","#2b90d9","#f26d5b")
M3237 <- specGene(dds_t,"MSTRG.3237.3","#2b90d9","#f26d5b")
M6104 <- specGene(dds_t,"MSTRG.6104.1","#2b90d9","#f26d5b")
M761 <- specGene(dds_t,"MSTRG.761.1","#2b90d9","#f26d5b")
M2272 <- specGene(dds_t,"MSTRG.2272.1","#2b90d9","#f26d5b")

(M4853+M2379+M6419+M4191+M2145)
ggsave(file="lncFilter/top10/top5.pdf",width=210*0.9,height=297/3,units="mm")
(M4853+M2379)/(M6419+M4191)/(M2145+M5561)/(M3237+M6104)/(M761+M2272)
ggsave(file="lncFilter/top10/top10.pdf",width=210*0.6,height=297*0.6,units="mm")

#####################################
#output
output(lnc_tmp1_res,"lncRNA",dds_t)
write.csv(non_pA_lnc,"lncFilter/non_polyA_lnc.csv")

#####################################
#function
DElnc <- function(dds,fi,res,title){
  #results subset
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < 0.05)
  resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
  lnc_res <- resSigAll[rownames(resSigAll) %in% fi,]
  
  #normalise expr
  rldvst <- vst(dds,blind = F)
  expr_vst <- assay(rldvst)
  lnc_exp <- expr_vst[rownames(expr_vst) %in% rownames(lnc_res),]
  
  #heatmap
  z_score_matrix <- t(scale(t(lnc_exp)))
  annotation_col <- data.frame(Sample = condition)
  rownames(annotation_col) <- salist
  pdf(paste(title,"_rlts.pdf",sep=""),width=6,height=6)
  pheatmap(mat = z_score_matrix, cluster_rows = T, show_colnames = T, cluster_cols = T,
           annotation_col = annotation_col, 
           annotation_colors = list(Sample=c(yeast = 'red', mycelium = 'blue')), 
           show_rownames=F,  border_color=NA,
           legend = T, legend_breaks = c(-2, 0, 2), breaks = unique(seq(-2, 2, length=100)), 
           cutree_cols = 2, cutree_rows = 2, lagend_labels = c('≤2','0','≥2'),
           main=paste('heat map of DE ',title))
  dev.off()
  lnc_res
}

output <- function(res,title,dds){
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < 0.05)
  resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
  resSigUp <- subset(resSig,log2FoldChange > 1)
  resSigDown <- subset(resSig,log2FoldChange < (-1))
  resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
  write.csv(resdata,file= paste(title,"_rlt_total.csv",sep=""),row.names = T)
  write.csv(resSigAll,file= paste(title,"_rlt_sig.csv",sep=""),row.names = T)
  write.csv(resSigUp,file= paste(title,"_rlt_up.csv",sep=""),row.names = T)
  write.csv(resSigDown,file= paste(title,"_rlt_down.csv",sep=""),row.names = T)
}




