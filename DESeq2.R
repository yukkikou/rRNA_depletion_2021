#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("tximportData")
#BiocManager::install("ashr")  
library(zoo)
library(sm)
library(pheatmap)
library(vioplot)
library(dplyr)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(parallel)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(IDPmisc)


################################################read#data###############################################
#transcript matrix: countData
countData_t <- read.csv("comparison/rD_tA/transcript_count_matrix.csv")
dim(countData_t)
#[1] 46914    13
countData_t <- countData_t[,c(1:6,10)]
head(countData_t)


#gene matrix
countData_g <- read.csv("comparison/rD_tA/gene_count_matrix.csv")
dim(countData_g)
#[1] 8127   13
countData_g <- countData_g[,c(1:6,10)]
head(countData_g)


#phenotype
salist <- c("C1","C2","EC3","EY3","EY4","RY3")
phase <- c(rep("mycelium",3),rep("yeast",3))
condition <- factor(c("mycelium","mycelium","mycelium","yeast","yeast","yeast"),levels = c("mycelium","yeast"))

#############################################DESeqobject
dds_g <- DESobject(countData_g)
resultsNames(dds_g)
dds_t <- DESobject(countData_t)
resultsNames(dds_t)
#threshold
#M2Y
res_g <- results(dds_g, contrast=c("condition","yeast","mycelium"), alpha=0.05)
res_t <- results(dds_t, contrast=c("condition","yeast","mycelium"), alpha=0.05)
#without foldchange
summary(res_g)
statUD(res_g,0.05,1)
statUD(res_g,0.01,3)
head(res_g[order(res_g$padj),])

plot.new()
specGene(dds_g,"MSTRG.6412","#2b90d9","#f26d5b")


summary(res_t)
statUD(res_t,0.05,1)
statUD(res_t,0.01,3)
head(res_t[order(res_t$padj),])
specGene(dds_t,"MSTRG.2600.1","#2b90d9","#f26d5b")

##############################################visual
#PCA
rld <- rlog(dds_g)
pdf("gene_rlts.pdf",width=6,height=6)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotPCA(rld, intgroup="condition")
countsboxplot(dds_g,"(gene)")
volcano(res_g,"(gene)")
fallrank(dds_g,"(gene)")
fallplot(res_g,"(gene)")
heatmap(res_g,dds_g,"(gene)")
hdenscatter(dds_g, "(gene)")
dev.off()

rld_t <- rlog(dds_t)
pdf("transcript_rlts.pdf",width=6,height=6)
plotCounts(dds_t, gene=which.min(res_t$padj), intgroup="condition")
plotPCA(rld_t, intgroup="condition")
countsboxplot(dds_t,"(transcript)")
volcano(res_t,"(transcript)")
fallrank(dds_t,"(transcript)")
fallplot(res_t,"(transcript)")
heatmap(res_t,dds_t,"(transcript)")
hdenscatter(dds_t, "(transcript)")
dev.off()

#############################################output to excel
output(res_g,"gene",dds_g)
output(res_t,"transcript",dds_t)


######################################################################################################
#function
DESobject <- function(countData){
  #condition
  rownames(countData) <- countData[,1]
  countData <- countData[,-1]
  condition <- factor(phase,levels = c("yeast","mycelium"))
  colData <- data.frame(row.names=colnames(countData), condition)
  #check same
  #all(rownames(colData) %in% colnames(countData))
  #all(rownames(colData) == colnames(countData))
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design=~ condition)
  #filteringt row full of zero
  dds <- dds[ rowSums(counts(dds)) > 1, ] #47360
  #or filter gene only in one samples
  #dds <- dds[ rowVars(texpr(dds)) > 1, ]
  #View(counts(dds))
  #normalization:  median of ratios method
  dds <- DESeq(dds)
  dds
}

volcano <- function(res,title){
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < 0.05)
  resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
  resSigUp <- subset(resSig,log2FoldChange > 1)
  resSigDown <- subset(resSig,log2FoldChange < (-1))
  
  
  #volcano plot
  #ggplot2
  for_volcano <- data.frame('log2FoldChange'=res$log2FoldChange, 
                            'padj'=res$padj, 
                            'trend' = rep('no', length(res$log2FoldChange))) 
  up_sig_indices <- intersect(which(for_volcano$log2FoldChange > 1), which(for_volcano$padj < 0.05)) 
  down_sig_indices <- intersect(which(for_volcano$log2FoldChange < -1), which(for_volcano$padj < 0.05)) 
  for_volcano[up_sig_indices,'trend'] <- 'up'
  for_volcano[down_sig_indices, 'trend'] <- 'down'
  for_volcano$trend <- as.factor(for_volcano$trend)
  for_volcano$padj <- -log10(for_volcano$padj)
  p <- ggplot(for_volcano, aes(x=log2FoldChange, y=padj, colour=trend))+ 
    geom_point(size=I(1.0))+ 
    scale_color_manual(values = c('no'='#6E7783', 'up'="darkred", 'down'="darkblue"))+ 
    geom_vline(xintercept = c(1, -1), lty=2, size=I(0.4), colour = 'grey11')+ 
    geom_hline(yintercept = c(-log(x=0.05, base=10)), lty=2, size=I(0.4), colour = 'grey11')+ 
    theme_bw()+ 
    theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
          panel.grid = element_blank())+ 
    labs(x='log2FoldChange', y='-log10Pvalue',title = paste("Mycelium to Yeast",title))+
    theme(plot.title = element_text(hjust = 0.5))
  p
}

fallrank <- function(dds,title){
  #fall rank
  pre_ranked_all_genes <- as.data.frame(results(dds, alpha = 0.9999,
                                                contrast = c("condition","mycelium","yeast")))
  pre_ranked_all_genes <- pre_ranked_all_genes[order(pre_ranked_all_genes$log2FoldChange, decreasing = F),]
  trend_all <- sapply(pre_ranked_all_genes$log2FoldChange, function(x){if(x>0) 'up' else 'down'})
  pre_ranked_all_genes <- data.frame('gene_names'=rownames(pre_ranked_all_genes), pre_ranked_all_genes,
                                     'trend'=trend_all,'rank'=1:nrow(pre_ranked_all_genes),
                                     'hits'='medium',stringsAsFactors = F)
  pre_ranked_all_genes[1:5,'hits'] <- 'bottom5'
  pre_ranked_all_genes[(nrow(pre_ranked_all_genes)-4):nrow(pre_ranked_all_genes), 'hits'] <- 'top5'
  pre_ranked_all_genes$padj <- -log10(pre_ranked_all_genes$padj)
  set.seed(12)
  to_be_pointed_out_all <- pre_ranked_all_genes[c(1:5, (nrow(pre_ranked_all_genes)-4):nrow(pre_ranked_all_genes)),]
  fall_plot <- ggplot(pre_ranked_all_genes, aes(x=rank, y=log2FoldChange, color=hits))+
    geom_point(aes(size=padj))+ 
    geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+ 
    geom_hline(yintercept = c(0), linetype=1, size=0.5)+ 
    geom_vline(xintercept = sum(pre_ranked_all_genes$trend == 'down')+0.5, 
               linetype=2, size=0.25)+ 
    ggrepel::geom_text_repel(inherit.aes = F, data = to_be_pointed_out_all, aes(x=rank, y=log2FoldChange, label=gene_names, color=hits), 
                             size=3, direction = 'y', xlim = c(1000, 15000))+ 
    scale_color_manual(values=c('bottom5'='#274c5e', 'top5'='red3', 'medium'='#a3a1a1'))+ 
    scale_size_continuous(name='-log10(FDR)')+ 
    scale_y_continuous(breaks=c(-30, -10, -5, -2, 0, 2, 5, 10, 30))+ 
    xlab('rank of differentially expressed genes')+ 
    theme_bw()+ theme(panel.grid = element_line(color='white'), 
                      legend.title.align = 0.5)+
    labs(title = paste('Yeast to Mycelium',title))+
    theme(plot.title = element_text(hjust = 0.5))
  fall_plot
}

fallplot <- function(res,title){
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < 0.05)
  resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
  #fall_plot
  #fall only DEG
  #draw the foldchange plot with specific annotations.
  pre_ranked_sig_genes <- as.data.frame(resSigAll)[,c(2,6)] %>% data.frame('gene_names'=resSigAll@rownames, stringsAsFactors = F)
  pre_ranked_sig_genes <- pre_ranked_sig_genes[order(pre_ranked_sig_genes$log2FoldChange, decreasing = T),]
  trend <- sapply(pre_ranked_sig_genes$log2FoldChange, function(x){if(x>0) 'up' else 'down'})
  pre_ranked_sig_genes <- data.frame(pre_ranked_sig_genes, 'trend'=trend, 'rank'=1:nrow(pre_ranked_sig_genes), stringsAsFactors = F)
  #randomly select 5 genes to be annotated in the plot.
  #to_be_pointed_out <- pre_ranked_sig_genes[sample(1:nrow(pre_ranked_sig_genes), 5),]
  set.seed(12)
  to_be_pointed_out <- pre_ranked_sig_genes[c(1:3, (nrow(pre_ranked_sig_genes)-2):nrow(pre_ranked_sig_genes)),]
  fall_sig_plot <- ggplot(pre_ranked_sig_genes, aes(x=rank, y=log2FoldChange, color=log2FoldChange))+ 
    geom_point(size=1)+ geom_hline(yintercept = c(2,-2), linetype=2, size=0.25)+ 
    geom_hline(yintercept = c(0), linetype=1, size=0.5)+ 
    #geom_vline(xintercept = 350, linetype=2, size=0.25)+ 
    scale_color_gradient2(low="navy", high="firebrick3", mid="white", midpoint=0)+ 
    geom_point(inherit.aes = F, data=to_be_pointed_out, aes(x=rank, y=log2FoldChange), 
               size = 3, color='black')+ 
    geom_point(inherit.aes = F, data=to_be_pointed_out, aes(x=rank, y=log2FoldChange), 
               size = 2, color='#71226e')+ 
    ggrepel::geom_text_repel(inherit.aes = F, data = to_be_pointed_out, 
                             aes(x=rank, y=log2FoldChange, label=gene_names), size=3)+
    xlab('rank of differentially expressed genes')+ 
    theme_bw()+ theme(panel.grid = element_line(color='white'), legend.title.align = 0.5)+
    labs(title = paste('Yeast to Mycelium',title))+
    theme(plot.title = element_text(hjust = 0.5))
  fall_sig_plot
}

countsboxplot <- function(dds,title){
  #reads counts boxplot
  #show normalized reads counts distribution
  rldvst <- vst(dds,blind = F)
  expr_vst <- assay(rldvst)
  par(cex = 1.0)
  #margin
  par(mar=c(3,3,3,3))
  n.sample <- ncol(expr_vst)
  cols <- rainbow(n.sample*1.2)
  boxplot(expr_vst,col=cols, main=paste('expression value',title),las=2)  
}

heatmap <- function(res,dds,title){
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < 0.05)
  resSigAll <- subset(resSig, (log2FoldChange < (-1)|log2FoldChange > 1))
  rldvst <- vst(dds,blind = F)
  expr_vst <- assay(rldvst)
  #heat map
  DEG_expr_matr <- expr_vst[resSigAll@rownames,]
  #scale() will calculate the z-scores of each [column, that is gene] (that's why we use t function twice!).
  z_score_matrix <- t(scale(t(DEG_expr_matr)))
  annotation_col <- data.frame(Sample = condition)
  rownames(annotation_col) <- salist
  pheatmap(mat = z_score_matrix, cluster_rows = T, show_colnames = F, cluster_cols = T,
           annotation_col = annotation_col, 
           annotation_colors = list(Sample=c(yeast = 'red', mycelium = 'blue')), 
           show_rownames=F, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
           legend = T, legend_breaks = c(-2, 0, 2), breaks = unique(seq(-2, 2, length=100)), 
           cutree_cols = 2, cutree_rows = 2, lagend_labels = c('≤2','0','≥2'),
           main=paste('heat map of DE genes',title))
}

specGene <- function(dds,id,colorm,colory){
  rldvst <- vst(dds,blind = F)
  expr_vst <- assay(rldvst)
  #for interseted genes
  intended_gene <- id
  yeast_gene <- data.frame('values'=expr_vst[intended_gene, condition == 'yeast'],
                           'group'= rep('yeast', sum(condition == 'yeast')))
  mycelium_gene <- data.frame('values'=expr_vst[intended_gene, condition == 'mycelium'],
                              'group'= rep('mycelium', sum(condition == 'mycelium')))
  for_violin <- rbind(yeast_gene, mycelium_gene)
  #To reorder the rank of groups. Otherwise cancer group will be visualized first by alphabetical order.
  for_violin$group <- as.factor(for_violin$group) %>% relevel(ref = 'mycelium')
  pvo <- ggplot(for_violin, aes(x=group, y=values, fill=group))+ 
    geom_violin(position = position_dodge(width = 1), scale = 'width')+ 
    scale_fill_manual(values = c('mycelium'=colorm, 'yeast'=colory))+
    geom_boxplot(position = position_dodge(width = 1), 
                 outlier.size = 0.7, width = 0.2, show.legend = FALSE)+
    labs(x = NULL, y = 'normalized read counts',title = id)+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5),
          text=element_text(size=8),
          axis.text.x = element_text( hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position="none")
}


hdenscatter <- function(dds,title){
  rldvst <- vst(dds,blind = F)
  expr_vst <- assay(rldvst)
  #Draw the high-density scatterplot.
  intended_pixel_num = 70
  expr_for_visualization <- data.frame('mycelium'=rowMeans(expr_vst[,1:3]), 
                                       'yeast'=rowMeans(expr_vst[,4:6]), 
                                       stringsAsFactors = F)
  mycelium_notches <- seq(from = min(expr_for_visualization$mycelium-0.1),
                          to = max(expr_for_visualization$mycelium+0.1), 
                          length.out = (intended_pixel_num + 1)) 
  yeast_notches <- seq(from = min(expr_for_visualization$yeast-0.1),
                       to = max(expr_for_visualization$yeast+0.1), 
                       length.out = (intended_pixel_num + 1)) 
  mycelium_groups <- as.character(cut(expr_for_visualization$mycelium,
                                      breaks = mycelium_notches)) 
  yeast_groups <- as.character(cut(expr_for_visualization$yeast, 
                                   breaks = yeast_notches)) 
  expr_for_visualization <- data.frame('genes'=rownames(expr_for_visualization),
                                       expr_for_visualization, 
                                       'groups'=paste(mycelium_groups,yeast_groups, sep = '_'), stringsAsFactors = F)
  frequency <- as.data.frame(table(expr_for_visualization$groups))
  frequency$Var1 <- as.character(frequency$Var1)
  colnames(frequency) <- c('groups','freq')
  expr_for_visualization <- merge(expr_for_visualization, frequency, by='groups', all.x = T)
  high_density_scatterplot <- ggplot(expr_for_visualization, aes(x=mycelium, y=yeast, colour=freq))+ 
    geom_point(size=0.1)+ 
    scale_color_gradientn(colours=c('black','blue','yellow','red'))+ 
    ggtitle(paste('correlation scattorplot with density',title))+ 
    theme_bw()+theme(panel.grid = element_line(colour = 'white'), 
                     plot.title = element_text(hjust = 0.5)) 
  high_density_scatterplot
}

statUD <- function(res,p,trd){
  resOrdered <- res[order(res$padj),]
  resSig <- subset(resOrdered,padj < p)
  resSigAll <- subset(resSig, (log2FoldChange < (-trd)|log2FoldChange > trd))
  resSigUp <- subset(resSig,log2FoldChange > trd)
  resSigDown <- subset(resSig,log2FoldChange < (-trd))
  statAll <- data.frame("All"=dim(resSigAll)[1],"Up"=dim(resSigUp)[1],"Down"=dim(resSigDown)[1])
  statAll
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
