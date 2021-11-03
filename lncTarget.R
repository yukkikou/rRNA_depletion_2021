##after DESeq2.R & lnc.R
#trans prediction
library(limma)
library(tidyverse)
#LncTar

corFilter=0.4
pvalueFilter=0.001


#gene matrix
#countData_g <- read.csv("gene_count_matrix.csv")
dim(countData_g) #[1] 7843    7
#filtering: > 1 in all samples
contg_fil <- filter_if(countData_g, is.numeric, all_vars((.) != 0))
rownames(contg_fil) <- contg_fil[,1]
contg_fil <- contg_fil[,-1]
dim(contg_fil) #[1]6677    6

#transcript matrix
#countData_t <- read.csv("transcript_count_matrix.csv")
dim(countData_t) #[1] 46914    13

#DE lnc matrix: different expression lnc
lncDE_t <- read.csv("lncRNA_rlt_sig.csv")
lncDE_count <- countData_t[countData_t$transcript_id %in% lncDE_t$X,]
rownames(lncDE_count) <- lncDE_count$transcript_id
lncDE_count <- lncDE_count[,-1] #158  6
dim(lncDE_count)

#transvert to log2
co_log2_g <- t(log2Exp(contg_fil))
co_log2_lnc <- t(log2Exp(lncDE_count))


Score <- cor(co_log2_lnc,co_log2_g,method="pearson")
pheatmap(Score)


######################################################################################
#cis prediction: FEELnc claasification
cis_lnc <- read.delim("lncRNA_classes.txt")
dim(cis_lnc)
cis_besdt_lnc <- subset(cis_lnc, cis_lnc$isBest==1)
dim(cis_besdt_lnc)
length(unique(cis_besdt_lnc$lncRNA_gene)) #932
length(cis_besdt_lnc$lncRNA_transcript) #1744 (all lncRNA)

cis_besdt_de_lnc <- subset(cis_besdt_lnc, cis_besdt_lnc$lncRNA_gene %in% diff_gene_deseq2)
dim(cis_besdt_de_lnc) #92 10
length(unique(cis_besdt_de_lnc$lncRNA_gene))

pie_labels <- names(table(cis_besdt_de_lnc$location))
pie_x <- table(cis_besdt_de_lnc$location)
pie_percent <- paste(round(100*pie_x/sum(pie_x)), "%")
pie_col <- c('#FFEEE4','#dddfe6','#d3e0f7','#fec9c9')
pie(pie_x, labels=pie_percent, radius = 1.0,
    main='location of lncRNA',col=pie_col,border = NA)
legend("topright", pie_labels, cex=0.8, fill=pie_col)


##################
write.csv(cis_besdt_de_lnc,file= "cis_best_de_lnc.csv",row.names = T)

######################################################
#function
log2Exp <- function(exp_matrix){
  dimnames <- list(rownames(exp_matrix),colnames(exp_matrix))
  data <- matrix(as.numeric(as.matrix(exp_matrix)),nrow=nrow(exp_matrix),dimnames=dimnames)
  data <- avereps(data)
  data <- log2(data+1)
}

