#based on david results txt
library(ggplot2)
library(patchwork)
library(stringr)

#pvalue

BP <- read.csv("GO/DOWN_GOTERM_BP_DIRECT.csv")
CC <- read.csv("GO/DOWN_GOTERM_CC_DIRECT.csv")
MF <- read.csv("GO/DOWN_GOTERM_MF_DIRECT.csv")
KEGG <- read.csv("GO/DOWN_KEGG_PATHWAY.csv")

colnames(BP)
View(BP)
table(sort(BP$Count))
table(sort(CC$Count))
table(sort(MF$Count))


BP[, c("GO_ID", "GO_Term")] <- str_split_fixed(BP$Term, "~", 2)
CC[, c("GO_ID", "GO_Term")] <- str_split_fixed(CC$Term, "~", 2)
MF[, c("GO_ID", "GO_Term")] <- str_split_fixed(MF$Term, "~", 2)
#UP
#MF$GO_Term[4] <-  "transcriptional activator activity"
KEGG[, c("KEGG_ID", "KEGG_Term")] <- str_split_fixed(KEGG$Term, ":", 2)

bubbleGO <- function(res,num=10,title="GO term BP"){
  res <- res[c(1:num),]
  p = ggplot(res,aes(Fold.Enrichment,GO_Term))
  p = p + geom_point()
  p= p + geom_point(aes(size=Count))
  pbubble = p + geom_point(aes(size=Count,color=-1*log10(PValue)))
  pbubble = pbubble + scale_colour_gradient(low="blue",high="red")
  pr = pbubble + scale_colour_gradient(low="blue",high="red") + 
    labs(color=expression(-log[10](PValue)),size="Gene counts",
         x="Fold Enrichment",y="GO Term",title=title)
  pr = pr + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  pr
}

bubbleKEGG <- function(res,num=10,title="KEGG pathway enrichment"){
  res <- res[c(1:num),]
  p = ggplot(res,aes(Fold.Enrichment,KEGG_Term))
  p = p + geom_point()
  p= p + geom_point(aes(size=Count))
  pbubble = p + geom_point(aes(size=Count,color=-1*log10(PValue)))
  pbubble = pbubble + scale_colour_gradient(low="blue",high="red")
  pr = pbubble + scale_colour_gradient(low="blue",high="red") + 
    labs(color=expression(-log[10](PValue)),size="Gene counts",
         x="Fold Enrichment",y="KEGG Term",title=title)
  pr = pr + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  pr
}

bp <- bubbleGO(BP,num=10,title="BP")
cc <-bubbleGO(CC,num=length(CC$Category),title="CC")
mf <-bubbleGO(MF,num=10,title="MF")
kegg <-bubbleKEGG(KEGG,num=length(KEGG$Category),title="KEGG")

#####################OUTPUT#############################
pdf(file = "DOWN_EnrichmentBubble.pdf")
bp / cc 
mf / kegg
dev.off()
