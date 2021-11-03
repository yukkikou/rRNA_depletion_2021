#setwd("C:/Users/16103/Desktop/trival/scripts/R/PM1_lncRNA/lncRNA_id")

library(tidyverse)
library(grid)
library(futile.logger)
library(VennDiagram)
library(RColorBrewer)
#install.packages("UpSetR")
#install.packages("venneuler")
library(UpSetR)
library(venneuler)

display.brewer.all()
cols <- brewer.pal(4, "Pastel2")

cpc2_rd <- read_lines("CPC2_id.txt") #930
cnci_rd <- read_lines("CNCI_id.txt") #1549
pfam_rd <- read_lines("Pfam_id.txt") #834
feelnc_rd <- read_lines("FEElnc_id.txt") #727

check(cpc2_rd)
check(cnci_rd)
check(pfam_rd)
check(feelnc_rd)

cc_rd <- intersect(cpc2_rd,cnci_rd) #698
pf_rd <- intersect(pfam_rd,feelnc_rd) #558
length(intersect(cpc2_rd,pfam_rd)) #382
length(intersect(cnci_rd,pfam_rd)) #646
length(intersect(cpc2_rd,feelnc_rd)) #400
length(intersect(cnci_rd,feelnc_rd)) #610

length(intersect(intersect(cpc2_rd,cnci_rd),pfam_rd)) #345


fi_rd <- intersect(intersect(cpc2_rd,cnci_rd),intersect(pfam_rd,feelnc_rd)) #246
length(fi_rd)

fi_all <- unique(c(cpc2_rd,cnci_rd,pfam_rd,feelnc_rd)) #1884
length(fi_all)


venn.plot <- draw.quad.venn(
  area1 = length(cpc2_rd),
  area2 = length(cnci_rd),
  area3 = length(pfam_rd),
  area4 = length(feelnc_rd),
  n12 = length(intersect(cpc2_rd,cnci_rd)),
  n13 = length(intersect(cpc2_rd,pfam_rd)),
  n14 = length(intersect(cpc2_rd,feelnc_rd)),
  n23 = length(intersect(cnci_rd,pfam_rd)),
  n24 = length(intersect(cnci_rd,feelnc_rd)),
  n34 = length(intersect(pfam_rd,feelnc_rd)),
  n123 = length(intersect(intersect(cpc2_rd,cnci_rd),pfam_rd)),
  n124 = length(intersect(intersect(cpc2_rd,cnci_rd),feelnc_rd)),
  n134 = length(intersect(intersect(cpc2_rd,pfam_rd),feelnc_rd)),
  n234 = length(intersect(intersect(cnci_rd,pfam_rd),feelnc_rd)),
  n1234 = length(fi_rd),
  category = c("CPC2", "CNCI", "PfamScan", "FEELnc"),
  fill = cols,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = cols
)
pdf(file = "../figure/final/lncVenn.pdf", width = 4.5, height = 4.5)
grid.draw(venn.plot)
dev.off()
#############################################################
#output
write(fi_rd,"final_intersect_lncRNA.id")
write(fi_all,"final_all_lncRNA.id")

pdf(file = "lnc_YC_venn.pdf",width = 10,height = 10)
grid.draw(venn.plot)
dev.off()

#############################################################
#funtion
check <- function(cpc2_pa){
  #return uniq length and true length
  un <- length(unique(cpc2_pa))
  l <- length(cpc2_pa)
  c(un,l)
}

############################################################
?upset()
# input <- c(
#   "CPC2"=1610,
#   "CNCI"=2702, 
#   "PfamScan"=1560, 
#   "FEELnc"=1489, 
#   "CPC2&CNCI"=1379, 
#   "CPC2&PfamScan"=734, 
#   "CPC2&FEELnc"=891, 
#   "CNCI&PfamScan"=1219, 
#   "CNCI&FEELnc"=1283, 
#   "PfamScan&FEELnc"=1124, 
#   "CPC2&CNCI&PfamScan"=674, 
#   "CPC2&CNCI&FEELnc"=833,
#   "CPC2&PfamScan&FEELnc"=556,
#   "CNCI&PfamScan&FEELnc"=944,
#   "CPC2&CNCI&PfamScan&FEELnc"=514
# )
# data <- fromExperssion(input)
# upset(data)

lnc_all <- read.csv("lncFilter/final_all_lncRNA.id",header = F)
cpc_idx <- lnc_all$V1 %in% cpc2_rd
cnci_idx <- lnc_all$V1 %in% cnci_rd
pf_idx <- lnc_all$V1 %in% pfam_rd
fl_idx <- lnc_all$V1 %in% feelnc_rd

dataplot <- data.frame("id"=lnc_all$V1,
                       "CPC2"=cpc_idx,
                       "CNCI"=cnci_idx,
                       "PfamScan"=pf_idx,
                       "FEELnc"=fl_idx)
head(dataplot)
dataplot_0 <- dataplot
dataplot_0[dataplot_0=="TRUE"]<-1
dataplot_0[dataplot_0=="FALSE"]<-0
head(dataplot_0)

plot1 <- function(mydata,x,y){
  myplot <- (ggplot(mydata,aes_string(x = x,y = y))
             + theme(plot.title = element_text(hjust = 0.5),
                     text=element_text(size=10),
                     axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
                     axis.text.y = element_text(size = 8,color = "black"),
                     legend.position="none"))
}
attributeplots <- list(gridrows = 100,
                       plots = list(list(plot = plot1)),
                       ncols = 1,nrows = 1)

upset(dataplot_0, nsets = 4,
      main.bar.color = "#D9AB42",sets.bar.color = "#D9AB42",
      mainbar.y.label = "Intersection of LncRNA",
      set_size.numbers_size = 6,sets.x.label = NULL,
      mb.ratio = c(0.7, 0.3),
      order.by = c("degree"), decreasing = F,line.size = 0.5,point.size = 1.5,
      attribute.plots = attributeplots)

#ggsave(file="lncRNA_id/vennbarpoint.pdf",width=210*0.9,height=297*0.9,units="mm")

  