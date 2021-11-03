library(RColorBrewer)
display.brewer.all()
palcol = brewer.pal(11, "Spectral")

head(rD_gene_tpm_log10)
head(hk_gene_tpm_log10)

pcaData = cbind(rD_gene_tpm_log10[,1:3], hk_gene_tpm_log10[,1:4])
head(pcaData)

df_pca <- prcomp(pcaData, scale. = T, retx = T, center = T)
df_pca_sam <- data.frame(df_pca$rotation, gruops = rownames(df_pca$rotation))  

# gene
# df_pcs <- data.frame(df_pca$x, gruops = rownames(df_pca$x))  

ggplot(df_pca_sam, aes(x=PC1, y=PC2, color = gruops))+
  geom_point(size = 4, alpha = 0.5) + 
  theme_classic()

library(FactoMineR)
res.pca <- PCA(pcaData, graph = TRUE)
eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
var <- get_pca_var(res.pca)
var
fviz_pca_var(res.pca, col.var = "black")

library(corrplot)
corrplot(var$cos2, is.corr=FALSE)
fviz_pca_var(res.pca, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
ind <- get_pca_ind(res.pca)
ind

fviz_pca_ind(res.pca, geom.ind = "point", 
             pointshape = 21, fill.ind = palcol[1], 
             alpha.ind = 0.3,
             repel = FALSE)


fviz_pca_ind(res.pca,geom.ind = "point",
             col.ind = rownames(pcaData), # color by groups
             addEllipses = TRUE, 
             legend.title = "Groups"
)

### Extended Data Fig
# setwd("C:/Users/16103/Desktop/lncRNA/_supplymentary")
# rna ratio
rro = read_csv("rRNA_ratio_rD_hk.csv")
rbar = rro[7:8,] %>% 
  column_to_rownames(var = "ratio") %>% t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(group = c(rep("rProbe", 6), rep("Kit", 4)))

rbar = rbar %>% 
  pivot_longer(cols = c(rRNA, subrRNA),
               names_to = "lib",
               values_to = "r3")

rbar$sample[rbar$sample == "Y3"] = "RY3"
rbar$sample[rbar$sample == "SRR5028789"] = "PY.3"
rbar$sample[rbar$sample == "SRR5028790"] = "PY.4"
rbar$sample[rbar$sample == "SRR6435879"] = "PY.2"
rbar$sample[rbar$sample == "SRR6435878"] = "PY.1"

rbar$sample = factor(rbar$sample, level = c("C1", "C2", "EC3", "EY3", "EY4", "RY3",
                                         "PY.1", "PY.2", "PY.3", "PY.4"))

ggplot(rbar, aes(x = sample, y = r3, fill = lib))+
  geom_bar(stat = "identity",
           position = "dodge",
           width = 0.6)+
  scale_fill_manual(values = palcol[c(2,10)]) +
  geom_hline(yintercept = 1,colour = "grey11",linetype="dashed")+
  labs(x=NULL, y="rRNA residual ratio (%)")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position=c(0.1, 0.85),
        legend.text = element_text(size=6,color = "black"),
        legend.title = element_text(size=8),
        legend.key.size = unit(11, "pt"))
ggsave(file="figure/revise/EDFig_bar.pdf",width=210/1.8,height=297/4,units="mm")

# transcriptome
subrD = read_csv("revise/subExp/rD_gene_tpm_matrix.csv") %>% column_to_rownames(var = "gene_id")
subrD_log10 = log10(subrD + 1)
subhk = read_csv("revise/subExp/hk_gene_tpm_matrix.csv") %>% column_to_rownames(var = "gene_id")
subhk_log10 = log10(subhk + 1)

subrD_log10$mean = apply(subrD_log10, 1, mean)
subhk_log10$mean.1 = apply(subhk_log10, 1, mean)


sub_log10 <- cbind(subrD_log10, subhk_log10) %>%
  filter(mean > log10(2), mean.1 > log10(2))


sub_log10$trend <- sapply((sub_log10$mean-sub_log10$mean.1), 
                             function(x){if(x>1) 'up' else if(x<(-1)) 'down' else 'none'})
#labeled point
up <- filter(sub_log10, mean-mean.1>2)
up$gene <- rownames(up)
down <- filter(sub_log10, mean-mean.1<(-2))
down$gene <- rownames(down)
if (nrow(down)>3){
  down$minus <- down$mean-down$mean.1
  down <- down[order(down$minus),]
  down$gene[4:nrow(down)] <- ""    
}
#lab <- c("Probe","Ribo-Zero Kit","poly(A) selection")
my.formula <- y ~ x
lm.scattor <- lm(sub_log10$mean.1 ~ sub_log10$mean)
#summary(lm.scattor)
ggplot(sub_log10,aes(x=mean,y=mean.1,color=trend))+
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
  labs(y=paste("log10(TPM+1) of Kit sublibrary"),x=paste("log10(TPM+1) of rProbe sublibrary "))+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=10),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        legend.position="none")

ggsave(file="figure/revise/sub_rD_hk_lm.pdf",width=210/2,height=210/2.3,units="mm")
