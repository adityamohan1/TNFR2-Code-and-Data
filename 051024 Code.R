setwd("/hpc/group/feccilab/aam131/TNFR2 KO Project RNASeq/031824 Analysis")
datatable <- read.csv("/hpc/group/feccilab/aam131/TNFR2 KO Project RNASeq/031824 Analysis/converted_counts_data (1)_Xannev2.csv")

###Setting DESeq----
library("org.Mm.eg.db")
library("clusterProfiler")
library("dplyr")
datatable$gene_id <- mapIds(org.Mm.eg.db, keys = gsub("\\..*","",datatable$gene_id),
                            column = "SYMBOL", keytype = "ENSEMBL")
datatable = datatable %>%
  group_by(gene_id) %>%
  summarise_all(sum) %>% data.frame()
countsmatrix<-datatable
countsmatrix = countsmatrix %>% mutate_if(is.numeric, round, digits=0)

ID <- c("WT_1","WT_2","WT_3","WT_4","WT_7","WT_8","KO_1","KO_2","KO_3","KO_4","KO_5","KO_6")
Group <- c("WT","WT","WT","WT","WT","WT","KO","KO","KO","KO","KO","KO")
directory <- data.frame(ID, Group)

var1 = "WT" # Define Groups
var1list = (directory %>% filter(Group == 'WT'))$ID

var2 = "KO" # Define Groups
var2list = (directory %>% filter(Group == 'KO'))$ID

#Mapping Table

NameList <- append(var1list, var2list)
idx <- match(NameList, names(countsmatrix))
NewDF <- countsmatrix[,idx]
NewDF <- na.omit(NewDF)
countsmatrix <- na.omit(countsmatrix)
rownames(NewDF) <-countsmatrix$Gene 

NameList <- append(var1list, var2list)
idx <- match(NameList, names(countsmatrix))
NewDF <- countsmatrix[,idx]
NewDF <- na.omit(NewDF)
countsmatrix <- na.omit(countsmatrix)
rownames(NewDF) <-countsmatrix$Gene 

# Mapping Table
MapTable = as.data.frame(matrix(nrow=length(NameList),ncol=2))
rownames(MapTable) = NameList
colnames(MapTable)[1] = "Group"
colnames(MapTable)[2] = "SampleID"
MapTable[2] = NameList
MapTable[1:length(var1list),1] = var1
MapTable[(length(var1list)+1):length(NameList),1] = var2

### Running DESEq

library(DESeq2)
NewDF <- NewDF[,unique(rownames(MapTable))]
rownames(NewDF) <-countsmatrix$gene_id 
all(colnames(NewDF) == rownames(MapTable))
MapTable$Group = factor(MapTable$Group,levels=c(var1,var2))
deseq2Data <- DESeqDataSetFromMatrix(countData=NewDF, colData=MapTable, design= ~ Group)
deseq2Data <- deseq2Data[rowSums(counts(deseq2Data)) > 2, ]
#Differential GeneExpression and Results
deseq2Data <- DESeq(deseq2Data)
deseq2Results <- results(deseq2Data)
deseq2Results
rm(ENSResultsSummary)
ENSResultsSummary<-as.data.frame(matrix(nrow=deseq2Results@nrows,ncol=7))
colnames(ENSResultsSummary)<-c("Gene","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
ENSResultsSummary$Gene<-deseq2Results@rownames
ENSResultsSummary$baseMean<-deseq2Results@listData$baseMean
ENSResultsSummary$log2FoldChange<-deseq2Results@listData$log2FoldChange
ENSResultsSummary$lfcSE<-deseq2Results@listData$lfcSE
ENSResultsSummary$stat<-deseq2Results@listData$stat
ENSResultsSummary$pvalue<-deseq2Results@listData$pvalue
ENSResultsSummary$padj<-deseq2Results@listData$padj 
GeneResultsSummary<- ENSResultsSummary

### Labels ----
GeneResultsSummary <- GeneResultsSummary %>% 
  mutate(
    Expression = case_when(log2FoldChange >= 0.5 & padj <= 0.01 ~ "Up-regulated",
                           log2FoldChange <= -0.5 & padj <= 0.01 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

top <- 20
top_genes <- bind_rows(
  GeneResultsSummary %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top),
  GeneResultsSummary %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(log2FoldChange))) %>% 
    head(top)
)

genes <- c("Pdcd1", "Slamf6", "Havcr2", "Tox", "Tox2", "Tox4", "Tnfrsf1b", "Tcf7", "Cd160",
           "Prdm1", "Ctla4", "Cd96", "Cd226", "Vsir", "Cd101", "Cd38", "Cd86", "Cd27",
           "Tnfrsf14", "Tnfsf9", "Tnfsf14", "Yy1", "Klrg1", "Tbx21", "Ptger4", "Stat1",
           "Eomes", "Foxo1", "Ets1", "Ep300", "Crebbp", "Irf4", "Rgs16", "Myc", "Satb1",
           "Mapk8", "Il7r", "Tcf4", "Atf3", "JunD", "JunB", "Fosb", "Fos", "Akt1", "Batf3",
           "Mafg", "Atf2", "Maf", "Mafb", "Fosl2", "Nfatc2", "Nfatc3", "Nfatc1", "Atf4",
           "Batf", "Mafa", "Klf2", "Gata3", "Cx3cr1", "Entpd1", "Klrk1", "Ccl5", "Stat5a",
           "Jun", "Myb", "Foxo3", "Klf4", "Il2ra", "Il1rl1", "Tnf", "Il12rb", "Tnfsf4",
           "Spp1", "Gzmf", "Tnfrsf25", "Ccr1", "Notch3", "Gzmd", "Cd70", "Cxcr6",
           "Tgfbr2", "Slamf7", "Nfat5")
genes_df <- data.frame(genes)


### Volcano Plot ----
library(ggplot2)
library(ggrepel)

library(ggpubr)
library(plyr) 
library(gridExtra) 
library(cowplot) 
theme_set(theme_cowplot())
theme(plot.title = element_text(hjust = 0.5,face = "plain"))
library(forcats) 
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library(ggbreak)
g2<- ggplot(GeneResultsSummary, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color = Expression), size = 1, alpha=0.8) +
  xlab("Log 2 Fold Change") +
  ylab("-Log10 BH Adjusted P-Value")+
  scale_color_manual(values = c("dodgerblue3", "gray50", "sienna2")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  geom_hline(yintercept=1.3, linetype="dashed",
             color = "black", size=0.5)+
  geom_vline(xintercept = -0.5, linetype="dotdash",
             color = "black", size=0.5)+
  geom_vline(xintercept = 0.5, linetype="dotdash",
             color = "black", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"))+
  theme(legend.position = "none")+
  geom_label_repel(data = top_genes,
                   mapping = aes(log2FoldChange, -log10(padj), label = top_genes$Gene),
                   size = 3)

### MAPlot ----
subset_data <- ENSResultsSummary[ENSResultsSummary$Gene %in% genes, ]


significance_threshold <- 0.05
ENSResultsSummary$Color <- ifelse(
  ENSResultsSummary$padj < significance_threshold & ENSResultsSummary$log2FoldChange > 0, "sienna2", 
  ifelse(ENSResultsSummary$padj < significance_threshold & ENSResultsSummary$log2FoldChange < 0, "dodgerblue3", 
         "gray50")
)

subset_data <- subset_data[subset_data$Color != "gray50", ]

ggplot(ENSResultsSummary, aes(x = baseMean, y = log2FoldChange, color = ENSResultsSummary$Color)) +
  geom_point(alpha = 0.6) +  # Adjust transparency
  scale_color_manual(values = c("dodgerblue3", "grey", "sienna2")) +
  geom_label_repel(data = subset_data,
                   mapping = aes(baseMean, log2FoldChange, label = subset_data$Gene),
                   size = 3)+
  labs(x = "Mean Counts in WT",
       y = "Log Fold Change")+
  scale_x_log10() +
  xlim(0, 1.5e5) +  # Limit x-axis to a maximum of 150,000
  ylim(-10, 10)+
  theme(legend.position = "none")

library(ggplot2)
library(ggrepel)

# Assuming ENSResultsSummary and subset_data are your data frames
ggplot(ENSResultsSummary, aes(x = baseMean, y = log2FoldChange, color = Color)) +
  geom_point(alpha = 0.6) +  # Adjust transparency
  scale_color_manual(values = c("dodgerblue3", "grey", "sienna2")) +
  geom_label_repel(data = subset_data,
                   mapping = aes(baseMean, log2FoldChange, label = Gene),
                   size = 3) +
  labs(x = "Mean Counts in WT",
       y = "Log Fold Change") +
  scale_x_log10() +
  xlim(0, 1.5e5) +  # Limit x-axis to a maximum of 150,000
  ylim(-10, 10) +
  theme(legend.position = "none")


# Create separate data frames for blue and red labels

highlighted_genes <- c("Tox", "Tox2", "Tox4", "Nfatc1", "Nfatc2", "Nfatc3", "Stat5a", "Batf3", "Jun", "Fox", "Akt1")

blue_labels <- subset_data[subset_data$Color == "dodgerblue3" & !subset_data$Gene %in% highlighted_genes, ]
red_labels <- subset_data[subset_data$Color == "sienna2" & !subset_data$Gene %in% highlighted_genes, ]
highlight_blue_labels <- subset_data[subset_data$Color == "dodgerblue3" & subset_data$Gene %in% highlighted_genes, ]
highlight_red_labels <- subset_data[subset_data$Color == "sienna2" & subset_data$Gene %in% highlighted_genes, ]

# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Create the plot with a logarithmic x-axis and specified limits
ggplot(ENSResultsSummary, aes(x = baseMean, y = log2FoldChange, color = ENSResultsSummary$Color)) +
  geom_point(alpha = 0.4) +  # Adjust transparency
  scale_color_manual(values = c("dodgerblue3", "gray50", "sienna2")) +
  geom_text_repel(data = blue_labels,
                  mapping = aes(baseMean, log2FoldChange, label = Gene),
                  size = 4,  # Increase label size
                  box.padding = 0.35,  # Add padding around labels
                  point.padding = 0.3,  # Add padding between points and labels
                  segment.color = 'dodgerblue3',  # Lighter color for line segments
                  segment.alpha = 0.5,  # Increase transparency of line segments
                  #direction = "y",  # Repel in the vertical direction
                  nudge_y = -5) +  # Nudge labels downwards
  geom_text_repel(data = red_labels,
                  mapping = aes(baseMean, log2FoldChange, label = Gene),
                  size = 4,  # Increase label size
                  box.padding = 0.35,  # Add padding around labels
                  point.padding = 0.3,  # Add padding between points and labels
                  segment.color = 'sienna2',  # Lighter color for line segments
                  segment.alpha = 0.5,  # Increase transparency of line segments
                  #direction = "y",  # Repel in the vertical direction
                  nudge_y = 5) +  # Nudge labels upwards
  geom_label_repel(data = highlight_blue_labels,
                   mapping = aes(baseMean, log2FoldChange, label = Gene),
                   size = 4,  # Increase label size
                   box.padding = 0.35,  # Add padding around labels
                   point.padding = 0.3,  # Add padding between points and labels
                   segment.color = 'dodgerblue3',  # Lighter color for line segments
                   segment.alpha = 0.5,  # Increase transparency of line segments
                   #direction = "y",  # Repel in the vertical direction
                   nudge_y = -2,  # Nudge labels downwards
                   fill = "white",  # Transparent fill for labels
                   color = 'black') +  # Blue border for highlighted blue labels
  geom_label_repel(data = highlight_red_labels,
                   mapping = aes(baseMean, log2FoldChange, label = Gene),
                   size = 4,  # Increase label size
                   box.padding = 0.35,  # Add padding around labels
                   point.padding = 0.3,  # Add padding between points and labels
                   segment.color = 'sienna2',  # Lighter color for line segments
                   segment.alpha = 0.5,  # Increase transparency of line segments
                   #direction = "y",  # Repel in the vertical direction
                   nudge_y = 2,  # Nudge labels upwards
                   fill = "white",  # Transparent fill for labels
                   color = 'black') +  # Red border for highlighted red labels
  labs(x = "Mean Counts in WT",
       y = "Log Fold Change") +
  scale_x_log10(limits = c(10^0.3, 10^5.5)) +  # Use a logarithmic scale for the x-axis with specified limits
  ylim(-10, 10) +
  theme(legend.position = "none")


ggplot(ENSResultsSummary, aes(x = baseMean, y = log2FoldChange, color = Color)) +
  geom_point(alpha = 0.4) +  # Adjust transparency
  scale_color_manual(values = c("dodgerblue3", "gray50", "sienna2")) +
  geom_text_repel(data = blue_labels,
                  mapping = aes(baseMean, log2FoldChange, label = Gene),
                  size = 4,  # Increase label size
                  box.padding = 0.35,  # Add padding around labels
                  point.padding = 0.3,  # Add padding between points and labels
                  segment.color = 'dodgerblue3',  # Lighter color for line segments
                  segment.alpha = 0.5,  # Increase transparency of line segments
                  # direction = "y",  # Repel in the vertical direction
                  nudge_y = -5) +  # Nudge labels downwards
  geom_text_repel(data = red_labels,
                  mapping = aes(baseMean, log2FoldChange, label = Gene),
                  size = 4,  # Increase label size
                  box.padding = 0.35,  # Add padding around labels
                  point.padding = 0.3,  # Add padding between points and labels
                  segment.color = 'sienna2',  # Lighter color for line segments
                  segment.alpha = 0.5,  # Increase transparency of line segments
                  # direction = "y",  # Repel in the vertical direction
                  nudge_y = 5) +  # Nudge labels upwards
  geom_label_repel(data = highlight_blue_labels,
                   mapping = aes(baseMean, log2FoldChange, label = Gene),
                   size = 4,  # Increase label size
                   box.padding = 0.35,  # Add padding around labels
                   point.padding = 0.3,  # Add padding between points and labels
                   segment.color = 'dodgerblue3',  # Lighter color for line segments
                   segment.alpha = 0.5,  # Increase transparency of line segments
                   # direction = "y",  # Repel in the vertical direction
                   nudge_y = -2,  # Nudge labels downwards
                   fill = "white",  # Transparent fill for labels
                   color = 'dodgerblue3',  # Text color for highlighted blue labels
                   fontface = "bold") +  # Bold text for highlighted blue labels
  geom_label_repel(data = highlight_red_labels,
                   mapping = aes(baseMean, log2FoldChange, label = Gene),
                   size = 4,  # Increase label size
                   box.padding = 0.35,  # Add padding around labels
                   point.padding = 0.3,  # Add padding between points and labels
                   segment.color = 'sienna2',  # Lighter color for line segments
                   segment.alpha = 0.5,  # Increase transparency of line segments
                   # direction = "y",  # Repel in the vertical direction
                   nudge_y = 2,  # Nudge labels upwards
                   fill = "white",  # Transparent fill for labels
                   color = 'sienna2',  # Text color for highlighted red labels
                   fontface = "bold") +  # Bold text for highlighted red labels
  labs(x = "Mean Counts in WT",
       y = "Log Fold Change") +
  scale_x_log10(limits = c(10^0.3, 10^5.5)) +  # Use a logarithmic scale for the x-axis with specified limits
  ylim(-10, 10) +
  theme(legend.position = "none")



#Generating a PCA Plot----
library(ComplexHeatmap)

vsd<-vst(deseq2Data,blind=FALSE)
data_for_pca <- assay(vsd)
pca_result <- prcomp(t(data_for_pca))
pca_data <- as.data.frame(pca_result$x)
pca_data$sample <- rownames(pca_data)
new_row_names <- c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5", "WT_6", "KO_1", "KO_2", "KO_3", "KO_4", "KO_5","KO_6")
rownames(pca_data) <- new_row_names
pca_data$sample <- rownames(pca_data)
pca_data$Color <- ifelse(pca_data$sample %in% c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5", "WT_6"), "WT", "KO")


#new_row_names <- c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5", "WT_6", "KO_1", "KO_2", "KO_3", "KO_4", "KO_5","KO_6")
#rownames(pca_data) <- new_row_names
#pca_data$sample <- rownames(pca_data)

g1<-ggplot(pca_data, aes(x = PC1, y = PC2, label = sample, colour = Color)) +
  geom_point() +
  geom_text_repel() +
  scale_color_manual(values = c("WT" = "maroon", "KO" = "blue")) +
  xlab(paste("PC1: ", round(pca_result$sdev[1], 2))) +
  ylab(paste("PC2: ", round(pca_result$sdev[2], 2)))+
  theme(legend.position = "none")
g1

#Fgsea ----


library(readr)
library(dplyr)
library(DESeq2)
library(readxl)
library(ggplot2)
library(tidyr) 
library(ggpubr)
library(plyr) 
library(gridExtra) 
library(cowplot) 
theme_set(theme_cowplot())
theme(plot.title = element_text(hjust = 0.5,face = "plain"))
library(forcats) 
library(ggrepel)
options(ggrepel.max.overlaps = Inf)
library("org.Hs.eg.db")
library(ReactomePA)
library(DOSE)
library("enrichplot")
library("clusterProfiler")
library("fgsea")
library("tibble")
library("msigdbr")
library(reshape2)
library(ComplexHeatmap)
library(edgeR)

pathways <- gmtPathways("MouseExhaustGMTx.gmt")
# If Error:
## gmt_file_path <- "FingerPrintModules.gmt"
## lines <- readLines(gmt_file_path, warn = FALSE)  # Set warn = FALSE to temporarily suppress the warning

### Check if the last line is an empty line, if not, add one
## if (length(lines) > 0 && nzchar(lines[length(lines)])) {lines <- c(lines, "")  # Add an empty line at the end }

### Write the lines back to the file
## writeLines(lines, gmt_file_path)

rankings <- sign(ENSResultsSummary$log2FoldChange)*(-log10(ENSResultsSummary$padj))
names(rankings) <- ENSResultsSummary$Gene 
rankings <- sort(rankings, decreasing = TRUE)
plot(rankings)
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 10)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 10)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

rankings2 <- sign(ENSResultsSummary$log2FoldChange)
names(rankings2) <- ENSResultsSummary$Gene 
rankings2 <- sort(rankings2, decreasing = TRUE) # sort genes by ranking


ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fgseaRes <- fgsea(pathways, rankings, minSize=4, maxSize=500)

#write.xlsx(fgseaRes, "FGSEA_Results.xlsx", sheetName = "Scores", row.names = FALSE)

fgseaRes2 <- fgseaRes
fgseaRes2 <- fgseaRes2[!pathway %in% c("Anergy", "Cytotoxicity", "DAP12 Signaling", "NK Activity", "NK Exhaustion", "TCR Signaling")]
library(openxlsx)

#write.xlsx(fgseaRes2, "FGSEA_Results2.xlsx", sheetName = "Scores", row.names = FALSE)

fgseaRes2$Color <- ifelse(fgseaRes2$padj < 0.1, 
                          ifelse(fgseaRes2$NES > 0, "orange", "blue"), 
                          "grey")


write.xlsx(fgseaRes2, "FGSEA_Results2.xlsx", sheetName = "Scores", row.names = FALSE)

library(readxl)
FGSEA_Results2 <- read_excel("FGSEA_Results.xlsx")
FGSEA_Results2 <- FGSEA_Results2[!(FGSEA_Results2$pathway %in% c("Anergy", "Cytotoxicity", "DAP12 Signaling", "NK Activity", "NK Exhaustion", "TCR Signaling")), ]
FGSEA_Results2$Color <- ifelse(FGSEA_Results2$padj < 0.1, 
                          ifelse(FGSEA_Results2$NES > 0, "orange", "blue"), 
                          "grey")

ggplot(FGSEA_Results2, aes(x = NES, y = reorder(pathway, NES), fill = Color)) +
  geom_bar(stat = "identity", width = 0.6) +  # Bars with width
  scale_fill_manual(values = c("dodgerblue3", "grey50")) +  # Define custom colors
  labs(x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme(legend.position = "none")+
  theme_minimal()+
  theme(legend.position = "none")


fgseaRes$GeneRatio <- 1

for(i in 1:(nrow(fgseaRes))) {        
  fgseaRes[i,9]<-(length(unlist(fgseaRes[i,8])))/fgseaRes[i,7]
}


fgseaRes <- fgseaRes %>% 
  mutate(
    Expression = case_when(NES >= 1 & padj <= 0.05 ~ "Up-regulated",
                           NES <= -1 & padj <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
top <- 6
top_Pathways <- bind_rows(
  fgseaRes %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(padj, desc(abs(NES))) %>% 
    head(top),
  fgseaRes %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(padj, desc(abs(NES))) %>% 
    head(top)
)


topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], rankings, fgseaRes, 
              gseaParam=0.5)



###Pathplot ----
pathplot = ggplot(fgseaRes, aes(x=NES, y=-log10(padj))) +
  geom_point(aes(color = Expression), size = (fgseaRes$GeneRatio*4)^2, alpha=0.4)+
  xlab("NES") +
  ylab("-Log10 BH Adjusted P-Value")+
  scale_color_manual(values = c( "gray50", "sienna2")) +
  guides(colour = guide_legend(override.aes = list(size=1.5)))+
  geom_hline(yintercept=1.3, linetype="dashed",
             color = "black", size=0.5)+
  geom_vline(xintercept = -1, linetype="dotdash",
             color = "black", size=0.5)+
  geom_vline(xintercept = 1, linetype="dotdash",
             color = "black", size=0.5)+
  theme(plot.title = element_text(hjust = 0.5,face = "plain"))+
  theme(legend.position = "none")+
  geom_label_repel(data = top_Pathways,
                   mapping = aes(NES, -log10(padj), label = top_Pathways$pathway),
                   size = 4.5)
pathplot




library(fgsea)
gmt_file <- "MouseExhaustGMT.gmt"
install.packages("GSA")
library(GSA)


gmt.file <- system.file("extdata", "MouseExhaustGMTx.gmt", package="fgsea")

gmt_file <-GSA.read.gmt(paste("MouseExhaustGMTx.gmt",sep = "/"))
pathways <- gmtPathways(gmt_file)
str(head(pathways))
fgseaRes <- fgsea(pathways, rankings, minSize=5, maxSize=500)



























#Exhaustion Heatmap
library(pheatmap)
Exhaustion<-c("Cd160","Prdm1","Ctla4","Entpd1","Cd96","Cd226","Vsir","Cd101","Cd38","Cd86","Cd27","Tnfrsf14","Tnfsf9","Tnfsf14","Ctla4","Dgka","Dgkz","Eomes","Fasl","Havcr2","Irf4","Pirb","Ptger4","Tox3","Yy1","Klrg1","Tbx21","Klf3","Ptger4","Tcf7","Irf1","Stat1","Eomes","Foxo1","Ets1","Ep300","Crebbp","Tox","Tox2","Tox4","Rgs16","Satb1","Il7r","Tcf4")
Checkpoints<-c("Havcr2","Icos","Icosl","Pdcd1","Pdcd1lg2","Pecam1","Pvrig","Tnfsf18","Tnfrsf4","Tnfsf4","Cd86","Fyn","Lck","Ppp2ca","Ppp2cb","Ppp2r1a","Ppp2r5a","Ppp2r5b","Ppp2r5c","Ptpn11","Yes1","Cd3g","H2-Aa","H2-D1","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-T23")
Anergy <- c("Dgka","Fasl","Gbp2b","Ldha","Ptprk","Ptprs","Rgs2","Slc7a5","Tle4","Atm","Btrc","Calm1","Ccnb1","Cdc16","Cdc25a","Cebpb","E2f2","Foxo1","Foxo3","Hipk1","Itpr1","Itpr2","Klf2","Phc3","Ppp1cb","Rbl2","Rheb","Rras","Sesn1","Sesn3","Smad2","Smad3","Tgfb1","Ube2d1","Vdac1")
AP1NFAT<-c("Nfatc1","Nfatc2","Nfatc3","Atf4","Batf","Cebpb","Fos","Jun","Junb","Jund","Mapk8","Akt1","Atf3","Fos","Fosb","Batf3","Mafg","Atf2","Mafg","Mafb","Fosl2")
Epigenetic <- c("Atm","Camk2d","Ccnb1","Cdyl","Crebbp","Dnmt3a","Gata3","Gnas","Gsk3b","Hdac7","Hdac9","Irf4","Jak2","Lef1","Lif","Naa50","Pml","Prkcb","Smad4","Taf6l","Tgfb1","Uhrf1","Usp15","Yy1")
NFKB <- c("Atm","Card11","Cflar","Cyld","Gadd45a","Icam1","Irak4","Lck","Malt1","Nfkb1","Nfkb2","Nfkbia","Plcg1","Prkcb","Syk","Ticam1","Tnfrsf1a","Xiap")

data2<-as.data.frame(data_for_pca)
colnames(data2) <- c("WT_1", "WT_2", "WT_3", "WT_4", "WT_5", "WT_6", "KO_1", "KO_2", "KO_3", "KO_4", "KO_5", "KO_6")

#data2 <- data2[, !(names(data2) %in% c("WT4", "KO4"))]
#rownames(data2) <- toupper(rownames(data2))
#Quiescence, Epigenetic, Exhaust, TF, Coinhibit, Costim, TF1, moreTF1, Effectors, TFNRF22TOX

filtered_Q_data <- data2[rownames(data2) %in% NFKB, ]
normalized_data <- t(scale(t(filtered_Q_data)))

# Create heatmap
Heatmap(normalized_data, 
        name = "Z-score", 
        cluster_columns = FALSE,
        show_row_names = TRUE, 
        show_column_names = TRUE)

