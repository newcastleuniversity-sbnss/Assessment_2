#Idenitfy current working directory
getwd()
#Download required data into desired working directory, then set the wokring directory to this location
setwd()

#Install and load packages
install.packages('BiocManager')
library(BiocManager)
install(c('tximport', 'DESeq2', 'biomaRt', 'pheatmap','tidyverse'))
library(tximport)
library(DESeq2)
library(biomaRt)
library(pheatmap)
library(tidyverse)
install.packages('RColorBrewer')
library(RColorBrewer)
install.packages('ggrepel')
library(ggrepel)

#Import the count data
sample_table = read.csv('https://raw.githubusercontent.com/sjcockell/mmb8052/main/practicals/practical_08/data/sample_table.csv')
files = pull(sample_table, Run)
#Merge data
files = paste0('counts/', files, '/quant.sf')
names(files) = pull(sample_table, Run)
#Create a vector for the tximport function
gene_map = read_csv('https://github.com/sjcockell/mmb8052/raw/main/practicals/practical_08/extdata/gene_map.csv')
txi = tximport(files, 
               type='salmon',
               tx2gene=gene_map,
               ignoreTxVersion=TRUE)


#Run the DESeq2 command for finding size factors and data normalisation
dds = DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ Group)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds)
dds = nbinomWaldTest(dds)
#Plot the estimates on a graph 
plotDispEsts(dds)
#Reduce the impact of outliers
rld = rlog(dds)
rld = vst(dds, blind=FALSE)
mat = assay(rld)
pca_plot = plotPCA(rld, intgroup=c('Group'))
show(pca_plot)

#Make data frame of PCA data and get group column to display only the experimental conditions
df_pca = plotPCA(rld, intgroup=c('Sample_Name'), returnData=TRUE)
df_pca$group = substr(df_pca$group, 8, nchar(df_pca$group))
#Add colour blind friendly colours
display.brewer.all(n=3, type="all", select=NULL, exact.n=TRUE, colorblindFriendly = TRUE)
custom_colors = brewer.pal(3, "Paired")
custom_fill_colors = brewer.pal(3, "Paired")
#Plot the PCA data onto a ggplot with aesthetics
gg_pca = ggplot(df_pca, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  ggtitle("PCA Plot of the distance between the three experimental groups") +
  xlab("PC1") +
  ylab("PC2") +  
  theme_bw() +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(aes(label = name), vjust = -0.5, hjust = -0.5, color = "black", 
                  box.padding = 0.5, 
                  point.padding = 0.3, 
                  force = 0.5,
                  size=3.5)
#Add ellipse around clustered samples, add colour blind friendly palette, and modify text size of axis and title
gg_pca + stat_ellipse(level = 0.95, geom = "polygon", aes(fill = group), alpha = 0.2) +
  scale_fill_manual(values = custom_fill_colors) +
  theme(plot.title = element_text(size = 15),
        axis.title = element_text(size = 15),  
        axis.text = element_text(size = 12))

#Cluster dendrogram  
#Colour blind friendly colours
cluster_colors = list(group= c(Allo24h="#A6CEE3", Allo2h= "#1F78B4", Naive="#B2DF8A"))
num_clusters = 3
cluster_colors=brewer.pal(3, 'Paired')
#Plot cluster dendrogram
gene_dist=dist(sample_distance_matrix_vst, method="euclidian")
gene_clust=hclust(gene_dist, method="complete")
plot(gene_clust, hang=-1, xlab = "Genes", ylab = "Euclidian distance")
rect.hclust(gene_clust,k=num_clusters, border = cluster_colors)

#Variance stabilising transformation of dds
vst = varianceStabilizingTransformation(dds)
#Create new colour blind friendly palette using 'purples' colours
palette= brewer.pal(3, "Purples")
color_low = "#EFEDF5"
color_medium = "#BCBDDC"
color_high = "#756BB1"
#Set up a color vector of the three colours
groupedcolors = c(color_low, color_medium, color_high)
#Set up another vector for the colour blind friendly palette 'Paired' used for the 3 experimental groups
group_colors = brewer.pal(3, "Paired")
group_colors = list(group= c(Allo24h="#A6CEE3", Allo2h= "#1F78B4", Naive="#B2DF8A"))
#Euclidian distance heatmap 
sample_distance_vst = dist(t(assay(vst)), method='euclidian')
sample_distance_matrix_vst = as.matrix(sample_distance_vst)
heatmap_annotation = data.frame(group=colData(vst)[,c('Group')], row.names=rownames(colData(vst)))
pheatmap(sample_distance_matrix_vst,
         clustering_distance_rows=sample_distance_vst,
         clustering_distance_cols=sample_distance_vst,
         annotation_col = heatmap_annotation,
         color = groupedcolors,
         annotation_colors = group_colors)


#Volcano plot for allo24h vs naive
#Differential gene expression: Allo24h vs naive
results_table_24h_vs_naive = results(dds, contrast= c('Group', 'Allo24h', 'Naive'))
summary(results_table_24h_vs_naive)
results_tibble_24h_vs_naive = as_tibble(results_table_24h_vs_naive, rownames='ensembl_gene_id')
#Retain only rows with complete cells
filtered_results_24h_vs_naive = filter(results_tibble_24h_vs_naive, complete.cases(results_tibble_24h_vs_naive))
#Merge required pvalues to vector
filtered_results_24h_vs_naive = mutate(filtered_results_24h_vs_naive, logPVal = -log10(padj))
filtered_results_24h_vs_naive = mutate(filtered_results_24h_vs_naive, significant=padj<0.05)
#Create new colour blind friendly palette 'Dark2'
palette= brewer.pal(3, "Dark2")
#Create the parameters required for volcano plot
filtered_results_24h_vs_naive = mutate(filtered_results_24h_vs_naive, 
  significant=ifelse(padj<0.05&abs(log2FoldChange)>1,'both', 
              ifelse(padj<0.05&abs(log2FoldChange)<1, 'pval only', 'neither')))
#Produce the volcano plot of allo24h vs naive
ggplot(filtered_results_24h_vs_naive, aes(x=log2FoldChange, y=logPVal, color = significant)) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(alpha = 0.3, size =1) + 
  scale_y_continuous(trans="log10") +
  theme_bw() +
  scale_color_manual(values =c(palette)) +
   labs(title = "Allo24h against Naive",
       x = "log2 Fold Change",
       y = "logPVal") +
        theme(plot.title = element_text(size = 20),  
         axis.title.x = element_text(size = 16),
         axis.title.y = element_text(size = 16),
         axis.text = element_text(size = 14))
#Calculate number of differentially expressed genes in 'both'
number_of_both_24hvsnaive <- sum(filtered_results_24h_vs_naive$significant == 'both')
cat("Number of points labeled as 'both':", number_of_both_24hvsnaive, "\n")


#Volcano plot for allo2h vs naive
#Differential gene expression: Allo2h vs naive
results_table_2h_vs_naive = results(dds, contrast= c('Group', 'Allo2h', 'Naive'))
summary(results_table_2h_vs_naive)
results_tibble_2h_vs_naive = as_tibble(results_table_2h_vs_naive, rownames='ensembl_gene_id')
#Retain only rows with complete cells
filtered_results_2h_vs_naive = filter(results_tibble_2h_vs_naive, complete.cases(results_tibble_2h_vs_naive))
#Merge required pvalues to vector
filtered_results_2h_vs_naive = mutate(filtered_results_2h_vs_naive, logPVal = -log10(padj))
filtered_results_2h_vs_naive = mutate(filtered_results_2h_vs_naive, significant=padj<0.05)
#Create the parameters required for volcano plot
filtered_results_2h_vs_naive = mutate(filtered_results_2h_vs_naive,
      significant=ifelse(padj<0.05&abs(log2FoldChange)>1,'both', 
                  ifelse(padj<0.05&abs(log2FoldChange)<1, 'pval only', 'neither')))
#Produce the volcano plot of allo2h vs naive
ggplot(filtered_results_2h_vs_naive, aes(x=log2FoldChange, y=logPVal, color = significant)) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(alpha = 0.3, size =1) + 
  scale_y_continuous(trans="log10") +
  theme_bw() +
  scale_color_manual(values =c(palette)) +
  labs(title = "Allo2h against Naive",
       x = "log2 Fold Change",
       y = "logPVal") +
  theme(plot.title = element_text(size = 20),  
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14))
#Calculate number of differentially expressed genes in 'both'
number_of_both_2hvsnaive <- sum(filtered_results_2h_vs_naive$significant == 'both')
cat("Number of points labeled as 'both':", number_of_both_2hvsnaive, "\n")


#Volcano plot for Allo2h vs allo24h
#Differential gene expression: Allo2h vs allo24h
results_table_2h_vs_24h = results(dds, contrast= c('Group', 'Allo2h', 'Allo24h'))
summary(results_table_2h_vs_24h)
results_tibble_2h_vs_24h = as_tibble(results_table_2h_vs_24h, rownames='ensembl_gene_id')
#Retain only rows with complete cells
filtered_results_2h_vs_24h = filter(results_tibble_2h_vs_24h, complete.cases(results_tibble_2h_vs_24h))
filtered_results_2h_vs_24h = mutate(filtered_results_2h_vs_24h, logPVal = -log10(padj))
#Merge required pvalues to vector
filtered_results_2h_vs_24h = mutate(filtered_results_2h_vs_24h, significant=padj<0.05)
#Create the parameters required for volcano plot
filtered_results_2h_vs_24h = mutate(filtered_results_2h_vs_24h,
     significant=ifelse(padj<0.05&abs(log2FoldChange)>1,'both', 
                 ifelse(padj<0.05&abs(log2FoldChange)<1, 'pval only', 'neither')))
#Produce the volcano plot of allo2h vs allo24h
ggplot(filtered_results_2h_vs_24h, aes(x=log2FoldChange, y=logPVal, color = significant)) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_point(alpha = 0.3, size =1) + 
  scale_y_continuous(trans="log10") +
  theme_bw() +
  scale_color_manual(values =c(palette)) +
  labs(title = "Allo2h against Allo24h",
       x = "log2 Fold Change",
       y = "logPVal") +
  theme(plot.title = element_text(size = 20),  
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 14))
#Calculate number of differentially expressed genes in 'both'
number_of_both_2hvs24h <- sum(filtered_results_2h_vs_24h$significant == 'both')
cat("Number of points labeled as 'both':", number_of_both_2hvs24h, "\n")


#Annotating the table and listing top differentially expressed genes 
ensembl108 = useEnsembl(biomart="ensembl", version=108)
ensembl108 = useDataset("mmusculus_gene_ensembl", mart=ensembl108)

#Annotation for allo2h-naive
annotation_2h_vs_naive = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                'start_position', 'end_position', 
                                'strand', 'gene_biotype', 'external_gene_name', 
                                'description'), 
                   filters = 'ensembl_gene_id', values = filtered_results_2h_vs_naive$ensembl_gene_id, 
                   mart = ensembl108)
annot_results_2h_vs_naive = left_join(filtered_results_2h_vs_naive, annotation_2h_vs_naive)
annot_results_2h_vs_naive = arrange(annot_results_2h_vs_naive, padj)
View(head(annot_results_2h_vs_naive, 10))
degs_2h_vs_naive = filter(annot_results_2h_vs_naive, abs(log2FoldChange) > 1 & padj < 0.05)
listAttributes(annotation_2h_vs_naive)
#View top 10 differentially expressed genes for allo2h-naive
View(head(degs_2h_vs_naive, 10))

#Annotation for allo24h-naive
annotation_24h_vs_naive = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                            'start_position', 'end_position', 
                                            'strand', 'gene_biotype', 'external_gene_name', 
                                            'description'), 
                               filters = 'ensembl_gene_id', values = filtered_results_24h_vs_naive$ensembl_gene_id, 
                               mart = ensembl108)
annot_results_24h_vs_naive = left_join(filtered_results_24h_vs_naive, annotation_24h_vs_naive)
annot_results_24h_vs_naive = arrange(annot_results_24h_vs_naive, padj)
View(head(annot_results_24h_vs_naive, 10))
degs_24h_vs_naive = filter(annot_results_24h_vs_naive, abs(log2FoldChange) > 1 & padj < 0.05)
listAttributes(annotation_24h_vs_naive)
#View top 10 differentially expressed genes for allo2h-naive
View(head(degs_24h_vs_naive, 10))



#Annotation for allo2h-allo24h
annotation_2h_vs_24h = getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 
                                             'start_position', 'end_position', 
                                             'strand', 'gene_biotype', 'external_gene_name', 
                                             'description'), 
                                filters = 'ensembl_gene_id', values = filtered_results_2h_vs_24h$ensembl_gene_id, 
                                mart = ensembl108)
annot_results_2h_vs_24h = left_join(filtered_results_2h_vs_24h, annotation_2h_vs_24h)
annot_results_2h_vs_24h = arrange(annot_results_2h_vs_24h, padj)
View(head(annot_results_2h_vs_24h, 10))
degs_2h_vs_24h = filter(annot_results_2h_vs_24h, abs(log2FoldChange) > 1 & padj < 0.05)
listAttributes(annotation_2h_vs_24h)

#View top 10 differentially expressed genes for allo2h-allo24h
View(head(degs_2h_vs_24h, 10))
