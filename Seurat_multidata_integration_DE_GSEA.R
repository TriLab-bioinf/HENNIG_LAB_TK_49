# Load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(HGNChelper)
library(viridis)
library(hrbrthemes)
library(ggpubr)


setwd("/Users/leeh15/Desktop/Stat5B/Results/scRNA-seq/GEX/")
### no folder except data set in working directory ###

### Note: Jump to line 189 if object "GEX_LN_seurat_analysis_HK.rds" is already saved in working directory 

# Load datasets
dirs <- list.dirs(path = "./", full.names = F, recursive = F)

object_name_list = c()
for (x in dirs){
  name = gsub('[0-9]+_GEX_','',x, perl = T)
  cts <- ReadMtx(mtx = paste0('./',x,'/filtered_feature_bc_matrix/matrix.mtx.gz'),
                 features = paste0('./',x,'/filtered_feature_bc_matrix/features.tsv.gz'),
                 cells = paste0('./',x,'/filtered_feature_bc_matrix/barcodes.tsv.gz')
  )
  assign(name,CreateSeuratObject(counts = cts))
  object_name_list <- c(object_name_list, name)
}
 
# Merge Seurat objects
liver_merged_seurat <- merge(x = Liver_1793_WT,
                             y=c(Liver_1812_WT,Liver_1797_WT,
                                 Liver_1794_Het,Liver_1798_Het,Liver_1811_Het,
                                 Liver_1795_Homo,Liver_1796_Homo),
                             add.cell.ids = c("Liver_1793_WT","Liver_1812_WT","Liver_1797_WT","Liver_1794_Het","Liver_1798_Het","Liver_1811_Het","Liver_1795_Homo","Liver_1796_Homo"),
                             project = 'GEX')

ln_merged_seurat <- merge(x = LN_1797_WT,
                          y=c(LN_1793_WT,LN_1812_WT,LN_1794_Het,LN_1798_Het,LN_1811_Het,LN_1795_Homo,LN_1796_Homo),
                          add.cell.ids = c("LN_1797_WT","LN_1793_WT","LN_1812_WT","LN_1794_Het","LN_1798_Het","LN_1811_Het","LN_1795_Homo","LN_1796_Homo"),
                          project = 'GEX')

# Add useful info to metadata table
ln_merged_seurat$sample <- rownames(ln_merged_seurat@meta.data)

ln_merged_seurat@meta.data <- separate(data = ln_merged_seurat@meta.data, col = "sample", into = c("tissue", "mouse", "genotype", "barcode"), sep = "_")

# Data quality control
# # Calculate mitochondrial percentage
ln_merged_seurat$mitoPercent <- PercentageFeatureSet(ln_merged_seurat, pattern = '^mt-')

# Explore QC
# Make Plot dir
dir.create(path = "./Plots", showWarnings = F)
# # Visualize QC metrics as a violin plot
p1 <- VlnPlot(ln_merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
ggsave2(filename = "./Plots/QC_feat_cnt_mt.pdf", plot = p1)

## Saving 7 x 5 in image
plot1 <- FeatureScatter(ln_merged_seurat, feature1 = "nCount_RNA", feature2 = "mitoPercent") + geom_hline(yintercept = 5, col = "gray") + geom_vline(xintercept = 500, col = "gray")
plot2 <- FeatureScatter(ln_merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept = c(200,4300), col = "gray") + geom_vline(xintercept = 500, col = "gray")
ggsave2(filename = "./Plots/QC_scatterplots.pdf",plot=grid.arrange(plot1,plot2, nrow = 2, ncol = 1), height = 11)

## Saving 7 x 11 in image
grid.arrange(plot1,plot2, nrow = 2, ncol = 1)

# Data filtering
ln_merged_seurat_filtered <- subset(ln_merged_seurat, subset = nCount_RNA > 500 &
                                      nFeature_RNA > 500 &
                                      nFeature_RNA < 4300 &
                                      mitoPercent < 5)

# Perform standard workflow to check if there is any batch effects
ln_merged_seurat_filtered <- NormalizeData(object = ln_merged_seurat_filtered, )
ln_merged_seurat_filtered <- FindVariableFeatures(object = ln_merged_seurat_filtered)
ln_merged_seurat_filtered <- ScaleData(object = ln_merged_seurat_filtered)
ln_merged_seurat_filtered <- RunPCA(object = ln_merged_seurat_filtered)

# Make PC selection
elbow.p1 <- ElbowPlot(object = ln_merged_seurat_filtered, ndims = 50)
ggsave2(filename = "./Plots/elbowplot.pdf", plot = elbow.p1)
elbow.p1

# Clustering and UMAP
ln_merged_seurat_filtered <- FindNeighbors(object = ln_merged_seurat_filtered, dims = 1:40)
ln_merged_seurat_filtered <- FindClusters(object = ln_merged_seurat_filtered, resolution = c(0.1, 0.3, 0.5, 0.7, 1), random.seed = 1234)
ln_merged_seurat_filtered <- RunUMAP(object = ln_merged_seurat_filtered, dims = 1:40)

# Dimension Plots
p1 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'mouse')
p2 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'genotype', cols = c('red','green','blue'))
p3 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'RNA_snn_res.0.1')
p4 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'RNA_snn_res.0.3')
p5 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'RNA_snn_res.0.5')
p6 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'RNA_snn_res.0.7')
p7 <- DimPlot(object = ln_merged_seurat_filtered, reduction = 'umap', group.by = 'RNA_snn_res.1')

ggsave2(filename = "./Plots/ln_umap_mouse_genotype.pdf", plot = grid.arrange(p1, p2,p3,p4,p5,p6,p7, ncol=2, nrow=4), width = 11, height = 15)
grid.arrange(p1, p2,p3, ncol = 2, nrow=2)

# Asign cluster number
ln_merged_seurat_filtered$seurat_clusters = ln_merged_seurat_filtered$RNA_snn_res.0.5


# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system"
# e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = ln_merged_seurat_filtered[["RNA"]]@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(ln_merged_seurat_filtered@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(ln_merged_seurat_filtered@meta.data[ln_merged_seurat_filtered@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(ln_merged_seurat_filtered@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

# Overlay the identified cell types on UMAP plot:
ln_merged_seurat_filtered@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  ln_merged_seurat_filtered@meta.data$customclassif[ln_merged_seurat_filtered@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

annot_clusters.p1 <- DimPlot(ln_merged_seurat_filtered, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')
annot_clusters_by_genotype.p1 <- DimPlot(ln_merged_seurat_filtered, reduction = "umap", label = FALSE, repel = TRUE, group.by = 'genotype')

ggsave2(filename = "./Plots/ln_annotated_clusters_umap.pdf", plot = annot_clusters.p1 + annot_clusters_by_genotype.p1, width = 21, height = 11)
annot_clusters.p1 + annot_clusters_by_genotype.p1


# ### Save Seurat analysis
# Save final classified pbmc dataset
saveRDS(ln_merged_seurat_filtered, file = "GEX_LN_seurat_analysis_03272023.rds")

### Batch correction
#### Perform integration of Seurat objects to correct for batch effect
obj.list <- SplitObject(ln_merged_seurat_filtered, split.by = 'mouse')
for (i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}

#### Select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)


#### Find integration anchors using cannonical correlation analysis (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list)


#### Integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


#### Scale data, run PCS and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:40)


#### Visualize integrated clusters
int.p1 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'mouse')
int.p2 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'genotype')
ggsave2(filename = "./Plots/ln_integrated_umap_mouse_genotype.pdf", plot = grid.arrange(int.p1, int.p2, ncol = 2), width = 11)
grid.arrange(int.p1, int.p2, ncol = 2)


### Read RDS object
ln_merged_seurat_filtered <- readRDS(file = "GEX_LN_seurat_analysis_HK.rds")


### Estimate proportion of each cell type per sample
ln_meta <- as_tibble(ln_merged_seurat_filtered@meta.data)

# Remove unwanted columns
keep <- grep(pattern = '^RNA_snn_res', colnames(ln_meta), invert = T)
ln_meta <- ln_meta[,keep]


#### Add columns with data summaries grouped by genotype, mouse ID and cell-type
cell_dist <- ln_meta %>% group_by(genotype, mouse) %>%
  count(genotype, customclassif) %>%
  mutate(percent = n / sum(n) * 100) %>%
  mutate(group = paste0(genotype,"_",customclassif))

colnames(cell_dist) <- c("Genotype","Mouse_ID","Cell_type","Cell_count","Percent","Group" )

write.csv(x = cell_dist, file = "Cell_distribution_summary.csv",col.names = F )
head(cell_dist)


#### Cell-type composition plots
cell_dist.p <- cell_dist %>% ggplot(aes(fill=Cell_type, y=Cell_count, x=Genotype)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Proportion of cell-types per genotype") +
  theme_ipsum(base_family = 'Helvetica', base_size = 10, axis_title_just = 0.5) +
  ylab('Percent') 

ggsave2(filename = "./Plots/cell_distrib_barchar.pdf", plot = cell_dist.p)
print(cell_dist.p)


cell_dist_perc <- cell_dist %>% ggplot(aes(fill=Cell_type, y=Percent, x=Genotype)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Proportion of cell-types per genotype") +
  theme_ipsum(axis_title_just = 0.5, base_family = 'Helvetica', strip_text_size=10) +
  facet_wrap(~Cell_type, labeller = label_wrap_gen())


ggsave2(filename = "./Plots/cell_distrib_barchar_facet.pdf", plot = cell_dist_perc, width = 13)
print(cell_dist_perc)


cell_dist.box.p <- cell_dist %>% ggplot(aes(x = Genotype, y = Percent, fill=Cell_type)) + 
  geom_boxplot() + 
  scale_fill_viridis(discrete = T) +
  theme_ipsum(axis_title_just = 0.5, 
              base_family = 'Helvetica', 
              strip_text_size=10) +
  facet_wrap(~Cell_type, labeller = label_wrap_gen())

ggsave2(filename = "./Plots/cell_distrib_boxplot_facet.pdf", plot = cell_dist.box.p, width = 13, height = 13)
print(cell_dist.box.p)

# sessionInfo()





library(stats4)
library(DESeq2)
library(ashr)
library(AnnotationDbi)
library(clusterProfiler)
library(msigdbr)
library(ReactomePA)
library(org.Mm.eg.db)
library(sp)


##########################################
### expression list and GSEA analysis ####
##########################################

### Aggregate raw read counts by sample and genotype
ln_merged_seurat_filtered@meta.data$samples <- paste(ln_merged_seurat_filtered@meta.data$mouse,ln_merged_seurat_filtered@meta.data$genotype, sep = '_')
#ln_merged_seurat_filtered@meta.data$samples <- paste(ln_merged_seurat_filtered@meta.data$genotype, sep = '_')

cts <- AggregateExpression(object = ln_merged_seurat_filtered,
                           group.by = c('customclassif','samples'),
                           assays = 'RNA', slot = 'counts', return.seurat = FALSE)

cts <- as.matrix(cts$RNA)

# transpose cts and convert to DF
cts.t <- as.data.frame(t(cts))

# create chr vector with factors used to split cts.t along rows bt cell type
splitRows <- gsub("_.*","", rownames(cts.t))

# Split cts.t dataframe along cell types creating a list of read count matrices, one per cell type.
cts.split <- split.data.frame(cts.t, f = factor(splitRows))

# Remove cell type from col names and transpose each list element 
cts.split.mod <- lapply(cts.split, function(x){
  rownames(x) <- gsub(".*_(.*_[Het|Homo|WT])","\\1", rownames(x))
  t(x)
})

## NOTE: it would be possible to generate a single data.frame with the read counts across all cell types and samples, but the current approach allows to work with simpler tables.

### Run DE analysis with DESeq2 on WT cells between CD8+ NKT-like and Naive CD4+ T cells

suppressMessages(library(DESeq2))

## NOTE: It would be possible to make the code below generalizable for any pair of cell types, and then create a function that uses as input the two cell types to be compared as well as the cts.split.mod object and that outputs the DE results.

# Extract counts for CD8+ NKT-like cells
counts_cd8nkt <- cts.split.mod$`CD8+ NKT-like cells`
#counts_cd4t <- cts.split.mod$`Naive CD4+ T cells`

# Add cell type flag to col names
#colnames(counts_cd8nkt) <- paste(colnames(counts_cd8nkt),"cd8nkt", sep = "_")
#colnames(counts_cd4t) <- paste(colnames(counts_cd4t),"cd4t", sep = "_")

# Merge both cell types
#counts_merged <- cbind(counts_cd4t, counts_cd8nkt)

# Create metadata table
metadata <- data.frame(samples = colnames(counts_cd8nkt))
metadata <- metadata %>% mutate(genotype = gsub(".*_(.*)","\\1", samples)) 

# Convert metadata columns to factors
metadata <- as.data.frame(unclass(metadata),stringsAsFactors=TRUE)

# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = counts_cd8nkt, colData = metadata, design = ~ genotype)

# Filter rows with less than 10 counts per gene across samples
keep <-  rowSums(counts(dds)) >= 10
dds <- dds[keep,] 

# Set WT as the reference
dds@colData$genotype <- relevel(x = dds@colData$genotype, ref = 'WT')

# run DESeq2
dds <- DESeq(object = dds)

# Print out coefficients for choosing comparison
resultsNames(dds)



# Get results for Het vs WT
coef_id <- "cd8nkt_Het_vs_WT"
coef = c("genotype","Het","WT")
DE_res <-results(dds, contrast=c(coef))
# Apply shrinkage to Log2FC
DE_res_shrink <- lfcShrink(dds, 
                           contrast = c(coef_id),
                           type = "ashr", 
                           res=DE_res )

# Sort results table based on Log2FC
DE_res_shrink_sorted <- DE_res_shrink[order(DE_res_shrink$log2FoldChange, decreasing = TRUE),]
write.table(x = DE_res_shrink_sorted, file = "cd8nkt_Het_vs_WT.txt", sep = "\t", col.names = NA)
write.table(x = DE_res_shrink_sorted, file = "cd8nkt_Het_vs_WT.csv", sep = ",", col.names = NA)

#####################################################################################################################
# Run GSEA analysis on DESeq2 result
#####################################################################################################################

# Load required libraries
library("clusterProfiler")
library("org.Mm.eg.db")
library("msigdbr")
library("ReactomePA")

# Function to replace gene symbols with Ensembl IDs or other types of IDs
replace_annotation_keys <- function(key_list, from, to){
  if(from %in% keytypes(org.Mm.eg.db) & to  %in% keytypes(org.Mm.eg.db)){
    if(typeof(key_list) == "S4"){
      ensemble_ids <- rownames(key_list)
      symbols <- mapIds(org.Hs.eg.db, keys = ensemble_ids, column = to, keytype = from, multiVals = "first")
      symbols <- symbols[!is.na(symbols)]
      to_name <- rownames(key_list) %in% names(symbols)
      key_list@rownames[to_name] <- as.vector(symbols)
      return(key_list)
    } else {
      to_ids <- mapIds(org.Mm.eg.db, keys = key_list, column = to, keytype = from, multiVals = "first")
      to_ids <- to_ids[!is.na(to_ids)]
      to_name <- key_list %in% names(to_ids)
      key_list[to_name] <- as.vector(to_ids)
      return(key_list)
    }
  } else {
    print("Annotation keys are not keytypes")
    stop()
  }
}

# Function to convert a character vector of strings of Ensembl IDs/Entrez IDs separated by "/" into string of Symbol IDs
get_symbols_from_gene_accessions <- function(accession_list, accession_type = "ENSEMBL"){
  #' Converts a string of gene accessions into a string of gene symbols
  #' @param accession_list character vector. List of gene accession (Ensembl or Entrez)
  #' @param accession_type character string. It should be either "ENSEMBL" or ENTREZID. Defaults to "ENSEMBL".
  #' @usage get_symbols_from_gene_accessions(accession_list=my_string_of_acc, accession_type="ENTREZID")
  #' @return The converted list of gene symbols as a string. Missing symbols will keep their original accession ID.
  
  
  symbol_list <- sapply(accession_list, function(x){
      x <- unlist(str_split(string = x, pattern = "/", simplify = FALSE))
      x <- replace_annotation_keys(key_list = x, from = accession_type, to = "SYMBOL")
      x <- str_flatten(x, collapse = "/")
    }
  )
  return(as.vector(symbol_list))
}


# Make output directory
dir.create(path = "./GSEA", showWarnings = FALSE)

# Prepare sorted gene list based on Log2FC results
gene_list <- DE_res_shrink_sorted$log2FoldChange
names(gene_list) <- replace_annotation_keys(key_list = rownames(DE_res_shrink_sorted), from = "SYMBOL", to = "ENSEMBL")

#########################################
# 1) Running GSEA on Gene Ontology DB (BP)
#########################################

gse <- gseGO(geneList=gene_list,
             ont ="BP",
             keyType = 'ENSEMBL',
             #nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Mm.eg.db,
             pAdjustMethod = "BH")

# Add gene symbols
gse.df <- as.data.frame(gse)
gse.df$core_enrichment_symbols <- get_symbols_from_gene_accessions(accession_list = gse.df$core_enrichment, 
                                                                   accession_type = 'ENSEMBL')

write.table(x = as.data.frame(gse.df), file = "./GSEA/GO_BP_cd8nkt_Het_vs_WT.csv", sep = ",", col.names = NA)


###################################
# 2) Using GSEA on MSigDB gene sets
###################################

# Fetch Hallmarks gene set for mouse genome
msig_h <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, ensembl_gene, ) 

# Fetch Reactome gene set for mouse genome
msig_c2_reactome <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:REACTOME") %>% dplyr::select(gs_name, ensembl_gene, )

# Run Gene Set Enrichment Analysis on Hallmarks gene set
msig_h.gsea <- GSEA(gene_list, TERM2GENE = msig_h, pvalueCutoff = 1, verbose = TRUE, eps = 0) 
msig_h.gsea.df <- as.data.frame(msig_h.gsea)
msig_h.gsea.df$core_enrichment_symbols <- get_symbols_from_gene_accessions(accession_list = msig_h.gsea.df$core_enrichment, 
                                                                   accession_type = 'ENSEMBL')

write.table(x = msig_h.gsea.df, file = "./GSEA/msig_hallmarks_cd8nkt_Het_vs_WT.csv", sep = ",", col.names = NA)

# Run Gene Set Enrichment Analysis on Reactome gene set
msig_reactome.gsea <- GSEA(gene_list, TERM2GENE = msig_c2_reactome, pvalueCutoff = 1, verbose = TRUE, eps = 0) 
msig_reactome.gsea.df <- as.data.frame(msig_reactome.gsea)
msig_reactome.gsea.df$core_enrichment_symbols <- get_symbols_from_gene_accessions(accession_list = msig_reactome.gsea.df$core_enrichment, 
                                                                   accession_type = 'ENSEMBL')

write.table(x = msig_reactome.gsea.df, file = "./GSEA/msig_reactome_cd8nkt_Het_vs_WT.csv", sep = ",", col.names = NA)

# Run Gene Set Enrichment Analysis on Reactome gene set using reactomePA function (it needs gene acc to be converted to Entrez IDs). This function outputs Reactome pathway IDs that can be used to visualize pathways in the Reactome website and also can be used to generate pathway networs with the results with the viewPathway function.

# 1) Convert symbols to Entrez IDs
gene_list.entrez <- gene_list
names(gene_list.entrez) <- replace_annotation_keys(key_list = rownames(DE_res_shrink_sorted), from = "SYMBOL", to = "ENTREZID")
# 2) Run GSEA with gsePathway function
msig_reactome.gsea2 <- gsePathway(geneList = gene_list.entrez, pvalueCutoff = 1, verbose = TRUE, eps = 0, organism = "mouse")
msig_reactome.gsea2.df <- as.data.frame(msig_reactome.gsea2)
msig_reactome.gsea2.df$core_enrichment_symbols <- get_symbols_from_gene_accessions(accession_list = msig_reactome.gsea2.df$core_enrichment, 
                                                                   accession_type = 'ENTREZID')

write.table(x = msig_reactome.gsea2.df, file = "./GSEA/msig_reactome2_cd8nkt_Het_vs_WT.csv", sep = ",", col.names = NA)

sessionInfo()
