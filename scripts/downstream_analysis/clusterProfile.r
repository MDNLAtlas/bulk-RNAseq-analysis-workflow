library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(tidyr)
library(dplyr)
library(gridExtra)
library(forcats)

### Define a function to convert entrezID to UCSC gene symbol
convert_entrez_to_symbol <- function(entrez_ids) {
    entrez_list <- unlist(strsplit(entrez_ids, "/"))
    
    # Use bitr to convert Entrez IDs to Gene Symbols
    gene_symbols <- bitr(entrez_list, 
                         fromType = "ENTREZID", 
                         toType = "SYMBOL", 
                         OrgDb = org.Hs.eg.db)
    
    return(paste(gene_symbols$SYMBOL, collapse = "/"))
}

### Convert geneRatio and bgRatio into numeric
ratioStr2ratioNum <- function(ratioStr) {
    num <- as.numeric(unlist(strsplit(ratioStr, "/")))
    return(num[1] / num[2])
}

PRJ_DIR <- "~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
DGE <- paste(PRJ_DIR, "DGE/edgeR", sep = "/")
GO <- paste(PRJ_DIR, "GO", sep = "/")

dir.create(GO, recursive = T, showWarnings = F)

DGE_res <- read.csv(paste(DGE, "expected_count.genes.results.DGE.csv", sep = "/"),
  header = T,
  na.strings = c("", ".", "NA"))
genes <- DGE_res[DGE_res$FDR < 0.05, "gene_symbol"]
entrez_ids <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
m <- match(entrez_ids$SYMBOL, DGE_res$gene_symbol)
entrez_ids$logFC <- DGE_res$logFC[m]
entrez_ids <- entrez_ids[!is.na(entrez_ids$logFC), ]

### GO BP pathway enrichment analysis
GO_BP <- enrichGO(gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH")

GO_BP <- as.data.frame(GO_BP)
GO_BP$geneSymbol <- sapply(GO_BP$geneID,
    function(entrezID_list) convert_entrez_to_symbol(entrezID_list))
geneRatio <- sapply(GO_BP$GeneRatio, ratioStr2ratioNum)
BgRatio <- sapply(GO_BP$BgRatio, ratioStr2ratioNum)
GO_BP$OR <- (geneRatio / (1-geneRatio)) / (BgRatio / (1-BgRatio))
write.csv(GO_BP, paste(GO, "GO_BP_pathways.csv", sep = "/"))

### KEGG pathway enrichment analysis
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                                organism = "hsa",  # for human
                                pvalueCutoff = 0.05)
kegg_results_df <- as.data.frame(kegg_results)
kegg_results_df$geneSymbol <- sapply(kegg_results_df$geneID,
    function(entrezID_list) convert_entrez_to_symbol(entrezID_list)) 
geneRatio <- sapply(kegg_results_df$GeneRatio, ratioStr2ratioNum)
BgRatio <- sapply(kegg_results_df$BgRatio, ratioStr2ratioNum)
kegg_results_df$OR <- (geneRatio / (1-geneRatio)) / (BgRatio / (1-BgRatio))
write.csv(kegg_results_df, paste(GO, "KEGG_pathways.csv", sep="/"), row.names = F)

# # Plot GO interaction
# # cnetplot(kegg_results)
# cnetplot(
#   kegg_results,
#   showCategory = 5,
#   foldChange = entrez_ids$logFC
# )

# Pathwiew for KEGG pathway
fc <- entrez_ids$logFC
names(fc) <- entrez_ids$ENTREZID

setwd(GO)
pathview(
  gene.data = fc,
  pathway.id = "hsa05010",
  species = "hsa",
  out.dir = GO
)

pathview(
  gene.data = fc,
  pathway.id = "hsa05012",
  species = "hsa",
  out.dir = GO
)

pathview(
  gene.data = fc,
  pathway.id = "hsa05016",
  species = "hsa",
  out.dir = GO
)

pathview(
  gene.data = fc,
  pathway.id = "hsa05022",
  species = "hsa",
  out.dir = GO
)

pathview(
  gene.data = fc,
  pathway.id = "hsa05020",
  species = "hsa",
  out.dir = GO
)

# hsa05014                     Amyotrophic lateral sclerosis  203/2655 368/9380
# hsa05012                                 Parkinson disease  151/2655 268/9380 # view
# hsa05010                                 Alzheimer disease  196/2655 388/9380 # view
# hsa05016                                Huntington disease  162/2655 308/9380 # view
# hsa04110                                        Cell cycle   97/2655 158/9380
# hsa05022 Pathways of neurodegeneration - multiple diseases  222/2655 480/9380 # view

### Draw GO pathways dotplot
GO_BP <- read.csv(paste(GO, "GO_BP_pathways.csv", sep = "/"), header = T)
relevant_pathways <- c("cytoplasmic translation",
  "cytoplasmic translational initiation",
  "formation of cytoplasmic translation initiation complex",
	"positive regulation of cytoplasmic translation",
	"regulation of translational fidelity",
	"aminoacyl-tRNA metabolism involved in translational fidelity",
	"stress granule assembly",
	"regulation of stress granule assembly",
  "response to amino acid starvation",
  "cellular response to amino acid starvation",
  "cellular response to glucose starvation",
	"IRE1-mediated unfolded protein response",
	"protein folding in endoplasmic reticulum",
	"regulation of IRE1-mediated unfolded protein response",
	"response to nitrosative stress",
	"mitochondrial depolarization")

GO_BP_subset <- GO_BP[GO_BP$Description %in% relevant_pathways, ]
# Prepare GO data (example: top 20 terms)
GO_BP_subset <- GO_BP_subset %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(Description = fct_reorder(Description, GeneRatio),
    GeneRatio = sapply(GeneRatio, function(x) {
        parts <- as.numeric(unlist(strsplit(x, "/")))
        parts[1] / parts[2]
    }),
    
    # Ensure Count is integer
    Count = as.integer(Count))  # reorder by GeneRatio

# Draw bubble plot
go_bp_plot <- ggplot(GO_BP_subset, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. P value") +
  scale_size(range = c(3, 10)) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    title = "GO BP",
    size = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
png(paste(GO, "GO_BP_subset_pathways.png", sep = "/"), width=10*300, height=6*300, res=300)
print(go_bp_plot)
dev.off()

### Draw KEGG pathway dotplot
kegg_results_df <- read.csv(paste(GO, "KEGG_pathways.csv", sep="/"), header = T)

kegg_results_subset <- kegg_results_df[1:10, ]
kegg_results_subset <- kegg_results_subset %>%
  arrange(p.adjust) %>%
  slice_head(n = 20) %>%
  mutate(Description = fct_reorder(Description, GeneRatio),
    GeneRatio = sapply(GeneRatio, function(x) {
        parts <- as.numeric(unlist(strsplit(x, "/")))
        parts[1] / parts[2]
    }),
    
    # Ensure Count is integer
    Count = as.integer(Count))  # reorder by GeneRatio

# Draw bubble plot
kegg_plot <- ggplot(kegg_results_subset, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", name = "Adj. P value") +
  scale_size(range = c(3, 10)) +
  labs(
    x = "Gene Ratio",
    y = NULL,
    title = "KEGG Pathways",
    size = "Gene Count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
png(paste(GO, "KEGG_pathways.png", sep = "/"), width=6*300, height=4*300, res=300)
print(kegg_plot)
dev.off()

# png(paste(GO, "KEGG_pathways.png", sep = "/"), width = 6*300, height = 4*300, res=300)
# kegg_plot <- dotplot(kegg_results, showCategory = 10) +
#   ggtitle("KEGG Pathways") +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold")
#   )
# dev.off()



