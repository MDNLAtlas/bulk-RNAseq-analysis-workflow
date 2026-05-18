library(synaptome.db)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(cowplot)
library(forcats)

PRJ_DIR <- "~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
DGE <- paste(PRJ_DIR, "DGE/edgeR", sep = "/")
GO <- paste(PRJ_DIR, "GO", sep = "/")

dir.create(GO, recursive = T, showWarnings = F)

DGE_res <- read.csv(paste(DGE, "expected_count.genes.results.DGE.csv", sep = "/"),
  header = T,
  na.strings = c("", ".", "NA"))

DGE_res <- DGE_res[!grepl("^ENSG", DGE_res$gene_symbol), ]

# Thresholds
logfc_cutoff <- 1
padj_cutoff  <- 0.05

# Convert p-value to -log10
DGE_res$negLogP <- -log10(DGE_res$FDR)

# Categorize significance
DGE_res$Category <- "Not significant"

# DGE_res$Category[
#   DGE_res$FDR < padj_cutoff &
#   DGE_res$logFC > logfc_cutoff
# ] <- "Upregulated"

# DGE_res$Category[
#   DGE_res$FDR < padj_cutoff &
#   DGE_res$logFC < -logfc_cutoff
# ] <- "Downregulated"

DGE_res$Category[
  DGE_res$FDR < padj_cutoff &
  DGE_res$logFC > 0
] <- "Upregulated"

DGE_res$Category[
  DGE_res$FDR < padj_cutoff &
  DGE_res$logFC < 0
] <- "Downregulated"

DGE_res$Category_2 <- "Not significant"
DGE_res$Category_2[
  DGE_res$FDR < padj_cutoff &
  DGE_res$logFC > 1
] <- "Upregulated"

DGE_res$Category_2[
  DGE_res$FDR < padj_cutoff &
  DGE_res$logFC < -1
] <- "Downregulated"

# Example gene lists
# dementia_risk_genes <- c(
#   "ABCA7", "APOE", "APP", "CHMP2B", "CSF1R",
#   "FUS", "GRN", "MAPT", "MT-ATP6", "MT-ATP8",
#   "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1",
#   "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5",
#   "MT-ND6", "MT-RNR1", "MT-RNR2", "MT-TA", "MT-TC",
#   "MT-TD", "MT-TE", "MT-TF", "MT-TG", "MT-TH",
#   "MT-TI", "MT-TK", "MT-TL1", "MT-TL2", "MT-TM",
#   "MT-TN", "MT-TP", "MT-TQ", "MT-TR", "MT-TS1",
#   "MT-TS2", "MT-TT", "MT-TV", "MT-TW", "MT-TY",
#   "PRNP", "PSEN1", "PSEN2", "RNF216", "SIGMAR1",
#   "SNCA", "SORL1", "TARDBP", "TREM2", "TUBA4A",
#   "UBE3A", "UBQLN2", "VCP"
# )

# Tijms, B. M. et al. Cerebrospinal fluid proteomics in patients with Alzheimer’s disease reveals 
# five molecular subtypes with distinct genetic risk profiles. Nat. Aging 4, 33–47 (2024).
# The role of genetics in neurodegenerative dementia: a large cohort study in South China.
# Genetic background of cognitive decline in Parkinson's disease Front in Cognition
dementia_risk_genes <- c(
  "PSEN1", "PSEN2", "APP", "APOE",
  "TREM2", "IDUA", "CLNK", "SCIMP",
  "BIN1", "PICALM", "IL-34", "CLNK", "ABCA7",  # AD, EOAD
  "MAPT", "GRN", "C9orf72", "CHCHD10", 
  "HTRA1", "TBK1", "OPTN", "SQSTM1",
  "VCP", "SIGMAR1", "HTT",  # FTD
  "GBA", "SNCA", # Dementia Lewy body
  "SNCA", "LRRK2", "GBA1", "PRKN",
  "PINK1", "DJ-1", "VPS35", "TMEM175", # PDD
  "CHMP2B", "FUS", "TARDBP", "PRNP" # AD neurodegenerative dementias
)

severe_tinnitus_genes <- c(
  "PRUNE2", "AKAP9", "SORBS1", "ITGAX", "ANK2",
  "KIF20B", "TSC2", "SPHK2", "SYNPO", "LRPPRC",
  "XYLT1", "ALCAM", "CDH13", "DOCK7", "BIN1",
  "FLII", "HSPA4L", "IQSEC1", "IQSEC3", "LLGL1", 
  "MADD", "MBP", "MPRIP", "NRCAM", "TRAP1",
  "VCAN", "MYO18A", "MYO5A", "PPP1R9A", "CCDC22",
  "EPX"
)

# hearing_loss_genes <- c(
#   "ABHD12", "ACTB", "ACTG1", "ADCY1", "ADGRV1",
#   "AIFM1", "ALMS1", "AMMECR1", "ANKH", "ATP2B2",
#   "ATP6V0A4", "ATP6V1B1", "ATP6V1B2", "BCS1L", "BDP1",
#   "BSND", "BTD", "CABP2", "CACNA1D", "CCDC50",
#   "CD164", "CDC14A", "CDH23", "CEACAM16", "CEP78",
#   "CHD7", "CHSY1", "CIB2", "CISD2", "CLDN14",
#   "CLDN9", "CLIC5", "CLPP", "CLRN1", "COCH",
#   "COL11A1", "COL11A2", "COL2A1", "COL4A3", "COL4A4",
#   "COL4A5", "COL4A6", "COL9A1", "COL9A2", "COL9A3",
#   "CRYM", "DCAF17", "DCDC2", "DIABLO", "DIAPH1",
#   "DIAPH3", "DLX5", "DMXL2", "DNMT1", "DSPP",
#   "EDN3", "EDNRB", "ELMOD3", "EPS8", "EPS8L2",
#   "ERAL1", "ESPN", "ESRRB", "EYA1", "EYA4",
#   "FDXR", "FGF3", "FGFR1", "FGFR2", "FGFR3",
#   "FITM2", "FOXI1", "GAB1", "GATA3", "GIPC3",
#   "GJB2", "GJB3", "GJB6", "GPRASP2", "GPSM2",
#   "GRAP", "GREB1L", "GRHL2", "GRXCR1", "GRXCR2",
#   "GSDME", "HARS2", "HGF", "HOMER2", "HOXA2",
#   "HOXB1", "HSD17B4", "IFNLR1", "ILDR1", "KARS1",
#   "KCNE1", "KCNJ10", "KCNQ1", "KCNQ4", "KITLG",
#   "KMT2D", "LARS2", "LHFPL5", "LHX3", "LMX1A", 
#   "LOXHD1", "LOXL3", "LRP2", "LRTOMT", "MAN2B1",
#   "MANBA", "MARVELD2", "MASP1", "MCM2", "MET",
#   "MGP", "MIR96", "MITF", "MPZL2", "MSRB3",
#   "MT-CO1", "MT-ND1", "MT-RNR1", "MT-TH", "MT-TI",
#   "MT-TK", "MT-TL1", "MT-TS1", "MT-TS2", "MYH14",
#   "MYH9", "MYO15A", "MYO3A", "MYO6", "MYO7A",
#   "NARS2", "NDP", "NEFL", "NF2", "NLRP3",
#   "NOG", "NR2F1", "OPA1", "OSBPL2", "OTOA",
#   "OTOF", "OTOG", "OTOGL", "P2RX2", "PAX1",
#   "PAX3", "PCDH15", "PDE1C", "PDZD7", "PEX1",
#   "PEX26", "PEX6", "PJVK", "PLS1", "PNPT1",
#   "POLR1B", "POLR1C", "POLR1D", "POU3F4", "POU4F3",
#   "PPIP5K2", "PRPS1", "PTPRQ", "RAI1", "RDX",
#   "REST", "RIPOR2", "ROR1", "S1PR2", "SEMA3E",
#   "SERPINB6", "SIX1", "SIX2", "SIX5", "SLC17A8",
#   "SLC19A2", "SLC22A4", "SLC26A4", "SLC26A5", "SLC33A1",
#   "SLC44A4", "SLC4A11", "SLC52A2", "SLC52A3", "SLITRK6",
#   "SMPX", "SNAI2", "SOX10", "SPATA5", "SPNS2",
#   "STRC", "SUCLA2", "SYNE4", "TBC1D24", "TBL1X",
#   "TBX1", "TCOF1", "TECTA", "TFAP2A", "TIMM8A",
#   "TJP2", "TMC1", "TMEM126A", "TMEM132E", "TMIE", 
#   "TMPRSS3", "TNC", "TPRN", "TRIOBP", "TRRAP",
#   "TSPEAR", "TUBB4B", "TWNK", "USH1C", "USH1G",
#   "USH2A", "WBP2", "WFS1", "WHRN"
# )

# Synaptic genes in Alzheimer's disease, Vascular, frontal temporal, and Lewy body
# Rare types (Parkinson's, Huntington's disease, Creutzfeldt-Jakob disease)
gcp <- findGeneByCompartmentPaperCnt(cnt = 1)
synap_genes <- unique(gcp$HumanName)
synap_genes_dt <- as.data.frame(getGeneDiseaseByName(synap_genes))
dementia_types <- c("dementia", "Alzheimer's_disease",
  "Frontaltemporal_dementia", "Vascular_dementia",
  "Parkinson's_disease", "Lewy_body_dementia")
pat <- paste(dementia_types, collapse = "|")
synap_dementia_dt <- synap_genes_dt[grep(pat, synap_genes_dt$Description), ]
synap_dementia_genes <- unique(synap_dementia_dt$HumanName)

# Synaptic genes with differentially expressed in dementia cases vs. control
# Synaptic markers of cognitive decline in neurodegenerative diseases: a proteomic approach
# de_synap_genes <- c("SV2C", "NRGN", "CBLN4", "BDNF", "GAP43", 
#   "GRIA3", "CAMK2A", "SYBU", "VDAC2", "ARC",
#   "RAB11A", "PDYN", "GRIA4", "SYT2", "CAMK2G", 
#   "CNIH2", "KCNIP2", "SNAP47", "TECR", "CACNG2",
#   "PVRL3", "LRFN2", "GRIK2", "CACNG3", "TNK2", "CAMKK1")



# Highlight groups
DGE_res$dementia_risk_gene <- ifelse(
  DGE_res$gene_symbol %in% dementia_risk_genes,
  "Yes",
  "No"
)

DGE_res$synap_dementia_gene <- ifelse(
  DGE_res$gene_symbol %in% synap_dementia_genes,
  "Yes",
  "No"
)

DGE_res$severe_tinnitus_gene <- ifelse(
  DGE_res$gene_symbol %in% severe_tinnitus_genes,
  "Yes",
  "No"
)

# DGE_res$Highlight <- NA

# DGE_res$Highlight[DGE_res$gene_symbol %in% dementia_risk_genes] <- "Dementia risk"

# DGE_res$Highlight[DGE_res$gene_symbol %in% synap_dementia_genes] <- "Synaptic dementia"

# DGE_res$Highlight[DGE_res$gene_symbol %in% severe_tinnitus_genes] <- "Severe tinnitus"

# dim(DGE_res[DGE_res$Category %in% c("Downregulated", "Upregulated"),])
# # [1] 5901   19
# dim(DGE_res[DGE_res$Category %in% c("Downregulated", "Upregulated") & DGE_res$synap_dementia_gene == "Yes",])
# # [1] 408  19
# dim(DGE_res[DGE_res$Category %in% c("Downregulated", "Upregulated") & DGE_res$synap_dementia_gene == "Yes",]) / dim(DGE_res[DGE_res$Category %in% c("Downregulated", "Upregulated"),])
# [1] 0.06914082


# Base volcano plot
phenotypes_map <- list(
  dementia_risk_gene = "Dementia risk genes",
  severe_tinnitus_gene = "Severe tinnitus genes",
  synap_dementia_gene = "Dementia-associated synaptic genes"
)

p <- list()

# Synaptic dementia genes
pnt <- "synap_dementia_gene"
res <- DGE_res[DGE_res[[pnt]] == "Yes", ]

# Highly down- and up-regulated proteins in dementia synaptic genes
hr_synap_dementia_genes <- res[res$Category_2 %in% c("Downregulated", "Upregulated"), "gene_symbol"]
hr_synap_dementia_dt <- synap_dementia_dt[synap_dementia_dt$HumanName %in% hr_synap_dementia_genes, ]

p[[pnt]] <- ggplot(res,
          aes(x = logFC,
              y = negLogP)) +

  geom_point(aes(color = Category),
              alpha = 0.7,
              size = 1.8) +

  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "Not significant" = "grey75"
  )) +

  # Threshold lines
  geom_hline(yintercept = -log10(padj_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  geom_vline(xintercept = c(-logfc_cutoff,
                            logfc_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  # Labels
  geom_text_repel(
    data = subset(res, Category_2 %in% c("Downregulated", "Upregulated")),
    aes(label = gene_symbol),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  ) +

  labs(
    title = phenotypes_map[[pnt]],
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P)
  ) +

  theme_bw(base_size = 14) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Dementia risk gene
pnt <- "dementia_risk_gene"
res <- DGE_res[DGE_res[[pnt]] == "Yes", ]

p[[pnt]] <- ggplot(res,
          aes(x = logFC,
              y = negLogP)) +

  geom_point(aes(color = Category),
              alpha = 0.7,
              size = 1.8) +

  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "Not significant" = "grey75"
  )) +

  # Threshold lines
  geom_hline(yintercept = -log10(padj_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  geom_vline(xintercept = c(-logfc_cutoff,
                            logfc_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  # Labels
  geom_text_repel(
    aes(label = gene_symbol),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  ) +

  labs(
    title = phenotypes_map[[pnt]],
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P)
  ) +

  theme_bw(base_size = 14) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.position = "none"
  )

# Severe tinnitus gene
pnt <- "severe_tinnitus_gene"
res <- DGE_res[DGE_res[[pnt]] == "Yes", ]
p[[pnt]] <- ggplot(res,
          aes(x = logFC,
              y = negLogP)) +

  geom_point(aes(color = Category),
              alpha = 0.7,
              size = 1.8) +

  scale_color_manual(values = c(
    "Upregulated" = "#D73027",
    "Downregulated" = "#4575B4",
    "Not significant" = "grey75"
  )) +

  # Threshold lines
  geom_hline(yintercept = -log10(padj_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  geom_vline(xintercept = c(-logfc_cutoff,
                            logfc_cutoff),
              linetype = "dashed",
              linewidth = 0.7) +

  # Labels
  geom_text_repel(
    aes(label = gene_symbol),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  ) +

  labs(
    title = phenotypes_map[[pnt]],
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~Adjusted~P)
  ) +

  theme_bw(base_size = 14) +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  )

### Part 2: GO and KEGG pathway dotplots
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

# First row
top_row <- plot_grid(
  p[["synap_dementia_gene"]],
  p[["dementia_risk_gene"]],
  p[["severe_tinnitus_gene"]],
  labels = c("A", "", ""),
  ncol = 3,
  rel_widths = c(1, 1, 1.3)
)

# Second row
bottom_row <- plot_grid(
  kegg_plot, go_bp_plot,
  labels = c("B", "C"),
  ncol = 2,
  rel_widths = c(1, 1)
)

# Combine
png(paste(GO, "volcano_plot_KEGG_pathways.png", sep = "/"), width = 20*300, height = 14*300, res=300)
pdf(paste(GO, "volcano_plot_KEGG_pathways.pdf", sep = "/"), width = 20, height = 14)
plot_grid(
  top_row,
  bottom_row,
  nrow = 2,
  rel_heights = c(1, 0.7)
)
dev.off()

