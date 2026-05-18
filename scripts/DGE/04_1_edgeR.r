library(edgeR)
library(ComplexHeatmap)


# nolint

PRJ_DIR <- "~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
QUANT <- paste(PRJ_DIR, "quantification/RSEM", sep = "/")
DGE <- paste(PRJ_DIR, "DGE/edgeR", sep = "/")
REF <- "~/myRDS/PRJ-MRefD/databases"
targetp <- paste(PRJ_DIR, "config/target.tsv", sep="/")
count_mat <- paste(DGE, "expected_count.genes.results.csv", sep = "/")
gene_info_p <- paste(REF, "Gencode/v48/gencode.v48.basic.annotation.genes.tsv", sep = "/")
FDR <- 0.05
logFC <- 1.5
case <- "LNCaP_PACT_KO"
control <- "LNCaP_WT"
pair <- "LNCaP_PACT_KO - LNCaP_WT"

filter_gene <- function(counts, gene_info_p) {
  message("Removing genes without gene symbol for real data ...")
  gene_info <- read.table(gene_info_p, header = T, na.strings = c(".", "NA", ""))
  m <- match(rownames(counts), gene_info$geneID)
  keep <- !is.na(m)
  counts <- counts[keep, ]
  annot <- gene_info[m[keep], ]
  rownames(annot) <- rownames(counts)
  counts <- data.matrix(counts)
  return(list("counts" = counts,
    "annot" = annot))
}

construct_model_matrix <- function(targetp) {
  message("Reading target ...")
  target <- read.table(targetp, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  group <- factor(target$group)
  samples <- target$sampleID
  
  message("Constructing model matrix ...")
  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  
  return(list("group" = group,
    "samples" = samples,
    "design" = design))
}

DGE_analysis <- function(counts, annot, group, design) {
  
  message("Constructing DGElist object ...")
  y <- DGEList(counts = counts, genes = annot, group = group)
  
  message("Filtering exons with low mapping reads ...")
  keep <- filterByExpr(y, group = group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  message("Normalizing lib sizes ...")
  y <- normLibSizes(y)
  
  message("Estimating dispersion ...")
  y <- estimateDisp(y, design, robust=TRUE)
  
  message("Fitting GLM-QL model for design matrix ...")
  fit <- glmQLFit(y, design, robust=TRUE)
  
  return(list("fit" = fit, "y" = y))
}

DGE_res <- function(fit, design, pair) {
  cmd <- paste0("makeContrasts(", pair, ",levels = design)")
  contr <- eval(parse(text=cmd))
  print(contr)
  qlf <- glmQLFTest(fit, contrast = contr)
  tt <- topTags(qlf, n=Inf, adjust.method = "BH")$table
  # add tt gene info
  return(tt)
}

# Only keep the basic genes
# Process gene counts and gene annotation
message("Processing gene counts and gene annotation")
counts <- read.csv(count_mat, row.names = 1)
counts <- round(counts) # used for edgeR output
flt_counts <- filter_gene(counts, gene_info_p)

# Construct model matrix
message("Constructing design matrix for contrast ...")
mat <- construct_model_matrix(targetp)

# Fit count model
message("Fitting count models ...")
counts <- flt_counts$counts[mat$samples] # match samples with design matrix
count_model_fit <- DGE_analysis(counts,
  flt_counts$annot,
  mat$group,
  mat$design)

# Running DEU
message("Running DEU ...")
res <- DGE_res(count_model_fit$fit, mat$design, pair)
write.csv(res, paste(DGE, "expected_count.genes.results.DGE.csv", sep = "/"), row.names = F)

# Volcano plot
pdf(paste(DGE, "volcano_plot_genes.pdf", sep = "/"), width = 3, height = 4, dpi = 300)
plot(tt$logFC, -log10(tt$PValue),
  type="n",
  xlab=paste0("logFC ", pair), ylab="-log10(p.value)")
text(tt$logFC, -log10(tt$PValue), labels = tt$gene_symbol, cex = 0.5)
abline(h = -log10(tt$PValue[sum(tt$FDR < FDR)]), col="red")
dev.off()

# MDS plot
logCPM <- cpm(count_model_fit$y, log = TRUE, prior.count = 2)
cols <- as.factor(mat$group)
plotMDS(logCPM, col = cols, pch = 19)
legend("topright",
  legend = levels(group),
  col = 1:length(levels(group)),
  pch = 19)

# Heatmap
vars <- apply(logCPM, 1, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:50]

mat <- logCPM[top_genes, ]
mat_z <- t(scale(t(mat)))


Heatmap(mat_z,
        name = "Z-score",
        column_split = group)

# Gene ontology

######################################################################################################
# # DE at gene level
# .libPaths("/scratch/lv15/lpl913/.R/3.6.1")
# library(edgeR)
# library(ggplot2)
# library(rtracklayer)
# library(data.table)
# library(goseq)
# library(GO.db)
# library(org.Hs.eg.db)
# library(optparse)

# # input <- "/g/data3/yo4/Cancer-Epigenetics/davban/thinh/projects/bpipe/RNAseq_hg38/EdgeR/input/expected_count.genes.results.csv"
# # output <- "/g/data3/yo4/Cancer-Epigenetics/davban/thinh/projects/bpipe/RNAseq_hg38/EdgeR/output/"
# # gtf <- "/g/data3/yo4/Cancer-Epigenetics/davban/thinh/data/RNASeq/GTF/GRCh38/gencode/gencode.v33.annotation.gtf"
# # # annotation file downloaded from http://feb2014.archive.ensembl.org/biomart/martview/; Ensembl 75: Feb 2014 (GRCh37.p13) - patched/updated gene set Sep 2013
# # # has following columns: ensembl.gene.ID, chr, start, end, strand, gene description, HGNC.symbol (external), entrez.gene.ID (external)
# # biomartp_g <- "/g/data3/yo4/Cancer-Epigenetics/davban/thinh/data/RNASeq/biomart/gene_hg38_p13.tsv"
# # designp <- "/g/data3/yo4/Cancer-Epigenetics/davban/thinh/projects/bpipe/RNAseq_hg38/config/design.tsv"
# # FWER <- 0.05
# # logFC <- 1.5
# # nread <- 3

# option_list <- list(
#   make_option("--input", type="character", default=NULL, help="dataset path of EdgeR input", metavar="character"),
#   make_option("--designp", type="character", default="out.txt", help="path to design document", metavar="character"),
#   make_option("--output", type="character", default=NULL, help="where to store output", metavar="character"),
#   make_option("--gtf", type="character", default=NULL, help="gencode gtf annotation", metavar="character"),
#   make_option("--biomartp_g", type="character", default=NULL, help="biomart annotation", metavar="character"),
#   make_option("--FWER", type="character", default=NULL, help="FWER", metavar="character"),
#   make_option("--logFC", type="character", default=NULL, help="logFC", metavar="character"),
#   make_option("--nread", type="character", default=NULL, help="nread", metavar="character")
# ) 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)

# input <- opt$input
# designp <- opt$designp
# output <- opt$output
# gtf <- opt$gtf
# biomartp_g <- opt$biomartp_g
# FWER <- as.numeric(opt$FWER)
# logFC <- as.numeric(opt$logFC)
# nread <- as.numeric(opt$nread)

# # read in various datasets
# counts <- read.csv(input, row.names=1)
# design <- read.table(designp, sep = "\t", header = T, stringsAsFactor = F)
# base <- design[design$Base==1, "Group"][1]
# notbase <- design[design$Base==0, "Group"][1]

# # round the rsem gene expected counts values to the nearest integer to input into edgeR
# counts <- round(counts)

# # matrix of counts with ENSGxxxxxxxx tags
# counts <- data.matrix(counts) 
# sample.names <- colnames(counts)

# rownames(design) <- design$Sample

# # class <- sapply(sample.names, function(x) {strsplit(x, "_")[[1]][1] } )
# # Filter out ENSGxxxx tags whose coverage is so low that any group differences aren't truly "real". 
# # filter out tags whose rowcount <= degrees of freedom.
# counts <- counts[rowSums(counts) >= nread,]

# # set up the design matrix. The given specification sets up LNCaP as the non-reference level.
# group <- factor(design[sample.names, "Group"])
# # The level which is chosen for the reference level or base level. This is the level which is 
# # contrasted against, and by default this is simply the first level alphabetically. 
# # Since we want PrEC to be the base level we use the relevel function to set this manually.
# group <- relevel(group, base)

# design_matrix <- model.matrix(~group) 
# y <- DGEList(counts=counts, group=group)

# # The calcNormFactors() function normalizes for RNA composition by finding a set of scaling factors 
# # for the library sizes that minimize the log-fold changes between the samples for most genes. 
# # The default method for computing these scale factors uses a trimmed mean of M-values 
# # (TMM) between each pair of samples. It is based on the hypothesis that most genes are not DE.
# y <- calcNormFactors(y)

# pdf(file.path(output, "MDS_plots_genes.pdf"))
# # produce a plot showing the sample relations based on multidimensional scaling.
# plotMDS(y, method="bcv", col=as.numeric(y$samples$group), main=paste0("MDS plot of genes: ", paste0(unique(group), collapse="/", " samples")))
# dev.off()

# # GLM estimates of variance (dispersion)
# # Fitting a model in edgeR takes several steps. 
# # First, you must fit the common dispersion. Then you need to fit a trended model 
# # (if you do not fit a trend, the default is to use the common dispersion as a trend). 
# # Then you can fit the tagwise dispersion which is a function of this model.
# y <- estimateGLMCommonDisp(y, design_matrix, verbose=TRUE)
# y <- estimateGLMTrendedDisp(y, design_matrix)
# y <- estimateGLMTagwiseDisp(y, design_matrix)

# # plot the genewise biological coefficient of variation (BCV) against gene abundance 
# # (in log2 counts per million). 
# # It displays the common, trended and tagwise BCV estimates.
# pdf(file.path(output, "BCV_plot_genes.pdf"))
# plotBCV(y)
# dev.off()

# # GLM testing for differential expression

# # fitting linear model
# fit <- glmFit(y, design_matrix)
# # find the tags that are interesting by using a LRT (Likelihood Ratio Test)
# log.ratio.test <- glmLRT(fit, coef=2)
# # or log.ratio.test <- glmLRT(fit, contrast=c(0,1))

# # Top table
# # the default method used to adjust p-values for multiple testing is BH.
# # Using BH and applying a criteria of abs(logFC) > 4 and FDR < 0.05
# # 16.5% of the total genes were reported as differentially expressed. This is
# # quite a large percentage so applied Bonferroni method to adjust p-values
# # which is a more conservative approach. This reduced the percentage of genes
# # being reported as differentially expressed to be 10%.
# tt <- topTags(log.ratio.test, n=nrow(counts), adjust.method="bonferroni")$table

# # annotation FROM GTF FILE
# gtf.anno <- import(gtf)
# m <- match(rownames(tt), mcols(gtf.anno)['gene_id'][[1]])

# # assign the rownames(tt) as the gene_id as more specific with version number at end
# # as originates from original gtf
# tt$gene.id <- rownames(tt)
# tt$gene.symbol <- mcols(gtf.anno[m])['gene_name'][[1]] 
# tt$chr <- as.character(seqnames(gtf.anno[m]))
# tt$start <- start(gtf.anno[m])
# tt$end <- end(gtf.anno[m])
# tt$strand <- as.character(strand(gtf.anno)[m])
# tt$gene.type <- mcols(gtf.anno[m])['gene_type'][[1]] 

# # annotation file downloaded from http://dec2013.archive.ensembl.org/biomart/martview/
# # has following columns: ensembl.gene.ID, chr, start, end, strand, description, HGNC.symbol, entrez.gene.ID
# DT <- fread(biomartp_g)
# setnames(DT, gsub(" ", ".", colnames(DT)))
# setkey(DT, Ensembl.Gene.ID)
# m <- match(gsub("\\..*", "", rownames(tt)), DT$Ensembl.Gene.ID)
# tt$description <- DT$description[m]
# tt$entrez.gene.id <- DT$entrez.gene.ID[m]

# # Volcano plot
# pdf(file.path(output, "volcano_plot_genes.pdf"), width=15, height=20)
# plot(tt$logFC, -log10(tt$PValue), type="n", xlab=paste0(base, " <- -> ", notbase, " logFC", ylab="-log10(p.value)"), main=paste0("Volcano plot of genes: ", notbase, " vs ", base))
# text(tt$logFC, -log10(tt$PValue), labels=tt$gene.symbol, cex=0.5)
# abline(h=-log10(tt$PValue[sum(tt$FWER < FWER)]), col="red")
# dev.off()

# # The function plotSmear generates a plot of the tagwise log-fold-changes against 
# # log-cpm (analogous to an MA-plot for microarray data). 
# # DE tags are highlighted on the plot
# pdf(file.path(output, "smear_plot_genes.pdf"))
# de2 <- decideTestsDGE(log.ratio.test, adjust.method="bonferroni", p.value=0.05, lfc=logFC)
# de2tags <- rownames(y)[as.logical(de2)]
# plotSmear(log.ratio.test, de.tags=de2tags, main=paste0("smear plot of genes with p.value<", FWER, " and LFC=", logFC," cutoffs"))
# abline(h = c(-logFC, logFC), col = "blue")
# dev.off()

# # add DGE.status column <- UP, DOWN, NC
# tt$DGE.status <- "NC"
# # defining significant as FWER < 0.05 and abs(logFC) > 4 for UP/DOWN
# tt[((tt$FWER < FWER)&(tt$logFC > logFC)),]$DGE.status <- "UP"
# tt[((tt$FWER < FWER)&(tt$logFC < -logFC)),]$DGE.status <- "DOWN"
# write.table(tt, file.path(output, "ALL_GENE_DGE.tsv"), sep="\t", quote=F, row.names=F)

# # defining significant as FWER < 0.05 and abs(logFC) > 4
# sigtt <- tt[((tt$FWER < FWER)&(abs(tt$logFC) > logFC)),]
# write.table(sigtt, file.path(output, paste0("GENE_DGE_", FWER, "_", logFC, "_", nread, ".tsv")), sep="\t", quote=F, row.names=F)


