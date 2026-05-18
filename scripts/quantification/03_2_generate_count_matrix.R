
PRJ_DIR <- "~/myRDS/PRJ-UNITI/DATA_ANALYSIS_UNITI/BulkRNA_KO_PACT"
QUANT <- paste(PRJ_DIR, "quantification/RSEM", sep = "/")
DGE <- paste(PRJ_DIR, "DGE/edgeR", sep = "/")
patterns <- c(".genes.results:gene_id", ".isoforms.results:transcript_id")
colDatas <- c("expected_count", "TPM", "effective_length")

dir.create(DGE, recursive = TRUE, showWarnings = FALSE)

lapply(patterns, function(x) {
  print(x)
	pattern <- strsplit(x, ":")[[1]][1]
	colID <- strsplit(x, ":")[[1]][2]
	lapply(colDatas, function(colData) {
    print(colData)
		columns <- c(colID, colData)
		fns <- list.files(QUANT, pattern=paste0("*", pattern), full.names=TRUE, recursive=T)
		dt <- lapply(fns, function(x){
		  t <- read.table(x, header = TRUE, sep = '\t', stringsAsFactors=F)[columns]
		  colnames(t) <- c(colID, gsub(pattern, "", basename(x)))
		  return(t)
		})
		tb <- Reduce(function(x,y) merge(x, y, all=T, by=colID), dt)
		# save data to csv
		# tb[is.na(tb)] <- ""
		write.csv(tb, file=paste0(DGE, "/", colData, pattern, ".csv"), row.names=F)
	})
})
