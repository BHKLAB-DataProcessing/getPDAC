getHamidi <- function(){
  ## HAMIDI
  HAMIDI = getGEO("GSE77858", GSEMatrix = TRUE)
  HAMIDI = HAMIDI$GSE77858_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(HAMIDI)
  gene_table = HAMIDI@featureData@data
  gene_table = gene_table[which(gene_table$GENE_SYMBOL != ""),]
  
  gene_table = gene_table[,c("ID","GENE_SYMBOL")]
  # gene_table$gene_symbol = gsub(pattern = " ", replacement = "",
  # x = gene_table$gene_symbol)
  gene_table = gene_table[-which(duplicated(gene_table$GENE_SYMBOL)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$GENE_SYMBOL
  expr$ID = NULL
  expr$GENE_SYMBOL = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = HAMIDI@phenoData@data
  clinData = clinData[which(clinData$`morphology:ch2` == "Tumor"),]
  clinData = clinData[,c("geo_accession","morphology:ch2")]
  colnames(clinData) = c("sample_id","sample_type")
  rownames(clinData) = NULL
  
  ## geneRanges
  edb <- EnsDb.Hsapiens.v86
  gene_ranges = genes(edb,
                      columns = c(listColumns(edb, "gene"), "entrezid"),
                      filter = GeneNameFilter(rownames(expr)),
                      return.type = "GRanges")
  gene_ranges = gene_ranges[-which(duplicated(gene_ranges$gene_name))]
  
  ## remove expression of genes which are deprecated
  ### NOTE: Do we want to keep deprecated genes with NAs in description?
  expr = expr[which(rownames(expr) %in% gene_ranges$gene_name),
              which(colnames(expr) %in% clinData$sample_id)]
  
  expr = expr %>% arrange(factor(rownames(expr), levels = gene_ranges$gene_name))
  
  ## metadata
  experimentData = MIAME(
    name = "HAMIDI",
    lab = "",
    contact = "",
    title = "Hamidi H et al",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77858",
    other = list(
      ExperimentalStrategy = "Agilent-012097 Human 1A Microarray (V2) G4110B",
      Units = "Normalized-Expression",
      normalisation = "Lowess Normalization"
    )
  )
  
  ## create summarized experiments
  HAMIDI_SE = SummarizedExperiment(
    assay = expr,
    rowData = rownames(expr),
    colData = clinData,
    metadata = list(
      experimentData,
      # replicating MetaGxPancreas TCGA format
      annotation = "",
      protocolData = AnnotatedDataFrame()
    )
  )
  rowRanges(HAMIDI_SE) = gene_ranges
  return(HAMIDI_SE)
}
