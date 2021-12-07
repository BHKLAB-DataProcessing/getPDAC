getLunardi <- function(){
  ## LUNARDI
  LUNARDI = getGEO("GSE55643", GSEMatrix = TRUE)
  LUNARDI = LUNARDI$GSE55643_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(LUNARDI)
  
  gene_table = LUNARDI@featureData@data
  gene_table = gene_table[which(gene_table$GENE_SYMBOL != ""),]
  gene_table = gene_table[,c("ID","GENE_SYMBOL")]
  gene_table = gene_table[-which(duplicated(gene_table$GENE_SYMBOL)),]
  
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$GENE_SYMBOL
  expr$ID = NULL
  expr$GENE_SYMBOL = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = LUNARDI@phenoData@data
  clinData = clinData[which(clinData$`tissue:ch1`== "Pancreas tumour"),]
  clinData = clinData[,c("geo_accession", "tissue:ch1", "gender:ch1")]
  
  colnames(clinData) = c("sample_id","sample_type","gender")
  clinData$unique_patient_id = clinData$sample_id
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
    name = "LUNARDI",
    lab = "",
    contact = "",
    title = "Lunardi S et al, 2014, Oncotarget",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55643",
    other = list(
      ExperimentalStrategy = "Agilent-014850 Whole Human Genome Microarray 4x44K G4112F",
      Units = "Normalized-Expression",
      normalisation = "Percentile normalization 75th"
    )
  )
  
  ## create summarized experiments
  LUNARDI_SE = SummarizedExperiment(
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
  rowRanges(LUNARDI_SE) = gene_ranges
  return(LUNARDI_SE)
}
