getZhang <- function(){
  ## ZHANG
  ZHANG = getGEO("GSE28735", GSEMatrix = TRUE)
  ZHANG = ZHANG$GSE28735_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(ZHANG)
  gene_table = ZHANG@featureData@data
  gene_table = gene_table[which(gene_table$gene_assignment != "---"),]
  
  gene_table$gene_symbol = str_split_fixed(gene_table$gene_assignment, "//", 3)[,2]
  gene_table = gene_table[,c("ID","gene_symbol")]
  gene_table = gene_table[which(gene_table$gene_symbol != ""),]
  gene_table$gene_symbol = gsub(pattern = " ", replacement = "",
                                x = gene_table$gene_symbol)
  gene_table = gene_table[-which(duplicated(gene_table$gene_symbol)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$gene_symbol
  expr$ID = NULL
  expr$gene_symbol = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = ZHANG@phenoData@data
  clinData = clinData[which(clinData$`tissue:ch1` == "T"),]
  clinData = clinData[,c("geo_accession","cancer_death:ch1",
                         "survival_month:ch1","tissue:ch1")]
  clinData$days_to_death = 30*as.numeric(clinData$`survival_month:ch1`)
  
  clinData$vital_status = NA
  clinData$vital_status[which(clinData$`cancer_death:ch1` == 1)] = "deceased"
  clinData$vital_status[which(clinData$`cancer_death:ch1` == 0)] = "living"
  
  clinData$sample_type = "Tumor"
  
  clinData = clinData[,c("geo_accession","days_to_death",
                         "vital_status","sample_type")]
  colnames(clinData) = c("sample_id","days_to_death",
                         "vital_status","sample_type")
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
    name = "ZHANG",
    lab = "",
    contact = "",
    title = "Zhang et al, PLoS One 2012",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735",
    other = list(
      ExperimentalStrategy = "Affymetrix GeneChip Human Gene 1.0 ST arrays",
      Units = "Normalized-Expression",
      normalisation = "RMA"
    )
  )
  
  ## create summarized experiments
  ZHANG_SE = SummarizedExperiment(
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
  rowRanges(ZHANG_SE) = gene_ranges
  return(ZHANG_SE)
}