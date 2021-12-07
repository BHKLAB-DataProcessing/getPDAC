getJanky <- function(){
  ## JANKY
  JANKY = getGEO("GSE62165", GSEMatrix = TRUE)
  JANKY = JANKY$GSE62165_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(JANKY)
  gene_table = JANKY@featureData@data
  gene_table = gene_table[,c("ID","Gene Symbol")]
  
  gene_table$gene_symbol = str_split_fixed(gene_table$`Gene Symbol`, "///",2)[,1]
  gene_table = gene_table[,c("ID","gene_symbol")]
  gene_table = gene_table[which(gene_table$gene_symbol != ""),]
  gene_table = gene_table[which(gene_table$gene_symbol != "---"),]
  gene_table = gene_table[-which(duplicated(gene_table$gene_symbol)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$gene_symbol
  expr$ID = NULL
  expr$gene_symbol = NULL
  
  ## clinical data
  clinData = JANKY@phenoData@data
  clinData = clinData[which(clinData$`tissue:ch1` == "pancreatic tumor"),]
  clinData = clinData[,c("geo_accession", "tissue:ch1","Stage:ch1")]
  clinData$sample_type = "Primary_tumor"
  
  clinData = clinData[,c("geo_accession","Stage:ch1","sample_type")]
  colnames(clinData) = c("sample_id","tumor_stage","sample_type")
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
    name = "JANKY",
    lab = "",
    contact = "",
    title = "Janky et al, BMC Cancer 2016",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62165",
    other = list(
      ExperimentalStrategy = "Affymetrix Human Genome U219 Array",
      Units = "Normalized-Expression",
      normalisation = "RMA"))
  
  ## create summarized experiments
  JANKY_SE = SummarizedExperiment(
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
  rowRanges(JANKY_SE) = gene_ranges
  return(JANKY_SE)
}
