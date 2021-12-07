getChen <- function(){
  ## chen
  chen = getGEO("GSE57495", GSEMatrix = TRUE)
  chen = chen$GSE57495_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(chen)
  gene_table = chen@featureData@data
  gene_table = gene_table[,c("ID","GeneSymbol")]
  gene_table = gene_table[-which(gene_table$GeneSymbol == ""),]
  gene_table = gene_table[-which(duplicated(gene_table$GeneSymbol)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$GeneSymbol
  expr$ID = NULL
  expr$GeneSymbol = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = chen@phenoData@data
  clinData = clinData[,c("geo_accession", "overall survival (month):ch1",
                         "Stage:ch1","vital.status:ch1")]
  clinData$days_to_death = 30*as.numeric(clinData$`overall survival (month):ch1`)
  
  clinData$vital_status = NA
  clinData$vital_status[which(clinData$`vital.status:ch1` == "DEAD")] = "deceased"
  clinData$vital_status[which(clinData$`vital.status:ch1` == "ALIVE")] = "living"
  
  clinData$sample_type = "Primary_tumor"
  
  clinData = clinData[,c("geo_accession","days_to_death",
                         "vital_status","sample_type",
                         "Stage:ch1")]
  colnames(clinData) = c("sample_id","days_to_death",
                         "vital_status","sample_type",
                         "tumor_stage")
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
    name = "CHEN",
    lab = "",
    contact = "",
    title = "Chen et al, PLoS One 2015",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57495",
    other = list(
      ExperimentalStrategy = "Affymetrix,Rosetta/Merck RSTA Custom 2.0",
      Units = "Normalized-Expression",
      normalisation = "RMA"
    )
  )
  
  ## create summarized experiments
  CHEN_SE = SummarizedExperiment(
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
  rowRanges(CHEN_SE) = gene_ranges
  return(CHEN_SE)
}
