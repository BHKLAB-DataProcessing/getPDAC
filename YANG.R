getYang <- function(){
  ## yang
  yang = getGEO("GSE62452", GSEMatrix = TRUE)
  yang = yang$GSE62452_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(yang)
  gene_table = yang@featureData@data
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
  clinData = yang@phenoData@data
  clinData = clinData[which(!is.na(clinData$`survival months:ch1`)),]
  clinData = clinData[,c("geo_accession", "grading:ch1","Stage:ch1",
                         "survival months:ch1","survival status:ch1",
                         "tissue:ch1")]
  clinData$days_to_death = 30*as.numeric(clinData$`survival months:ch1`)
  
  clinData$vital_status = NA
  clinData$vital_status[which(clinData$`survival status:ch1` == 1)] = "deceased"
  clinData$vital_status[which(clinData$`survival status:ch1` == 0)] = "living"
  
  clinData$sample_type = "Primary_tumor"
  
  clinData = clinData[,c("geo_accession","days_to_death",
                         "vital_status","sample_type",
                         "Stage:ch1","grading:ch1")]
  colnames(clinData) = c("sample_id","days_to_death",
                         "vital_status","sample_type",
                         "tumor_stage","tumor_grade")
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
    name = "YANG",
    lab = "",
    contact = "",
    title = "Yang et al, 2016, Cancer Research",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452",
    other = list(
      ExperimentalStrategy = "Affymetrix GeneChip Human Gene 1.0 ST",
      Units = "Normalized-Expression",
      normalisation = "RMA"
    )
  )
  
  ## create summarized experiments
  YANG_SE = SummarizedExperiment(
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
  rowRanges(YANG_SE) = gene_ranges
  return(YANG_SE)
}
