getHaider <- function(){
  ## HAIDER
  HAIDER = getGEO("GSE56560", GSEMatrix = TRUE)
  HAIDER = HAIDER$GSE56560_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(HAIDER)
  gene_table = HAIDER@featureData@data
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
  
  ## clinical data
  clinData = HAIDER@phenoData@data
  clinData = clinData[which(clinData$source_name_ch1 == "PDAC"),]
  clinData = clinData[,c("geo_accession", "cellularity:ch1","age:ch1",
                         "grading:ch1","Sex:ch1","pt:ch1")]
  
  clinData$sample_type = "Tumor"
  colnames(clinData) = c("sample_id","cellularity",
                         "age","tumor_grade",
                         "gender","T","sample_type")
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
    name = "HAIDER",
    lab = "",
    contact = "",
    title = "haider et al, Genome medicine, 2014",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56560",
    other = list(
      ExperimentalStrategy = "Affymetrix Human Exon 1.0 ST Array",
      Units = "Normalized-Expression",
      normalisation = "RMA and quantile normalisation"
    )
  )
  
  ## create summarized experiments
  HAIDER_SE = SummarizedExperiment(
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
  rowRanges(HAIDER_SE) = gene_ranges
  return(HAIDER_SE)
}
