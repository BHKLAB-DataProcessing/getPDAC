getBalagurnath <- function(){
  ## BALAGURANATH
  BALAGURANATH = getGEO("GSE11838", GSEMatrix = TRUE)
  BALAGURANATH = BALAGURANATH$GSE11838_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(BALAGURANATH)
  gene_table = BALAGURANATH@featureData@data
  gene_table = gene_table[which(gene_table$ORF != "--"),]
  gene_table = gene_table[which(gene_table$ORF != ""),]
  
  gene_table$gene_symbol = gene_table$ORF
  gene_table = gene_table[,c("ID","gene_symbol")]
  gene_table$gene_symbol = gsub(pattern = " ", replacement = "",
                                x = gene_table$gene_symbol)
  gene_table = gene_table[-which(duplicated(gene_table$gene_symbol)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$gene_symbol
  expr$ID = NULL
  expr$gene_symbol = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = BALAGURANATH@phenoData@data
  clinData = clinData[which(clinData$source_name_ch1 == "primary pancreatic tumor"),]
  clinData$sample_type = "tumor"
  clinData = clinData[,c("geo_accession", "sample_type")]
  colnames(clinData)[1] = c("sample_id")
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
    name = "BALAGURANATH",
    lab = "",
    contact = "",
    title = "Balagurunathan et al, Mol Cancer Ther 2008",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11838",
    other = list(
      ExperimentalStrategy = "Human 1A Microarray G4110A/G4110B",
      Units = "Normalized-Expression",
      normalisation = "Median normalization"
    )
  )
  
  ## create summarized experiments
  BALAGURANATH_SE = SummarizedExperiment(
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
  rowRanges(BALAGURANATH_SE) = gene_ranges
  return(BALAGURANATH_SE)
}
