getPei <- function(){
  ## PEI
  PEI = getGEO("GSE16515", GSEMatrix = TRUE)
  PEI = PEI$GSE16515_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(PEI)
  gene_table = PEI@featureData@data
  
  gene_table$gene_symbol = str_split_fixed(gene_table$`Gene Symbol`, "///", 2)[,1]
  gene_table = gene_table[,c("ID","gene_symbol")]
  gene_table$gene_symbol = gsub(pattern = " ", replacement = "",
                                x = gene_table$gene_symbol)
  gene_table = gene_table[which(gene_table$gene_symbol != ""),]
  gene_table = gene_table[-which(duplicated(gene_table$gene_symbol)),]
  expr = merge(gene_table, expr, by = "row.names")
  
  rownames(expr) = expr$gene_symbol
  expr$ID = NULL
  expr$gene_symbol = NULL
  expr$Row.names = NULL
  
  ## clinical data
  clinData = PEI@phenoData@data
  clinData = clinData[which(clinData$`tissue:ch1` == "Tumor Tissue in Pancreatic Cancer Sample"),]
  clinData = clinData[,c("geo_accession", "tissue:ch1",
                         "sex:ch1","age:ch1")]
  
  clinData$sample_type = "Tumor"
  
  clinData = clinData[,c("geo_accession","sample_type",
                         "sex:ch1","age:ch1")]
  colnames(clinData) = c("sample_id","sample_type",
                         "gender","age")
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
    name = "PEI",
    lab = "",
    contact = "",
    title = "Pei et al, Cancer Cell 2009",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse16515",
    other = list(
      ExperimentalStrategy = "Affymetrix Human Genome U133 Plus 2.0 Array",
      Units = "RNA Counts",
      normalisation = ""
    )
  )
  
  ## create summarized experiments
  PEI_SE = SummarizedExperiment(
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
  rowRanges(PEI_SE) = gene_ranges
  return(PEI_SE)
}
