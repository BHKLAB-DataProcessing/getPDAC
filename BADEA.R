getBadea <- function(){
  ## badea
  badea = getGEO("gse15471", GSEMatrix = TRUE)
  badea = badea$GSE15471_series_matrix.txt.gz
  
  ## expr data
  expr = exprs(badea)
  gene_table = badea@featureData@data
  gene_table = gene_table[which(gene_table$`Gene Symbol` != ""),]
  
  
  gene_table$gene_symbol = str_split_fixed(gene_table$`Gene Symbol`, "///", 3)[,1]
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
  clinData = badea@phenoData@data
  clinData = clinData[which(clinData$characteristics_ch1.1 == "sample: tumor"),]
  clinData = clinData[,c("geo_accession", "characteristics_ch1.1")]
  clinData$sample_type = "tumor"
  
  clinData = clinData[,c("geo_accession","sample_type")]
  colnames(clinData)[1] = c("sample_id")
  rownames(clinData) = NULL
  clinData$characteristics_ch1.1 = NULL
  
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
    name = "BADEA",
    lab = "",
    contact = "",
    title = "Badea et al, Hepatogastroenterology 2008",
    abstract = "",
    pubMedIds = "",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse15471",
    other = list(
      ExperimentalStrategy = "Affymetrix Human Genome U133 Plus 2.0 Array",
      Units = "Normalized-Expression",
      normalisation = "RMA"
    )
  )
  
  ## create summarized experiments
  badea_SE = SummarizedExperiment(
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
  rowRanges(badea_SE) = gene_ranges
  return(badea_SE)
}
