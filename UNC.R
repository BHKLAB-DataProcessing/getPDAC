getUNC <- function(){
  ## UNC
  unc = getGEO("GSE71729", GSEMatrix = TRUE)
  unc = unc$GSE71729_series_matrix.txt.gz
  
  ## expr data
  expr = as.data.frame(exprs(unc))
  
  ## clinical data
  clinData = unc@phenoData@data
  clinData = clinData[which(!is.na(clinData$`survival_months:ch2`)),]
  clinData = clinData[,c("geo_accession","tumor_subtype_0na_1classical_2basal:ch2",
                         "survival_months:ch2",
                         "death_event_1death_0censor:ch2")]
  clinData$days_to_death = 30*as.numeric(clinData$`survival_months:ch2`)
  
  clinData$vital_status = NA
  clinData$vital_status[which(clinData$`death_event_1death_0censor:ch2` == 1)] = "deceased"
  clinData$vital_status[which(clinData$`death_event_1death_0censor:ch2` == 0)] = "living"
  
  clinData$sample_type = "Primary_tumor"
  
  clinData$tumor_subtype = NA
  clinData$tumor_subtype[which(clinData$`tumor_subtype_0na_1classical_2basal:ch2` == 1)] = "classical"
  clinData$tumor_subtype[which(clinData$`tumor_subtype_0na_1classical_2basal:ch2` == 2)] = "basal"
  
  clinData = clinData[,c("geo_accession","days_to_death",
                         "vital_status","sample_type","tumor_subtype")]
  colnames(clinData)[1] = "sample_id"
  clinData$unique_patient_id = clinData$sample_id
  
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
    name = "Richard,,Moffitt",
    lab = "UNC",
    contact = "",
    title = "Moffitt et al, Nat Genet 2015",
    abstract = "",
    pubMedIds = "26343385",
    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71729",
    other = list(
      ExperimentalStrategy = "micro-array",
      Units = "Normalized-Expression",
      normalisation = "Quantile"
    )
  )
  
  ## create summarized experiments
  UNC_SE = SummarizedExperiment(
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
  rowRanges(UNC_SE) = gene_ranges
  return(UNC_SE)
}
