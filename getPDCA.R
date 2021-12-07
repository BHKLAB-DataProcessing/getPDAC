library(SummarizedExperiment)
library(stringr)
library(dplyr)

if(!require("GEOquery", quietly=TRUE)){
  BiocManager::install("GEOquery", ask=FALSE)
}
library(GEOquery)

if(!require("EnsDb.Hsapiens.v86", quietly=TRUE)){
  BiocManager::install("EnsDb.Hsapiens.v86", ask=FALSE)
}
library(EnsDb.Hsapiens.v86)

# script_dir <- "~/Documents/work/ORCESTRA/MetaGxPancreas-Update/code/"
script_dir <- "/pfs/getPDCA/"

source(paste0(script_dir, "BADEA.R"))
source(paste0(script_dir, "BALAGURANATH.R"))
source(paste0(script_dir, "CHEN.R"))
source(paste0(script_dir, "HAIDER.R"))
source(paste0(script_dir, "HAMIDI.R"))
source(paste0(script_dir, "JANKY.R"))
source(paste0(script_dir, "LUNARDI.R"))
source(paste0(script_dir, "PEI.R"))
source(paste0(script_dir, "UNC.R"))
source(paste0(script_dir, "YANG.R"))
source(paste0(script_dir, "ZHANG.R"))

# out_dir <- "~/Documents/work/ORCESTRA/MetaGxPancreas-Update/output/"
out_dir <- "/pfs/out/"

filenames <- c()

badea_SE <- getBadea()
saveRDS(badea_SE, paste0(out_dir, "PDAC_Badea.rds"))
filenames <- c(filenames, "PDAC_Badea.rds")
rm(badea_SE)

BALAGURANATH_SE <- getBalagurnath()
saveRDS(BALAGURANATH_SE, paste0(out_dir, "PDAC_Balaguranath.rds"))
filenames <- c(filenames, "PDAC_Balaguranath.rds")
rm(BALAGURANATH_SE)

CHEN_SE <- getChen()
saveRDS(CHEN_SE, paste0(out_dir, "PDAC_Chen.rds"))
filenames <- c(filenames, "PDAC_Chen.rds")
rm(CHEN_SE)

HAIDER_SE <- getHaider()
saveRDS(HAIDER_SE, paste0(out_dir, "PDAC_Haider.rds"))
filenames <- c(filenames, "PDAC_Haider.rds")
rm(HAIDER_SE)

HAMIDI_SE <- getHamidi()
saveRDS(HAMIDI_SE, paste0(out_dir, "PDAC_Hamidi.rds"))
filenames <- c(filenames, "PDAC_Hamidi.rds")
rm(HAMIDI_SE)

JANKY_SE <- getJanky()
saveRDS(JANKY_SE, paste0(out_dir, "PDAC_Janky.rds"))
filenames <- c(filenames, "PDAC_Janky.rds")
rm(JANKY_SE)

LUNARDI_SE <- getLunardi()
saveRDS(LUNARDI_SE, paste0(out_dir, "PDAC_Lunardi.rds"))
filenames <- c(filenames, "PDAC_Lunardi.rds")
rm(LUNARDI_SE)

PEI_SE <- getPei()
saveRDS(PEI_SE, paste0(out_dir, "PDAC_Pei.rds"))
filenames <- c(filenames, "PDAC_Pei.rds")
rm(PEI_SE)

UNC_SE <- getUNC()
saveRDS(UNC_SE, paste0(out_dir, "PDAC_UNC.rds"))
filenames <- c(filenames, "PDAC_UNC.rds")
rm(UNC_SE)

YANG_SE <- getYang()
saveRDS(YANG_SE, paste0(out_dir, "PDAC_Yang.rds"))
filenames <- c(filenames, "PDAC_Yang.rds")
rm(YANG_SE)

ZHANG_SE <- getZhang()
saveRDS(ZHANG_SE, paste0(out_dir, "PDAC_Zhang.rds"))
filenames <- c(filenames, "PDAC_Zhang.rds")
rm(ZHANG_SE)


filenames_df <- data.frame(matrix(nrow=length(filenames), ncol=3))
colnames(filenames_df) <- c("filename", "doi", "download_link")
filenames_df$filename <- filenames
write.csv(filenames_df, paste0(out_dir, "data_list.csv"))
