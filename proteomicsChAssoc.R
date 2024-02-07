# Content: Proteomics and CH association analysis

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

indelFeatures <- grep("indel_", colnames(clinical), value = TRUE)
cnvFeatures <- grep("cnv_", colnames(clinical), value = TRUE)
olinkFeatures <- colnames(data)[2:ncol(data)]

# Add CH class
geneGroup <- list(
  EpiM = c("DNMT3A", "TET2", "ASXL1", "IDH1", "IDH2", "EZH2", "BCOR", "BCORL1"),
  DnaDRP = c("PPM1D", "TP53", "CHEK1", "CHEK2", "ATM"),
  SF = c("SF3B1", "SRSF2", "U2AF1", "ZRSR2"),
  GFS = c("JAK2", "CBL", "GNB1", "GNAS")
)
clinical <- clinical %>%
  mutate(
    indel_EpiM = rowSums(clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$EpiM))) == "Yes") %>%
      {
        ifelse(. == 0, "No", "Yes")
      },
    indel_EpiM_maxVaf = clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$EpiM, "_vaf"))) %>% apply(1, max),
    indel_DnaDRP = rowSums(clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$DnaDRP))) == "Yes") %>%
      {
        ifelse(. == 0, "No", "Yes")
      },
    indel_DnaDRP_maxVaf = clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$DnaDRP, "_vaf"))) %>% apply(1, max),
    indel_SF = rowSums(clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$SF))) == "Yes") %>%
      {
        ifelse(. == 0, "No", "Yes")
      },
    indel_SF_maxVaf = clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$SF, "_vaf"))) %>% apply(1, max),
    indel_GFS = rowSums(clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$GFS))) == "Yes") %>%
      {
        ifelse(. == 0, "No", "Yes")
      },
    indel_GFS_maxVaf = clinical %>% dplyr::select(any_of(paste0("indel_", geneGroup$GFS, "_vaf"))) %>% apply(1, max)
  )

clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
indelFeatures <- grep("indel_", colnames(clinical), value = TRUE)
cnvFeatures <- grep("cnv_", colnames(clinical), value = TRUE)
chFeatures <- c(indelFeatures, cnvFeatures)

# full features analysis
fullFeatureAnalysis <- c("indel_chPD", "cnv_allCnv", "indel_or_cnv", "indel_maxVaf")
fullFeatureAnalysis <- c(fullFeatureAnalysis, paste0("indel_", names(geneGroup)), paste0("indel_", names(geneGroup), "_maxVaf"))
if (!chFeatures[1] %in% fullFeatureAnalysis) {
  # keep only important features
  clinicalFeatures <- setdiff(
    clinicalFeatures,
    c(
      grep("^genetic_|^collection", clinicalFeatures, value = TRUE),
      "smokingCat"
    )
  )
}

# Cox regression
allResult <- list()
for (chFeature in chFeatures) {
  result <- pbmclapply(olinkFeatures,
    mc.cores = 8, ignore.interactive = TRUE,
    function(col) {
      set.seed(1)
      df <- data %>% dplyr::select(eid, all_of(col))
      df <- clinical %>%
        dplyr::select(eid, all_of(chFeature), all_of(clinicalFeatures)) %>%
        left_join(df,
          by = "eid"
        ) %>%
        dplyr::select(-eid)
      if (!grepl("_vaf|_maxVaf|_frac", chFeature)) {
        df <- df %>% mutate(!!chFeature := ifelse(!!sym(chFeature) == "Yes", 1, 0))
      }

      if (all(df[[chFeature]] %in% c(0, 1))) {
        family <- binomial(link = "logit")
      } else {
        df[[chFeature]] <- ifelse(df[[chFeature]] < 0, 0, df[[chFeature]])
        df[[chFeature]] <- df[[chFeature]] * 100
        family <- gaussian(link = "identity")
      }

      form <- as.formula(str_glue("{chFeature} ~ ."))
      model <- glm(form, data = df, family = family)
      tmp <- summary(model)

      df <- df %>% drop_na()

      list(
        pVal = tmp$coefficients[nrow(tmp$coefficients), ncol(tmp$coefficients)],
        coef = tmp$coefficients,
        numCh = sum(df[[chFeature]] != 0),
        numNonCh = sum(df[[chFeature]] == 0)
      )
    }
  )
  allResult <- result
  names(allResult) <- olinkFeatures
  saveRDS(allResult, str_glue("Result/proChAssoc/proChAssoc_{chFeature}.rds"))
}
