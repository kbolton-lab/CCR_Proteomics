# Content: Association between protein biomarkers and blood count

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
library(argparser)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

# Only use the clinical
clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
olinkFeatures <- colnames(data)[2:ncol(data)]

chFeatures <- grep("indel_", colnames(clinical), value = TRUE)

bloodCountFeatures <- grep("_count", colnames(clinical), value = TRUE)
bloodCountFeatures <- c(bloodCountFeatures, "mean_corpuscular_volume", "haemoglobin_concentration", "red_blood_cell_erythrocyte_distribution_width")
clinicalFeatures <- setdiff(clinicalFeatures, bloodCountFeatures)

allResult <- list()
for (bc in bloodCountFeatures) {
  result <- pbmclapply(olinkFeatures,
    mc.cores = 8, ignore.interactive = TRUE,
    function(col) {
      set.seed(1)
      df <- data %>% dplyr::select(eid, all_of(col))
      df <- df %>% left_join(
        clinical %>%
          dplyr::select(eid, all_of(clinicalFeatures), all_of(chFeatures), all_of(bc)),
        by = "eid"
      )

      df <- df %>% drop_na()

      family <- gaussian()

      form <- as.formula(str_glue("{col} ~ ."))
      model <- try(glm(form, data = df, family = family))
      modelSummary <- summary(model)

      list(
        pVal = modelSummary$coefficients[bc, "Pr(>|t|)"],
        coef = modelSummary$coefficients,
        correlation = cor(df[[col]], df[[bc]], use = "complete.obs")
      )
    }
  )
  names(result) <- olinkFeatures
  allResult[[bc]] <- result
}

saveRDS(allResult, file = str_glue("Result/proBcAssoc.rds"))
