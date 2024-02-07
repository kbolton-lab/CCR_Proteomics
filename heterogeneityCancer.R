# Content: Heterogeneity analysis for each Olink feature

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

diagnosis <- diagnosis %>%
  dplyr::filter(MPN_all == 1 | MDS_AML_all == 1) %>%
  dplyr::filter(!(MPN_all == 1 & MDS_AML_all == 1)) %>%
  mutate(cancer = ifelse(MPN_all == 1, "MPN", "MDS_AML")) %>%
  mutate(cancer = factor(cancer, levels = c("MDS_AML", "MPN")))

clinical <- clinical %>%
  dplyr::filter(eid %in% diagnosis$eid)

data <- data %>%
  dplyr::filter(eid %in% diagnosis$eid)

# Only use the clinical
clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
olinkFeatures <- colnames(data)[2:ncol(data)]

for (chCorrection in c("", "ChCorrection")) {
  result <- pbmclapply(olinkFeatures,
    mc.cores = 8, ignore.interactive = TRUE,
    function(col) {
      set.seed(1)
      df <- data %>% dplyr::select(eid, all_of(col))
      df <- df %>%
        left_join(clinical %>% dplyr::select(eid, all_of(clinicalFeatures)),
          by = "eid"
        ) %>%
        left_join(diagnosis %>% dplyr::select(eid, cancer), by = "eid")

      if (chCorrection == "ChCorrection") {
        df <- df %>%
          left_join(
            clinical %>% dplyr::select(eid, cnv_allCnv, indel_maxVaf),
            by = "eid"
          ) %>%
          relocate(cancer, .after = last_col())
      }

      df <- df %>% dplyr::select(-eid)

      df <- df %>% drop_na()

      form <- as.formula(str_glue("{col} ~ ."))
      model <- glm(form, data = df)
      modelSummary <- summary(model)

      list(
        pVal = modelSummary$coefficients[, 4] %>% tail(1),
        coef = modelSummary$coefficients,
        numMPN = sum(df$cancer == "MPN"),
        numMDS_AML = sum(df$cancer == "MDS_AML")
      )
    }
  )
  names(result) <- olinkFeatures

  saveRDS(result, str_glue("Result/heterogeneityCancer{chCorrection}.rds"))
}
