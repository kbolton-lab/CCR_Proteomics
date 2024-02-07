# Content: SNP score and CH association analysis

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("snpScoreData.RData")


clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
indelFeatures <- grep("indel_", colnames(clinical), value = TRUE)
cnvFeatures <- grep("cnv_", colnames(clinical), value = TRUE)
chFeatures <- c(indelFeatures, cnvFeatures)
proteomicFeatures <- colnames(data)[2:ncol(data)]

# Cox regression
for (chFeature in chFeatures) {
  allResult <- list()
  result <- pbmclapply(proteomicFeatures,
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
  names(allResult) <- proteomicFeatures
  saveRDS(allResult, str_glue("Result/snpScoreChAssoc/snpScoreChAssoc_{chFeature}.rds"))
}
