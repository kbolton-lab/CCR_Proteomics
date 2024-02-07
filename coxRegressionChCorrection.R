# Content: Cox regression for each Olink feature with CH correction (CNV and indel)

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

# Only use the clinical
clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
allEvents <- c("Myeloid")
olinkFeatures <- colnames(data)[2:ncol(data)]

# Cox regression
allResult <- list()
for (event in allEvents) {
  col <- olinkFeatures[1]
  result <- pbmclapply(olinkFeatures,
    mc.cores = 8, ignore.interactive = TRUE,
    function(col) {
      set.seed(1)
      df <- data %>% dplyr::select(eid, all_of(col))
      df <- df %>% left_join(
        clinical %>%
          dplyr::select(eid, all_of(clinicalFeatures), cnv_allCnv, indel_maxVaf),
        by = "eid"
      )

      tmp <- c(event, paste0("survtime_", event))
      df <- df %>%
        left_join(diagnosis %>% dplyr::select(eid, all_of(tmp)), by = "eid") %>%
        dplyr::select(-eid)

      df[[event]] <- replace_na(df[[event]], 0)

      df <- df %>% drop_na()

      survObj <- df %>% dplyr::select(all_of(tmp))
      colnames(survObj) <- c("event", "time")
      survObj <- Surv(survObj$time, survObj$event)

      dfClinical <- df %>% dplyr::select(-all_of(tmp))

      form <- as.formula("survObj ~ .")
      model <- coxph(form, data = dfClinical)
      modelSummary <- summary(model)

      cIndex <- concordance(model)$concordance

      list(
        pVal = modelSummary$coefficients[col, "Pr(>|z|)"],
        cIndex = cIndex,
        coef = modelSummary$coefficients,
        numCase = sum(df[[event]] == 1),
        numControl = sum(df[[event]] == 0)
      )
    }
  )
  allResult[[event]] <- result
  names(allResult[[event]]) <- olinkFeatures
}

saveRDS(allResult, str_glue("Result/coxRegChCorrection_clinical.rds"))
