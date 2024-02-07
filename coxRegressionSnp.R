# Content: Cox regression for SNP score of each Olink feature

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("snpScoreData.RData")

allEvents <- c("Myeloid")
clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
indelFeatures <- grep("indel_", colnames(clinical), value = TRUE)
cnvFeatures <- grep("cnv_", colnames(clinical), value = TRUE)
chFeatures <- c(indelFeatures, cnvFeatures)
proteomicFeatures <- colnames(data)[2:ncol(data)]

featureSetList <- c(list(clinicalFeatures), list(c(clinicalFeatures, "indel_maxVaf", "cnv_allCnv")))
names(featureSetList) <- c("clinical", "cliVafCnv")
featureSetList <- featureSetList[sort(names(featureSetList))]

# Cox regression
for (featureSet in names(featureSetList)) {
  allResult <- list()
  for (event in allEvents) {
    result <- pbmclapply(proteomicFeatures,
      mc.cores = nCores, ignore.interactive = TRUE,
      function(col) {
        set.seed(1)
        df <- data %>% dplyr::select(eid, all_of(col))
        df <- df %>% left_join(clinical %>% dplyr::select(eid, any_of(featureSetList[[featureSet]])),
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
    names(allResult[[event]]) <- proteomicFeatures
  }
  saveRDS(allResult, str_glue("Result/coxRegSnp_{featureSet}.rds"))
}
