# Content: Proteomics and indel VAF association analysis

library(tidyverse)
library(data.table)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
chFeature <- c("indel_maxVaf")
olinkFeatures <- colnames(data)[2:ncol(data)]

clinical <- clinical %>%
  dplyr::select(eid, any_of(c(clinicalFeatures, chFeature)), indel_JAK2) %>%
  dplyr::filter(!is.na(!!sym(chFeature))) %>%
  mutate(!!chFeature := case_when(
    !!sym(chFeature) == 0 ~ "No CH",
    !!sym(chFeature) >= 0.1 ~ "High Vaf",
    !!sym(chFeature) < 0.1 ~ "Low Vaf"
  )) %>%
  mutate(!!chFeature := factor(!!sym(chFeature), levels = c("No CH", "Low Vaf", "High Vaf")))

result <- pbmclapply(olinkFeatures,
  mc.cores = 8, ignore.interactive = TRUE,
  function(col) {
    set.seed(1)
    df <- data %>% dplyr::select(eid, all_of(col))
    df <- clinical %>%
      dplyr::select(eid, all_of(clinicalFeatures), indel_JAK2, all_of(chFeature)) %>%
      left_join(df,
        by = "eid"
      ) %>%
      dplyr::select(-eid)

    family <- gaussian(link = "identity")
    form <- as.formula(str_glue("{col} ~ ."))
    model <- glm(form, data = df, family = family)
    tmp <- summary(model)

    df <- df %>% drop_na()

    list(
      pVal = tmp$coefficients[nrow(tmp$coefficients), ncol(tmp$coefficients)],
      coef = tmp$coefficients,
      numHighVaf = sum(df[[chFeature]] == "High Vaf"),
      numLowVaf = sum(df[[chFeature]] == "Low Vaf"),
      numNoCh = sum(df[[chFeature]] == "No CH")
    )
  }
)
allResult <- result
names(allResult) <- olinkFeatures
saveRDS(allResult, str_glue("Result/proVafAssoc.rds"))

# Repeat without normal (non-CH) samples
load("data.RData")

clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
chFeature <- c("indel_maxVaf")
olinkFeatures <- colnames(data)[2:ncol(data)]

clinical <- clinical %>%
  dplyr::select(eid, all_of(c(clinicalFeatures, chFeature)), indel_JAK2) %>%
  dplyr::filter(!is.na(!!sym(chFeature)) & !!sym(chFeature) != 0) %>%
  mutate(!!chFeature := ifelse(!!sym(chFeature) >= 0.1, "High Vaf", "Low Vaf")) %>%
  mutate(!!chFeature := factor(!!sym(chFeature), levels = c("Low Vaf", "High Vaf")))

result <- pbmclapply(olinkFeatures,
  mc.cores = 8, ignore.interactive = TRUE,
  function(col) {
    set.seed(1)
    df <- data %>% dplyr::select(eid, all_of(col))
    df <- clinical %>%
      dplyr::select(eid, all_of(clinicalFeatures), indel_JAK2, all_of(chFeature)) %>%
      left_join(df,
        by = "eid"
      ) %>%
      dplyr::select(-eid)

    family <- gaussian(link = "identity")
    form <- as.formula(str_glue("{col} ~ ."))
    model <- glm(form, data = df, family = family)
    tmp <- summary(model)

    df <- df %>% drop_na()

    list(
      pVal = tmp$coefficients[nrow(tmp$coefficients), ncol(tmp$coefficients)],
      coef = tmp$coefficients,
      numHighVaf = sum(df[[chFeature]] == "High Vaf"),
      numLowVaf = sum(df[[chFeature]] == "Low Vaf"),
      numNoCh = sum(df[[chFeature]] == "No CH")
    )
  }
)
allResult <- result
names(allResult) <- olinkFeatures
saveRDS(allResult, str_glue("Result/proVafAssocWoNormal.rds"))
