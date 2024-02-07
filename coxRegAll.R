# Content: Cox regression using all features, cross validation with 100 folds

library(tidyverse)
library(data.table)
library(glmnet)
library(parallel)
library(pbmcapply)
library(survival)
library(argparser)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")

allEvents <- c("Myeloid")
clinicalFeatures <- grep("indel_|cnv_|eid", colnames(clinical), value = TRUE, invert = TRUE)
indelFeatures <- grep("indel_", colnames(clinical), value = TRUE)
cnvFeatures <- grep("cnv_", colnames(clinical), value = TRUE)
chFeatures <- c(indelFeatures, cnvFeatures)
proFeatures <- grep("eid", colnames(data), value = TRUE, invert = TRUE)

# impute NA to min
data <- data %>%
  mutate(across(all_of(colnames(data)), \(x) replace_na(x, min(x, na.rm = T))))

data <- data %>% left_join(clinical, by = "eid")

featureSetList <- list(
  clinical = clinicalFeatures,
  cliCh = c(clinicalFeatures, chFeatures),
  cliProCh = c(clinicalFeatures, proFeatures, chFeatures)
)

featureSetList <- featureSetList[sort(names(featureSetList))]

# Cox regression
for (featureSet in names(featureSetList)) {
  allResult <- list()
  for (event in allEvents) {
    print(featureSet)
    features <- featureSetList[[featureSet]]
    # check if proteomics feature available in the set, remove non sig proteomics
    if (any(features %in% proFeatures)) {
      features <- c(setdiff(features, proFeatures), intersect(features, allSigProteomics$Myeloid))
    }

    if (any(features %in% chFeatures)) {
      tmp <- intersect(features, allSigCh$Myeloid)
      features <- c(setdiff(features, chFeatures), tmp)
    }

    df <- data
    df <- df %>% dplyr::select(eid, all_of(features))

    tmp <- c(event, paste0("survtime_", event))
    df <- df %>%
      left_join(diagnosis %>% dplyr::select(eid, all_of(tmp)), by = "eid") %>%
      dplyr::select(-eid)

    df[[event]] <- replace_na(df[[event]], 0)

    df <- df %>% drop_na()

    K <- 5
    folds <- NULL
    # 20 different seeds, 5 folds each => 100 folds
    for (seed in 1:20) {
      set.seed(seed)
      folds <- c(folds, split(sample(seq_len(nrow(df))), cut(seq_len(nrow(df)), K, labels = FALSE)))
    }

    i <- 1
    allResult[[event]] <- pbmclapply(seq_along(folds),
      mc.cores = 8, ignore.interactive = TRUE,
      mc.preschedule = FALSE,
      \(i) {
        testIdx <- folds[[i]]
        tmp <- c(event, paste0("survtime_", event))

        survObjTest <- df[testIdx, ] %>% dplyr::select(all_of(tmp))
        colnames(survObjTest) <- c("status", "time")

        survObjTrain <- df[-testIdx, ] %>% dplyr::select(all_of(tmp))
        colnames(survObjTrain) <- c("status", "time")
        survObjTrain <- Surv(survObjTrain$time, survObjTrain$status)

        dfTest <- df[testIdx, ] %>%
          dplyr::select(-all_of(tmp)) %>%
          data.matrix() %>%
          apply(., 2, \(x) {
            if (all(x == x[1])) {
              return(x)
            } else {
              return((x - min(x)) / (max(x) - min(x)))
            }
          })

        dfTrain <- df[-testIdx, ] %>%
          dplyr::select(-all_of(tmp)) %>%
          data.matrix() %>%
          apply(., 2, \(x) {
            if (all(x == x[1])) {
              return(x)
            } else {
              return((x - min(x)) / (max(x) - min(x)))
            }
          })

        doParallel::registerDoParallel(cores = 5)
        # alpha = 0 : ridge, 0.5 : elastic net, 1 : lasso
        cvfit <- cv.glmnet(dfTrain, survObjTrain,
          thresh = 1e-5,
          maxit = 1e4,
          standardize = FALSE,
          family = "cox", alpha = 1, nfolds = 5, parallel = TRUE
        )
        doParallel::stopImplicitCluster()

        print("Predict CV test")
        pred <- predict(cvfit, newx = dfTest, s = cvfit$lambda.min, family = "coxph", type = "response")

        cutoff <- c(5, 10) * 365
        modRF <- timeROC::timeROC(
          T = survObjTest$time,
          delta = survObjTest$status,
          marker = pred,
          cause = 1,
          weighting = "cox",
          times = cutoff
        )
        modRF <- modRF[c("TP", "FP", "AUC")]

        smoothAuc <- function(df) {
          df %>%
            as.data.frame() %>%
            mutate(sub = cut(seq(n()), 500, labels = FALSE)) %>%
            group_by(sub) %>%
            group_split(.keep = FALSE) %>%
            lapply(colMeans) %>%
            bind_rows() %>%
            as.data.frame()
        }
        modRF$TP <- smoothAuc(modRF$TP)
        modRF$FP <- smoothAuc(modRF$FP)

        result <- list(
          coef = coef(cvfit, s = cvfit$lambda.min),
          modRF = modRF,
          numCaseTrain = sum(df[-testIdx, ][[event]] == 1),
          numControlTrain = sum(df[-testIdx, ][[event]] == 0),
          numCaseTest = sum(df[testIdx, ][[event]] == 1),
          numControlTest = sum(df[testIdx, ][[event]] == 0)
        )
        result
      }
    )
    names(allResult[[event]]) <- as.character(seq_along(folds))
  }
  saveRDS(allResult, str_glue("Result/coxRegAll_{featureSet}.rds"))
}
