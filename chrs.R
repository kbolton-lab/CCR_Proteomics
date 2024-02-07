# Content: Compare the performance of Cox regression with different feature sets (chrs vs chrs + proteomics) in different data subsets (all, nonChAndLowRisk, medHighRisk)

library(tidyverse)
library(data.table)
library(glmnet)
library(parallel)
library(pbmcapply)
library(survival)
RhpcBLASctl::blas_set_num_threads(2)

load("data.RData")
allEvents <- c("Myeloid")

# Calculate chrs
tmp <- clinical
tmp <- tmp %>%
    mutate(
        chrs_singleDNMT3A = case_when(
            indel_chGene == "" ~ 0,
            indel_chGene == "DNMT3A" ~ 0.5,
            TRUE ~ 1
        ),
        chrs_highRiskMut = case_when(
            indel_chGene == "" ~ 0,
            Reduce(\(x, y) x | y, lapply(c("SF3B1", "SRSF2", "ZRSR2", "JAK2", "TP53", "RUNX1", "FLT3", "IDH1", "IDH2"), \(x) grepl(x, indel_chGene))) ~ 2.5,
            TRUE ~ 1
        ),
        chrs_mutNumber = case_when(
            indel_chGene == "" ~ 0,
            sapply(indel_chGene, \(x) length(str_split(x, ",")[[1]])) >= 2 ~ 2,
            TRUE ~ 1
        ),
        chrs_vaf = case_when(
            indel_chGene == "" ~ 0,
            indel_maxVaf >= 0.2 ~ 2,
            TRUE ~ 1
        ),
        chrs_rdw = case_when(
            # indel_chGene == "" ~ 0,
            red_blood_cell_erythrocyte_distribution_width >= 15 ~ 2.5,
            TRUE ~ 1
        ),
        chrs_mcv = case_when(
            # indel_chGene == "" ~ 0,
            mean_corpuscular_volume > 100 ~ 2.5,
            TRUE ~ 1
        ),
        chrs_age = case_when(
            # indel_chGene == "" ~ 0,
            ageBaseline >= 65 ~ 1.5,
            TRUE ~ 1
        ),
        ccus = case_when(
            indel_chPD == 0 ~ "",
            (sex == "Female" & haemoglobin_concentration < 12) | (sex == "Male" & haemoglobin_concentration < 13) | (platelet_count < 150) | (neutrophill_count < 1.8) ~ "CCUS",
            TRUE ~ "CHIP"
        ),
        chrs_cytopenia = case_when(
            ccus == "CCUS" ~ 1.5,
            ccus == "CHIP" ~ 1,
            TRUE ~ 0
        )
    )
tmp$chrs <- rowSums(tmp[, grep("^chrs_", colnames(tmp))], na.rm = TRUE)

clinical <- clinical %>%
    left_join(tmp %>% dplyr::select(eid, chrs), by = "eid")
clinical$one <- 1

chrsFeature <- c("one", "chrs")
proFeatures <- grep("eid", colnames(data), value = TRUE, invert = TRUE)

data <- data %>%
    mutate(across(all_of(colnames(data)), \(x) replace_na(x, min(x, na.rm = T))))

allData <- data %>%
    left_join(clinical, by = "eid")

featureSetList <- list(
    chrs = chrsFeature,
    chrsPro = c(chrsFeature, proFeatures)
)
featureSetList <- featureSetList[sort(names(featureSetList))]

combi <- expand.grid(
    featureSet = names(featureSetList),
    dataSubset = c("all", "nonChAndLowRisk", "medHighRisk")
)

for (batch in seq_len(nrow(combi))) {
    featureSet <- combi$featureSet[batch] %>% as.character()
    dataSubset <- combi$dataSubset[batch] %>% as.character()

    data <- switch(dataSubset,
        all = allData,
        nonChAndLowRisk = allData %>% dplyr::filter(indel_chGene == "" | chrs <= 9.5),
        medHighRisk = allData %>% dplyr::filter(indel_chGene != "" & chrs > 9.5)
    )

    # Cox regression
    allResult <- list()
    for (event in allEvents) {
        print(str_glue("Processing {event}: {featureSet} {dataSubset}"))
        features <- featureSetList[[featureSet]]
        # check if proteomics feature available in the set, remove non sig proteomics
        if (any(features %in% proFeatures)) {
            features <- c(setdiff(features, proFeatures), allSigProteomics$Myeloid)
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

                smoothAuc <- function(df, nBin = 500) {
                    if (nrow(df) < nBin) {
                        df
                    } else {
                        df %>%
                            as.data.frame() %>%
                            mutate(sub = cut(seq_len(n()), nBin, labels = FALSE)) %>%
                            group_by(sub) %>%
                            group_split(.keep = FALSE) %>%
                            lapply(colMeans) %>%
                            bind_rows() %>%
                            as.data.frame()
                    }
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

    saveRDS(allResult, str_glue("Result/chrs/coxRegChrs_{featureSet}_{dataSubset}.rds"))
}
