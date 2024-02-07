# Content: Perform Mendelian randomization analysis for cis-pQTLs

library(data.table)
library(tidyverse)
library(MendelianRandomization)
library(pbmcapply)
data.table::setDTthreads(2)

event <- "Myeloid"

{
  allAssoc <- readRDS("Result/proChAssoc/proChAssoc_indel_or_cnv.rds")
  allAssoc <- data.frame(
    type = "chAssoc",
    assay = names(allAssoc),
    pVal = sapply(allAssoc, function(x) x$pVal),
    pFDR = p.adjust(sapply(allAssoc, function(x) x$pVal), method = "fdr"),
    coef = sapply(allAssoc, \(x) x$coef[nrow(x$coef), 1])
  )

  allSigProteomics <- readRDS("Result/coxReg_clinical.rds")
  allSigProteomics <- allSigProteomics$Myeloid
  olinkFeatures <- names(allSigProteomics)
  allSigProteomics <- lapply(olinkFeatures, function(feature) {
    data.frame(
      assay = feature,
      pVal = allSigProteomics[[feature]][["coef"]][feature, "Pr(>|z|)"],
      coef = allSigProteomics[[feature]][["coef"]][feature, 1]
    )
  }) %>%
    do.call(what = rbind) %>%
    mutate(pFDR = p.adjust(pVal, method = "fdr"))

  allLevels <- c(str_glue("Associated with {tolower(event)} and CH"), str_glue("Associated with {tolower(event)}"), "Associated with CH", "Not significant")
  allSigProteomics <- allSigProteomics %>%
    dplyr::mutate(
      type = "myAssoc"
    ) %>%
    dplyr::select(type, assay, pFDR, coef) %>%
    rbind(allAssoc %>% dplyr::select(type, assay, pFDR, coef)) %>%
    group_by(assay) %>%
    mutate(typeSig = paste(type[pFDR <= 0.05], collapse = "|")) %>%
    ungroup() %>%
    mutate(
      typeSig = case_when(
        typeSig == "myAssoc|chAssoc" ~ allLevels[1],
        typeSig == "chAssoc" ~ allLevels[3],
        typeSig == "myAssoc" ~ allLevels[2],
        TRUE ~ allLevels[4]
      ),
      typeSig = factor(typeSig, levels = allLevels)
    ) %>%
    dplyr::arrange(desc(as.numeric(typeSig))) %>%
    dplyr::filter(type == "myAssoc" & typeSig != "Not significant") %>%
    mutate(typeSig = typeSig %>% droplevels()) %>%
    dplyr::select(-type) %>%
    mutate(assay = toupper(assay))
}


prefix <- c("myeloid")
names(prefix) <- c("Myeloid")
prefix <- prefix[event]

allResult <- fread(str_glue("Result/MR/allSnpResult_{event}.csv"), data.table = FALSE) %>%
  dplyr::filter(str_split(snpId, "_", simplify = TRUE)[, 1] != "chrX")
allRsqResult <- readRDS(str_glue("Result/MR/allRsqResult_{event}.rds"))
allRsqResult <- allRsqResult[allRsqResult > 0.1]

allResult <- allResult %>%
  dplyr::filter(targetGene %in% names(allRsqResult))

allMrResult <- pbmclapply(allResult$targetGene %>% unique(),
  mc.cores = 8,
  \(target) {
    res <- allResult %>%
      dplyr::filter(targetGene == !!target)

    if (nrow(res) < 3) {
      return(NULL)
    }

    mrInput <- mr_input(
      bx = res$olinkCoef,
      bxse = res$olinkSe,
      by = res[[str_glue("{prefix}Coef")]],
      byse = res[[str_glue("{prefix}Se")]],
      exposure = target,
      outcome = event
    )

    mrObject <- mr_median(mrInput,
      weighting = "weighted",
      distribution = "normal",
      alpha = 0.05,
      iterations = 10000,
      seed = 1
    )

    mrObjectAll <- mr_allmethods(mrInput, seed = 1)

    list(
      mrInput = mrInput,
      mrRes = data.frame(
        target = target,
        mrCoef = mrObject@Estimate,
        mrSe = mrObject@StdError,
        mrPvalue = mrObject@Pvalue,
        rSquared = allRsqResult[target],
        numGenotype = nrow(res),
        stringsAsFactors = FALSE
      ),
      fullMrRes = mrObjectAll
    )
  }
)
names(allMrResult) <- allResult$targetGene %>% unique()
allMrResult <- allMrResult[!sapply(allMrResult, is.null)]

saveRDS(allMrResult, str_glue("Result/MR/allMrResult_{event}.rds"))

geneMap <- fread("Data/olinkGeneRegion.csv", data.table = FALSE)

allMrResultTable <- lapply(allMrResult, \(x) x$mrRes) %>% bind_rows()
allMrResultTable <- allMrResultTable %>%
  left_join(geneMap %>% dplyr::select(V1, V5) %>%
    dplyr::rename(target = V1, assay = V5), by = "target") %>%
  left_join(allSigProteomics, by = c("assay" = "assay")) %>%
  arrange(mrPvalue) %>%
  dplyr::filter(typeSig != "Associated with CH") %>%
  mutate(sameSign = sign(coef) == sign(mrCoef)) %>%
  dplyr::rename(
    pFDRProteomics = pFDR,
    coefProteomics = coef,
    typeSigProteomics = typeSig
  )
allMrResultTable <- allMrResultTable %>% arrange(desc(sameSign), mrPvalue)

fwrite(allMrResultTable, str_glue("Result/mrGene_{event}.csv"))
