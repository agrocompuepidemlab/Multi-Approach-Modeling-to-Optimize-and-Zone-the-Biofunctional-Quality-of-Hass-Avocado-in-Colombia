## 05_tables.R — tables for the manuscript and supplementary material
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
O <- readRDS("models2.rds"); L <- readRDS("labels.rds"); d <- readRDS("clean.rds"); D <- readRDS("design.rds")
res <- O$res; cb <- O$combos; best <- O$best; M <- O$M
mv <- function(m) if(inherits(m,"ranger")) m$forest$independent.variable.names else all.vars(formula(m))[-1]
nice <- function(v) paste(ifelse(is.na(L$vs[v]), v, L$vs[v]), collapse="; ")
T1 <- best %>% rowwise() %>% mutate(Trait=L$tl[resp], `Best model (selection)`=paste0(model, " (", method, ")"),
  Predictors=nice(mv(res[[resp]]$full[[which(cb$id==id)]])),
  `RMSE training`=signif(RMSE_train,3), `RMSE LOFO`=signif(RMSE_test,3), `Q² LOFO`=sprintf("%.2f",Q2_test),
  `y-rand. P`=ifelse(is.na(P_yrand),"–",sprintf("%.2f",P_yrand)), Valid=ifelse(Valid,"Yes","No")) %>% ungroup() %>%
  arrange(match(resp, names(L$tl))) %>% select(Trait, `Best model (selection)`, Predictors, `RMSE training`, `RMSE LOFO`, `Q² LOFO`, `y-rand. P`, Valid)
dir.create("out/tables", showWarnings=FALSE); write.csv(T1, "out/tables/Table1_best_models.csv", row.names=FALSE)
S2 <- M %>% transmute(Trait=L$tl[resp], Model=id, Algorithm=model, Selection=method, `RMSE training`=signif(RMSE_train,3), `RMSE LOFO`=signif(RMSE_test,3), `MAE LOFO`=signif(MAE_test,3), `Q² LOFO`=sprintf("%.2f",Q2_test)) %>% arrange(match(Trait, L$tl), Model)
write.csv(S2, "out/tables/TableS2_all_models.csv", row.names=FALSE)
LOG <- read.csv("out/data/correction_log.csv"); names(LOG) <- c("Code","Variable","Rows","Correction"); write.csv(LOG, "out/tables/TableS3_corrections.csv", row.names=FALSE)
## Table S1: variables, units, observed range (corrected data)
rng <- function(x) paste0(format(signif(min(x),3), big.mark=",", trim=TRUE), "–", format(signif(max(x),3), big.mark=",", trim=TRUE))
resp <- names(L$tl); preds <- setdiff(names(D$X), c("Harvest_off","Sub_East","Sub_SW"))
S1 <- bind_rows(
  data.frame(Category="Response", Variable=sub(" \\(.*","",L$tu[resp]), Units=ifelse(grepl("\\(", L$tu[resp]), sub(".*\\((.*)\\)$","\\1",L$tu[resp]), "–"), Range=sapply(resp, function(r) rng(d[[r]]))),
  data.frame(Category=ifelse(grepl("_[12]$",preds),"Soil",ifelse(preds %in% c("Altitude","DM"),"Covariate","Climate")),
    Variable=sub(" \\(.*$","",L$vl[preds]), Units=ifelse(grepl("\\(", L$vl[preds]), sub("^[^(]*\\((.*)\\)$","\\1",L$vl[preds]), "–"), Range=sapply(preds, function(p) rng(D$X[[p]]))),
  data.frame(Category="Covariate", Variable=c("Harvest","Subregion","Farm","Tree"), Units="factor", Range=c("Main, off-season","North, East, Southwest","8 farms","2 trees per farm")))
write.csv(S1, "out/tables/TableS1_variables.csv", row.names=FALSE)
A <- read.csv(list.files("out/data", pattern="zoning_area_by_department", full.names=TRUE)[1], check.names=FALSE)
A$Department <- ifelse(is.na(A$Department), "Not assigned (border cells)", A$Department); A$Total <- rowSums(A[,-1])
A <- A %>% arrange(desc(Total)); write.csv(A, "out/tables/TableS4_zoning_area.csv", row.names=FALSE)
print(T1, width=200); print(head(A,8))
