## 01_clean.R — Audit and correction of the original spreadsheet
suppressPackageStartupMessages({library(readxl);library(dplyr);library(tidyr)})
raw <- read_excel("data/Consolidado_broma_fisi_ambas_cosechas_15032019.xlsx", sheet=1, .name_repair="minimal")
names(raw)[60] <- "porc_relativo_ac_palmitoleico"      # duplicated name in original (relative palmitoleic)
log <- list(); L <- function(...) log[[length(log)+1]] <<- data.frame(...)
d <- raw
## C1 altitude in km -> m
d$altitud <- d$altitud*1000
L(id="C1", variable="altitud", rows="all", action="Converted from km to m a.s.l. (values 1.753–2.453 were stored in km)")
## C2 fatty acids: absolute (mg/100 g) and relative (%)
abs_nm <- c("ac_palmitico","ac_palmitoleico","ac_estearico","ac_oleico","ac_linoleico","ac_linolenico","ac_araquidonico")
rel_nm <- c("porc_relativo_ac_palmitico","porc_relativo_ac_palmitoleico","porc_relativo_ac_estearico","proc_relativo_ac_oleico",
            "porc_relativo_ac_linoleico","porc_relativo_ac_linolenico","porc_relativo_ac_araquidonico")
A <- as.matrix(d[,abs_nm]); R <- as.matrix(d[,rel_nm]); Rc <- A/rowSums(A)*100
cp <- d$cosecha=="CP"
dev_cp <- apply(abs(R-Rc),1,max)
## C2a main-harvest rows whose relative profile does not equal abs/sum (rule satisfied exactly by 14 of 16 rows)
bad_cp <- which(cp & dev_cp > 0.01)
for(i in bad_cp){ R[i,] <- Rc[i,]
  L(id="C2a", variable="relative fatty acids", rows=paste(d$finca[i],d$cosecha[i],d$bloque_arbol[i]),
    action=sprintf("Relative profile recomputed as absolute/total x 100 (max deviation %.2f points); all other main-harvest rows satisfy this rule exactly", dev_cp[i])) }
## C2b off-season relative stearic acid = rounded copy of main-harvest values -> recompute from absolute stearic/total
ct <- which(!cp); old <- R[ct,3]; new <- Rc[ct,3]
R[ct,-3] <- R[ct,-3] * (100-new)/(100-old)   # keep the profile closed at 100%
R[ct,3] <- new
L(id="C2b", variable="porc_relativo_ac_estearico", rows="16 off-season rows",
  action="Values were main-harvest values rounded to 2 decimals; recomputed as absolute stearic/total absolute x 100; other acids rescaled so the profile sums 100%")
## C2c label rotation in relative columns: 'linoleico' (3–8%) = unassigned peak; 'linolenico' (7–15%) = linoleic C18:2; 'araquidonico' (0.5–1.2%) = linolenic C18:3
colnames(R) <- c("palmitic","palmitoleic","stearic","oleic","fa_unassigned","linoleic","linolenic")
L(id="C2c", variable="relative linoleic/linolenic/arachidonic", rows="all",
  action="Column labels rotated: values labelled linolenic (6.8–14.8%) are linoleic C18:2; values labelled arachidonic (0.5–1.2%) are linolenic C18:3 (typical Hass ranges; confirmed by the off-season absolute columns, stored in the expected order); values labelled linoleic (3–8%) kept as an unassigned fatty acid")
## sums
FA <- as.data.frame(R)
FA$sfa <- FA$palmitic + FA$stearic
FA$mufa <- FA$palmitoleic + FA$oleic
FA$pufa <- FA$linoleic + FA$linolenic
FA$sfa_ufa <- FA$sfa/(100-FA$sfa)
L(id="C2d", variable="SFA, MUFA, PUFA, SFA/UFA", rows="all", action="Recomputed: SFA = palmitic + stearic; MUFA = palmitoleic + oleic; PUFA = linoleic + linolenic (unassigned peak excluded); SFA/UFA = SFA/(100 - SFA)")
## C3 soil: values follow the sampling zone (row order), tree codes are swapped in the off-season rows of CS, EB and EG
soil_cols <- grep("\\(p[12]a[123]\\)", names(d), value=TRUE)
for(f in unique(d$finca)){
  icp <- which(d$finca==f & d$cosecha=="CP"); ict <- which(d$finca==f & d$cosecha=="CT")
  for(i in ict){ j <- icp[d$bloque_arbol[icp]==d$bloque_arbol[i]]
    if(any(abs(unlist(d[i,soil_cols])-unlist(d[j,soil_cols]))>1e-9, na.rm=TRUE)){
      d[i,soil_cols] <- d[j,soil_cols]
      L(id="C3", variable="soil (all years and depths)", rows=paste(f,"CT",d$bloque_arbol[i]), action="Soil record re-assigned by tree code (it followed row order and belonged to the other tree)") } } }
## C4 EG 2017 soil (a3) identical to CS 2017 for all variables -> set to NA for EG
a3 <- grep("a3\\)", soil_cols, value=TRUE)
d[d$finca=="EG", a3] <- NA
L(id="C4", variable="soil, third sampling (a3)", rows="EG (4 rows)", action="All 34 values identical to farm CS (copy error); set to missing; EG soil = mean of the first two samplings")
## C5 soil predictors = mean of the three annual samplings per depth (López-Hoyos et al., 2024: no year effect)
props <- c("Ph"="pH","acidez"="Acidity","Acidez"="Acidity","M_O"="OM","P"="P","S"="S","Al"="Al","Ca"="Ca","Mg"="Mg","K"="K","Na"="Na","CICE"="ECEC","C_E"="EC","Fe"="Fe","Cu"="Cu","Mn"="Mn","Zn"="Zn","B"="B")
soil <- data.frame(row=1:nrow(d))
for(dep in 1:2) for(p in unique(props)){
  keys <- names(props)[props==p]
  cols <- soil_cols[sapply(soil_cols, function(cn){ base <- sub("_*\\(p[12]a[123]\\)_*$","",cn); base %in% keys }) & grepl(paste0("p",dep,"a"), soil_cols)]
  stopifnot(length(cols)==3)
  soil[[paste0(p,"_",dep)]] <- rowMeans(as.matrix(d[,cols]), na.rm=TRUE) }
soil$row <- NULL
L(id="C5", variable="soil predictors", rows="all", action="Each soil property = mean of the three annual samplings (2015–2017) per depth (1: 0–30 cm; 2: 31–60 cm), as in López-Hoyos et al. (2024)")
## assemble clean dataset (English names)
clim <- d %>% transmute(Tmean=`promedio-tem_amb`, Tmax=max_tem_amb, Tmin=min_tem_amb, Tsd=desvest_tem_amb, RHmean=`promedio-hr_amb`, RHsd=desvest_hr_amb,
  PARmean=promedio_par, PARsd=desvest_par, SRmean=promedio_rad_solar, SRsd=desvest_rad_solar, Rain=suma_pliv)
clean <- data.frame(Harvest=ifelse(d$cosecha=="CP","Main","Off-season"), Subregion=recode(d$region, Norte="North", Oriente="East", Suroeste="Southwest"),
  Farm=d$finca, Tree=d$bloque_arbol, Altitude=d$altitud, DM=d$M_S_fv, Days_FH=d$`dias-flor-cosecha`,
  VitE=d$`VitaE_gr/100gr`, EE=d$`extracto_etereo_g/100gr`, FA[,c("palmitic","palmitoleic","stearic","oleic","linoleic","linolenic","fa_unassigned","sfa","mufa","pufa","sfa_ufa")],
  clim, soil, check.names=FALSE)
dir.create("out/data", showWarnings=FALSE, recursive=TRUE)
write.csv(clean, "out/data/avocado_quality_edaphoclimatic_corrected.csv", row.names=FALSE)
LOG <- bind_rows(log); write.csv(LOG, "out/data/correction_log.csv", row.names=FALSE)
saveRDS(clean, "clean.rds")
print(LOG[,c("id","rows")]); cat("\nprofile sums:", range(rowSums(FA[,1:7])), "\n")
print(summary(clean[,c("stearic","linoleic","linolenic","sfa","mufa","pufa","sfa_ufa","Altitude")]))
