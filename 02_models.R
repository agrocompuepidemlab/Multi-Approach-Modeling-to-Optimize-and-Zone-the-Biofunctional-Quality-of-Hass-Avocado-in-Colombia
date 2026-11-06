## 02_models.R — 9 model x selection-method combinations (as in the original design), validated by leave-one-farm-out
suppressPackageStartupMessages({library(dplyr);library(mgcv);library(ranger);library(caret);library(randomForest);library(parallel)})
RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
d <- readRDS("clean.rds")
resp <- c("VitE","EE","stearic","palmitic","oleic","linoleic","sfa","mufa","pufa","sfa_ufa")
soil <- grep("_[12]$", names(d), value=TRUE)
clim <- c("Tmean","Tmax","Tmin","Tsd","RHmean","RHsd","PARmean","PARsd","SRmean","SRsd","Rain")
X <- d[,c(clim, soil)]
X$Altitude <- d$Altitude; X$DM <- d$DM
X$Harvest_off <- as.numeric(d$Harvest=="Off-season"); X$Sub_East <- as.numeric(d$Subregion=="East"); X$Sub_SW <- as.numeric(d$Subregion=="Southwest")
X <- X[, apply(X,2,sd) > 0]
preds <- names(X); saveRDS(list(X=X, preds=preds, resp=resp), "design.rds")
## --- selection methods ---
sel_step <- function(x, y){ df <- data.frame(y=y, x); f0 <- lm(y~1, df)
  f <- step(f0, scope=as.formula(paste("~", paste(names(x), collapse="+"))), direction="both", k=log(nrow(df)), steps=6, trace=0)
  v <- attr(terms(f),"term.labels"); if(length(v)==0) v <- names(x)[which.max(abs(cor(x,y)))]; v }
sbf_funcs <- rfSBF; sbf_funcs$score <- function(x, y) { p <- summary(lm(y ~ x))$coefficients; if(nrow(p) < 2) 1 else p[2,4] }
sbf_funcs$filter <- function(score, x, y) { k <- score < 0.05; if(sum(k)==0) k[which.min(score)] <- TRUE; k }
sel_sbf <- function(x, y){ s <- sbf(x, y, sbfControl=sbfControl(functions=sbf_funcs, method="repeatedcv", number=5, repeats=5, verbose=FALSE))
  v <- s$optVariables; sc <- sapply(v, function(z) sbf_funcs$score(x[[z]], y)); v[order(sc)] }      # ordered by p-value
sel_rfe <- function(x, y){ r <- rfe(x, y, sizes=1:8, rfeControl=rfeControl(functions=rfFuncs, method="repeatedcv", number=5, repeats=5, verbose=FALSE)); r$optVariables }
## --- model fits (cap terms to keep models estimable with n = 28-32) ---
decor <- function(v, df, thr=0.8){ keep <- c(); for(z in v){ if(!length(keep) || all(abs(cor(df[[z]], df[,keep,drop=FALSE])) < thr)) keep <- c(keep, z) }; keep }
fit_lm <- function(df, v){ v <- head(decor(v, df), 6); m <- lm(as.formula(paste("y ~", paste(v, collapse="+"))), df)
  al <- names(which(is.na(coef(m)))); if(length(al)) { v <- setdiff(v, al); m <- lm(as.formula(paste("y ~", paste(v, collapse="+"))), df) }; m }
gterm <- function(v, df) if(length(unique(df[[v]])) >= 6) paste0("s(", v, ",k=4,bs='cr')") else v
fit_gam <- function(df, v){ v <- head(decor(v, df), 4); gam(as.formula(paste("y ~", paste(sapply(v, gterm, df=df), collapse="+"))), data=df, method="REML") }
fit_rf <- function(df, v){ ranger(as.formula(paste("y ~", paste(v, collapse="+"))), data=df, num.trees=1000, min.node.size=3, importance="permutation", seed=1) }
combos <- data.frame(model=c("MLR","MLR","MLR","GAM","GAM","RF","RF","RF","RF"), method=c("Stepwise","SBF","RFE","SBF","RFE","Saturated","Stepwise","SBF","RFE"))
combos$id <- c("lm2","lm3","lm4","gam3","gam4","rf1","rf2","rf3","rf4")
select_all <- function(x, y) list(Saturated=names(x), Stepwise=sel_step(x,y), SBF=sel_sbf(x,y), RFE=sel_rfe(x,y))
fit_one <- function(model, df, v) switch(model, MLR=fit_lm(df,v), GAM=fit_gam(df,v), RF=fit_rf(df,v))
pred_one <- function(m, nd) if(inherits(m,"ranger")) predict(m, nd)$predictions else as.numeric(predict(m, nd))
farms <- unique(d$Farm)
run_resp <- function(r){
  y <- d[[r]]; df <- data.frame(y=y, X)
  # full-data selection and fits (training RMSE and final models)
  S <- select_all(X, y); full <- list(); trainP <- matrix(NA, nrow(d), nrow(combos))
  for(k in seq_len(nrow(combos))){ m <- tryCatch(fit_one(combos$model[k], df, S[[combos$method[k]]]), error=function(e) NULL); full[[k]] <- m
    if(!is.null(m)) trainP[,k] <- pred_one(m, df) }
  # leave-one-farm-out with nested selection
  testP <- matrix(NA, nrow(d), nrow(combos))
  for(f in farms){ te <- which(d$Farm==f); tr <- setdiff(seq_len(nrow(d)), te)
    xs <- X[tr,]; xs <- xs[, apply(xs,2,sd)>0]; Sf <- select_all(xs, y[tr])
    for(k in seq_len(nrow(combos))){ m <- tryCatch(fit_one(combos$model[k], df[tr,], Sf[[combos$method[k]]]), error=function(e) NULL)
      if(!is.null(m)) testP[te,k] <- tryCatch(pred_one(m, df[te,]), error=function(e) NA) } }
  list(resp=r, S=S, full=full, trainP=trainP, testP=testP) }
res <- mclapply(resp, run_resp, mc.cores=2, mc.preschedule=FALSE)
names(res) <- resp
saveRDS(list(res=res, combos=combos), "models.rds")
## metrics
M <- bind_rows(lapply(res, function(o){ y <- d[[o$resp]]; sst <- sum((y-mean(y))^2)
  bind_rows(lapply(seq_len(nrow(combos)), function(k) data.frame(resp=o$resp, combos[k,],
    RMSE_train=sqrt(mean((y-o$trainP[,k])^2)), RMSE_test=sqrt(mean((y-o$testP[,k])^2)),
    Q2_test=1-sum((y-o$testP[,k])^2)/sst, R2_train=1-sum((y-o$trainP[,k])^2)/sst))) }))
write.csv(M, "out/data/model_performance.csv", row.names=FALSE)
print(M %>% mutate(across(where(is.numeric), ~round(.,3))) %>% select(resp,id,RMSE_train,RMSE_test,Q2_test), n=100)
