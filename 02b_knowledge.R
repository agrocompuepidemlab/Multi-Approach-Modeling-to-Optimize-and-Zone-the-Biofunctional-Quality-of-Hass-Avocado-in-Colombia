## 02b_knowledge.R — knowledge-based selection: variables of greatest weight in López-Hoyos et al. (2024)
## (mean air temperature ~ altitude, subsoil pH, harvest season) + fruit dry matter; MLR, GAM and RF; LOFO validation; y-randomization
suppressPackageStartupMessages({library(dplyr);library(mgcv);library(ranger);library(parallel)})
RNGkind("L'Ecuyer-CMRG"); set.seed(2025)
d <- readRDS("clean.rds"); D <- readRDS("design.rds"); X <- D$X; resp <- D$resp
O <- readRDS("models.rds"); res <- O$res; combos <- O$combos
KB <- c("Tmean","pH_2","DM","Harvest_off")
gt <- function(v, df) if(length(unique(df[[v]])) >= 6) paste0("s(", v, ",k=4,bs='ts')") else v
fitK <- function(model, df) switch(model,
  MLR = lm(as.formula(paste("y ~", paste(KB, collapse="+"))), df),
  GAM = gam(as.formula(paste("y ~", paste(sapply(KB, gt, df=df), collapse="+"))), data=df, method="REML"),
  RF  = ranger(as.formula(paste("y ~", paste(KB, collapse="+"))), data=df, num.trees=1000, min.node.size=3, importance="permutation", seed=1))
pr <- function(m, nd) if(inherits(m,"ranger")) predict(m, nd)$predictions else as.numeric(predict(m, nd))
newc <- data.frame(model=c("MLR","GAM","RF"), method="Knowledge-based", id=c("lm5","gam5","rf5"))
farms <- unique(d$Farm)
lofo <- function(y, model){ df <- data.frame(y=y, X); p <- numeric(nrow(d))
  for(f in farms){ te <- which(d$Farm==f); p[te] <- pr(fitK(model, df[-te,]), df[te,]) }; p }
for(r in resp){ y <- d[[r]]; df <- data.frame(y=y, X); res[[r]]$S$`Knowledge-based` <- KB
  for(k in 1:3){ m <- fitK(newc$model[k], df); res[[r]]$full[[length(res[[r]]$full)+1]] <- m
    res[[r]]$trainP <- cbind(res[[r]]$trainP, pr(m, df)); res[[r]]$testP <- cbind(res[[r]]$testP, lofo(y, newc$model[k])) } }
combos <- rbind(combos, newc)
M <- bind_rows(lapply(res, function(o){ y <- d[[o$resp]]; sst <- sum((y-mean(y))^2)
  bind_rows(lapply(seq_len(nrow(combos)), function(k) data.frame(resp=o$resp, combos[k,],
    RMSE_train=sqrt(mean((y-o$trainP[,k])^2)), RMSE_test=sqrt(mean((y-o$testP[,k])^2)),
    MAE_test=mean(abs(y-o$testP[,k])), Q2_test=1-sum((y-o$testP[,k])^2)/sst, R2_train=1-sum((y-o$trainP[,k])^2)/sst))) }))
## best model per trait = lowest LOFO RMSE; y-randomization (99 permutations) of its LOFO Q2
best <- M %>% group_by(resp) %>% slice_min(RMSE_test, n=1, with_ties=FALSE) %>% ungroup()
yrand <- function(r, k){ y <- d[[r]]; sst <- sum((y-mean(y))^2); mod <- combos$model[k]; meth <- combos$method[k]
  if(meth != "Knowledge-based") return(NA)
  nul <- unlist(mclapply(1:99, function(i){ yy <- sample(y); p <- lofo(yy, mod); 1-sum((yy-p)^2)/sum((yy-mean(yy))^2) }, mc.cores=2))
  q <- 1-sum((y-res[[r]]$testP[,k])^2)/sst; (1+sum(nul>=q))/100 }
best$P_yrand <- mapply(function(r,id) yrand(r, which(combos$id==id)), best$resp, best$id)
## validity: (i) LOFO Q2 > 0 and y-randomization P < 0.05; (ii) for MLR/GAM, Shapiro-Wilk and Breusch-Pagan-type checks on full-data residuals
valid <- function(r, id){ k <- which(combos$id==id); m <- res[[r]]$full[[k]]
  if(inherits(m,"ranger")) return(c(sw=NA, bp=NA)); e <- residuals(m); f <- fitted(m)
  c(sw=shapiro.test(e)$p.value, bp=summary(lm(e^2 ~ f))$coefficients[2,4]) }
V <- t(mapply(valid, best$resp, best$id)); best$p_normality <- V[,1]; best$p_homosced <- V[,2]
best$Valid <- with(best, Q2_test > 0 & (is.na(P_yrand) | P_yrand < 0.05) & (is.na(p_normality) | p_normality > 0.05) & (is.na(p_homosced) | p_homosced > 0.05))
best$Variables <- mapply(function(r, m) paste(res[[r]]$S[[m]], collapse=", "), best$resp, best$method)
saveRDS(list(res=res, combos=combos, M=M, best=best, KB=KB), "models2.rds")
write.csv(M, "out/data/model_performance.csv", row.names=FALSE); write.csv(best, "out/data/best_models.csv", row.names=FALSE)
print(as.data.frame(best %>% mutate(across(where(is.numeric), ~round(.,3))) %>% select(resp,id,RMSE_train,RMSE_test,Q2_test,P_yrand,p_normality,p_homosced,Valid)))
print(as.data.frame(M %>% group_by(method, model) %>% summarise(medQ2=round(median(Q2_test),2), .groups="drop")))
