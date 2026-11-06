## 03_figures.R — Figures 2–6 and supplementary figures in the style of the original manuscript
suppressPackageStartupMessages({library(dplyr);library(tidyr);library(ggplot2);library(patchwork);library(mgcv);library(ranger);library(viridis)})
d <- readRDS("clean.rds"); D <- readRDS("design.rds"); X <- D$X; O <- readRDS("models2.rds")
res <- O$res; combos <- O$combos; best <- O$best; M <- O$M
dir.create("out/figs", showWarnings=FALSE)
## ---------- labels ----------
tl <- c(VitE="Vitamin E", EE="Ether extract", stearic="Stearic acid", palmitic="Palmitic acid", oleic="Oleic acid", linoleic="Linoleic acid",
        sfa="Sum of saturated fatty acids", mufa="Sum of monounsaturated fatty acids", pufa="Sum of polyunsaturated fatty acids", sfa_ufa="Saturated/unsaturated ratio")
tu <- c(VitE="Vitamin E (mg kg⁻¹)", EE="Ether extract (g 100 g⁻¹)", stearic="Stearic acid (%)", palmitic="Palmitic acid (%)", oleic="Oleic acid (%)",
        linoleic="Linoleic acid (%)", sfa="SFA (%)", mufa="MUFA (%)", pufa="PUFA (%)", sfa_ufa="SFA/UFA ratio")
cm <- "cmol(+) kg⁻¹"; mk <- "mg kg⁻¹"
sl <- c(pH="pH", Acidity=paste0("Exchangeable acidity (",cm,")"), OM="Organic matter (%)", P=paste0("P (",mk,")"), S=paste0("S (",mk,")"),
        Al=paste0("Exchangeable Al (",cm,")"), Ca=paste0("Exchangeable Ca (",cm,")"), Mg=paste0("Exchangeable Mg (",cm,")"), K=paste0("Exchangeable K (",cm,")"),
        Na=paste0("Exchangeable Na (",cm,")"), ECEC=paste0("ECEC (",cm,")"), EC="Electrical conductivity (dS m⁻¹)", Fe=paste0("Fe (",mk,")"),
        Cu=paste0("Cu (",mk,")"), Mn=paste0("Mn (",mk,")"), Zn=paste0("Zn (",mk,")"), B=paste0("B (",mk,")"))
vl <- c(Tmean="Mean temperature (°C)", Tmax="Maximum temperature (°C)", Tmin="Minimum temperature (°C)", Tsd="SD of temperature (°C)",
        RHmean="Mean relative humidity (%)", RHsd="SD of relative humidity (%)", PARmean="Mean PAR (µmol m⁻² s⁻¹)", PARsd="SD of PAR (µmol m⁻² s⁻¹)",
        SRmean="Mean solar radiation (W m⁻²)", SRsd="SD of solar radiation (W m⁻²)", Rain="Precipitation (mm)", Altitude="Altitude (m a.s.l.)",
        DM="Dry matter (%)", Harvest_off="Off-season harvest", Sub_East="East subregion", Sub_SW="Southwest subregion")
for(p in names(sl)) for(k in 1:2){ dep <- ifelse(k==1,"0–30 cm","31–60 cm"); u <- regmatches(sl[p], regexpr("\\(.*\\)$", sl[p])); b <- trimws(sub("\\(.*\\)$","",sl[p]))
  vl[paste0(p,"_",k)] <- if(length(u)) paste0(b, ", ", dep, " ", u) else paste0(b, ", ", dep) }
vs <- trimws(sub(" \\((?!0).*$", "", vl, perl=TRUE)); vs <- sub("Exchangeable", "Exch.", vs)   # short labels (no units)
saveRDS(list(tl=tl,tu=tu,vl=vl,vs=vs), "labels.rds")
binvars <- c("Harvest_off","Sub_East","Sub_SW")
## ---------- helpers ----------
getm <- function(r){ id <- best$id[best$resp==r]; k <- which(combos$id==id); list(m=res[[r]]$full[[k]], id=id, model=combos$model[k], method=combos$method[k]) }
mvars <- function(m) if(inherits(m,"ranger")) m$forest$independent.variable.names else all.vars(formula(m))[-1]
typical <- function(){ t <- as.data.frame(lapply(X, median)); t$Harvest_off <- 0; t$Sub_East <- 1; t$Sub_SW <- 0; t }
predb <- function(m, nd){ # mean and bands (inner ~68%, outer 95%)
  if(inherits(m,"ranger")){ P <- predict(m, nd, predict.all=TRUE)$predictions
    data.frame(fit=rowMeans(P), lo1=apply(P,1,quantile,.25), hi1=apply(P,1,quantile,.75), lo2=apply(P,1,quantile,.05), hi2=apply(P,1,quantile,.95)) }
  else { p <- predict(m, nd, se.fit=TRUE); f <- as.numeric(p$fit); s <- as.numeric(p$se.fit)
    data.frame(fit=f, lo1=f-s, hi1=f+s, lo2=f-1.96*s, hi2=f+1.96*s) } }
importance <- function(r){ g <- getm(r); m <- g$m; v <- mvars(m)
  if(inherits(m,"ranger")) return(data.frame(var=names(m$variable.importance), imp=as.numeric(m$variable.importance), type="Importance"))
  t0 <- typical(); out <- lapply(v, function(z){ xs <- if(z %in% binvars) c(0,1) else seq(min(X[[z]]),max(X[[z]]),length.out=50)
    nd <- t0[rep(1,length(xs)),]; nd[[z]] <- xs; p <- as.numeric(predict(m, nd)); data.frame(var=z, imp=diff(range(p)), type=ifelse(cor(xs,p) >= 0, "Positive","Negative")) })
  bind_rows(out) }
topvars <- function(r, n=3){ I <- importance(r) %>% filter(!var %in% binvars) %>% arrange(desc(imp)); head(I$var, n) }
th <- theme_minimal(base_size=11) + theme(plot.title=element_text(face="bold", hjust=0.5))
## ---------- Figure 2: exploratory response curves (lm, GAM df=2, GAM df=3) ----------
f2 <- lapply(names(tl), function(r){ v <- topvars(r, 3)
  dd <- bind_rows(lapply(v, function(z) data.frame(var=vs[z], x=X[[z]], y=d[[r]])))
  fits <- bind_rows(lapply(v, function(z){ df <- data.frame(x=X[[z]], y=d[[r]]); xs <- seq(min(df$x),max(df$x),length.out=100)
    k2 <- min(3, length(unique(df$x))); k3 <- min(4, length(unique(df$x)))
    bind_rows(data.frame(var=vs[z], x=xs, y=predict(lm(y~x,df), data.frame(x=xs)), Model="lm"),
      data.frame(var=vs[z], x=xs, y=as.numeric(predict(gam(y~s(x,k=k2,fx=TRUE),data=df), data.frame(x=xs))), Model="gam df=2"),
      data.frame(var=vs[z], x=xs, y=as.numeric(predict(gam(y~s(x,k=k3,fx=TRUE),data=df), data.frame(x=xs))), Model="gam df=3")) }))
  fits$Model <- factor(fits$Model, levels=c("gam df=2","gam df=3","lm")); dd$var <- factor(dd$var, levels=vs[v]); fits$var <- factor(fits$var, levels=vs[v])
  ggplot(dd, aes(x, y)) + geom_point(size=1) + geom_line(data=fits, aes(color=Model), linewidth=0.7) + facet_wrap(~var, scales="free_x", nrow=1) +
    scale_color_manual(values=c("gam df=2"="blue","gam df=3"="green3","lm"="red")) + labs(title=tl[r], x="Value", y=tu[r]) + th + theme(strip.text=element_text(size=9)) })
F2 <- wrap_plots(f2, ncol=2) + plot_layout(guides="collect")
ggsave("out/figs/Fig2_response_curves_exploratory.png", F2, width=16, height=17, dpi=300, bg="white", device=ragg::agg_png)
## ---------- Figure 3: variable importance of the best model ----------
f3 <- lapply(names(tl), function(r){ g <- getm(r); I <- importance(r); I$lab <- vs[I$var]; I$lab[is.na(I$lab)] <- I$var[is.na(I$lab)]
  I <- I %>% arrange(imp) %>% mutate(lab=factor(lab, levels=unique(lab)))
  p <- ggplot(I, aes(imp, lab)) + labs(title=paste0(tl[r]," (", g$model, ")"), x=ifelse(g$model=="RF","Permutation importance","Effect size (response units)"), y="Variables") + th
  if(g$model=="RF") p + geom_col(fill="grey35") else p + geom_col(aes(fill=type)) + scale_fill_manual(values=c(Negative="#F8766D", Positive="#00BFC4"), name="Type", drop=FALSE) })
F3 <- wrap_plots(f3, ncol=2)
ggsave("out/figs/Fig3_variable_importance.png", F3, width=13, height=16, dpi=300, bg="white", device=ragg::agg_png)
## ---------- response curve and surface (style of original Figs 4–7) ----------
tc <- theme_classic(base_size=12) + theme(axis.line=element_line(color="grey70", linewidth=1.2), panel.grid.major=element_line(color="grey92"))
curve_plot <- function(r, z){ g <- getm(r); t0 <- typical(); xs <- seq(min(X[[z]]),max(X[[z]]),length.out=120); nd <- t0[rep(1,120),]; nd[[z]] <- xs
  b <- predb(g$m, nd); b$x <- xs; if(r %in% c('VitE','EE','stearic','linoleic','pufa')) b[,1:5] <- pmax(as.matrix(b[,1:5]), 0)
  ggplot(b, aes(x)) + geom_ribbon(aes(ymin=lo2, ymax=hi2), fill="#CCFF99") + geom_ribbon(aes(ymin=lo1, ymax=hi1), fill="#66CC44") +
    geom_line(aes(y=fit), color="red3", linewidth=0.7) + labs(x=vl[z], y=tu[r]) + tc }
inhull <- function(px, py, hx, hy){ h <- chull(hx, hy); P <- cbind(hx[h], hy[h]); n <- nrow(P); ins <- rep(FALSE, length(px)); j <- n
  for(i in 1:n){ ins <- xor(ins, ((P[i,2]>py)!=(P[j,2]>py)) & (px < (P[j,1]-P[i,1])*(py-P[i,2])/(P[j,2]-P[i,2])+P[i,1])); j <- i }; ins }
surf_plot <- function(r, z1, z2){ if(abs(cor(X[[z1]], X[[z2]])) > 0.8){ alt <- setdiff(topvars(r, 6), z1); alt <- alt[sapply(alt, function(a) abs(cor(X[[z1]], X[[a]])) < 0.8)]; if(length(alt)) z2 <- alt[1] }; g <- getm(r); t0 <- typical(); gx <- seq(min(X[[z1]]),max(X[[z1]]),length.out=80); gy <- seq(min(X[[z2]]),max(X[[z2]]),length.out=80)
  G <- expand.grid(a=gx, b=gy); nd <- t0[rep(1,nrow(G)),]; nd[[z1]] <- G$a; nd[[z2]] <- G$b
  G$z <- predb(g$m, nd)$fit; if(r %in% c('VitE','EE','stearic','linoleic','pufa')) G$z <- pmax(G$z, 0); ok <- inhull(G$a, G$b, jitter(X[[z1]],1e-6), jitter(X[[z2]],1e-6)); G$z[!ok] <- NA
  pts <- unique(data.frame(a=X[[z1]], b=X[[z2]]))
  ggplot(G, aes(a, b, fill=z)) + geom_raster() + scale_fill_viridis_c(name=tu[r], na.value="white") +
    geom_point(data=pts, aes(a, b), inherit.aes=FALSE, shape=21, fill="white", size=1.4) + labs(x=vl[z1], y=vl[z2]) + tc + theme(legend.key.height=unit(0.9,"cm")) }
make_fig <- function(rs, file, h){
  rows <- lapply(seq_along(rs), function(i){ r <- rs[i]; v <- topvars(r, 3); a <- letters[2*i-1]; bb <- letters[2*i]
    cp <- lapply(seq_along(v), function(j){ p <- curve_plot(r, v[j]); if(j==1) p <- p + labs(tag=paste0("(",a,")")); p })
    sp <- surf_plot(r, v[1], v[2]) + labs(tag=paste0("(",bb,")"))
    wrap_plots(list(wrap_plots(cp, nrow=1), sp), nrow=1, widths=c(3,1.25)) })
  F <- wrap_plots(rows, ncol=1) & theme(plot.tag=element_text(size=14, face="plain"))
  ggsave(file, F, width=20, height=h, dpi=300, bg="white", device=ragg::agg_png); setNames(lapply(rs, topvars, n=3), rs) }
v4 <- make_fig(c("VitE","EE"), "out/figs/Fig4_vitE_EE.png", 9)
v5 <- make_fig(c("stearic","palmitic","oleic","linoleic"), "out/figs/Fig5_fatty_acids.png", 17)
v6 <- make_fig(c("sfa","mufa","pufa","sfa_ufa"), "out/figs/Fig6_fatty_acid_sums.png", 17)
saveRDS(c(v4,v5,v6), "fig_vars.rds"); print(c(v4,v5,v6))
## ---------- Fig S1: RMSE train/test (LOFO) for every model and method ----------
M$lab <- factor(M$id, levels=c("gam3","gam4","gam5","lm2","lm3","lm4","lm5","rf1","rf2","rf3","rf4","rf5"))
S1 <- M %>% select(resp, lab, Train=RMSE_train, Test=RMSE_test) %>% pivot_longer(c(Train,Test), names_to="Data", values_to="RMSE")
fs1 <- lapply(names(tl), function(r) ggplot(S1 %>% filter(resp==r), aes(lab, RMSE, fill=Data)) + geom_col(position="dodge") +
  scale_fill_manual(values=c(Test="#F8766D", Train="#00BFC4")) + labs(title=tl[r], x="Model", y="RMSE") + th + theme(axis.text.x=element_text(angle=45, hjust=1)))
ggsave("out/figs/FigS1_RMSE_models.png", wrap_plots(fs1, ncol=2) + plot_layout(guides="collect"), width=14, height=17, dpi=300, bg="white", device=ragg::agg_png)
## ---------- Fig S2: observed vs predicted (leave-one-farm-out) for the best models ----------
OV <- bind_rows(lapply(names(tl), function(r){ k <- which(combos$id==best$id[best$resp==r]); data.frame(trait=tl[r], obs=d[[r]], pred=res[[r]]$testP[,k], Harvest=d$Harvest,
  lab=sprintf("%s: Q² = %.2f", best$id[best$resp==r], best$Q2_test[best$resp==r])) }))
OV$trait <- factor(OV$trait, levels=tl); L <- OV %>% group_by(trait) %>% summarise(lab=first(lab), x=min(obs), y=max(pred))
fs2 <- ggplot(OV, aes(obs, pred)) + geom_abline(linetype="dashed", color="grey40") + geom_point(aes(color=Harvest), size=1.8) + facet_wrap(~trait, scales="free", ncol=4) +
  geom_text(data=L, aes(x, y, label=lab), hjust=0, vjust=1, size=3) + scale_color_manual(values=c(Main="#4C72B0", `Off-season`="#DD5555")) +
  labs(x="Observed", y="Predicted (leave-one-farm-out)") + theme_bw(base_size=10) + theme(legend.position="bottom")
ggsave("out/figs/FigS2_observed_predicted.png", fs2, width=12, height=9, dpi=300, bg="white", device=ragg::agg_png)
## ---------- Fig S3: LOFO Q2 by model, method and trait ----------
M$Trait <- factor(tl[M$resp], levels=tl)
fs3 <- ggplot(M, aes(Trait, lab, fill=Q2_test)) + geom_tile(color="white") + geom_text(aes(label=sprintf("%.2f",Q2_test)), size=2.8) +
  scale_fill_gradient2(low="#B2182B", mid="white", high="#2166AC", midpoint=0, limits=c(-1,0.6), oob=scales::squish, name=expression(Q^2~"(LOFO)")) +
  labs(x=NULL, y="Model") + theme_minimal(base_size=10) + theme(axis.text.x=element_text(angle=35, hjust=1))
ggsave("out/figs/FigS3_Q2_heatmap.png", fs3, width=10, height=6, dpi=300, bg="white", device=ragg::agg_png)
cat("figures done\n")
