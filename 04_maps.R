## 04_maps.R — Spatial extrapolation of lipid quality to Colombian avocado areas
## Uses the high-resolution layers in gis/final/ when available (WorldClim 2.1 30", SoilGrids 2.0 pH 30–60 cm, UPRA Hass aptitude);
## otherwise a provisional version with WorldClim 10' (temperature only) and the thermal envelope of the sampled farms.
suppressPackageStartupMessages({library(terra);library(sf);library(dplyr);library(ggplot2);library(patchwork);library(mgcv);library(viridis)})
sf_use_s2(FALSE)
d <- readRDS("clean.rds"); lab <- readRDS("labels.rds")
load("gis/states10.rda"); col <- st_as_sf(states10) %>% filter(admin=="Colombia") %>% select(name)
col$name <- iconv(col$name, "UTF-8", "UTF-8")
FINAL <- file.exists("gis/final/tavg_col.tif") && file.exists("gis/final/phh2o_30_60_col.tif") && length(list.files("gis/final/upra", pattern="\\.(shp|gpkg)$"))>0
if(FINAL){
  tmean <- mean(rast("gis/final/tavg_col.tif")); names(tmean) <- "Tmean"
  ph <- resample(rast("gis/final/phh2o_30_60_col.tif")/10, tmean); names(ph) <- "pH_2"
  upra <- st_read(list.files("gis/final/upra", pattern="\\.(shp|gpkg)$", full.names=TRUE)[1], quiet=TRUE) %>% st_transform(4326)
  zones <- vect(upra); env <- c(tmean, ph); version <- "final"
} else {
  b1 <- rast("gis/wc10/bio1.bil")/10; tmean <- crop(b1, ext(-80,-66,-4.5,13)); names(tmean) <- "Tmean"
  env <- tmean; version <- "provisional"
}
env <- mask(crop(env, vect(col)), vect(col))
## mappable models: harvest fixed (main crop), fruit dry matter at the commercial median; predictors available as free layers
vars <- if(version=="final") c("Tmean","pH_2") else c("Tmean")
traits <- c(oleic="Oleic acid (%)", sfa="SFA (%)", mufa="MUFA (%)", sfa_ufa="SFA/UFA ratio")
dd <- d; dd$Harvest_off <- as.numeric(d$Harvest=="Off-season")
farms <- unique(d$Farm); Q <- c(); fits <- list()
for(r in names(traits)){ f <- as.formula(paste(r, "~ Harvest_off + DM +", paste0("s(", vars, ",k=4,bs='ts')", collapse="+")))
  p <- numeric(nrow(dd)); for(fm in farms){ te <- dd$Farm==fm; p[te] <- predict(gam(f, data=dd[!te,], method="REML"), dd[te,]) }
  Q[r] <- 1 - sum((dd[[r]]-p)^2)/sum((dd[[r]]-mean(dd[[r]]))^2); fits[[r]] <- gam(f, data=dd, method="REML") }
print(round(Q,2))
## area of applicability: predictor values within the range sampled by the farms (+-2.5% of the range)
aoa <- env[[1]] >= min(dd$Tmean) - 0.025*diff(range(dd$Tmean)) & env[[1]] <= max(dd$Tmean) + 0.025*diff(range(dd$Tmean))
if(version=="final"){ aoa <- aoa & env[["pH_2"]] >= min(dd$pH_2)-0.025*diff(range(dd$pH_2)) & env[["pH_2"]] <= max(dd$pH_2)+0.025*diff(range(dd$pH_2))
  aoa <- mask(aoa, rasterize(zones, aoa)) }
aoa <- classify(aoa, cbind(0, NA))
nd <- as.data.frame(env, xy=TRUE, na.rm=TRUE); nd$Harvest_off <- 0; nd$DM <- median(dd$DM)
for(r in names(traits)) nd[[r]] <- as.numeric(predict(fits[[r]], nd))
A <- as.data.frame(aoa, xy=TRUE, na.rm=TRUE); nd <- semi_join(nd, A, by=c("x","y"))
nd$class <- cut(nd$oleic, c(-Inf, 58, 65, Inf), labels=c("Lower (oleic < 58%)","Intermediate (58–65%)","Higher (oleic ≥ 65%)"))
write.csv(nd, sprintf("out/data/map_predictions_%s.csv", version), row.names=FALSE)
## approximate locations of the sampled farms (municipal seats)
sites <- data.frame(farm=c("LA","LE","EG/EC","IM","BV","CS/EB"), mun=c("El Peñol","Rionegro","El Retiro","Amagá","Jardín","San Pedro de los Milagros"),
  x=c(-75.243,-75.374,-75.502,-75.704,-75.819,-75.556), y=c(6.219,6.155,6.058,6.039,5.598,6.461))
ant <- col %>% filter(grepl("Antioquia", name))
base <- function() list(geom_sf(data=col, fill="grey95", color="grey60", linewidth=0.2), theme_minimal(base_size=10), theme(axis.title=element_blank(), panel.grid=element_line(color="grey90")))
cz <- function() list(coord_sf(xlim=c(-79.2,-71.3), ylim=c(0.3,11.6), expand=FALSE), scale_x_continuous(breaks=c(-78,-76,-74,-72), labels=function(x) paste0(abs(x),"\u00b0W")), scale_y_continuous(breaks=c(2,4,6,8,10), labels=function(y) paste0(y,"\u00b0N")))
pa <- ggplot() + base() + geom_raster(data=A, aes(x, y), fill="#7FBF7F") + geom_sf(data=ant, fill=NA, color="black", linewidth=0.5) +
  geom_point(data=sites, aes(x, y), shape=21, fill="red", size=1.8) +
  labs(title=ifelse(version=="final","UPRA Hass aptitude within the sampled envelope","Areas within the thermal envelope of the farms")) + cz()
mp <- function(v, t){ ggplot() + base() + geom_raster(data=nd, aes(x, y, fill=.data[[v]])) + scale_fill_viridis_c(name=t) +
  geom_sf(data=ant, fill=NA, color="black", linewidth=0.4) + geom_point(data=sites, aes(x, y), shape=21, fill="white", size=1.4) + labs(title=paste("Predicted", sub(" \\(.*","",t))) + cz() }
pd <- ggplot() + base() + geom_raster(data=nd, aes(x, y, fill=class)) + scale_fill_manual(values=c("#440154","#21908C","#FDE725"), name="Lipid quality zone") +
  geom_sf(data=ant, fill=NA, color="black", linewidth=0.4) + geom_point(data=sites, aes(x, y), shape=21, fill="white", size=1.4) + labs(title="Quality zones (main harvest)") + cz()
F7 <- (pa + mp("oleic", traits["oleic"]) + mp("sfa_ufa", traits["sfa_ufa"]) + pd) + plot_layout(ncol=2) +
  plot_annotation(tag_levels="a", tag_prefix="(", tag_suffix=")") & theme(plot.title=element_text(size=10, face="bold"))
ggsave(sprintf("out/figs/Fig7_zoning_maps_%s.png", version), F7, width=10, height=12, dpi=300, bg="white", device=ragg::agg_png)
## summary by department (area in km2 of each class)
cell_km2 <- cellSize(env[[1]], unit="km")
nd$km2 <- extract(cell_km2, as.matrix(nd[,c("x","y")]))[,1]
pts <- st_as_sf(nd, coords=c("x","y"), crs=4326); pts <- st_join(pts, col)
S <- pts %>% st_drop_geometry() %>% group_by(Department=name, class) %>% summarise(km2=round(sum(km2)), .groups="drop") %>%
  tidyr::pivot_wider(names_from=class, values_from=km2, values_fill=0)
write.csv(S, sprintf("out/data/zoning_area_by_department_%s.csv", version), row.names=FALSE); print(as.data.frame(S))
saveRDS(list(Q=Q, version=version), "maps_meta.rds")
