## download_gis_layers_Colombia.R
## Run on a computer with internet access. It downloads free layers, crops them to Colombia
## and writes a small zip (gis_final.zip) to upload back for the final zoning maps (Fig. 7).
## Requires: install.packages(c("geodata", "terra"))
library(geodata); library(terra)
out <- "gis_final"; dir.create(out, showWarnings = FALSE)
col <- gadm("COL", level = 1, path = out)                                   # departments (GADM)
tavg <- worldclim_country("COL", var = "tavg", res = 0.5, path = out)       # WorldClim 2.1, 30 arc-s, 12 months
writeRaster(mean(tavg), file.path(out, "tavg_col.tif"), overwrite = TRUE)   # annual mean temperature (°C)
prec <- worldclim_country("COL", var = "prec", res = 0.5, path = out)
writeRaster(sum(prec), file.path(out, "prec_col.tif"), overwrite = TRUE)    # annual precipitation (mm)
ph <- soil_world(var = "phh2o", depth = 60, path = out)                     # SoilGrids 2.0, pH(H2O) x10, 30-60 cm
writeRaster(crop(ph, col, mask = TRUE), file.path(out, "phh2o_30_60_col.tif"), overwrite = TRUE)
## UPRA layer: download manually from https://sipra.upra.gov.co (Zonificación de aptitud para el cultivo
## comercial de aguacate Hass, escala 1:100.000) and copy the shapefile (.shp, .shx, .dbf, .prj) into gis_final/upra/
dir.create(file.path(out, "upra"), showWarnings = FALSE)
unlink(list.files(out, pattern = "^wc2|soil_world|gadm", full.names = TRUE, recursive = TRUE))  # keep only cropped layers
zip("gis_final.zip", list.files(out, recursive = TRUE, full.names = TRUE))
