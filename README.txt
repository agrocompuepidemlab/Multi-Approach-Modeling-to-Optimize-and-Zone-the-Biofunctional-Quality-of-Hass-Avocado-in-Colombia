Supplementary code — Multi-approach modeling and zoning of Hass avocado quality (Colombia)
Run from this folder in R >= 4.3 (UTF-8 locale), in order:
  R/01_clean.R        audit and correction of the original spreadsheet (writes out/data/*corrected.csv, correction_log.csv)
  R/02_models.R       MLR, GAM, RF x selection methods 1-4 (saturated, stepwise, SBF, RFE), leave-one-farm-out validation
  R/02b_knowledge.R   selection method 5 (knowledge-based: Tmean, subsoil pH, dry matter, harvest), best model, y-randomization, validity
  R/03_figures.R      Figures 2-6 and S1-S3 (style of the original manuscript)
  R/04_maps.R         Figure 7 and zoning tables. Uses gis/final/ (WorldClim 30", SoilGrids pH 30-60 cm, UPRA layer) if present;
                      otherwise the provisional version with gis/wc10 (WorldClim 1.4, 10', from github.com/GBIF-Europe/nordic_oikos_2018_r)
                      and gis/states10.rda (Natural Earth admin-1, from github.com/ropensci/rnaturalearthhires).
  R/05_tables.R       Tables 1, S1-S4
download_gis_layers_Colombia.R: run on a computer with internet to build gis_final.zip; unzip into gis/final/ (add the UPRA shapefile in gis/final/upra/) and rerun R/04_maps.R.
Packages: readxl dplyr tidyr caret randomForest ranger mgcv ggplot2 patchwork viridis terra sf ragg. Seeds fixed (2025; L'Ecuyer-CMRG).
