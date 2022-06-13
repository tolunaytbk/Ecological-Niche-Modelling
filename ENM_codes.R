setwd("C:/Users/tolun/Desktop/term_project_new/")
getwd()

##Library installation
#install.packages("tibble")
#install.packages("dplyr")
#install.packages("biomod2")
#install.packages("raster")
#install.packages("rasterVis")
#install.packages("gridExtra")

##Library loading
library("dplyr")
library("rgdal")
library("ggplot2")
library("rnaturalearth")
library("tibble")
library("readr")
library("biomod2")
library("rasterVis")
library("gridExtra")

# This packages are necessary for the markdown layout probably,
#   they should be installed before opening of markdown.
library("yaml")
library("rmarkdown")
library("rticles")

# Data Set
gbif_lynx <- read.delim("data/gbif/0075557-210914110416597.csv")
iucn_map_lynx <- readOGR("data/iucn/data_0.shp")

# Data Cleaning
gbif_lynx<- select(gbif_lynx, c("countryCode", "decimalLatitude", "decimalLongitude"))

gbif_lynx <- filter(gbif_lynx, +
                      decimalLatitude != is.na(1), +
                      decimalLongitude!= is.na(1))
                    


#  Clade naming 
## Clade A: CN
## Clade B: TR, ME, RS, MK
## Clade C: PL, LV, RV, SK, RO, NO, MN, CH

# Clade A
cn_data <- filter(gbif_lynx, countryCode == c("CN")) # 5 obs.

clade_a_data <- cn_data

clade_a_data <- clade_a_data %>%
  add_column(lineage = c("clade_a"))
head(clade_a_data)

# Clade B

tr_data <- filter(gbif_lynx, countryCode == c("TR")) # 1 obs.
#me_data <- filter(gbif_lynx, countryCode == c("ME")) # 0 obs.
rs_data <- filter(gbif_lynx, countryCode == c("RS")) # 2 obs.

clade_b_data <- full_join(tr_data, rs_data)

clade_b_data <- clade_b_data %>%
  add_column(lineage = c("clade_b"))
head(clade_b_data)


# Clade C

no_data <- filter(gbif_lynx, countryCode == c("NO"))  #25813 obs.
new_no_data <- sample_n(no_data, 200, replace = FALSE) 

ch_data <- filter(gbif_lynx, countryCode == c("CH"))  #29602 obs.
new_ch_data <- sample_n(ch_data, 75, replace = FALSE)

pl_data <- filter(gbif_lynx, countryCode == c("PL"))  #48 obs.
new_pl_data <- sample_n(pl_data, 7, replace = FALSE)

ro_data <- filter(gbif_lynx, countryCode == c("RO"))  #02 obs.
new_ro_data <- sample_n(ro_data, 1, replace = FALSE)

ru_data <- filter(gbif_lynx, countryCode == c("RU"))  #49 obs.
new_ru_data <- sample_n(ru_data, 7, replace = FALSE)

#lv_data <- filter(gbif_lynx, countryCode == c("LV")) #00 obs.
#new_lv_data <- sample_n(lv_data, 50, replace = FALSE)

sk_data <- filter(gbif_lynx, countryCode == c("SK"))  #12 obs.
new_sk_data <- sample_n(sk_data, 2, replace = FALSE)

mn_data <- filter(gbif_lynx, countryCode == c("MN"))  #07 obs.
new_mn_data <- sample_n(mn_data, 1, replace = FALSE)


#### Clade C Lineage Data 

clade_c_data <- full_join(new_no_data, new_ch_data)
clade_c_data <- full_join(clade_c_data, pl_data)
clade_c_data <- full_join(clade_c_data, ro_data)
clade_c_data <- full_join(clade_c_data, ru_data)
clade_c_data <- full_join(clade_c_data, sk_data)
clade_c_data <- full_join(clade_c_data, mn_data)

clade_c_data <- clade_c_data %>%
  add_column(lineage = c("clade_c"))
head(clade_c_data)

# Data frame of all clades 

clades_data <- full_join(clade_a_data, clade_b_data)
clades_data <- full_join(clades_data, clade_c_data)

clades_data <- clades_data %>%
  add_column(occurrenceStatus = c(1))
head(clades_data)

clade_b_data <- clade_b_data %>%
  add_column(occurrenceStatus = c(1))
head(clade_b_data)

# Base Map downloading
Land50 <- ne_download(scale = 50, type = "land", category = "physical")
Coast50 <- ne_download(scale = 50, type = "coastline", category = "physical")

#Plotting the graph

ggplot()+
  geom_polygon(data = Land50, aes(x=long, y=lat, group=group), fill="#BC9354") +
  geom_polygon(data = iucn_map_lynx,aes(x=long, y=lat, group=group), fill="orange") +
  geom_point(data = clades_data, aes(x=decimalLongitude, y=decimalLatitude, fill= lineage), shape=21, size=2)+
  theme(panel.background = element_rect(fill="#2C8598"), element_line(color ="#2C8598")) +
  coord_fixed(xlim = c(0,155), ylim= c(35,70), 1.3)#for the cutting the map
  

# Saving the file
ggsave("distribution_map3.png", type= "cairo-png", dpi = 1200, width = 8, height = 4)


#WorldClim Data

library(raster)

bioclim_data <- raster :: stack (
  c(
    bio1='data/wc2.1_5m_bio_1.tif',
    bio6='data/wc2.1_5m_bio_6.tif',
    bio11='data/wc2.1_5m_bio_11.tif',
    bio15='data/wc2.1_5m_bio_15.tif'
  )
)
plot(bioclim_data)

#Data Formating for Entire Lynx Lynx

lynx_data <- BIOMOD_FormatingData(
  resp.var = clades_data$occurrenceStatus,
  expl.var = bioclim_data,
  resp.xy = clades_data[, c("decimalLongitude", "decimalLatitude")],
  resp.name = "Lynx.lynx",
  PA.nb.rep = 1,
  PA.nb.absences = 500,
  PA.strategy = 'random',
  na.rm=TRUE
) 





# Modeling for Entire Lynx Lynx

lynx_opt <- BIOMOD_ModelingOptions(
  GLM = list(type='quadratic', interaction.level=1),
  GAM = list(algo ='GAM_mgcv')
)

lynx_models <- BIOMOD_Modeling(
  data = lynx_data,
  models = c("GLM","GAM"),
  models.options = lynx_opt,
  NbRunEval = 1,
  DataSplit = 80,
  VarImport = 3,
  do.full.models = FALSE,
  modeling.id = "demo1",
)

lynx_models_scores <- get_evaluations(lynx_models)
clade_b_models_scores <-get_evaluations(clade_b_models)


dim(lynx_models_scores)
dimnames(lynx_models_scores)

dim(clade_b_models_scores)
dimnames(clade_b_models_scores)

# Model Score Graph for Entire Lynx Lynx

models_scores_graph(
  lynx_models,
  by ="models",
  metrics = c("ROC", "TSS"),
  xlim = c(0.5,1), 
  ylim = c(0.5,1)
)

models_scores_graph(
  lynx_models,
  by ="data_set",
  metrics = c("ROC", "TSS"),
  xlim = c(0.5,1), 
  ylim = c(0.5,1)
)


# Model Score Graph for Clade B

models_scores_graph(
  clade_b_models,
  by ="models",
  metrics = c("ROC", "TSS"),
  xlim = c(0,0.5), 
  ylim = c(0,0.5)
)

models_scores_graph(
  clade_b_models,
  by ="data_set",
  metrics = c("ROC", "TSS"),
  xlim = c(0,0.5), 
  ylim = c(0,0.5)
)


(lynx_models_var_import <- get_variables_importance(lynx_models))
apply(lynx_models_var_import, c(1, 2), mean)

lynx_glm <- BIOMOD_LoadModels(lynx_models, models = 'GLM')
lynx_gam <- BIOMOD_LoadModels(lynx_models, models = 'GAM')


(clade_b_models_var_import <- get_variables_importance(clade_b_models))
apply(clade_b_models_var_import, c(1, 2), mean)

clade_b_glm <- BIOMOD_LoadModels(clade_b_models, models = 'GLM')
clade_b_gam <- BIOMOD_LoadModels(clade_b_models, models = 'GAM')

# Response Plot

glm_eval_strip <- biomod2::response.plot2(
  models = lynx_glm,
  Data = get_formal_data(lynx_models, 'expl.var'),
  show.variables = get_formal_data(lynx_models, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species =get_formal_data(lynx_models, 'resp.var')
)


gam_eval_strip <- biomod2::response.plot2(
  models = lynx_gam,
  Data = get_formal_data(lynx_models, 'expl.var'),
  show.variables = get_formal_data(lynx_models, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species =get_formal_data(lynx_models, 'resp.var')
)

# Response Plot for Clade B

b_glm_eval_strip <- biomod2::response.plot2(
  models = clade_b_glm,
  Data = get_formal_data(clade_b_models, 'expl.var'),
  show.variables = get_formal_data(clade_b_models, 'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  legend = FALSE,
  display_title = FALSE,
  data_species =get_formal_data(clade_b_models, 'resp.var')
)


# Ensemble Models for Entire Lynx Lynx

lynx_ensemble_models <- BIOMOD_EnsembleModeling(
  modeling.output = lynx_models,
  em.by = 'all',
  eval.metric = 'ROC',
  eval.metric.quality.threshold = 0.8,
  models.eval.meth = c('ROC'),
  prob.mean = FALSE,
  prob.cv = TRUE,
  committee.averaging = TRUE, 
  prob.mean.weight = TRUE,
  VarImport = 0,
)


(lynx_ensemble_models_scores <- get_evaluations(lynx_ensemble_models))

# Projection of Models for Entire Lynx Lynx

lynx_models_project_current <- BIOMOD_Projection(
  modeling.output = lynx_models,
  new.env = bioclim_data,
  proj.name = "current",
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)


lynx_ensemble_models_current <- BIOMOD_EnsembleForecasting(
  EM.output = lynx_ensemble_models,
  projection.output = lynx_models_project_current,
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)

plot(lynx_ensemble_models_current, str.grep= "EMca|EMwmean")


# Future Data

bioclim_future <- raster::stack("data/share/spatial03/worldclim/cmip6/7_fut/5m/MIROC6/ssp370/wc2.1_5m_bioc_MIROC6_ssp370_2081-2100.tif")

summary(bioclim_future)


bio1_f <- subset(bioclim_future, "wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.1")
bio6_f <- subset(bioclim_future, "wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.6")
bio11_f <- subset(bioclim_future, "wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.11")
bio15_f <- subset(bioclim_future, "wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.15")

bioclim_future_data <- stack(bio1_f, bio6_f, bio11_f, bio15_f)

summary(bioclim_future_data)
summary(bioclim_data)

names(bioclim_future_data)[names(bioclim_future_data) == 'wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.1'] <- "bio1"
names(bioclim_future_data)[names(bioclim_future_data) == 'wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.6'] <- "bio6"
names(bioclim_future_data)[names(bioclim_future_data) == 'wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.11'] <- "bio11"
names(bioclim_future_data)[names(bioclim_future_data) == 'wc2.1_5m_bioc_MIROC6_ssp370_2081.2100.15'] <- "bio15"

summary(bioclim_future_data)
summary(bioclim_data)

lynx_models_project_2081 <- BIOMOD_Projection(
  modeling.output = lynx_models,
  new.env = bioclim_future_data,
  proj.name = "future",
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)

lynx_ensemble_models_2081 <- BIOMOD_EnsembleForecasting(
  EM.output = lynx_ensemble_models,
  projection.output = lynx_models_project_2081,
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)

plot(lynx_ensemble_models_2081, str.grep= "EMca|EMwmean")

lynx_bin_proj_current <- stack(
  c(
    ca = "Lynx.lynx/proj_current/individual_projections/Lynx.lynx_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd",
    wm = "Lynx.lynx/proj_current/individual_projections/Lynx.lynx_EMwmeanByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd"
    )
  )


lynx_bin_proj_future <- stack(
  c(
    ca = "Lynx.lynx/proj_future/individual_projections/Lynx.lynx_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd",
    mw = "Lynx.lynx/proj_future/individual_projections/Lynx.lynx_EMwmeanByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd"
  )
)

SRC_lynx_current_2081_ssp370 <- BIOMOD_RangeSize(
  lynx_bin_proj_current, 
  lynx_bin_proj_future
)

SRC_lynx_current_2081_ssp370$Compt.By.Models

lynx_SRC_map <- stack(SRC_lynx_current_2081_ssp370$Diff.By.Pixel)


plot(lynx_SRC_map)

my.at <- seq(-2.5, 1.5, 1)
my.Colorkey <- list(
  at = my.at,
  labels =
    list(
      labels = c("lost", "presence", "absence", "gain"),
      at = my.at[-1]-0.5
    )
)

rasterVis::levelplot(
  lynx_SRC_map,
  main = "Lynx Lynx Range Change",
  colorkey = my.Colorkey,
  col.regions = c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1),
  dpi = 600,
)

#Data Formating for Clade B of Lynx Lynx

clade_b_biomod <- BIOMOD_FormatingData(
  resp.var = clade_b_data$occurrenceStatus,
  expl.var = bioclim_data,
  resp.xy = clade_b_data[, c("decimalLongitude", "decimalLatitude")],
  resp.name = "clade_b",
  PA.nb.rep = 1,
  PA.nb.absences = 5,
  PA.strategy = 'random',
  na.rm= TRUE
) 


# Modelling

clade_b_models <- BIOMOD_Modeling(
  data = clade_b_biomod,
  models = c("GLM"),
  NbRunEval = 1,
  DataSplit = 80,
  VarImport = 3,
  do.full.models = 'FALSE',
  modeling.id = 'demo2',
  Prevalence = 0.5,
  models.eval.meth = 'ROC', 
  SaveObj = TRUE,
  rescal.all.models = TRUE
)

# Ensemble Models for Clade B

clade_b_ensemble_models <- BIOMOD_EnsembleModeling(
  modeling.output = clade_b_models,
  em.by = 'all',
  eval.metric = 'ROC',
  eval.metric.quality.threshold = NULL,
  models.eval.meth = c('ROC'),
  prob.mean = FALSE,
  prob.cv = TRUE,
  committee.averaging = TRUE, 
  prob.mean.weight = TRUE,
  VarImport = 0
)

# Clade B Projection (Current) 

clade_b_proj_current <- BIOMOD_Projection(
  modeling.output = clade_b_models,
  new.env = bioclim_data,
  proj.name = "current",
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE,
  NbRunEval = 1,
  do.full.models = 'FALSE'
)

clade_b_ensemble_models_current <- BIOMOD_EnsembleForecasting(
  EM.output = clade_b_ensemble_models,
  projection.output = clade_b_proj_current,
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE,
  NbRunEval = 1,
  do.full.models = 'FALSE'
)


clade_b_models_project_2081 <- BIOMOD_Projection(
  modeling.output = clade_b_models,
  new.env = bioclim_future_data,
  proj.name = "future",
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)

clade_b_ensemble_models_2081 <- BIOMOD_EnsembleForecasting(
  EM.output = clade_b_ensemble_models,
  projection.output = clade_b_models_project_2081,
  binary.meth = "ROC",
  output_format = ".grd",
  do.stack = FALSE
)

plot(clade_b_ensemble_models_2081, str.grep= "EMca|EMwmean")



# Clade B Binary Data Stack

clade_b_bin_proj_current <- stack(
  c(
    ca = "clade.b/proj_current/individual_projections/Clade.b_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd",
    wm = "clade.b/proj_current/individual_projections/Clade.b_EMwmeanByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd"
  )
)


clade_b_bin_proj_future <- stack(
  c(
    ca = "clade.b/proj_Future/individual_projections/Clade.b_EMcaByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd",
    wm = "clade.b/proj_Future/individual_projections/Clade.b_EMwmeanByROC_mergedAlgo_mergedRun_mergedData_ROCbin.grd"
  )
)

# Clade B Mapping

SRC_clade_b_current_2081_ssp370 <- BIOMOD_RangeSize(
  clade_b_bin_proj_current, 
  clade_b_bin_proj_future
)

SRC_clade_b_current_2081_ssp370$Compt.By.Models

clade_b_SRC_map <- stack(SRC_clade_b_current_2081_ssp370$Diff.By.Pixel)

# Clade B Map plotting

plot(clade_b_SRC_map)
  my.at <- seq(-2.5, 1.5, 1)
  my.Colorkey <- list(
  at = my.at,
  labels =
    list(
      labels = c("lost", "presence", "absence", "gain"),
      at = my.at[-1]-0.5)) 

  rasterVis::levelplot(
    clade_b_SRC_map,
    main = "Clade B Range Change",
    colorkey = my.Colorkey,
    col.regions = c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
    layout = c(2,1),
    dpi = 600
    )



