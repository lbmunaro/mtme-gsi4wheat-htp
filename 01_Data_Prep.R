# Pheno data prep

rm(list = objects())

# Packages ----
library(tidyverse)
library(janitor)
source('Functions_HTP.R')

# Load data ----

## Pheno ----

ILYT_Pheno_HTP <- readRDS('Data/ILYT_Pheno.rds') |>
  filter(Year%in%c('2022','2023')) |>
  filter(Loc%in%c('Neoga, IL', 'St Peter, IL', 'Urbana, IL')) |>
  droplevels() |>
  group_by(Trait) |>
  mutate(
    Pheno_z = as.vector(scale(Pheno_SI, center = T)),
    Pheno_mean = mean(Pheno_SI, na.rm = TRUE),
    Pheno_sd = sd(Pheno_SI, na.rm = TRUE)
  ) |>
  ungroup() |>
  glimpse()

saveRDS(ILYT_Pheno_HTP, file = 'Data/ILYT_Pheno_HTP.rds')

## Env info ----
env.info <- data.frame(
  Env=unique(ILYT_Pheno_HTP$Env),
  PlantDate = c('10/22/2021',
                '09/29/2021',
                '09/28/2021',
                '10/10/2022',
                '10/18/2022',
                '09/30/2022'),
  lat = rep(c(39.23274, 38.8712128, 40.05833),2),
  long = rep(c(-88.38207, -88.8825766, -88.22937),2)
) |>
  glimpse()

## HTP ----
ILYT_Pheno_HTP_all <- read.csv('Data/HTP_22-24_2025.05.07.csv') |>
  remove_empty(which = c('cols')) |> # Remove empty columns.
  clean_names() |> # Clean names.
  # Replace location name.
  mutate(
    study_name = str_replace(study_name,'YT_Stj_22','YT_Addie_22'),
    location_name = str_replace(location_name,'St. Jacob Township, IL','Addieville, IL')
  ) |>
  select(study_name, study_year, location_name, germplasm_name,
         block_number, col_number, row_number, plot_number,
         canopy_reflectance_blue_reflectance_ratio_day_102_comp_0000079:canopy_reflectance_red_edge_reflectance_ratio_day_97_comp_0001188) |>
  # Simplify study name and convert it to a factor.
  mutate(study_name=as.factor(gsub('^YT_', '', study_name)),
         study_name=as.factor(gsub('^Addie_', 'Adv_', study_name)),
         study_name=str_replace(study_name, '(\\w+)_(\\d+)', '\\2-\\1'),
  ) |>
  rename(
    Env=study_name,
    Year=study_year,
    Loc=location_name,
    Gen=germplasm_name,
    Block=block_number,
    Col=col_number,
    Row=row_number,
    Plot=plot_number
  ) |> 
  rename_with(
    .cols = starts_with("canopy_reflectance_"),
    .fn   = ~ .x |>
      str_remove_all("^canopy_reflectance_") |>
      str_remove_all("reflectance_ratio_") |>
      str_remove_all("day_") |>
      str_remove_all("_comp_\\d+$") |>
      str_replace_all("red_edge","rededge")
  ) |> 
  pivot_longer(
    cols = matches("^(blue|green|red|nir|rededge)_"),
    names_to = c('MSband', 'FlyDoy'),
    names_sep = '_'
  ) |>
  mutate(
    FlyDoy = as.integer(FlyDoy),
    FlyDate = as.Date(FlyDoy-1, origin=paste0(Year,'-01-01'))
  ) |>
  pivot_wider(
    names_from = MSband,
    values_from = value
  ) |>
  mutate(IDEU = as.factor(paste(Env,Plot, sep='_'))) |>
  relocate(IDEU,.before = Env) |>
  mutate_at(vars(Env:FlyDoy), ~as.factor(.)) |>
  arrange(IDEU, FlyDoy) |> 
  group_by(IDEU,FlyDoy) |>
  mutate(NDVI = ((nir-red)/(nir+red)),
         NDRE = ((nir-rededge)/(nir+rededge)),
         GNDVI = ((nir-green)/(nir+green))
  ) |>
  ungroup() |>
  select(-Gen) |>
  group_by(Env, FlyDoy) |>
  filter(!is.na(mean(NDVI, na.rm = TRUE))) |>
  ungroup() |>
  filter(Year%in%c('2022','2023')) |>
  left_join(
    ILYT_Pheno_HTP |>
      select(IDEU, Gen, Gkeep, Gdrop, Trait, Pheno_SI) |>
      pivot_wider(names_from = Trait, values_from = Pheno_SI), 
    by = join_by(IDEU)
  ) |>
  left_join(
    env.info,
    by = join_by(Env)
  ) |>
  select(IDEU:Plot,Gen:TW,blue:GNDVI,FlyDoy,FlyDate,PlantDate:long) |>
  mutate(Env_FlyDoy = paste(Env,FlyDoy,sep='_')) |>
  filter(!Env_FlyDoy %in% c("22-Urb_122", "22-Stp_130", "23-Neo_164", '24-Neo_106')) |>
  mutate(PlantDate = as.Date(PlantDate, format = "%m/%d/%Y")) |>
  mutate(across(blue:GNDVI, ~ ifelse(IDEU == '24-Urb_958', NA, .))) |>
  droplevels() |>
  # calculate gdd for each Env and Flight
  group_by(Env, FlyDate) |>
  mutate(GDD = gdd(lat = unique(lat), long = unique(long),
                   plntDate = unique(PlantDate), date = unique(FlyDate)),
         GDD = round(GDD,0)) |>
  ungroup()

str(ILYT_Pheno_HTP_all)

## Save data ----
saveRDS(ILYT_Pheno_HTP_all, file = 'Data/ILYT_Pheno_HTP_all.rds')
