rm(list = objects())

library(tidyverse)
source('Functions_HTP.R')


set.seed(1)
n_folds <- 5

ILYT_Pheno_HTP <- readRDS('Data/ILYT_Pheno_HTP.rds') |>
  mutate(Index='none') |>
  glimpse()

ILYT_Pheno_HTP <- ILYT_Pheno_HTP |>
  # --- CV1: genotype-based fold assignment ---
  left_join(
    tibble(Gen = unique(ILYT_Pheno_HTP$Gen)) |>
      mutate(CV1_fold = sample(rep(1:n_folds, length.out = n()))),
    by = "Gen"
  ) |>
  # --- CV2: IDEU (GxE) unit-level fold assignment ---
  group_by(Gen) |>
  mutate(CV2_fold = sample(1:n_folds, size = n(), replace = n() > n_folds)) |>
  ungroup() |>
  glimpse()

# Check result
ILYT_Pheno_HTP |> count(CV1_fold)
ILYT_Pheno_HTP |> count(CV2_fold)

ILYT_Pheno_HTP_all <- readRDS('Data/ILYT_Pheno_HTP_all.rds') |>
  select(-c(blue:rededge)) |>
  glimpse()

str(ILYT_Pheno_HTP_all)

ggplot(ILYT_Pheno_HTP_all, aes(x = as.numeric(as.character(GDD)), y=NDRE, group = Plot, colour = Plot)) +
  geom_line(linewidth=0.5, alpha=0.1) +
  facet_wrap(~Env) +
  theme_minimal() +
  theme(legend.position = "none")

# test function
ILYT_Pheno_HTP_all |>
  filter(Env == "22-Neo", Plot == "321") |>
  extract_htp_derivatives(time_col = "GDD", value_col = "NDVI")

# Apply for each VI per plot within Env
HTP_Deriv_out <- purrr::map_dfr(index_vars, function(index) {
  ILYT_Pheno_HTP_all |>
    group_by(Env, Plot) |>
    group_modify(~ extract_htp_derivatives(.x, time_col = "GDD", value_col = index)) |>
    mutate(Index = index) |>
    ungroup()
})

HTP_Deriv <- HTP_Deriv_out |>
  mutate(Index=as.factor(Index)) |>
  pivot_longer(meanVI:GDD_.9MI, names_to = 'Trait', values_to = 'Pheno') |>
  mutate(Trait=as.factor(Trait)) |>
  group_by(Trait,Index) |>
  mutate(Pheno_SI = Pheno,
         Pheno_z = as.vector(scale(Pheno_SI, center = T)),
         Pheno_mean = mean(Pheno_SI, na.rm = TRUE),
         Pheno_sd = sd(Pheno_SI, na.rm = TRUE)) |>
  ungroup() |>
  mutate(TraitEnv=as.factor(paste(Trait,Env,sep = '-'))) |>
  left_join(
    ILYT_Pheno_HTP |>
      select(Env:Gdrop) |>
      group_by(IDEU) |>
      summarise_all(~unique(.))
  ) |>
  arrange(TraitEnv, Col, Row) |>
  glimpse()

ILYT_Pheno_Deriv <- bind_rows(ILYT_Pheno_HTP, HTP_Deriv) |>
  mutate(Index=as.factor(Index)) |>
  glimpse()

# Split dataset

# Define HTP and Conv traits

conv_traits <- as.character(unique(ILYT_Pheno_HTP$Trait))
htp_traits <- as.character(unique(HTP_Deriv$Trait))

ILYT_Pheno_Dsplit_NDVI <- map(htp_traits, function(htp_trait) {
  ILYT_Pheno_Deriv |> 
    filter(Index%in%c('none','NDVI')) |>
    filter(Trait %in% c(conv_traits, htp_trait)) |>
    droplevels()
})

# Name each element in the list with its HTP trait
names(ILYT_Pheno_Dsplit_NDVI) <- htp_traits

ILYT_Pheno_Dsplit_NDRE <- map(htp_traits, function(htp_trait) {
  ILYT_Pheno_Deriv |> 
    filter(Index%in%c('none','NDRE')) |>
    filter(Trait %in% c(conv_traits, htp_trait)) |>
    droplevels()
})

# Name each element in the list with its HTP trait
names(ILYT_Pheno_Dsplit_NDRE) <- htp_traits

ILYT_Pheno_Dsplit_GNDVI <- map(htp_traits, function(htp_trait) {
  ILYT_Pheno_Deriv |> 
    filter(Index%in%c('none','GNDVI')) |>
    filter(Trait %in% c(conv_traits, htp_trait)) |>
    droplevels()
})

# Name each element in the list with its HTP trait
names(ILYT_Pheno_Dsplit_GNDVI) <- htp_traits

ILYT_Pheno_HTP.Dsplit <- c()
ILYT_Pheno_HTP.Dsplit$Conv <- ILYT_Pheno_Deriv|>filter(Trait%in%conv_traits)
ILYT_Pheno_HTP.Dsplit$NDVI <- ILYT_Pheno_Dsplit_NDVI
ILYT_Pheno_HTP.Dsplit$NDRE <- ILYT_Pheno_Dsplit_NDRE
ILYT_Pheno_HTP.Dsplit$GNDVI <- ILYT_Pheno_Dsplit_GNDVI

glimpse(ILYT_Pheno_HTP.Dsplit)


saveRDS(ILYT_Pheno_HTP.Dsplit, file = 'Data/ILYT_Pheno_HTP.Dsplit.rds')

