#### Script for the analyses of the paper "Community assembly of 
#### serpentine plants is driven by interspecific differences in 
#### functional traits". 
#### Created by Guillaume Delhaye. g.delhaye@kew.org
#### Last update 14/02/24

library(tidyverse)
library(brms)
library(tidybayes)
library(FD)
library(rstan)
library(labdsv)
library(marginaleffects)
library(modelr)
source("Functions.R")

# load data
soil.b <- read.csv("./data/soil_clean.csv") |> select(-X)
traits_sp_site <- read.csv("./data/traits_sp_site.csv") |> select(-X)
traits_sp <- read.csv("./data/traits_sp.csv") |> select(-X)
abun <- read.csv("./data/abundance_clean.csv") |> select(-X)

######  Test for trait - environment relationships -----
## 1 Using MLM 3 of ter Braak 2019  ----
# Assemble dataset for the model 
ab <- data.frame(abun) |>
  # transform as integers using "ceiling" to avoid transforming small values into 0
  ceiling() |> 
  rownames_to_column("site") |>
  pivot_longer(cols = !site, names_to = "species", values_to = "abundance") |>
  left_join(traits_sp |> rename(species = Sp)) |>  
  left_join(soil.b |> select(id, Ni) |> rename(site = id))|>
  mutate_at(c("LS","LeafArea", "LDMC", "SLA", "LT", "LNC", "Ni"), scale2)
hist(ab$abundance, breaks = 20)

# To run analyses withou A.lesbiacum
# ab <- ab |> filter(species != "ALYLESB")

mod_SLA <- brm(abundance ~ Ni * SLA  + (1 + SLA | site) + (1 + Ni | species), 
                family = poisson(link = "log"),
                control=list(adapt_delta=0.95, 
                         max_treedepth=13),
                iter = 3000, warmup = 1000, chains = 4, cores = 4,
                data=ab)

mod_LT <- update(mod_SLA, abundance ~ Ni * LT  + (1 + LT|site) + (1 + Ni|species), newdata = ab)
mod_LA <- update(mod_SLA, abundance ~ Ni * LeafArea  + (1 + LeafArea|site) + (1 + Ni|species), newdata = ab)
mod_LDMC <- update(mod_SLA, .~Ni*LDMC+(1+LDMC|site)+(1+Ni|species), newdata = ab)
mod_LNC <- update(mod_SLA, .~Ni*LNC+(1+LNC|site)+(1+Ni|species), newdata = ab)
mod_LS <- update(mod_SLA, .~Ni*LS+(1+LS|site)+(1+Ni|species), newdata = ab)

saveRDS(mod_LT, file = "../output/trait_gradient_model/mod_LT_pois.rds")
saveRDS(mod_SLA, file = "../output/trait_gradient_model/mod_SLA_pois.rds")
saveRDS(mod_LDMC, file = "../output/trait_gradient_model/mod_LDMC_pois.rds")
saveRDS(mod_LNC, file = "../output/trait_gradient_model/mod_LNC_pois.rds")
saveRDS(mod_LA, file = "../output/trait_gradient_model/mod_LA_pois.rds")
saveRDS(mod_LS, file = "../output/trait_gradient_model/mod_LS_pois.rds")


plot_sla <- plot_predictions(mod_SLA, condition = list("Ni", SLA = "threenum"))
plot_lt <- plot_predictions(mod_LT, condition = list("Ni", LT = "threenum"))
plot_la <- plot_predictions(mod_LA, condition = list("Ni", LeafArea = "threenum"))
plot_ldmc <- plot_predictions(mod_LDMC, condition = list("Ni", LDMC = "threenum"))
plot_lnc <- plot_predictions(mod_LNC, condition = list("Ni", LNC = "threenum"))
plot_ls <- plot_predictions(mod_LS, condition = list("Ni", LS = "threenum"))

ggpubr::ggarrange(plot_la, plot_sla, plot_lt,
                  plot_ldmc, plot_lnc, plot_ls, 
                  common.legend = TRUE)

mod_draws <- rbind(
  mlt = m_ex(mod_LT, "LT"),
  mla = m_ex(mod_LA, "LA"),
  msla = m_ex(mod_SLA, "SLA"),
  mldmc = m_ex(mod_LDMC, "LDMC"),
  mlnc = m_ex(mod_LNC, "LNC"),
  mls <- m_ex(mod_LS, "LS"))|>
  mutate(trait = fct_relevel(c("LS", "LNC","LDMC","LT", "SLA","LA")))
#mod_draws$trait <- factor(mod_draws$trait, levels = c("LS", "LNC","LDMC","LT", "SLA","LA"))

mod_draws |>
  ggplot(aes(y = trait, x = var, fill = after_stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme(legend.position = "")

## Quantiles of the posterior
extract_quantile(mod_LA)
extract_quantile(mod_LT)
extract_quantile(mod_SLA)
extract_quantile(mod_LDMC)
extract_quantile(mod_LNC)
extract_quantile(mod_LS)

## 2. Community weighted mean and FDis change along the gradient -----
# Data preparation 
nrep = 500

## To run the analyses withou A. lesbiacum
# abun <- abun |> select(-ALYLESB)

ab_long <- abun |> 
  rename(id= site) |> 
  pivot_longer(where(is.numeric), names_to = "Sp", values_to = "abundance") |>
  separate("id", into = c("site"), remove = F) 

# Calculate CWM with median data   
CWM_sp <-  ab_long |>
  left_join(traits_sp) |>
  mutate(s_LA = abundance * LeafArea, 
         s_SLA = abundance * SLA,
         s_LT = abundance * LT,
         s_LDMC = abundance * LDMC,
         s_LNC = abundance * LNC,
         s_LS = abundance * LS) |>
  group_by(id) |>
  summarise(abundance = sum(abundance, na.rm = T), 
            sum_LA = sum(s_LA, na.rm = T), 
            sum_SLA = sum(s_SLA, na.rm = T), 
            sum_LT = sum(s_LT, na.rm = T), 
            sum_LDMC = sum(s_LDMC, na.rm = T), 
            sum_LNC = sum(s_LNC, na.rm = T), 
            sum_LS = sum(s_LS, na.rm = T)) |>
  mutate(CWM_LA_inter = sum_LA/abundance,
         CWM_SLA_inter = sum_SLA/abundance,
         CWM_LT_inter = sum_LT/abundance,
         CWM_LDMC_inter = sum_LDMC/abundance,
         CWM_LNC_inter = sum_LNC/abundance,
         CWM_LS_inter = sum_LS/abundance) |>
  select(id, starts_with("CWM_"))

# Calculate CWM with sp_site data
CWM_sp_site <- ab_long |>
  left_join(traits_sp_site) |>
  mutate(s_LA = abundance * LeafArea, 
         s_SLA = abundance * SLA,
         s_LT = abundance * LT,
         s_LDMC = abundance * LDMC,
         s_LNC = abundance * LNC,
         s_LS = abundance * LS) |>
  group_by(id) |>
  summarise(abundance = sum(abundance, na.rm = T), 
            sum_LA = sum(s_LA, na.rm = T), 
            sum_SLA = sum(s_SLA, na.rm = T), 
            sum_LT = sum(s_LT, na.rm = T), 
            sum_LDMC = sum(s_LDMC, na.rm = T), 
            sum_LNC = sum(s_LNC, na.rm = T), 
            sum_LS = sum(s_LS, na.rm = T)) |>
  mutate(CWM_LA_intra = sum_LA/abundance,
         CWM_SLA_intra = sum_SLA/abundance,
         CWM_LT_intra = sum_LT/abundance,
         CWM_LDMC_intra = sum_LDMC/abundance,
         CWM_LNC_intra = sum_LNC/abundance,
         CWM_LS_intra = sum_LS/abundance) |>
  select(id, starts_with("CWM_"))

# assemble the CWM dataset with Ni gradient
CWM <- CWM_sp |> 
  left_join(CWM_sp_site) |> 
  rename(site = id)

### FDis 
abun_site <- abun |> 
  separate("site", c("id", "no"), remove = F) |>
  column_to_rownames("site")

ab_amp <- abun_site |> filter(id == "Amp") |> select(-c(id, no))
ab_oly <- abun_site |> filter(id == "Oly")|> select(-c(id, no))
ab_vat <- abun_site |> filter(id == "Vat")|> select(-c(id, no))
ab_lout <- abun_site |> filter(id == "Lout")|> select(-c(id, no))

tr_amp <- traits_sp_site |> filter(site == "Amp") |> column_to_rownames("Sp") |> select(-site)
tr_oly <- traits_sp_site |> filter(site == "Oly")|> column_to_rownames("Sp") |> select(-site)
tr_vat <- traits_sp_site |> filter(site == "Vat")|> column_to_rownames("Sp") |> select(-site)
tr_lout <- traits_sp_site |> filter(site == "Lout")|> column_to_rownames("Sp") |> select(-site)

## Multivariate FDis
fdism_amp <- sesFDism(traits = tr_amp, comm = ab_amp, numberReps = nrep)
fdism_oly <- sesFDism(traits = tr_oly, comm = ab_oly, numberReps = nrep)
fdism_vat <- sesFDism(traits = tr_vat, comm = ab_vat, numberReps = nrep)
fdism_lout<- sesFDism(traits = tr_lout, comm = ab_lout, numberReps = nrep)
FDM_ITV <- data.frame(FDis = c(fdism_amp, fdism_oly, fdism_vat, fdism_lout)) |>
  rownames_to_column("site")

## Univariate FDis
fdisu_amp_LS <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'LS', numberReps = nrep)
fdisu_oly_LS <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'LS', numberReps = nrep)
fdisu_vat_LS <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'LS', numberReps = nrep)
fdisu_lout_LS<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'LS', numberReps = nrep)

fdisu_amp_LA <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'LeafArea', numberReps = nrep)
fdisu_oly_LA <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'LeafArea', numberReps = nrep)
fdisu_vat_LA <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'LeafArea', numberReps = nrep)
fdisu_lout_LA<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'LeafArea', numberReps = nrep)

fdisu_amp_SLA <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'SLA', numberReps = nrep)
fdisu_oly_SLA <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'SLA', numberReps = nrep)
fdisu_vat_SLA <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'SLA', numberReps = nrep)
fdisu_lout_SLA<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'SLA', numberReps = nrep)

fdisu_amp_LT <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'LT', numberReps = nrep)
fdisu_oly_LT <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'LT', numberReps = nrep)
fdisu_vat_LT <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'LT', numberReps = nrep)
fdisu_lout_LT<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'LT', numberReps = nrep)

fdisu_amp_LDMC <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'LDMC', numberReps = nrep)
fdisu_oly_LDMC <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'LDMC', numberReps = nrep)
fdisu_vat_LDMC <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'LDMC', numberReps = nrep)
fdisu_lout_LDMC<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'LDMC', numberReps = nrep)

fdisu_amp_LNC <- sesFDisu(traits = tr_amp,  comm = ab_amp,  i = 'LNC', numberReps = nrep)
fdisu_oly_LNC <- sesFDisu(traits = tr_oly,  comm = ab_oly,  i = 'LNC', numberReps = nrep)
fdisu_vat_LNC <- sesFDisu(traits = tr_vat,  comm = ab_vat,  i = 'LNC', numberReps = nrep)
fdisu_lout_LNC<- sesFDisu(traits = tr_lout, comm = ab_lout, i = 'LNC', numberReps = nrep)

FDisU_ITV <- data.frame(FDis.LeafArea = c(fdisu_amp_LA, fdisu_oly_LA,fdisu_vat_LA, fdisu_lout_LA), 
                        FDis.SLA = c(fdisu_amp_SLA, fdisu_oly_SLA,fdisu_vat_SLA, fdisu_lout_SLA) , 
                        FDis.LDMC = c(fdisu_amp_LDMC, fdisu_oly_LDMC,fdisu_vat_LDMC, fdisu_lout_LDMC), 
                        FDis.LT =  c(fdisu_amp_LT, fdisu_oly_LT,fdisu_vat_LT, fdisu_lout_LT),
                        FDis.LNC =  c(fdisu_amp_LNC, fdisu_oly_LNC,fdisu_vat_LNC, fdisu_lout_LNC),
                        FDis.LS = c(fdisu_amp_LS, fdisu_oly_LS,fdisu_vat_LS, fdisu_lout_LS)) |>
  rownames_to_column("site")
# 
# Result_FD_ITV <- CWM |> left_join(FDM_ITV)  |> left_join(FDisU_ITV) |>
#   left_join(soil.b |> select(site = id, Ni))
# write.csv(Result_FD_ITV, "./output/FD_ITV.csv")

# # without A. lesbiacum
# Result_FD_ITV_noal <- CWM |> left_join(FDM_ITV)  |> left_join(FDisU_ITV) |>
# left_join(soil.b |> select(site = id, Ni))
# 
# colnames(Result_FD_ITV) == colnames(Result_FD_ITV_noal)
# 
# FD_Final_ITV <- rbind(Result_FD_ITV, Result_FD_ITV_noal) |>
#   mutate(Com = rep(c("Tot", "No_Alesb"), each = 26))

write.csv(FD_Final_ITV, "./output/FD_Final_ITV.csv")

# --------- END ------------
