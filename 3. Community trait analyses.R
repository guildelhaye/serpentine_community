#### 3. Community level trait variation

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
soil <- read.csv("./data/soil.csv") |> select(-X)
traits_sp_site <- read.csv("./data/traits_sp_site.csv") |> select(-X)
traits_sp <- read.csv("./data/traits_sp.csv") |> select(-X)
abun <- read.csv("./data/abundance.csv") |> select(-X)

## Create directory for each result
# dir.create("./output/trait_gradient_model")
# dir.create("./output/CWM")
# dir.create("./output/FD")

######  Test for trait - environment relationships -----
## 1 Using MLM 3 of ter Braak 2019  ----
# Assemble dataset for the model 
# Transform abundance in long format
abun_long <- abun |> 
  rename(Plot = site) |> 
  pivot_longer(where(is.numeric), names_to = "Sp", values_to = "abundance") |>
  separate("Plot", into = c("site"), remove = F) 

# Data for MLM3 model using only species level traits
data_MLM3 <- abun_long |> 
  mutate(abundance = ceiling(abundance)) |> # transform abundance in integers
  left_join(traits_sp) |>  
  left_join(soil |> select(site, Plot, Ni)) |> 
  mutate_at(c("LS","LeafArea", "LDMC", "SLA", "LT", "LNC", "Ni"), scale2)


## Run models 
mod_SLA <- brm(abundance ~ Ni * SLA  + (1 + SLA | site) + (1 + Ni | Sp), 
               family = poisson(link = "log"),
               control=list(adapt_delta = 0.99, 
                            max_treedepth = 13),
               chains = 4, cores = 4,
               data = data_MLM3, 
               file = "./output/MLM3/mod_SLA.rds")

mod_LT    <- update(mod_SLA, 
                    abundance ~ Ni * LT  + (1 + LT|site) + (1 + Ni|Sp), 
                    newdata = data_MLM3, 
                    file = "./output/trait_gradient_model/mod_LT.rds")

mod_LA    <- update(mod_SLA, 
                    abundance ~ Ni * LeafArea  + (1 + LeafArea|site) + (1 + Ni|Sp), 
                    newdata = data_MLM3, 
                    file = "./output/trait_gradient_model/mod_LA.rds" )

mod_LDMC  <- update(mod_SLA, 
                    abundance ~ Ni * LDMC + (1 + LDMC|site) + (1  +Ni|Sp), 
                    newdata = data_MLM3, 
                    file = "./output/trait_gradient_model/mod_LDMC.rds")

mod_LNC   <- update(mod_SLA, 
                    abundance ~ Ni * LNC + (1 + LNC|site) + (1 + Ni|Sp), 
                    newdata = data_MLM3, 
                    file = "./output/trait_gradient_model/mod_LNC.rds")

mod_LS    <- update(mod_SLA, 
                    abundance ~ Ni * LS + (1 + LS|site) + (1 + Ni|Sp), 
                    newdata = data_MLM3, 
                    file = "./output/trait_gradient_model/mod_LS.rds")

# Without O. lesbiaca
data_MLM3_no_OL <- data_MLM3 |> filter(Sp != "ALYLESB")

mod_SLA_no_OL <- update(mod_SLA, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_SLA_no_OL.rds")
mod_LA_no_OL <- update(mod_LA, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_LA_no_OL.rds")
mod_LT_no_OL <- update(mod_LT, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_LT_no_OL.rds")
mod_LDMC_no_OL <- update(mod_LDMC, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_LDMC_no_OL.rds")
mod_LNC_no_OL <- update(mod_LNC, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_LNC_no_OL.rds")
mod_LS_no_OL <- update(mod_LS, newdata = data_MLM3_no_OL, file = "./output/trait_gradient_model/mod_LS_no_OL.rds")

mod_SLA <- readRDS("./output/trait_gradient_model/mod_SLA.rds")
mod_LT    <- readRDS("./output/trait_gradient_model/mod_LT.rds")
mod_LA    <- readRDS("./output/trait_gradient_model/mod_LA.rds" )
mod_LDMC  <- readRDS("./output/trait_gradient_model/mod_LDMC.rds")
mod_LNC   <- readRDS("./output/trait_gradient_model/mod_LNC.rds")
mod_LS    <- readRDS("./output/trait_gradient_model/mod_LS.rds")
mod_SLA_no_OL   <- readRDS("./output/trait_gradient_model/mod_SLA_no_OL.rds")
mod_LA_no_OL    <- readRDS("./output/trait_gradient_model/mod_LA_no_OL.rds")
mod_LT_no_OL    <- readRDS("./output/trait_gradient_model/mod_LT_no_OL.rds")
mod_LDMC_no_OL  <- readRDS("./output/trait_gradient_model/mod_LDMC_no_OL.rds")
mod_LNC_no_OL   <- readRDS("./output/trait_gradient_model/mod_LNC_no_OL.rds")
mod_LS_no_OL    <- readRDS("./output/trait_gradient_model/mod_LS_no_OL.rds")

## Extract posteriors for figure
mod_draws_mlm3 <- rbind(
  mlt = m_ex(mod_LT, "LT"),
  mla = m_ex(mod_LA, "LA"),
  msla = m_ex(mod_SLA, "SLA"),
  mldmc = m_ex(mod_LDMC, "LDMC"),
  mlnc = m_ex(mod_LNC, "LNC"),
  mls <- m_ex(mod_LS, "LS"))
mod_draws_mlm3$trait <- factor(mod_draws_mlm3$trait, levels = c("LS", "LNC","LDMC","LT", "SLA","LA"))

## figure
fig_fdis <- mod_draws_mlm3 |>
  ggplot(aes(y = trait, x = var, fill = after_stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme(legend.position = "")

## calculate proportions of posterieur >< 0
mod_draws_mlm3 |> 
  group_by(trait) |> 
  count(var > 0) |> 
  mutate(prop = n/4000) |> 
  select(-n) |>
  pivot_wider(names_from = `var > 0`, values_from = prop) |>
  select(trait, 'p > 0' = `TRUE`, 'p < 0' = `FALSE`)

## Extract posteriors for figure without O. lesbiaca
mod_draws_mlm3_OL <- rbind(
  mlt = m_ex(mod_LT_no_OL, "LT"),
  mla = m_ex(mod_LA_no_OL, "LA"),
  msla = m_ex(mod_SLA_no_OL, "SLA"),
  mldmc = m_ex(mod_LDMC_no_OL, "LDMC"),
  mlnc = m_ex(mod_LNC_no_OL, "LNC"),
  mls <- m_ex(mod_LS_no_OL, "LS"))
mod_draws_mlm3_OL$trait <- factor(mod_draws_mlm3_OL$trait, levels = c("LS", "LNC","LDMC","LT", "SLA","LA"))

## figure
fig_fdis_OL <- mod_draws_mlm3_OL |>
  ggplot(aes(y = trait, x = var, fill = after_stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme(legend.position = "")

## calculate proportions of posterieur >< 0
mod_draws_mlm3_OL |> 
  group_by(trait) |> 
  count(var >0) |> 
  mutate(prop = n/12000) |> 
  select(-n) |>
  pivot_wider(names_from = `var > 0`, values_from = prop) |>
  select(trait, 'p > 0' = `TRUE`, 'p < 0' = `FALSE`)

ggpubr::ggarrange(fig_fdis, fig_fdis_OL, ncol = 2)
# ggsave("./output/figure/MLM3.pdf", height = 4, width = 7)

####################################################################################################
############# 2. Community weighted mean and FDis change along the gradient -----
#### CWM -----
# Data preparation 
ab_long <- abun |> 
  rename(id= site) |> 
  pivot_longer(where(is.numeric), names_to = "Sp", values_to = "abundance") |>
  separate("id", into = c("site"), remove = F) 

##  2.1 CWM
# Calculate CWM   
CWM_sp <- cwm(ab_long, traits_sp) # with species median
CWM_sp_site <- cwm(ab_long, traits_sp_site) # with site specific trait values

# Without O. lesbiaca
ab_long_no_OL <- ab_long |> filter(Sp!= "ALYLESB")
traits_sp_no_OL <- traits_sp |> filter(Sp!= "ALYLESB")
traits_sp_site_no_OL <- traits_sp_site |> filter(Sp!= "ALYLESB")

CWM_sp_no_OL <- cwm(ab_long_no_OL, traits_sp_no_OL) # with species median
CWM_sp_site_no_OL <- cwm(ab_long_no_OL, traits_sp_site_no_OL) # with site specific trait values

write.csv(CWM_sp, "./output/CWM/CWM_sp.csv")
write.csv(CWM_sp_site, "./output/CWM/CWM_sp_site.csv")
write.csv(CWM_sp_no_OL, "./output/CWM/CWM_sp_OL.csv")
write.csv(CWM_sp_site_no_OL, "./output/CWM/CWM_sp_site_OL.csv")

# check for association between CWM intra and inter
cor(CWM_sp$CWM_LA,  CWM_sp_site$CWM_LA)
cor(CWM_sp$CWM_SLA,  CWM_sp_site$CWM_SLA)
cor(CWM_sp$CWM_LDMC,  CWM_sp_site$CWM_LDMC)
cor(CWM_sp$CWM_LT,  CWM_sp_site$CWM_LT)
cor(CWM_sp$CWM_LNC,  CWM_sp_site$CWM_LNC)
cor(CWM_sp$CWM_LS,  CWM_sp_site$CWM_LS)

### Model the response of CWM traits along the Ni gradient for the 4 types of CWM (
### (with and without ITV and with and without O. lesbiaca)
# transform dataset
CWM_sp <- CWM_sp |> 
  left_join(soil |> select(id = Plot, Ni)) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) |>
  separate("id", c("site", "no"))

CWM_sp_site <- CWM_sp_site |> 
  left_join(soil |> select(id = Plot, Ni)) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) |>
  separate("id", c("site", "no"))

CWM_sp_no_OL <- CWM_sp_no_OL |> 
  left_join(soil |> select(id = Plot, Ni)) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) |>
  separate("id", c("site", "no"))

CWM_sp_site_no_OL <- CWM_sp_site_no_OL |> 
  left_join(soil |> select(id = Plot, Ni)) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) |>
  separate("id", c("site", "no"))

## Leaf Area
mod_cwm_la <- brm(CWM_LA ~ Ni + (1 | site), 
                  data = CWM_sp, 
                  prior = c(prior(normal(0, 2), class = b, coef = Ni),
                            prior(normal(0, 2), class = Intercept)),
                  control = list(adapt_delta = 0.99), 
                  file = "./output/CWM/mod_cwm_la.rds")
c_cwm_la = summary(mod_cwm_la$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_la_itv <- update(mod_cwm_la, newdata = CWM_sp_site, 
                         file = "./output/CWM/mod_cwm_la_itv.rds")
c_cwm_la_itv = summary(mod_cwm_la_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_la_OL <- update(mod_cwm_la, newdata = CWM_sp_no_OL, 
                        file = "./output/CWM/mod_cwm_la_ol.rds")
c_cwm_la_OL = summary(mod_cwm_la_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_la_itv_OL <- update(mod_cwm_la, newdata = CWM_sp_site_no_OL, 
                            file = "./output/CWM/mod_cwm_la_itv_ol.rds")
c_cwm_la_itv_OL = summary(mod_cwm_la_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## SLA
mod_cwm_sla <- brm(CWM_SLA ~ Ni + (1 | site), 
                   data = CWM_sp, 
                   prior = c(prior(normal(0, 2), class = b, coef = Ni),
                             prior(normal(0, 2), class = Intercept)),
                   control = list(adapt_delta = 0.99), 
                   file = "./output/CWM/mod_cwm_sla.rds")
c_cwm_sla = summary(mod_cwm_sla$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_sla_itv <- update(mod_cwm_sla, newdata = CWM_sp_site, 
                          file = "./output/CWM/mod_cwm_sla_itv.rds")
c_cwm_sla_itv = summary(mod_cwm_sla_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_sla_OL <- update(mod_cwm_sla, newdata = CWM_sp_no_OL, 
                         file = "./output/CWM/mod_cwm_sla_ol.rds")
c_cwm_sla_OL = summary(mod_cwm_sla_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_sla_itv_OL <- update(mod_cwm_sla, newdata = CWM_sp_site_no_OL, 
                             file = "./output/CWM/mod_cwm_sla_itv_ol.rds")
c_cwm_sla_itv_OL = summary(mod_cwm_sla_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LDMC
mod_cwm_ldmc <- brm(CWM_LDMC ~ Ni + (1 | site), 
                    data = CWM_sp, 
                    prior = c(prior(normal(0, 2), class = b, coef = Ni),
                              prior(normal(0, 2), class = Intercept)),
                    control = list(adapt_delta = 0.99), 
                    file = "./output/CWM/mod_cwm_ldmc.rds")
c_cwm_ldmc = summary(mod_cwm_ldmc$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ldmc_itv <- update(mod_cwm_ldmc, newdata = CWM_sp_site, 
                           file = "./output/CWM/mod_cwm_ldmc_itv.rds")
c_cwm_ldmc_itv = summary(mod_cwm_ldmc_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ldmc_OL <- update(mod_cwm_ldmc, newdata = CWM_sp_no_OL, 
                          file = "./output/CWM/mod_cwm_ldmc_ol.rds")
c_cwm_ldmc_OL = summary(mod_cwm_ldmc_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ldmc_itv_OL <- update(mod_cwm_ldmc, newdata = CWM_sp_site_no_OL, 
                              file = "./output/CWM/mod_cwm_ldmc_itv_ol.rds")
c_cwm_ldmc_itv_OL = summary(mod_cwm_ldmc_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LT
mod_cwm_lt <- brm(CWM_LT ~ Ni + (1 | site), 
                  data = CWM_sp, 
                  prior = c(prior(normal(0, 2), class = b, coef = Ni),
                            prior(normal(0, 2), class = Intercept)),
                  control = list(adapt_delta = 0.99), 
                  file = "./output/CWM/mod_cwm_lt.rds")
c_cwm_lt = summary(mod_cwm_lt$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lt_itv <- update(mod_cwm_lt, newdata = CWM_sp_site, 
                         file = "./output/CWM/mod_cwm_lt_itv.rds")
c_cwm_lt_itv = summary(mod_cwm_lt_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lt_OL <- update(mod_cwm_lt, newdata = CWM_sp_no_OL, 
                        file = "./output/CWM/mod_cwm_lt_ol.rds")
c_cwm_lt_OL = summary(mod_cwm_lt_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lt_itv_OL <- update(mod_cwm_lt, newdata = CWM_sp_site_no_OL, 
                            file = "./output/CWM/mod_cwm_lt_itv_ol.rds")
c_cwm_lt_itv_OL = summary(mod_cwm_lt_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LNC
mod_cwm_lnc <- brm(CWM_LNC ~ Ni + (1 | site), 
                   data = CWM_sp, 
                   prior = c(prior(normal(0, 2), class = b, coef = Ni),
                             prior(normal(0, 2), class = Intercept)),
                   control = list(adapt_delta = 0.99), 
                   file = "./output/CWM/mod_cwm_lnc.rds")
c_cwm_lnc = summary(mod_cwm_lnc$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lnc_itv <- update(mod_cwm_lnc, newdata = CWM_sp_site, 
                          file = "./output/CWM/mod_cwm_lnc_itv.rds")
c_cwm_lnc_itv = summary(mod_cwm_lnc_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lnc_OL <- update(mod_cwm_lnc, newdata = CWM_sp_no_OL, 
                         file = "./output/CWM/mod_cwm_lnc_ol.rds")
c_cwm_lnc_OL = summary(mod_cwm_lnc_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_lnc_itv_OL <- update(mod_cwm_lnc, newdata = CWM_sp_site_no_OL, 
                             file = "./output/CWM/mod_cwm_lnc_itv_ol.rds")
c_cwm_lnc_itv_OL = summary(mod_cwm_lnc_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LS
mod_cwm_ls <- brm(CWM_LS ~ Ni + (1 | site), 
                  data = CWM_sp, 
                  prior = c(prior(normal(0, 2), class = b, coef = Ni),
                            prior(normal(0, 2), class = Intercept)),
                  control = list(adapt_delta = 0.99), 
                  file = "./output/CWM/mod_cwm_ls.rds")
c_cwm_ls = summary(mod_cwm_ls$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ls_itv <- update(mod_cwm_ls, newdata = CWM_sp_site, 
                         file = "./output/CWM/mod_cwm_ls_itv.rds")
c_cwm_ls_itv = summary(mod_cwm_ls_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ls_OL <- update(mod_cwm_ls, newdata = CWM_sp_no_OL, 
                        file = "./output/CWM/mod_cwm_ls_ol.rds")
c_cwm_ls_OL = summary(mod_cwm_ls_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_cwm_ls_itv_OL <- update(mod_cwm_ls, newdata = CWM_sp_site_no_OL, 
                            file = "./output/CWM/mod_cwm_ls_itv_ol.rds")
c_cwm_ls_itv_OL = summary(mod_cwm_ls_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

# table
CWM_table <- data.frame(
  trait = rep(c("LA",  "SLA", "LDMC",  "LT",  "LNC",  "LS"), each = 4),
  O.lesb = rep(c("present", "absent"), each = 2),
  itv = rep(c("yes", "no")),
  intercept = c(c_cwm_la[1, 1], c_cwm_la_itv[1, 1],c_cwm_la_OL[1, 1],c_cwm_la_itv_OL[1, 1],
                c_cwm_sla[1, 1], c_cwm_sla_itv[1, 1],c_cwm_sla_OL[1, 1],c_cwm_sla_itv_OL[1, 1],
                c_cwm_ldmc[1, 1], c_cwm_ldmc_itv[1, 1],c_cwm_ldmc_OL[1, 1],c_cwm_ldmc_itv_OL[1, 1],
                c_cwm_lt[1, 1], c_cwm_lt_itv[1, 1],c_cwm_lt_OL[1, 1],c_cwm_lt_itv_OL[1, 1],
                c_cwm_lnc[1, 1], c_cwm_lnc_itv[1, 1],c_cwm_lnc_OL[1, 1],c_cwm_lnc_itv_OL[1, 1],
                c_cwm_ls[1, 1], c_cwm_ls_itv[1, 1], c_cwm_ls_OL[1, 1], c_cwm_ls_itv_OL[1, 1]),
  slope_mean = c(c_cwm_la[2, 1], c_cwm_la_itv[2, 1],c_cwm_la_OL[2, 1],c_cwm_la_itv_OL[2, 1],
                 c_cwm_sla[2, 1], c_cwm_sla_itv[2, 1],c_cwm_sla_OL[2, 1],c_cwm_sla_itv_OL[2, 1],
                 c_cwm_ldmc[2, 1], c_cwm_ldmc_itv[2, 1],c_cwm_ldmc_OL[2, 1],c_cwm_ldmc_itv_OL[2, 1],
                 c_cwm_lt[2, 1], c_cwm_lt_itv[2, 1],c_cwm_lt_OL[2, 1],c_cwm_lt_itv_OL[2, 1],
                 c_cwm_lnc[2, 1], c_cwm_lnc_itv[2, 1],c_cwm_lnc_OL[2, 1],c_cwm_lnc_itv_OL[2, 1],
                 c_cwm_ls[2, 1], c_cwm_ls_itv[2, 1], c_cwm_ls_OL[2, 1], c_cwm_ls_itv_OL[2, 1]),
  slope_2.5 = c(c_cwm_la[2, 2], c_cwm_la_itv[2, 2],c_cwm_la_OL[2, 2],c_cwm_la_itv_OL[2, 2],
                c_cwm_sla[2, 2], c_cwm_sla_itv[2, 2],c_cwm_sla_OL[2, 2],c_cwm_sla_itv_OL[2, 2],
                c_cwm_ldmc[2, 2], c_cwm_ldmc_itv[2, 2],c_cwm_ldmc_OL[2, 2],c_cwm_ldmc_itv_OL[2, 2],
                c_cwm_lt[2, 2], c_cwm_lt_itv[2, 2],c_cwm_lt_OL[2, 2],c_cwm_lt_itv_OL[2, 2],
                c_cwm_lnc[2, 2], c_cwm_lnc_itv[2, 2],c_cwm_lnc_OL[2, 2],c_cwm_lnc_itv_OL[2, 2],
                c_cwm_ls[2, 2], c_cwm_ls_itv[2, 2], c_cwm_ls_OL[2, 2], c_cwm_ls_itv_OL[2, 2]),
  slope_97.5 = c(c_cwm_la[2, 3], c_cwm_la_itv[2, 3],c_cwm_la_OL[2, 3],c_cwm_la_itv_OL[2, 3],
                 c_cwm_sla[2, 3], c_cwm_sla_itv[2, 3],c_cwm_sla_OL[2, 3],c_cwm_sla_itv_OL[2, 3],
                 c_cwm_ldmc[2, 3], c_cwm_ldmc_itv[2, 3],c_cwm_ldmc_OL[2, 3],c_cwm_ldmc_itv_OL[2, 3],
                 c_cwm_lt[2, 3], c_cwm_lt_itv[2, 3],c_cwm_lt_OL[2, 3],c_cwm_lt_itv_OL[2, 3],
                 c_cwm_lnc[2, 3], c_cwm_lnc_itv[2, 3],c_cwm_lnc_OL[2, 3],c_cwm_lnc_itv_OL[2, 3],
                 c_cwm_ls[2, 3], c_cwm_ls_itv[2, 3], c_cwm_ls_OL[2, 3], c_cwm_ls_itv_OL[2, 3]))

# write.csv(CWM_table, "./output/CWM/cwm_reg_table.csv")


# Graph
mod_draws_cwm <- rbind(
  mlt = m_ex2(mod_cwm_la_itv, "LA"),
  mla = m_ex2(mod_cwm_sla_itv, "SLA"),
  msla = m_ex2(mod_cwm_lt_itv, "LT"),
  mldmc = m_ex2(mod_cwm_ldmc_itv, "LDMC"),
  mlnc = m_ex2(mod_cwm_lnc_itv, "LNC"),
  mls <- m_ex2(mod_cwm_ls_itv, "LS"))
mod_draws_cwm$trait <- factor(mod_draws_cwm$trait, levels = c("LS", "LNC","LDMC","LT", "SLA","LA"))

graph_cwm <- mod_draws_cwm |>
  ggplot(aes(y = trait, x = var, fill = after_stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme(legend.position = "")  +
  xlab("Slope") +
  ggtitle("A) ")

## calculate proportions of posterieur >< 0
mod_draws_cwm |> 
  group_by(trait) |> 
  count(var > 0) |> 
  mutate(prop = n/4000) |> 
  select(-n) |>
  pivot_wider(names_from = `var > 0`, values_from = prop) |>
  select(trait, 'p > 0' = `TRUE`, 'p < 0' = `FALSE`)

# without O. lesbiaca
mod_draws_cwm_noal <- rbind(
  mlt = m_ex2(mod_cwm_la_itv_OL, "LA"),
  mla = m_ex2(mod_cwm_sla_itv_OL, "SLA"),
  msla = m_ex2(mod_cwm_lt_itv_OL, "LT"),
  mldmc = m_ex2(mod_cwm_ldmc_itv_OL, "LDMC"),
  mlnc = m_ex2(mod_cwm_lnc_itv_OL, "LNC"),
  mls <- m_ex2(mod_cwm_ls_itv_OL, "LS"))
mod_draws_cwm_noal$trait <- factor(mod_draws_cwm$trait, levels = c("LS", "LNC","LDMC","LT", "SLA","LA"))

graph_cwm_noal <- mod_draws_cwm_noal |>
  ggplot(aes(y = trait, x = var, fill = after_stat(x > 0))) +
  stat_halfeye() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme(legend.position = "") +
  xlab("Slope") +
  ggtitle("B)")

## calculate proportions of posterieur >< 0
mod_draws_cwm_noal |> 
  group_by(trait) |> 
  count(var > 0) |> 
  mutate(prop = n/4000) |> 
  select(-n) |>
  pivot_wider(names_from = `var > 0`, values_from = prop) |>
  select(trait, 'p > 0' = `TRUE`, 'p < 0' = `FALSE`)

ggpubr::ggarrange(graph_cwm, graph_cwm_noal, ncol = 2)
#ggsave("./output/figure/cwm.pdf", height = 4, width = 7)

## 2.2 FDis
### FDis --------
nrep = 500

## Multivariate FDis
fdism_amp <- sesFDism(traits = traits_sp_site, si = "Amp", comm = abun, numberReps = nrep)
fdism_oly <- sesFDism(traits = traits_sp_site, si = "Oly", comm = abun, numberReps = nrep)
fdism_vat <- sesFDism(traits = traits_sp_site, si = "Vat", comm = abun, numberReps = nrep)
fdism_lout<- sesFDism(traits = traits_sp_site, si = "Lout", comm = abun, numberReps = nrep)
fdism_itv <- c(fdism_amp, fdism_lout, fdism_oly, fdism_vat)

## Univariate FDis
fdisu_amp_LS <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp", i = 'LS', numberReps = nrep)
fdisu_oly_LS <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly", i = 'LS', numberReps = nrep)
fdisu_vat_LS <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat", i = 'LS', numberReps = nrep)
fdisu_lout_LS<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'LS', numberReps = nrep)
FDis.LS = c(fdisu_amp_LS, fdisu_lout_LS,fdisu_oly_LS, fdisu_vat_LS)

fdisu_amp_LA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp",   i = 'LeafArea', numberReps = nrep)
fdisu_oly_LA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly",  i = 'LeafArea', numberReps = nrep)
fdisu_vat_LA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat",  i = 'LeafArea', numberReps = nrep)
fdisu_lout_LA<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'LeafArea', numberReps = nrep)
FDis.LeafArea = c(fdisu_amp_LA, fdisu_lout_LA,fdisu_oly_LA, fdisu_vat_LA)

fdisu_amp_SLA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp",  i = 'SLA', numberReps = nrep)
fdisu_oly_SLA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly",  i = 'SLA', numberReps = nrep)
fdisu_vat_SLA <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat",  i = 'SLA', numberReps = nrep)
fdisu_lout_SLA<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'SLA', numberReps = nrep)
FDis.SLA = c(fdisu_amp_SLA, fdisu_lout_SLA,fdisu_oly_SLA, fdisu_vat_SLA)

fdisu_amp_LT <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp",  i = 'LT', numberReps = nrep)
fdisu_oly_LT <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly",  i = 'LT', numberReps = nrep)
fdisu_vat_LT <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat",  i = 'LT', numberReps = nrep)
fdisu_lout_LT<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'LT', numberReps = nrep)
FDis.LT =  c(fdisu_amp_LT, fdisu_lout_LT,fdisu_oly_LT, fdisu_vat_LT)

fdisu_amp_LDMC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp",  i = 'LDMC', numberReps = nrep)
fdisu_oly_LDMC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly",  i = 'LDMC', numberReps = nrep)
fdisu_vat_LDMC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat",  i = 'LDMC', numberReps = nrep)
fdisu_lout_LDMC<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'LDMC', numberReps = nrep)
FDis.LDMC = c(fdisu_amp_LDMC, fdisu_lout_LDMC,fdisu_oly_LDMC, fdisu_vat_LDMC)

fdisu_amp_LNC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Amp",  i = 'LNC', numberReps = nrep)
fdisu_oly_LNC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Oly",  i = 'LNC', numberReps = nrep)
fdisu_vat_LNC <- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Vat",  i = 'LNC', numberReps = nrep)
fdisu_lout_LNC<- sesFDisu(traits = traits_sp_site,  comm = abun,  si = "Lout", i = 'LNC', numberReps = nrep)
FDis.LNC =  c(fdisu_amp_LNC, fdisu_lout_LNC,fdisu_oly_LNC, fdisu_vat_LNC)

fdis_ITV <- data.frame(Plot = abun$site, multi = fdism_itv, LA = FDis.LeafArea, SLA = FDis.SLA, 
                       LDMC = FDis.LDMC, LT = FDis.LT, LNC = FDis.LNC, LS = FDis.LS)

## FDis using average values 
fdism <- sesFDism(traits = traits_sp, si = "all", comm = abun, numberReps = nrep)

fdisu_LNC       <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'LNC', numberReps = nrep)
fdisu_med_LS    <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'LS', numberReps = nrep)
fdisu_med_LA    <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'LeafArea', numberReps = nrep)
fdisu_med_SLA   <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'SLA', numberReps = nrep)
fdisu_med_LDMC  <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'LDMC', numberReps = nrep)
fdisu_med_LT    <- sesFDisu(traits = traits_sp,  comm = abun, si = "all",  i = 'LT', numberReps = nrep)

fdis_median <- data.frame(Plot = abun$site, multi = fdism, LA = fdisu_med_LA, SLA = fdisu_med_SLA, 
                    LDMC = fdisu_med_LDMC, LT = fdisu_med_LT, LNC = fdisu_LNC, LS = fdisu_med_LS)

## Without O. lesbiaca --
abun_OL <- abun |> select(-"ALYLESB")
traits_sp_OL <- traits_sp |> filter(Sp != "ALYLESB")
traits_sp_site_OL <- traits_sp_site |> filter(Sp != "ALYLESB")

## Multivariate FDis
fdism_amp_OL <- sesFDism(traits = traits_sp_site_OL, si = "Amp", comm = abun_OL, numberReps = nrep)
fdism_oly_OL <- sesFDism(traits = traits_sp_site_OL, si = "Oly", comm = abun_OL, numberReps = nrep)
fdism_vat_OL <- sesFDism(traits = traits_sp_site_OL, si = "Vat", comm = abun_OL, numberReps = nrep)
fdism_lout_OL<- sesFDism(traits = traits_sp_site_OL, si = "Lout", comm = abun_OL, numberReps = nrep)
fdism_itv_OL <- c(fdism_amp_OL, fdism_lout_OL, fdism_oly_OL, fdism_vat_OL)

## Univariate FDis
fdisu_amp_LS_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp", i = 'LS', numberReps = nrep)
fdisu_oly_LS_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly", i = 'LS', numberReps = nrep)
fdisu_vat_LS_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat", i = 'LS', numberReps = nrep)
fdisu_lout_LS_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'LS', numberReps = nrep)
FDis.LS_OL = c(fdisu_amp_LS_OL, fdisu_lout_LS_OL,fdisu_oly_LS_OL, fdisu_vat_LS_OL)

fdisu_amp_LA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp",   i = 'LeafArea', numberReps = nrep)
fdisu_oly_LA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly",  i = 'LeafArea', numberReps = nrep)
fdisu_vat_LA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat",  i = 'LeafArea', numberReps = nrep)
fdisu_lout_LA_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'LeafArea', numberReps = nrep)
FDis.LeafArea_OL = c(fdisu_amp_LA_OL, fdisu_lout_LA_OL,fdisu_oly_LA_OL, fdisu_vat_LA_OL)

fdisu_amp_SLA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp",  i = 'SLA', numberReps = nrep)
fdisu_oly_SLA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly",  i = 'SLA', numberReps = nrep)
fdisu_vat_SLA_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat",  i = 'SLA', numberReps = nrep)
fdisu_lout_SLA_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'SLA', numberReps = nrep)
FDis.SLA_OL = c(fdisu_amp_SLA_OL, fdisu_lout_SLA_OL,fdisu_oly_SLA_OL, fdisu_vat_SLA_OL)

fdisu_amp_LT_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp",  i = 'LT', numberReps = nrep)
fdisu_oly_LT_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly",  i = 'LT', numberReps = nrep)
fdisu_vat_LT_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat",  i = 'LT', numberReps = nrep)
fdisu_lout_LT_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'LT', numberReps = nrep)
FDis.LT_OL =  c(fdisu_amp_LT_OL, fdisu_lout_LT_OL,fdisu_oly_LT_OL, fdisu_vat_LT_OL)

fdisu_amp_LDMC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp",  i = 'LDMC', numberReps = nrep)
fdisu_oly_LDMC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly",  i = 'LDMC', numberReps = nrep)
fdisu_vat_LDMC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat",  i = 'LDMC', numberReps = nrep)
fdisu_lout_LDMC_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'LDMC', numberReps = nrep)
FDis.LDMC_OL = c(fdisu_amp_LDMC_OL, fdisu_lout_LDMC_OL,fdisu_oly_LDMC_OL, fdisu_vat_LDMC_OL)

fdisu_amp_LNC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Amp",  i = 'LNC', numberReps = nrep)
fdisu_oly_LNC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Oly",  i = 'LNC', numberReps = nrep)
fdisu_vat_LNC_OL <- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Vat",  i = 'LNC', numberReps = nrep)
fdisu_lout_LNC_OL<- sesFDisu(traits = traits_sp_site_OL,  comm = abun_OL,  si = "Lout", i = 'LNC', numberReps = nrep)
FDis.LNC_OL =  c(fdisu_amp_LNC_OL, fdisu_lout_LNC_OL,fdisu_oly_LNC_OL, fdisu_vat_LNC_OL)

fdis_ITV_OL <- data.frame(Plot = abun$site, multi = fdism_itv_OL, LA = FDis.LeafArea_OL, SLA = FDis.SLA_OL, 
                          LDMC = FDis.LDMC_OL, LT = FDis.LT_OL, LNC = FDis.LNC_OL, LS = FDis.LS_OL)

## FDis using average values 
fdism_OL <- sesFDism(traits = traits_sp_OL, si = "all", comm = abun_OL, numberReps = nrep)

fdisu_LNC_OL       <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'LNC', numberReps = nrep)
fdisu_med_LS_OL    <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'LS', numberReps = nrep)
fdisu_med_LA_OL    <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'LeafArea', numberReps = nrep)
fdisu_med_SLA_OL   <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'SLA', numberReps = nrep)
fdisu_med_LDMC_OL  <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'LDMC', numberReps = nrep)
fdisu_med_LT_OL    <- sesFDisu(traits = traits_sp_OL,  comm = abun_OL, si = "all",  i = 'LT', numberReps = nrep)

fdis_median_OL <- data.frame(Plot = abun$site, multi = fdism_OL, LA = fdisu_med_LA_OL, SLA = fdisu_med_SLA_OL, 
                             LDMC = fdisu_med_LDMC_OL, LT = fdisu_med_LT_OL, LNC = fdisu_LNC_OL, LS = fdisu_med_LS_OL)

write.csv(fdis_median, "./output/FD/fdis_median.csv")
write.csv(fdis_ITV, "./output/FD/fdis_ITV.csv")
write.csv(fdis_median_OL, "./output/FD/fdis_median_OL.csv")
write.csv(fdis_ITV_OL, "./output/FD/fdis_ITV_OL.csv")

### Transform data for model 
data_fdis_ITV <- fdis_ITV |> 
  left_join(soil |> select(Plot, Ni)) |>
  separate(Plot, into = c("site", "num")) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) 

data_fdis <- fdis_median |> 
  left_join(soil |> select(Plot, Ni)) |>
  separate(Plot, into = c("site", "num")) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector))

data_fdis_ITV_OL <- fdis_ITV_OL |> 
  left_join(soil |> select(Plot, Ni)) |>
  separate(Plot, into = c("site", "num")) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector)) 

data_fdis_OL <- fdis_median_OL |> 
  left_join(soil |> select(Plot, Ni)) |>
  separate(Plot, into = c("site", "num")) |>
  mutate_if(is.numeric, ~(scale(.) %>% as.vector))

#### Model of change in FDis along the Ni gradient

## FDis total
mod_fdis <- brm(multi ~ Ni + (1 | site), 
                data = data_fdis, 
                prior = c(prior(normal(0, 1), class = b, coef = Ni),
                          prior(normal(0, 1), class = Intercept)),
                control = list(adapt_delta = 0.99), 
                file = "./output/FD/mod_fdis_multi.rds")
c_fdis = summary(mod_fdis$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_itv <- update(mod_fdis, newdata = data_fdis_ITV, 
                       file = "./output/FD/mod_fdis_multi_itv.rds")
c_fdis_itv = summary(mod_fdis_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_OL <- update(mod_fdis, newdata = data_fdis_OL, 
                      file = "./output/FD/mod_fdis_multi_ol.rds")
c_fdis_OL = summary(mod_fdis_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_itv_OL <- update(mod_fdis, newdata = data_fdis_ITV_OL, 
                          file = "./output/FD/mod_fdis_multi_itv_ol.rds")
c_fdis_itv_OL = summary(mod_fdis_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## Leaf Area
mod_fdis_la <- brm(LA ~ Ni + (1 | site), 
                    data = data_fdis, 
                    prior = c(prior(normal(0, 1), class = b, coef = Ni),
                              prior(normal(0, 1), class = Intercept)),
                    control = list(adapt_delta = 0.99), 
                   file = "./output/FD/mod_fdis_la.rds")
c_fdis_la = summary(mod_fdis_la$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_la_itv <- update(mod_fdis_la, newdata = data_fdis_ITV, 
                          file = "./output/FD/mod_fdis_la_itv.rds")
c_fdis_la_itv = summary(mod_fdis_la_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_la_OL <- update(mod_fdis_la, newdata = data_fdis_OL, 
                         file = "./output/FD/mod_fdis_la_ol.rds")
c_fdis_la_OL = summary(mod_fdis_la_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_la_itv_OL <- update(mod_fdis_la, newdata = data_fdis_ITV_OL, 
                             file = "./output/FD/mod_fdis_la_itv_ol.rds")
c_fdis_itv_la_OL = summary(mod_fdis_la_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## SLA
mod_fdis_sla <- brm(SLA ~ Ni + (1 | site), 
                       data = data_fdis, 
                       prior = c(prior(normal(0, 1), class = b, coef = Ni),
                                 prior(normal(0, 1), class = Intercept)),
                       control = list(adapt_delta = 0.99), 
                    file = "./output/FD/mod_fdis_sla.rds")
c_fdis_sla = summary(mod_fdis_sla$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_sla_itv <- update(mod_fdis_sla, newdata = data_fdis_ITV, 
                           file = "./output/FD/mod_fdis_sla_itv.rds")
c_fdis_sla_itv = summary(mod_fdis_sla_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_sla_OL <- update(mod_fdis_sla, newdata = data_fdis_OL, 
                          file = "./output/FD/mod_fdis_sla_ol.rds")
c_fdis_sla_OL = summary(mod_fdis_sla_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_sla_itv_OL <- update(mod_fdis_sla, newdata = data_fdis_ITV_OL, 
                              file = "./output/FD/mod_fdis_sla_itv_ol.rds")
c_fdis_itv_sla_OL = summary(mod_fdis_sla_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LDMC
mod_fdis_ldmc <- brm(LDMC ~ Ni + (1 | site), 
                        data = data_fdis, 
                        prior = c(prior(normal(0, 1), class = b, coef = Ni),
                                  prior(normal(0, 1), class = Intercept)),
                        control = list(adapt_delta = 0.99), 
                     file = "./output/FD/mod_fdis_ldmc.rds")
c_fdis_ldmc = summary(mod_fdis_ldmc$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ldmc_itv <- update(mod_fdis_ldmc, newdata = data_fdis_ITV, 
                            file = "./output/FD/mod_fdis_ldmc_itv.rds")
c_fdis_ldmc_itv = summary(mod_fdis_ldmc_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ldmc_OL <- update(mod_fdis_ldmc, newdata = data_fdis_OL, 
                           file = "./output/FD/mod_fdis_ldmc_ol.rds")
c_fdis_ldmc_OL = summary(mod_fdis_ldmc_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ldmc_itv_OL <- update(mod_fdis_ldmc, newdata = data_fdis_ITV_OL, 
                               file = "./output/FD/mod_fdis_ldmc_itv_ol.rds")
c_fdis_itv_ldmc_OL = summary(mod_fdis_ldmc_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LT
mod_fdis_lt <- brm(LT ~ Ni + (1 | site), 
                         data = data_fdis, 
                         prior = c(prior(normal(0, 1), class = b, coef = Ni),
                                   prior(normal(0, 1), class = Intercept)),
                         control = list(adapt_delta = 0.99), 
                   file = "./output/FD/mod_fdis_lt.rds")
c_fdis_lt = summary(mod_fdis_lt$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lt_itv <- update(mod_fdis_lt, newdata = data_fdis_ITV, 
                          file = "./output/FD/mod_fdis_lt_itv.rds")
c_fdis_lt_itv = summary(mod_fdis_lt_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lt_OL <- update(mod_fdis_lt, newdata = data_fdis_OL, 
                         file = "./output/FD/mod_fdis_lt_ol.rds")
c_fdis_lt_OL = summary(mod_fdis_lt_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lt_itv_OL <- update(mod_fdis_lt, newdata = data_fdis_ITV_OL, 
                             file = "./output/FD/mod_fdis_lt_itv_ol.rds")
c_fdis_itv_lt_OL = summary(mod_fdis_lt_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LNC
mod_fdis_lnc <- brm(LNC ~ Ni + (1 | site), 
                       data = data_fdis, 
                       prior = c(prior(normal(0, 1), class = b, coef = Ni),
                                 prior(normal(0, 1), class = Intercept)),
                       control = list(adapt_delta = 0.99), 
                    file = "./output/FD/mod_fdis_lnc.rds")
c_fdis_lnc = summary(mod_fdis_lnc$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lnc_itv <- update(mod_fdis_lnc, newdata = data_fdis_ITV, 
                           file = "./output/FD/mod_fdis_lnc_itv.rds")
c_fdis_lnc_itv = summary(mod_fdis_lnc_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lnc_OL <- update(mod_fdis_lnc, newdata = data_fdis_OL, 
                          file = "./output/FD/mod_fdis_lnc_ol.rds")
c_fdis_lnc_OL = summary(mod_fdis_lnc_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_lnc_itv_OL <- update(mod_fdis_lnc, newdata = data_fdis_ITV_OL, 
                              file = "./output/FD/mod_fdis_lnc_itv_ol.rds")
c_fdis_itv_lnc_OL = summary(mod_fdis_lnc_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]

## LS
mod_fdis_ls <- brm(LS ~ Ni + (1 | site), 
                        data = data_fdis, 
                        prior = c(prior(normal(0, 1), class = b, coef = Ni),
                                  prior(normal(0, 1), class = Intercept)),
                        control = list(adapt_delta = 0.99), 
                   file = "./output/FD/mod_fdis_ls.rds")
c_fdis_ls = summary(mod_fdis_ls$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ls_itv <- update(mod_fdis_ls, newdata = data_fdis_ITV, 
                          file = "./output/FD/mod_fdis_ls_itv.rds")
c_fdis_ls_itv = summary(mod_fdis_ls_itv$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ls_OL <- update(mod_fdis_ls, newdata = data_fdis_OL, 
                         file = "./output/FD/mod_fdis_ls_ol.rds")
c_fdis_ls_OL = summary(mod_fdis_ls_OL$fit)$summary[c(1:2), c(1, 4, 8)]
#
mod_fdis_ls_itv_OL <- update(mod_fdis_ls, newdata = data_fdis_ITV_OL, 
                             file = "./output/FD/mod_fdis_ls_itv_ol.rds")
c_fdis_itv_ls_OL = summary(mod_fdis_ls_itv_OL$fit)$summary[c(1:2), c(1, 4, 8)]


# table
fdis_table <- data.frame(
  trait = rep(c("Multivariate", "LA",  "SLA", "LDMC",  "LT",  "LNC",  "LS"), each = 4),
  O.lesb = rep(c("present", "absent"), each = 2),
  itv = rep(c("yes", "no")),
  intercept = c(c_fdis[1, 1], c_fdis_itv[1, 1],c_fdis_OL[1, 1],c_fdis_itv_OL[1, 1],
                c_fdis_la[1, 1], c_fdis_la_itv[1, 1],c_fdis_la_OL[1, 1],c_fdis_itv_la_OL[1, 1],
                c_fdis_sla[1, 1], c_fdis_sla_itv[1, 1],c_fdis_sla_OL[1, 1],c_fdis_itv_sla_OL[1, 1],
                c_fdis_ldmc[1, 1], c_fdis_ldmc_itv[1, 1],c_fdis_ldmc_OL[1, 1],c_fdis_itv_ldmc_OL[1, 1],
                c_fdis_lt[1, 1], c_fdis_lt_itv[1, 1],c_fdis_lt_OL[1, 1],c_fdis_itv_lt_OL[1, 1],
                c_fdis_lnc[1, 1], c_fdis_lnc_itv[1, 1],c_fdis_lnc_OL[1, 1],c_fdis_itv_lnc_OL[1, 1],
                c_fdis_ls[1, 1], c_fdis_ls_itv[1, 1], c_fdis_ls_OL[1, 1], c_fdis_itv_ls_OL[1, 1]),
  slope_mean = c(c_fdis[2, 1], c_fdis_itv[2, 1],c_fdis_OL[2, 1],c_fdis_itv_OL[2, 1],
                 c_fdis_la[2, 1], c_fdis_la_itv[2, 1],c_fdis_la_OL[2, 1],c_fdis_itv_la_OL[2, 1],
                 c_fdis_sla[2, 1], c_fdis_sla_itv[2, 1],c_fdis_sla_OL[2, 1],c_fdis_itv_sla_OL[2, 1],
                 c_fdis_ldmc[2, 1], c_fdis_ldmc_itv[2, 1],c_fdis_ldmc_OL[2, 1],c_fdis_itv_ldmc_OL[2, 1],
                 c_fdis_lt[2, 1], c_fdis_lt_itv[2, 1],c_fdis_lt_OL[2, 1],c_fdis_itv_lt_OL[2, 1],
                 c_fdis_lnc[2, 1], c_fdis_lnc_itv[2, 1],c_fdis_lnc_OL[2, 1],c_fdis_itv_lnc_OL[2, 1],
                 c_fdis_ls[2, 1], c_fdis_ls_itv[2, 1], c_fdis_ls_OL[2, 1], c_fdis_itv_ls_OL[2, 1]),
  slope_2.5 = c(c_fdis[2, 2], c_fdis_itv[2, 2],c_fdis_OL[2, 2],c_fdis_itv_OL[2, 2],
                c_fdis_la[2, 2], c_fdis_la_itv[2, 2],c_fdis_la_OL[2, 2],c_fdis_itv_la_OL[2, 2],
                c_fdis_sla[2, 2], c_fdis_sla_itv[2, 2],c_fdis_sla_OL[2, 2],c_fdis_itv_sla_OL[2, 2],
                c_fdis_ldmc[2, 2], c_fdis_ldmc_itv[2, 2],c_fdis_ldmc_OL[2, 2],c_fdis_itv_ldmc_OL[2, 2],
                c_fdis_lt[2, 2], c_fdis_lt_itv[2, 2],c_fdis_lt_OL[2, 2],c_fdis_itv_lt_OL[2, 2],
                c_fdis_lnc[2, 2], c_fdis_lnc_itv[2, 2],c_fdis_lnc_OL[2, 2],c_fdis_itv_lnc_OL[2, 2],
                c_fdis_ls[2, 2], c_fdis_ls_itv[2, 2], c_fdis_ls_OL[2, 2], c_fdis_itv_ls_OL[2, 2]),
  slope_97.5 = c(c_fdis[2, 3], c_fdis_itv[2, 3],c_fdis_OL[2, 3],c_fdis_itv_OL[2, 3],
                 c_fdis_la[2, 3], c_fdis_la_itv[2, 3],c_fdis_la_OL[2, 3],c_fdis_itv_la_OL[2, 3],
                 c_fdis_sla[2, 3], c_fdis_sla_itv[2, 3],c_fdis_sla_OL[2, 3],c_fdis_itv_sla_OL[2, 3],
                 c_fdis_ldmc[2, 3], c_fdis_ldmc_itv[2, 3],c_fdis_ldmc_OL[2, 3],c_fdis_itv_ldmc_OL[2, 3],
                 c_fdis_lt[2, 3], c_fdis_lt_itv[2, 3],c_fdis_lt_OL[2, 3],c_fdis_itv_lt_OL[2, 3],
                 c_fdis_lnc[2, 3], c_fdis_lnc_itv[2, 3],c_fdis_lnc_OL[2, 3],c_fdis_itv_lnc_OL[2, 3],
                 c_fdis_ls[2, 3], c_fdis_ls_itv[2, 3], c_fdis_ls_OL[2, 3], c_fdis_itv_ls_OL[2, 3]))

write.csv(fdis_table, "./output/FD/FDis_table.csv")

#### Figure
## including itv
plot_multi_fdis_itv <- plot_fdis(mod_fdis_itv, mod_fdis_itv_OL, "FDis multivariate")

# ggsave("./output/figure/fdis_multi_itv.pdf", height = 4, width = 4)

plot_la_fdis_itv <- plot_fdis(mod_fdis_la_itv, mod_fdis_la_itv_OL, "Leaf area")
plot_sla_fdis_itv <- plot_fdis(mod_fdis_sla_itv, mod_fdis_sla_itv_OL, "SLA")
plot_lt_fdis_itv <- plot_fdis(mod_fdis_lt_itv, mod_fdis_lt_itv_OL, "Leaf thickness")
plot_ldmc_fdis_itv <- plot_fdis(mod_fdis_ldmc_itv, mod_fdis_ldmc_itv_OL, "LDMC")
plot_lnc_fdis_itv <- plot_fdis(mod_fdis_lnc_itv, mod_fdis_lnc_itv_OL, "Leaf N content")
plot_ls_fdis_itv <- plot_fdis(mod_fdis_ls_itv, mod_fdis_ls_itv_OL, "Leaf length:width")

ggpubr::ggarrange(plot_la_fdis_itv, plot_sla_fdis_itv, plot_lt_fdis_itv,
                  plot_ldmc_fdis_itv, plot_lnc_fdis_itv, plot_ls_fdis_itv)
# ggsave("./output/figure/fdis_uni_itv.pdf", height = 4, width = 7)

# with median trait value
(plot_multi_fdis_med <- plot_fdis(mod_fdis, mod_fdis_OL, "FDis multivariate"))

plot_la_fdis_med <- plot_fdis(mod_fdis_la, mod_fdis_la_OL, "Leaf area")
plot_sla_fdis_med <- plot_fdis(mod_fdis_sla, mod_fdis_sla_OL, "SLA")
plot_lt_fdis_med <- plot_fdis(mod_fdis_lt, mod_fdis_lt_OL, "Leaf thickness")
plot_ldmc_fdis_med <- plot_fdis(mod_fdis_ldmc, mod_fdis_ldmc_OL, "LDMC")
plot_lnc_fdis_med <- plot_fdis(mod_fdis_lnc, mod_fdis_lnc_OL, "Leaf N content")
plot_ls_fdis_med <- plot_fdis(mod_fdis_ls, mod_fdis_ls_OL, "Leaf length:width")

ggpubr::ggarrange(plot_la_fdis_med, plot_sla_fdis_med, plot_lt_fdis_med,
                  plot_ldmc_fdis_med, plot_lnc_fdis_med, plot_ls_fdis_med)
# ggsave("./output/figure/fdis_uni_med.pdf", height = 4, width = 7)
