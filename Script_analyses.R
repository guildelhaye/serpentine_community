#### Script for the analyses of the paper "Community assembly of 
#### serpentine plants is driven by interspecific differences in 
#### functional traits". 
#### Created by Guillaume Delhaye. g.delhaye@kew.org
#### Last update 19/09/23

library(tidyverse)
library(ggplot2)
library(ade4)
library(vegan)
library(openxlsx)
library(ggfortify)
library(brms)
library(tidybayes)
library(lmerTest)
library(FD)
source("Script_FD.R") #script containing FD indices function and null model

### Data --------------------------------------------------------
setwd("./data")
soil <- read.xlsx("soil.xlsx")
traits.itv <- read.xlsx("traitsITV.xlsx")
abundance <- read.xlsx("abundance.xlsx", rowNames = T) 

# apply(abundance>0, 1, sum) #species richness in each plot

traits.itv <- traits.itv %>%
  mutate(LS = LL/as.numeric(LW), #add leaf shape
         Serp = case_when(site %in% c("Amp", "Oly") ~ "Serpentine",
                          site %in% c("Lout", "Vat") ~ "Non_Serpentine")) 


##### Soil analyses #--------------------------------------------------------------
soil.b <- soil %>% 
  mutate(MgCaRatio = Mg/Ca) %>%
  select(Site, Soil, Nsite, Ni, Cr, Co, MgCaRatio, Mn, Zn, P, K) %>%
  mutate(Ni = log(Ni), Cr = log(Cr), Co = log(Co), Mn = log(Mn), 
         MgCaRatio = log(MgCaRatio), P = log(P),  K = log(K)) 

## correlation begtween variables
plot(soil.b[,4:11])
round(cor(soil.b[,4:11]), 2)

## PCA
PC <- dudi.pca(soil.b[,4:11], center = TRUE, scale = TRUE, 
               scannf = F, nf = 3)
PC2 <- prcomp(soil.b[,4:11], center = TRUE, scale = TRUE)#check

summary(PC)
par(mfrow = c(1,2))
biplot(PC) ## Mostly 1 axis of soil condition variation 
PC$c1 # loadings
soil.pc <-PC$li #principal components stored in soil.pc
plot(soil.pc$Axis1, soil.pc$Axis2, pch = soil$ID); abline(h = 0, v = 0)

#Plot
pdf(file = './output/my_plot.pdf')
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
autoplot(PC2, data = soil.b, colour = "Site", type = "Site", loadings = TRUE,
         loadings.label = TRUE, loadings.label.size = 4, loadings.colour = 'black', 
         loadings.label.colour = 'black', frame = TRUE) +
    theme_minimal() +
  theme(legend.position = "none")
dev.off()

##### Variation partitionning # -----------------------------------------------------
setwd("./Analyses")

## Helper function to extract the proportion of variance
varpart_fit <- function(model, trait){
  model %>% 
    tidy_draws() %>% 
    select(starts_with("sd_"), sigma) %>% 
    transmute_all(.funs = list(sq = ~(. ^ 2))) %>% 
    mutate(total_var = rowSums(.)) %>% 
    mutate_at(.vars = vars(-total_var), 
              .funs = list(pct = ~(. / total_var))) %>% 
    map_df(.f = ~ median_hdci(., .width = 0.95), .id = "component") %>% 
    mutate(
      component = str_remove(component, "sd_"),
      component = str_remove(component, "__Intercept_sq"),
      component = str_replace(component, "sigma_sq", "residual")) %>%
    mutate(trait = paste0(trait))
}

## Data preparation
traits.norm <- traits.itv  %>%
  unite("Plot", site, soil, sep = "_", remove = FALSE) %>%
  select(-c(year, LL, LW)) %>%
  mutate(LDMC = scale(sqrt(LDMC)),
         LeafArea = scale(log10(LeafArea)),
         SLA = scale(sqrt(SLA)),
         LT = scale(log10(LT)),
         C = scale(log10(C)),
         N = scale(log10(N)), 
         LS = scale(log10(LS)))

#### Test variance partition between intra and inter
#### inter -> 1|Sp and intra is 1|Sp:soil
fit_N <- brm(formula = N ~ 1 + (1 | Sp) + (1 | Sp:Serp),
             data = traits.norm,
             control = list(adapt_delta = 0.95, max_treedepth = 10),
             iter = 4000, warmup = 1000, file = "output/N")

fit_SLA <- update(fit_N, SLA ~ ., newdata = traits.norm, file = "./output/SLA")
fit_LA <- update(fit_N, LeafArea ~ ., newdata = traits.norm, file = "./output/LA")
fit_LT <- update(fit_N, LT ~ ., newdata = traits.norm, file = "./output/LT")
fit_LDMC <- update(fit_N, LDMC ~ ., newdata = traits.norm, file = "./output/LDMC")
fit_LS <- update(fit_N, LS ~ ., newdata = traits.norm, file = "./output/LS")

fit_SLA <- readRDS("./output/SLA.RDS")
fit_LA <- readRDS("./output/LA.RDS")
varpart_result <-rbind(
  varpart_fit(fit_SLA, "SLA"),
  varpart_fit(fit_N, "N"), 
  varpart_fit(fit_LA, "LA"), 
  varpart_fit(fit_LT, "LT"), 
  varpart_fit(fit_LDMC, "LDMC"), 
  varpart_fit(fit_LS, "LS")
  )

# write.csv(varpart_result, "./output/varpart.csv")

##### Traits Var intra ##------------------------------------------------------------------
# itv.median = traits.itv %>%
#   group_by(Sp, species, Serp) %>%
#   summarise(
#     LS = median(LS, na.rm=TRUE), 
#     LeafArea = median(LeafArea, na.rm=TRUE),
#     LDMC = median(LDMC, na.rm=TRUE),
#     SLA = median(SLA, na.rm=TRUE),
#     LT = median(LT, na.rm=TRUE),
#     LNC = median(N,  na.rm=TRUE)) %>%
#   ungroup
# write.csv(itv.median, "./output/itv.median.csv")
# 
# ### Remove three species present only on one type of soil
# itv.median.b = itv.median %>% filter(!Sp %in% c("BROMCOM", "HIRSINCA", "TRIFCAMP"))
# sp.list = itv.median.b %>% select(Sp) %>% unique

MOD = data.frame(Species = NA, Est = NA, Tval = NA, Pval = NA)

for(i in 1:length(sp.list$Sp)){
  name = sp.list$Sp[i]
mod = summary(lmer(LS ~ Serp + (1 | site) , data = traits.itv %>% filter(Sp == name)))
#summary(lm(LeafArea ~ Serp , data = traits.itv %>% filter(Sp == "FILAERIO")))
MOD[i,] = c(name, round(mod$coefficients[2,1], 2), 
  round(mod$coefficients[2,4], 2), 
  round(mod$coefficients[2,5], 4))
}

write.xlsx(MOD, "/output/LeafShape_test.xlsx")

## Different species list for N because not N for all species
itv.median.c = itv.median %>% filter(!Sp %in% c("BROMCOM", "HIRSINCA", "TRIFCAMP", "PHLEPRAT"))
sp.list = itv.median.c %>% select(Sp) %>% unique
for(i in 1:length(sp.list$Sp)){
  name = sp.list$Sp[i]
  mod = summary(lmer(N ~ Serp + (1 | site) , data = traits.itv %>% filter(Sp == name)))
  #summary(lm(LeafArea ~ Serp , data = traits.itv %>% filter(Sp == "FILAERIO")))
  MOD[i,] = c(name, round(mod$coefficients[2,1], 2), 
              round(mod$coefficients[2,4], 2), 
              round(mod$coefficients[2,5], 4))
}
write.xlsx(MOD, "Results/LNC_test.xlsx")

##### Community assembly including intraspecific variation #----------------------------
## Data preparation ####
nrep = 1000
# transform abundance using hellinger
abundance <- decostand(abundance, method = "hellinger")

traits_itv_soil <- traits.itv |> 
  select(Sp, species, LeafArea, LDMC, SLA, LT, N, LS, Serp) |>
  mutate(LeafArea = LeafArea*100) %>%
  group_by(Sp, species, Serp) |>
  summarise_all(median, na.rm = T) |>
  reframe(
    Serp = Serp,
    LS = log(LS), 
    LDMC = LDMC, 
    SLA = SLA,
    LeafArea = log(LeafArea),
    LT = log(LT),
    LNC = log(N))

traits_itv_low <- traits_itv_soil |> 
  filter(Serp == "Non_Serpentine") |> 
  select(-c(species, Serp)) |> 
  column_to_rownames("Sp")
traits_itv_high <- traits_itv_soil |> 
  filter(Serp == "Serpentine") |> 
  select(-c(species, Serp)) |> 
  column_to_rownames("Sp")

abundance_high <- abundance[c(1:5, 15:19),colnames(abundance) %in% rownames(traits_itv_high) == TRUE]
traits_itv_high <- traits_itv_high[rownames(traits_itv_high)%in% colnames(abundance_high),]

abundance_low  <- abundance[c(6:14, 20:26),colnames(abundance) %in% rownames(traits_itv_low) == TRUE]
traits_itv_low <- traits_itv_low[rownames(traits_itv_low)%in% colnames(abundance_low),]

## Trifolium arvense occurs in none of the plots -> remove
abundance_low <- abundance_low |> 
  select(-TRIFARVE)
traits_itv_low <- traits_itv_low |> 
  filter(rownames(traits_itv_low) != "TRIFARVE")

## Community total ----------------------------
## CWM 
FD_low <- dbFD(traits_itv_low, abundance_low, corr = "lingoes")
FD_high <- dbFD(traits_itv_high, abundance_high, corr = "lingoes")
CWM_ITV <- data.frame(CWM = rbind(FD_low$CWM, FD_high$CWM))

## backtransform in normal units
CWM_ITV <- CWM_ITV %>%
  mutate(CWM.LeafArea= exp(CWM.LeafArea),
         CWM.LT= exp(CWM.LT),
         CWM.LNC = exp(CWM.LNC),
         CWM.LS = exp(CWM.LS)) |>
  # rename(CWM.LeafArea_ITV = CWM.LeafArea, CWM.LDMC_ITV = CWM.LDMC, CWM.SLA_ITV = CWM.SLA,
  #        CWM.LNC_ITV = CWM.LNC, CWM.LT_ITV = CWM.LT, CWM.LS_ITV = CWM.LS)|> 
  rownames_to_column("site")

## Multivariate FDis
fdism_low <- sesFDism(traits = traits_itv_low, 
                      comm = abundance_low, 
                      numberReps = nrep)
fdism_high <- sesFDism(traits = traits_itv_high, 
                       comm = abundance_high, 
                       numberReps = nrep)

FDM_ITV <- data.frame(FDis.Tot = c(fdism_low, fdism_high)) |>
  #rename(FDis.Tot_ITV = fdism) |>
  rownames_to_column("site")

## Univariate FDis
fdisu.LS_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'LS', numberReps = nrep)
fdisu.LS_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'LS', numberReps = nrep)

fdisu.LA_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'LeafArea', numberReps = nrep)
fdisu.LA_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'LeafArea', numberReps = nrep)

fdisu.SLA_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'SLA', numberReps = nrep)
fdisu.SLA_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'SLA', numberReps = nrep)

fdisu.LT_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'LT', numberReps = nrep)
fdisu.LT_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'LT', numberReps = nrep)

fdisu.LDMC_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'LDMC', numberReps = nrep)
fdisu.LDMC_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'LDMC', numberReps = nrep)

fdisu.LNC_l <- sesFDisu(traits = traits_itv_low, comm = abundance_low, i = 'LNC', numberReps = nrep)
fdisu.LNC_h <- sesFDisu(traits = traits_itv_high, comm = abundance_high, i = 'LNC', numberReps = nrep)


FDisU_ITV <- data.frame(FDis.LS = c(fdisu.LS_l, fdisu.LS_h), FDis.LeafArea = c(fdisu.LA_l, fdisu.LA_h), 
                        FDis.SLA = c(fdisu.SLA_l, fdisu.SLA_h) , FDis.LDMC = c(fdisu.LDMC_l, fdisu.LDMC_h), 
                        FDis.LT = c(fdisu.LT_l, fdisu.LT_h), FDis.LNC = c(fdisu.LNC_l, fdisu.LNC_h)) |>
  rownames_to_column("site")

Result_FD_ITV <- CWM_ITV |> left_join(FDM_ITV)  |> left_join(FDisU_ITV) |>
  left_join(soil |> unite(site, c(Site, Soil, Nsite), sep = "_") |> select(site, Ni))

## Without A.lesbiacum -----------------------
abundance_high_noal <- abundance_high[,-2]
traits_itv_high_noal <- traits_itv_high[-2,]

abundance_low_noal  <- abundance_low[,-2]
traits_itv_low_noal <- traits_itv_low[-2,]

FD_low_noal <- dbFD(traits_itv_low_noal, abundance_low_noal, corr = "lingoes")
FD_high_noal <- dbFD(traits_itv_high_noal, abundance_high_noal, corr = "lingoes")
CWM_ITV_noal <- data.frame(CWM = rbind(FD_low_noal$CWM, FD_high_noal$CWM))

## backtransform in normal units
CWM_ITV_noal <- CWM_ITV_noal %>%
  mutate(CWM.LeafArea= exp(CWM.LeafArea),
         CWM.LT= exp(CWM.LT),
         CWM.LNC = exp(CWM.LNC),
         CWM.LS = exp(CWM.LS)) |>
  # rename(CWM.LeafArea_ITV = CWM.LeafArea, CWM.LDMC_ITV = CWM.LDMC, CWM.SLA_ITV = CWM.SLA,
  #        CWM.LNC_ITV = CWM.LNC, CWM.LT_ITV = CWM.LT, CWM.LS_ITV = CWM.LS)|> 
  rownames_to_column("site")

## Multivariate FDis
fdism_low_noal <- sesFDism(traits = traits_itv_low_noal, 
                           comm = abundance_low_noal, 
                           numberReps = nrep)
fdism_high_noal <- sesFDism(traits = traits_itv_high_noal, 
                            comm = abundance_high_noal, 
                            numberReps = nrep)

FDM_ITV_noal <- data.frame(FDis.Tot = c(fdism_low_noal, fdism_high_noal)) |>
  #rename(FDis.Tot_ITV = fdism) |>
  rownames_to_column("site")

## Univariate FDis
fdisu.LS_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'LS', numberReps = nrep)
fdisu.LS_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'LS', numberReps = nrep)

fdisu.LA_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'LeafArea', numberReps = nrep)
fdisu.LA_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'LeafArea', numberReps = nrep)

fdisu.SLA_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'SLA', numberReps = nrep)
fdisu.SLA_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'SLA', numberReps = nrep)

fdisu.LT_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'LT', numberReps = nrep)
fdisu.LT_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'LT', numberReps = nrep)

fdisu.LDMC_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'LDMC', numberReps = nrep)
fdisu.LDMC_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'LDMC', numberReps = nrep)

fdisu.LNC_l_noal <- sesFDisu(traits = traits_itv_low_noal, comm = abundance_low_noal, i = 'LNC', numberReps = nrep)
fdisu.LNC_h_noal <- sesFDisu(traits = traits_itv_high_noal, comm = abundance_high_noal, i = 'LNC', numberReps = nrep)


FDisU_ITV_noal <- data.frame(FDis.LS = c(fdisu.LS_l_noal, fdisu.LS_h_noal), FDis.LeafArea = c(fdisu.LA_l_noal, fdisu.LA_h_noal), 
                             FDis.SLA = c(fdisu.SLA_l_noal, fdisu.SLA_h_noal) , FDis.LDMC = c(fdisu.LDMC_l_noal, fdisu.LDMC_h_noal), 
                             FDis.LT = c(fdisu.LT_l_noal, fdisu.LT_h_noal), FDis.LNC = c(fdisu.LNC_l_noal, fdisu.LNC_h_noal)) |>
  rownames_to_column("site")

Result_FD_ITV_noal <- CWM_ITV_noal |> left_join(FDM_ITV_noal)  |> left_join(FDisU_ITV_noal) |>
  left_join(soil |> unite(site, c(Site, Soil, Nsite), sep = "_") |> select(site, Ni))

colnames(Result_FD_ITV) == colnames(Result_FD_ITV_noal)

Results_FD_Final_ITV <- rbind(Result_FD_ITV, Result_FD_ITV_noal) |> 
  mutate(Com = rep(c("Tot", "No_Alesb"), each = 26))

write.csv(Results_FD_Final_ITV, 
          "../output/Results_FD_Final_ITV_hel.csv")

#### Community assembly with median trait value #####
## Data preparation ####
traits.median = traits.itv %>%
  group_by(Sp, species) %>%
  mutate(LeafArea = LeafArea*100) %>% #transform in mm2
  summarise(
    LS = log(median(LS, na.rm=TRUE)), 
    LeafArea = log(median(LeafArea, na.rm=TRUE)),
    LDMC = median(LDMC, na.rm=TRUE),
    SLA = median(SLA, na.rm=TRUE),
    LT = log(median(LT, na.rm=TRUE)),
    LNC = log(median(N,  na.rm=TRUE))) %>%
  ungroup %>%
  select(-species) %>%
  column_to_rownames("Sp") %>%
  filter(!row_number() %in% c(16)) # Remove a species present in none of the sites

# Community total ----------------------------
## CWM 
FD <- dbFD(traits.median, abundance, corr = "lingoes")
CWM <- data.frame(CWM = FD$CWM)

## backtransform in normal units
CWM <- CWM %>%
  mutate(CWM.LeafArea = exp(CWM.LeafArea),
         CWM.LT = exp(CWM.LT),
         CWM.LNC = exp(CWM.LNC),
         CWM.LS = exp(CWM.LS))|> 
  rownames_to_column("site")

## Multivariate FDis
fdism <- sesFDism(traits = traits.median, comm = abundance, numberReps = nrep)
FDM <- data.frame(FDis.Tot = fdism) |> rownames_to_column("site")

## Univariate FDis
FDis.LS <- sesFDisu(traits = traits.median, comm = abundance, 
                     i = 'LS', numberReps = nrep)
FDis.LeafArea <- sesFDisu(traits = traits.median, comm = abundance, 
                     i = 'LeafArea', numberReps = nrep)
FDis.SLA <- sesFDisu(traits = traits.median, comm = abundance, 
                     i = 'SLA', numberReps = nrep)
FDis.LT <- sesFDisu(traits = traits.median, comm = abundance, 
                      i = 'LT', numberReps = nrep)
FDis.LDMC <- sesFDisu(traits = traits.median, comm = abundance, 
                     i = 'LDMC', numberReps = nrep)
FDis.LNC <- sesFDisu(traits = traits.median, comm = abundance, 
                      i = 'LNC', numberReps = nrep) 

FDisU <- data.frame(FDis.LS, FDis.LeafArea, FDis.SLA, 
                    FDis.LDMC, FDis.LT, FDis.LNC) |> 
  rownames_to_column("site")

Result_FD <- CWM |> left_join(FDM)  |> left_join(FDisU) |>
  left_join(soil |> unite(site, c(Site, Soil, Nsite), sep = "_") |> select(site, Ni))

### Without A.lesbiacum -----------------------
traits.median_noal <- traits.median[-2,]
abundance_noal <- abundance[,-2]

## CWM 
FD_noal <- dbFD(traits.median_noal, abundance_noal, corr = "lingoes")
CWM_noal <- data.frame(CWM = FD_noal$CWM)

## backtransform in normal units
CWM_noal <- CWM_noal %>%
  mutate(CWM.LeafArea = exp(CWM.LeafArea),
         CWM.LT = exp(CWM.LT),
         CWM.LNC = exp(CWM.LNC),
         CWM.LS = exp(CWM.LS)) |> 
  rownames_to_column("site")

## Multivariate FDis
FDis.Tot_noal <- sesFDism(traits = traits.median_noal, comm = abundance_noal, numberReps = nrep)
FDM_noal <- data.frame(FDis.Tot = FDis.Tot_noal) |> 
  rownames_to_column("site")

## Univariate FDis
FDis.LS_noal <- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                     i = 'LS', numberReps = nrep)
FDis.LeafArea_noal <- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                     i = 'LeafArea', numberReps = nrep)
FDis.SLA_noal <- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                      i = 'SLA', numberReps = nrep)
FDis.LT_noal <- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                     i = 'LT', numberReps = nrep)
FDis.LDMC_noal <- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                       i = 'LDMC', numberReps = nrep)
FDis.LNC_noal<- sesFDisu(traits = traits.median_noal, comm = abundance_noal, 
                     i = 'LNC', numberReps = nrep) 

FDisU_noal <- data.frame(FDis.LS=FDis.LS_noal, FDis.LeafArea=FDis.LeafArea_noal, FDis.SLA=FDis.SLA_noal, 
                         FDis.LDMC=FDis.LDMC_noal, FDis.LT=FDis.LT_noal, FDis.LNC=FDis.LNC_noal) |> 
  rownames_to_column("site")

Result_FD_noal <- CWM_noal |> left_join(FDM_noal)  |> left_join(FDisU_noal) |>
  left_join(soil |> unite(site, c(Site, Soil, Nsite), sep = "_") |> select(site, Ni))

colnames(Result_FD) == colnames(Result_FD_noal)

Results_FD_Final_Aver <- rbind(Result_FD, Result_FD_noal) |> 
  mutate(Com = rep(c("Tot", "No_Alesb"), each = 26))

## Results final 
write.csv(Results_FD_Final_Aver, 
          "../output/Results_FD_Final_Aver_hel.csv")

### END --------------------------------------------------------------------------------