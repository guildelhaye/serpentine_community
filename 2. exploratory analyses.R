### Exploratory analyses 
library(tidyverse)
library(vegan)
library(ggfortify)

## load data
soil.b <- read.csv("./data/soil_clean.csv") |> select(-X)
traits_sp_site <- read.csv("./data/traits_sp_site.csv") |> select(-X)
abun <- read.csv("./data/abundance_clean.csv") |> select(-X)


###### 1. Correlation structures in soil and trait data -----
## Soil analyses ----------------------------------------------
# correlation begtween variables
plot(soil.b[,3:10])
round(cor(soil.b[,3:10]), 2)

## PCA
PC2 <- prcomp(soil.b[,3:10], center = TRUE, scale = TRUE)

#Plot
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
autoplot(PC2, data = soil.b, colour = "site", type = "site", loadings = TRUE,
         loadings.label = TRUE, loadings.label.size = 4, loadings.colour = 'black', 
         loadings.label.colour = 'black', frame = TRUE) +
  theme_minimal() +
  theme(legend.position = "none")

##Relationship between traits ----------------------------
traits <- traits_sp_site |> 
  na.omit()|>
  left_join(soil.b |> 
              select(site, Ni) |> 
              summarise(.by = site, Ni = log(mean(Ni))))
# PCA
PC2 <- prcomp(na.omit(traits[,c(3:8)]), center = TRUE, scale = TRUE)#check

#Plot
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
autoplot(PC2, data = traits, colour = "Sp", 
         # type = "species", 
         loadings = TRUE,
         loadings.label = TRUE, 
         loadings.label.size = 4, 
         loadings.colour = 'black', 
         loadings.label.colour = 'black', 
         frame = TRUE) +
  theme_minimal() 
# theme(legend.position = "none")

ggsave("./output/trait_pca.pdf")
round(cor(traits[,c(3, 4, 5, 7, 8)]), 2)

##  RDA to investigate the effect of environment on community composition 
rda <- rda(abun ~ Ni + Cr + Co + MgCaRatio + Mn + Zn + P + K, 
           data = soil.b)

summary(rda)$biplot 
# Ni and Co positively and K negatively correlated to PC1

par(mfrow = c(1,1))
plot(rda)
