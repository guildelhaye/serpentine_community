#### Null models 
## sesFRic uses a independant swap randomisation that shuffle species identities between sites
## keeping the number of species within sites identical and the prevalence of species identical
## within the dataset. Guillaume Delhaye 2022. g.delhaye@kew.org
library(tidyverse)
library(FD)
library(picante)

### Multivariate FRic ####
sesFRicm <- function(traits, comm, ndim, numberReps){#ndim is the number of axes to keep in the PCoA
  FDobs <- dbFD(traits, comm, m = ndim, calc.FGR = F, calc.CWM = F, calc.FDiv = F, corr = c("lingoes"))
  #Lets create a matrix to store results from each iteration (one column per iteration) 
  resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                        dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  for(rep in 1:numberReps){ 
    commrandom <- randomizeMatrix(
      samp = comm, null.model = "independentswap") 
    FDnull <- dbFD(traits, commrandom, m = ndim, calc.FGR = F, 
                   calc.CWM = F, calc.FDiv = F, message = F, corr = c("lingoes"))
    resultsRandom[,rep] <- FDnull$FRic
  }

  obsFric <- FDobs$FRic
  meanNull <- rowMeans(resultsRandom) 
  ES <- obsFric - meanNull
  sdNull <- apply(resultsRandom, 1, sd) 
  SESFRic <- ES / sdNull
  return(SESFRic)
}

### Multivariate FDis ####
sesFDism <- function(traits, comm, ndim, numberReps){#ndim is the number of axes to keep in the PCoA
  FDobs <- dbFD(traits, comm, m = ndim, calc.FGR = F, 
                calc.CWM = F, calc.FDiv = F, calc.FRic = F, corr = c("lingoes"))
  #Lets create a matrix to store results from each iteration (one column per iteration) 
  resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                          dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  for(rep in 1:numberReps){ 
    traitsRand <- traits[sample(1:nrow(traits)),] 
    rownames(traitsRand) <- rownames(traits) 
    FDnull <- dbFD(traitsRand, comm, m = ndim, calc.FGR = F, 
                   calc.CWM = F, calc.FDiv = F, calc.FRic = F, message = F, corr = c("lingoes"))
    resultsRandom[,rep] <- FDnull$FDis
  }
  
  obsFDis <- FDobs$FDis
  meanNullFDis <- rowMeans(resultsRandom) 
  ES <- obsFDis - meanNullFDis
  sdNull <- apply(resultsRandom, 1, sd) 
  SESFDis <- ES / sdNull
  return(SESFDis)
}

### Multivariate FTI ####
sesFTI<- function(traits, comm, numberReps){#ndim is the number of axes to keep in the PCoA

  # select numeric traits because PCA
  n_sites <- nrow(comm)  
  FTI.obs = NULL
  # matrix to store results from each iteration (one column per iteration) 
  FTI.Random <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                          dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  
  for (a in 1:n_sites){ 

    # Select the communities
    com.obs=comm[a,comm[a,]>0] # Keeping species with abundance >0
    com.obs = com.obs/sum(com.obs)
    traits.obs=traits[names(com.obs),] #Keep species present 
    
    # Calculate observed index
    cor.mat <- cor(traits.obs, use = "complete.obs", method = "pearson")
    eigenvalues <- eigen(cor.mat)$values/sum(eigen(cor.mat)$values) #relative eigenvalues
    FTI.obs[a]<-sd(eigenvalues)

  ## randomize traits values per species but keeping traits relationships
  for(rep in 1:numberReps){ 
    #Random communities
    traitsRand <- traits[sample(1:nrow(traits)),] 
    rownames(traitsRand) <- rownames(traits) 
    (traits.obs.rand=traitsRand[names(com.obs),]) #Keep species present 
      # FTI for random communities
      cor.mat.rand <- cor(traits.obs.rand, use = "complete.obs", method = "pearson")
      eigenvalues.rand <- eigen(cor.mat.rand)$values/sum(eigen(cor.mat.rand)$values) #relative eigenval
      FTI.Random[a,rep]<-sd(eigenvalues.rand)
    }
  }

  meanNullFTI <- rowMeans(FTI.Random) 
  ES <- FTI.obs - meanNullFTI
  sdNull <- apply(FTI.Random, 1, sd) 
  SESFTI <- ES / sdNull
  return(SESFTI)
}

### traits range ####
sesRange <- function(traits, comm, i, numberReps){# i is trait name as character
  
  CommBin <- (comm >0)*1
  TrMat <- (t(CommBin) * na.omit(traits[,i]))
  ObsRange <- apply(TrMat, 2, max) - apply(TrMat, 2, min)

  #Lets create a matrix to store results from each iteration (one column per iteration) 
  resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                        dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  for(rep in 1:numberReps){ 
    commrandom <- randomizeMatrix(
      samp = comm, null.model = "independentswap") 
    CommBinNull <- (commrandom >0)*1
    TrMatNull <- (t(CommBinNull) * na.omit(traits[,i]))
    rangeNull <- apply(TrMatNull, 2, max) - apply(TrMatNull, 2, min)
    resultsRandom[,rep] <- rangeNull
  }

  meanNullRange <- rowMeans(resultsRandom) 
  ES <- ObsRange - meanNullRange
  sdNullRange <- apply(resultsRandom, 1, sd) 
  SESRange <- ES / sdNullRange
  return(SESRange)
}

## Univariate FDis
sesFDisu <- function(traits, comm, i, numberReps){ # i is the trait name in ' ', ex: 'SLA'

Tr <- traits %>% select(i) %>% filter(is.na(.) == F)
comm <- comm %>% select(any_of(rownames(Tr)))

FDobs <- dbFD(Tr, comm, calc.FGR = F, 
              calc.CWM = F, calc.FDiv = F, calc.FRic = F)
#Lets create a matrix to store results from each iteration (one column per iteration) 
resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                        dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
for(rep in 1:numberReps){ 
  traitsRand <- data.frame(Tr[sample(1:nrow(Tr)),])
  rownames(traitsRand) <- rownames(Tr) 
  FDnull <- dbFD(traitsRand, comm, m = ndim, calc.FGR = F, 
              calc.CWM = F, calc.FDiv = F, calc.FRic = F, message = F)
  resultsRandom[,rep] <- FDnull$FDis
}
obsFDis <- FDobs$FDis
meanNullFDis <- rowMeans(resultsRandom) 
ES <- obsFDis - meanNullFDis
sdNull <- apply(resultsRandom, 1, sd) 
SESFDis <- ES / sdNull
return(SESFDis)
}
