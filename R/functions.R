# do CCA for each depth and season and determine total inertia explained by each temperature ===========
doCCA <- function(region, indx){
  datRegion <- coretopSpp[[which(names(coretopSpp) == region)]]
  spp <- datRegion$species[indx, ]
  env <- datRegion$env[indx, ]
  map_dfr(env, function(x){
    # remove samples where no environmental data available (shallow sites)
    naIndx <- which(!is.na(x))
    e <- x[naIndx]
    spp <- spp[naIndx, ]
    mod <- cca(spp ~ e)
    eigenvals(mod)[1]/sum(eigenvals(mod))
  }, .id = 'SeasonDepth') %>%
    rename(EX = CCA1) %>%
    separate(SeasonDepth, c('season', 'depth'))
}

# predict temperatures from training set  ===========================================================
# note that temperature data are not available for each depth for shallow core top samples
# so I calculate a dissimilarity matrix each time (there may be a more efficient way of doing this)

# spp: core top species
# fos: fossil species
# env: environmental variables
# k: # of analogues

doMAT <- function(spp, fos = NULL, env, k = 10){
  require(vegan)
  require(rioja)
  # check which env is NA and remove env and species for which no env available
  naIdx <- which(!is.na(env))
  envSub <- env[naIdx]
  sppIdx <- which(spp$meta$forcensID %in% names(envSub))
  samples <- spp$meta$forcensID[sppIdx]
  modSpp <- lapply(spp, function(x) x[sppIdx,])
  
  # dissimilarity
  # for core top
  disModern <- function(modSpp){
    diss <- as.matrix(paldist(modSpp$species, dist.method = 'sq.chord'))
    ind <- apply(diss, 2, order)
    dist.n <- t(apply(diss, 2, sort)[2:(2 + k - 1), , drop = FALSE]) # take first k, leave out identical sample
    rownames(dist.n) <- rownames(modSpp$species) # sample names
    colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep = "")
    if (any(dist.n < 1e-06)) {
      dist.n[dist.n < 1e-06] <- 1e-06
    }
    nms <- t(matrix(as.integer(rownames(modSpp$species)[ind[2:(2 + k - 1), , drop = FALSE]]), nrow = k))
    rownames(nms) <- rownames(modSpp$species)
    colnames(nms) <- paste("N", sprintf("%02d", 1:k), sep = "")
    uIDs <- t(matrix(modSpp$meta$forcensID[ind[2:(2 + k - 1), , drop = FALSE]], nrow = k))
    rownames(uIDs) <- rownames(modSpp$species)
    colnames(uIDs) <- paste("N", sprintf("%02d", 1:k), sep = "")
    list(dist = dist.n, ID = uIDs, nms = nms)
  }
  
  # for fossil
  disFossil <- function(modSpp, fosSpp){
    diss <- paldist2(modSpp$species, fosSpp$species, dist.method = 'sq.chord')
    ind <- apply(diss, 2, order)
    dist.n <- t(apply(diss, 2, sort)[1:(1 + k - 1), , drop = FALSE])
    if (any(dist.n < 1e-06)) {
      dist.n[dist.n < 1e-06] <- 1e-06
    }
    rownames(dist.n) <- rownames(fosSpp)
    colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep = "")
    nms <- t(matrix(as.integer(rownames(modSpp$species)[ind[1:(1 + k - 1), , drop = FALSE]]), nrow = k))
    rownames(nms) <- rownames(fosSpp$species)
    colnames(nms) <- paste("N", sprintf("%02d", 1:k), sep = "")
    uIDs <- t(matrix(modSpp$meta$forcensID[ind[1:(1 + k - 1), , drop = FALSE]], nrow = k))
    rownames(uIDs) <- rownames(fosSpp$species)
    colnames(uIDs) <- paste("N", sprintf("%02d", 1:k), sep = "")
    list(dist = dist.n, ID = uIDs, nms = nms)
  }
  
  dis <- if(is.null(fos)){
    disModern(modSpp)
  } else {
    disFossil(modSpp, fos)
  }
  
  # predict environmental variable
  predictEnv <- function(env, dis){
    x.n <- t(matrix(sapply(dis$ID, function(x) env[names(env) == x]), nrow = k, byrow = TRUE)) # get values from env of k closest analogues
    predT <- rowSums(x.n/dis$dist)/rowSums(1/dis$dist) # weighted mean
    predT
  }
  out <- tibble(predT = predictEnv(envSub, dis))
  out$latitude <- if(is.null(fos)){modSpp$meta$Latitude} else {fos$meta$Latitude}
  out$samples <- if(is.null(fos)){samples}
  out$sppIdx <- if(is.null(fos)){sppIdx}
  return(out)
}

# variance explained by predicted temperatures ==================================================================
# spp: species, either modern or fossil
# pred: output from doRandomMAT
# indx: index for subsetting spp (tropical | extratropical)

varEx <- function(spp, pred, indx = NULL){
  require(vegan)
  spp <- if(!is.null(pred$sppIdx)){spp$species[pred$sppIdx,]} else {
    spp$species
  }
  Spp <- if(is.null(indx)){spp} else {spp[indx,]}
  predT <- if(is.null(indx)){pred$predT} else {pred$predT[indx]}
  
  rdaspp <- rda(Spp ~ predT)
  rdaspp$CCA$tot.chi/rdaspp$tot.chi
}

# sort taxonomy ==================================================================
# needed when merging regional data for similarity decay
replaceRuber <- function(x){
  x$Globigerinoides_ruber_Globigerinoides_white <- x$Globigerinoides_ruber + x$Globigerinoides_white
  select(x , -c(Globigerinoides_ruber, Globigerinoides_white))
}

replaceMenardii <- function(x){
  x$Globorotalia_menardii_Globorotalia_tumida <- x$Globorotalia_menardii + x$Globorotalia_tumida
  select(x , -c(Globorotalia_menardii, Globorotalia_tumida))
}

# make df with environmental and community dissimilarity =========================================
makeDist <- function(x){ # x = modernSppDF
  require(vegan)
  require(reshape2)
  # similarity matrix for species
  modSppDist <- as.matrix(vegdist(select(x, Dentigloborotalia_anfracta:Globorotalia_menardii_Globorotalia_tumida), method = "bray"))
  # convert to similarity
  modSppDist <- 1- modSppDist
  modSppDist <- melt(modSppDist)[melt(upper.tri(modSppDist, diag = FALSE))$value,]
  names(modSppDist)[3] <- 'commdist'
  # for environmental values
  modEnvDist <- as.matrix(dist(select(x, annual_50)))
  modEnvDist <- melt(modEnvDist)[melt(upper.tri(modEnvDist, diag = FALSE))$value,]
  names(modEnvDist)[3] <- 'envdist'
  # combine
  Dist <-bind_cols(modSppDist, select(modEnvDist, envdist))
  tibble(Dist)
}

# get environmental distances from model output =========================================
getEnvDist_avg <- function(x){ # x = temperature, y = community dissimilarity
  LGMEnvDist <- as.matrix(dist(unname(x)))
  LGMEnvDist <- melt(LGMEnvDist)[melt(upper.tri(LGMEnvDist, diag = FALSE))$value,]
  names(LGMEnvDist)[3] <- 'envdist'
  LGMEnvDist
  left_join(LGMEnvDist,
            lgmSppDist_avg,
            by = c('Var1', 'Var2'))
}


# function to grid similarity decay data =========================================
gridFunction <- function(x, resEnv, resSim){
  
  maxEnv <- ceiling(max(x$envdist, na.rm = TRUE)/resEnv)
  
  envBins <- seq(0, maxEnv*resEnv, resEnv)
  envLabels <- envBins[-1] - 0.5*resEnv
  
  simBins <- seq(0, 1, resSim)
  simLabels <- simBins[-1] - 0.5*resSim
  
  x %>%
    select(envdist, commdist) %>%
    drop_na() %>%
    mutate(envBin = cut(envdist, breaks = envBins, include.lowest = TRUE, right = FALSE, labels = envLabels),
           simBin = cut(commdist, breaks = simBins, include.lowest = TRUE, right = FALSE, labels = simLabels)) %>%
    group_by(envBin, simBin) %>%
    summarise(n = n(), .groups = 'drop') %>%
    ungroup() %>%
    mutate(envBin = as.numeric(levels(envBin))[envBin],
           simBin = as.numeric(levels(simBin))[simBin],
           logn = log10(n),
           nnorm = n/max(n),
           lognnorm = logn/max(logn))
}


# make named list  =========================================
named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = "_")))
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}

# define assemblage limits =========================================
limfunction <- function(x, resSim, quant){
  
  simBins <- seq(0, 1, resSim)
  simLabels <- simBins[-1] - 0.5*resSim
  
  res <- x %>%
    select(envdist, commdist) %>%
    drop_na() %>%
    mutate(simBin = cut(commdist, breaks = simBins, include.lowest = TRUE, right = FALSE, labels = simLabels)) %>%
    group_by(simBin) %>%
    summarise(Tlim = quantile(envdist, quant)) %>%
    mutate(simBinPlot = as.numeric(levels(simBin))[simBin],
           simBinPlot = simBinPlot - resSim/2
    )
  bind_rows(res, tibble(simBin = NA, Tlim  = tail(res$Tlim, 1), simBinPlot = 1))  
}


