
library(terra)
library(tidyverse)

select <- dplyr::select

upscale <- 1000


# MERGE TILES ===========

## merge ####
cdl_merge <- function(x){
      fout <- paste0("big_data/merged/", x, ".tif")
      if(file.exists(fout)) return("skipped")
      message(fout)
      list.files("big_data/tiles", pattern = paste0(x, "-"), full.names = T) %>%
            map(rast) %>%
            map(function(r) if(x == "confidence"){mean(r)}else{r}) %>%
            c(list(filename = fout)) %>% 
            do.call("merge", .)
}
f <- c("crops", "groups",
       "crop_abandon_moving_8e_7r", "group_abandon_moving_8e_7r",
       "crop_introduce_moving_8e_7r", "group_introduce_moving_8e_7r",
       "crop_diversity_2008_2022", "crop_diversity_2008_2015", "crop_diversity_2015_2022",
       "group_diversity_2008_2022", "group_diversity_2008_2015", "group_diversity_2015_2022",
       "reclamation_abandon_moving_8e_7r", "reclamation_introduce_moving_8e_7r",
       "reclassified")
f %>% map(cdl_merge)

## crop to common extent ####
f <- list.files("big_data/merged", full.names = T)
x <- f %>%
      map(rast) %>%
      map(ext) %>%
      map(as.vector) %>%
      do.call("rbind", .)
x <- ext(max(x[,1]), min(x[,2]), max(x[,3]), min(x[,4]))
f %>% map(function(ff) ff %>% rast() %>% crop(x) %>% writeRaster(ff, overwrite = T))











# SWITCHING ==========

agg_switch <- function(level = "crop", 
                       phenom = "introduce",
                       window = "moving", 
                       outcome = "switch", 
                       timeframe = "8e_7r"){
      
      # input admin
      f <- list.files("big_data/merged", full.names = T,
                      pattern = paste0(paste(level, phenom, window, timeframe, sep = "_"), "\\.tif"))
      if(length(f) != 1) return("input file problem")
      
      # output admin
      outfile <- paste0("data/aggregated/", 
                        paste(level, phenom, window, outcome, timeframe, sep = "_"), ".tif")
      if(file.exists(outfile)) return("skipped")
      message(outfile)
      
      # compute
      out <- ifelse(outcome == "switch", 2, 1)
      z <- rast(f) %>% 
            aggregate(fun = function(x, ...) sum(x == out), 
                      fact = upscale, filename = outfile)
      return(f)
}

expand_grid(level = c("crop", "group"),
            phenom = c("introduce", "abandon"), # note, "abandon" is called "discontinuation" in manuscript
            window = "moving",
            outcome = c("switch", "stasis"),
            timeframe = c("8e_7r")) %>%
      pmap(agg_switch)


# FALLOW ======###########

agg_fallow <- function(phenom = "introduce", 
                       window = "moving", 
                       timeframe = "8e_7r"){
      
      # input admin
      f <- list.files("big_data/merged", full.names = T,
                      pattern = paste0(paste(
                            "group", phenom, window, timeframe, sep = "_"), "\\.tif"))
      f2 <- list.files("big_data/merged", full.names = T,
                      pattern = paste0(paste(
                            "reclamation", phenom, window, timeframe, sep = "_"), "\\.tif"))
      if(length(f) != 1) return("input file problem")
      
      # output admin
      outfile <- paste0("data/aggregated/", 
                        paste("fallow", phenom, window, "switch", timeframe, sep = "_"), ".tif")
      if(file.exists(outfile)) return("skipped")
      message(outfile)
      
      # compute
      n <- nlyr(rast(f))
      i <- switch(phenom,
                  "introduce" = sort(16 - 1:n),
                  "abandon" = 1:n)
      z <- ((rast(f) == 2 & rast("big_data/merged/groups.tif")[[i]] == 305) | # either switched to/from fallow
                  (rast(f2) == 2)) %>% # or entire ref period is fallow
            aggregate(fun = sum, fact = upscale, na.rm = T, filename = outfile, overwrite = T)
      return(outfile)
}

expand_grid(phenom = c("abandon", "introduce"),
            window = "moving",
            timeframe = c("8e_7r")) %>%
      pmap(agg_fallow)



# DIVERSITY =====

## temporal ======
agg_diversity <- function(x){
      outfile <- paste0("data/aggregated/", x, ".tif")
      if(file.exists(outfile)) return("skipped")
      r <- list.files("big_data/merged", pattern = x, full.names = T) %>%
            rast() %>%
            aggregate(fun = sum, fact = upscale, na.rm = T, 
                      filename = outfile, overwrite = T)
      return(outfile)
}
c("crop_diversity_2008_2022", "crop_diversity_2008_2015", "crop_diversity_2015_2022", 
  "group_diversity_2008_2022", "group_diversity_2008_2015", "group_diversity_2015_2022") %>%
      map(agg_diversity)


## spatial =====
shannon <- function(x){
      f <- as.matrix(table(na.omit(x)))
      p <- f / sum(f)
      - sum(p * log(p))
}

crops <- rast("big_data/merged/crops.tif")
NAflag(crops) <- 0
div <- aggregate(crops, fun = shannon, fact = upscale,
                 filename = "data/aggregated/crop_spdiversity.tif")


## spatiotemporal =======
div <- aggregate(crops, fun = shannon, fact = c(upscale, upscale, nlyr(crops)),
                 filename = "data/aggregated/crop_stdiversity.tif")
div <- aggregate(crops[[1:8]], fun = shannon, fact = c(upscale, upscale, nlyr(crops)),
                 filename = "data/aggregated/crop_stdiversity_2008_2015.tif")
div <- aggregate(crops[[8:15]], fun = shannon, fact = c(upscale, upscale, nlyr(crops)),
                 filename = "data/aggregated/crop_stdiversity_2015_2022.tif")




# DIVERGENCE ==============================

intro <- rast("big_data/merged/crop_introduce_moving_8e_7r.tif")
aband <- rast("big_data/merged/crop_abandon_moving_8e_7r.tif")
crops <- rast("big_data/merged/crops.tif")
NAflag(crops) <- 0

convergence <- function(x, stat = "diff", refyr = 1, ...){
      if(length(x) != 3e6) return(NA)
      
      m <- matrix(x, ncol = 3)
      m <- m[is.finite(m[, 1]),] # cropland only
      if(length(dim(m)) == 0) m <- matrix(m, 1)
      if(nrow(m) == 0) return(NA)
      
      freq <- table(m[, refyr])
      freq <- freq / sum(freq)
      i <- as.integer(names(freq))
      
      s <- m[m[, 3] == 2, 1:2] # switched cells
      s[is.na(s)] <- 0
      if(length(dim(s)) == 0) s <- matrix(s, 1)
      s[, 1] <- freq[match(s[, 1], i)]
      s[, 2] <- freq[match(s[, 2], i)]
      s[is.na(s)] <- 0
      if(stat == "diff") y <- mean(s[, 2] - s[, 1]) # positive values indicate convergence
      if(stat == "sign") y <- mean(s[, 2] > s[, 1])
      y
}

conv <- function(year = 2015, phenom = "intro", stat = "sign", refyr = 1){
      message(year)
      if(phenom == "intro"){
            cs <- c(crops[[year - 2015 + c(7, 8)]], # crops switched from, to
                    intro[[year - 2014]]) # switched pixels
      }else{
            cs <- c(crops[[year - 2008 + c(1, 2)]], # crops switched from, to
                    aband[[year - 2007]]) # switched pixels
      }
      terra::aggregate(cs, fun = convergence, fact = c(upscale, upscale, nlyr(cs)), 
                       stat = stat, refyr = refyr)
}

cis <- map(2015:2022, conv, phenom = "intro", stat = "sign") %>% rast() %>% 
      writeRaster("data/aggregated/crop_convergence_introduce.tif")
cis <- map(2015:2022, conv, phenom = "intro", stat = "sign", refyr = 2) %>% rast() %>% 
      writeRaster("data/aggregated/crop_convergence2_introduce.tif")

cis <- map(2008:2015, conv, phenom = "abandon", stat = "sign") %>% rast() %>% 
      writeRaster("data/aggregated/crop_convergence_abandon.tif")
cis <- map(2008:2015, conv, phenom = "abandon", stat = "sign", refyr = 2) %>% rast() %>% 
      writeRaster("data/aggregated/crop_convergence2_abandon.tif")





# PARIWISE ===========

backbone <- expand_grid(b = c(0, 301:305),
                        r = c(0, 301:305)) %>%
      mutate(br = b*1000 + r) %>%
      select(br) %>%
      distinct() %>%
      pull(br) %>%
      as.character()

intro <- rast("big_data/merged/crop_introduce_moving_8e_7r.tif")
aband <- rast("big_data/merged/crop_abandon_moving_8e_7r.tif")

fx <- function(x, ...){
      f <- table(x)
      z <- rep(0, 36)
      z[match(names(f), backbone)] <- f
      z
}

pairwise_sw <- function(y, ref, phenom = "introduce"){
      
      # admin
      outfile <- paste0("data/aggregated/pairwise_", phenom, "_", y, "_", ref, ".tif")
      message(outfile)
      if(file.exists(outfile)) return(NULL)
      
      # switching data
      s <- rast(paste0("big_data/merged/group_", phenom, "_moving_8e_7r.tif"))
      i <- ifelse(phenom == "introduce", y-2015, y-2008)
      s <- s[[grepl(i, names(s))]] == 2
      
      # group type data
      g <- rast("big_data/merged/groups.tif")
      NAflag(g) <- 0
      i <- ifelse(phenom == "introduce", -ref, ref)
      g <- g[[grepl(paste0(y , "|", y+i, collapse = ""), 
                    names(g))]]
      
      g0 <- g[[1]]
      g1 <- g[[2]]
      z <- (g0*1000 + g1) * s
      z %>% 
            aggregate(upscale, fun = fx) %>%
            setNames(str_pad(backbone, 6, "left", 0)) %>%
            writeRaster(outfile, overwrite = T)
}

bind_rows(expand_grid(y = 2008:2015, ref = 1:7, phenom = "abandon"),
          expand_grid(y = 2015:2022, ref = 1:7, phenom = "introduce")) %>%
      pmap(pairwise_sw)




# COMBOS ==========

convergence <- function(x, out = 1, ...){
      if(length(x) != 3e6) return(NA)
      
      m <- matrix(x, ncol = 3)
      m <- m[is.finite(m[, 1]),] # cropland only
      if(length(dim(m)) == 0) m <- matrix(m, 1)
      # if(nrow(m) > 0) browser()
      freq <- table(m[, 1])
      freq <- freq / sum(freq)
      i <- as.integer(names(freq))
      
      s <- m[m[, 3] > 0, 1:2] # switched cells
      s[is.na(s)] <- 0
      if(length(dim(s)) == 0) s <- matrix(s, 1)
      s[, 1] <- freq[match(s[, 1], i)]
      s[, 2] <- freq[match(s[, 2], i)]
      s[is.na(s)] <- 0
      
      conv <- as.integer(s[, 2] > s[, 1])
      y <- table(m[m[, 3] > 0, 3] * 10 + conv)
      z <- rep(0, 4)
      z[match(names(y), c("10", "11", "20", "21"))] <- y
      z[out]
}

combo_intro <- function(year = 2015){
      
      message(year)
      outfile <- paste0("data/aggregated/tidc_combo_introduce_", year, ".tif")
      if(file.exists(outfile)) return("skipped")
      
      phenom = "introduce"
      window = "moving"
      timeframe = "8e_7r"
      
      # transformative vs incremental switching
      fc <- list.files("big_data/merged", full.names = T,
                       pattern = paste0(paste("crop", phenom, window, timeframe, sep = "_"), "\\.tif"))
      fg <- list.files("big_data/merged", full.names = T,
                       pattern = paste0(paste("group", phenom, window, timeframe, sep = "_"), "\\.tif"))
      cs <- rast(fc)[[year - 2014]] == 2 # switched crops
      gs <- rast(fg)[[year - 2014]] == 2 # switched groups
      ti <- cs + gs # 2 = transformative, 1 = incremental, 0 = none
      
      # remove switches to fallow
      not_fallow <- rast("big_data/merged/groups.tif")[[year - 2015 + 8]] != 305
      ti <- ti * not_fallow
      
      # divergent vs convergent switching
      crops <- rast("big_data/merged/crops.tif")
      NAflag(crops) <- 0
      cti <- c(crops[[year - 2015 + c(7, 8)]], # crops switched from, to
               ti) # switched pixels
      
      cmb <- c(terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 1),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 2),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 3),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 4))
      names(cmb) <- c("incr_div", "incr_conv", "trans_div", "trans_conv")
      
      writeRaster(cmb, outfile)
}

combo_aband <- function(year = 2015){
      
      message(year)
      outfile <- paste0("data/aggregated/tidc_combo_abandon_", year, ".tif")
      if(file.exists(outfile)) return("skipped")
      
      phenom = "abandon"
      window = "moving"
      timeframe = "8e_7r"
      
      # transformative vs incremental switching
      fc <- list.files("big_data/merged", full.names = T,
                       pattern = paste0(paste("crop", phenom, window, timeframe, sep = "_"), "\\.tif"))
      fg <- list.files("big_data/merged", full.names = T,
                       pattern = paste0(paste("group", phenom, window, timeframe, sep = "_"), "\\.tif"))
      cs <- rast(fc)[[year - 2007]] == 2 # switched crops
      gs <- rast(fg)[[year - 2007]] == 2 # switched groups
      ti <- cs + gs # 2 = transformative, 1 = incremental, 0 = none
      
      # remove switches from fallow
      not_fallow <- rast("big_data/merged/groups.tif")[[year - 2015 + 8]] != 305
      ti <- ti * not_fallow
      
      # divergent vs convergent switching
      crops <- rast("big_data/merged/crops.tif")
      NAflag(crops) <- 0
      cti <- c(crops[[year - 2015 + c(8, 9)]], # crops switched from, to
               ti) # switched pixels
      
      cmb <- c(terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 1),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 2),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 3),
               terra::aggregate(cti, fun = convergence, fact = c(upscale, upscale, nlyr(cti)), out = 4))
      names(cmb) <- c("incr_div", "incr_conv", "trans_div", "trans_conv")
      
      writeRaster(cmb, outfile)
}

map(2015:2022, combo_intro)
map(2015:2005, combo_aband)



# OTHER ==========================================

# average number of years switched by pixels that switched at least once
for(phenom in c("abandon", "introduce")){
      
      # crop switches
      cs <- rast(paste0("big_data/merged/crop_", phenom, "_moving_8e_7r.tif")) == 2
      
      # remove fallow switches
      n <- nlyr(cs)
      i <- switch(phenom,
                  "introduce" = sort(16 - 1:n),
                  "abandon" = 1:n)
      fs <- ((rast(paste0("big_data/merged/group_", phenom, "_moving_8e_7r.tif")) == 2 & 
              rast("big_data/merged/groups.tif")[[i]] == 305) | # either switched to/from fallow
                  (rast(paste0("big_data/merged/reclamation_", phenom, "_moving_8e_7r.tif")) == 2)) # or entire ref period is fallow
      cs <- cs & !fs
      
      mcs <- max(cs) # switched at least once
      scs <- sum(cs) # number of switches
      y <- aggregate(mcs, fun = sum, fact = upscale,
                     filename = paste0("data/aggregated/crop_", phenom, "_moving_switchONCE_8e_7r.tif"),
                     overwrite = T)
      y <- aggregate(scs * mcs, fun = sum, fact = upscale,
                     filename = paste0("data/aggregated/crop_", phenom, "_moving_switchMEAN_8e_7r.tif"),
                     overwrite = T)
}



# mean number of years reclassified for mean pixel in each cell
r <- rast("big_data/merged/reclassified.tif")
NAflag(r) <- 0
mr <- (aggregate(r, "mean", fact = 1000, na.rm = T) - 100) / 15
writeRaster(mr, "data/aggregated/reclassified.tif")
