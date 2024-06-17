
library(tidyverse)
library(terra)
library(furrr)
library(sf)
plan(multisession, workers = 5)
options(future.globals.maxSize = 4e9)

from_scratch <- F

# spatiotemporal references
states <- tigris::fips_codes %>%
      select(state, state_code) %>%
      distinct() %>%
      filter(! state %in% c("AK", "HI")) %>%
      slice(1:49)
years <- 2008:2022

# load CDL confusion matrices
matrices <- years %>%
      map(function(year){
            dirs <- list.dirs("big_data/cdl_confusion")
            dir <- dirs[str_detect(dirs, as.character(year))]
            
            states$state %>%
                  future_map(function(state){
                        if(state == "DC") state <- "MD"
                        file <- list.files(dir, full.names = T,
                                           pattern = paste0("_", state))
                        cm <- readxl::read_xlsx(file[1], "All - Matrix") %>% 
                              select("1":"256") %>% 
                              slice(1:256) %>%
                              as.matrix()
                        cm[!is.finite(cm)] <- 0
                        return(cm)
                  }) %>%
                  setNames(states$state)
      }) %>%
      setNames(years)
matrices_national <- lapply(matrices, function(x) Reduce("+", x))


# sample crop time series (using regular sample for comp speed)
if(from_scratch){
      crops <- rast("data/merged/crops.tif")
      s <- spatSample(crops, c(4000, 8000), method = "regular", as.raster = T)
      writeRaster(s, "data/crop_samples_4k_8k.tif")
}else{
      s <- rast("data/crop_samples_4k_8k.tif")
}

# map grid cells to states
state_shapes <- tigris::states(cb = TRUE)
fips_rast <- state_shapes %>% mutate(fips = as.integer(STATEFP)) %>% vect() %>% project(s) %>% rasterize(s, "fips")

# map grid cells to ERS regions
region <- st_read("data/regions.shp") %>% vect() %>% project(s) %>% rasterize(s, "region")

# crop data to matrix
crp <- c(fips_rast, region, s)
c1 <- crds(crp[[2]], na.rm = F)
c2 <- values(crp, na.rm = F)
crp <- na.omit(cbind(c1, c2))
crp <- crp[crp[, 5] != 0, ]
rm(c1, c2); gc()

# function to convert confusion counts to probabilities
conf_to_prob <- function(conf, margin = 1){
      cm <- conf %>% apply(margin, function(x) x/sum(x)) %>% t()
      cm[!is.finite(cm)] <- 0
      cm
}

# probability of a switch happening
switch_prob <- function(crops, # 8-column matrix of crops over space x time 
                        cm){# confusion matrices (array)
      ps <- function(x){
            m <- sapply(1:8, function(i) cm[x[i],,i]) %>% t()
            y <- sum(apply(m, 2, function(p) p[8] * prod(1 - p[1:7])))
            f <- m[8, 61L] * prod(1 - m[1:7, 61L]) # prob switched to fallow (code 61) 
            r <- (1 - m[8, 61L]) * prod(m[1:7, 61L]) # prob all ref yrs were fallow and eval wasn't
            y - f - r
      }
      apply(crops[, 1:8], 1, ps)
}

# if margin = "user", return prob of a switch truly having happened given measured crop data
# if margin = "prod", return prob of a switch being measured assuming crop data are true
swp <- function(state, state_code, start_year, margin = "user", introduce = T){
      
      message(paste(state, ifelse(introduce, start_year + 7, start_year), 
                    margin, ifelse(introduce, "introduction", "discontinuation")))
      
      # crop data
      d <- crp[crp[,3] == as.integer(state_code), (5:12) + start_year - 2008]
      if(nrow(d) == 0) return(tibble())
      if(!introduce) d <- d[, 8:1]
      
      # state-level confusion data
      m <- matrices[as.character(start_year:(start_year + 7))] %>% map(state)
      a <- array(dim = c(dim(m[[1]]), length(m)))
      for(i in 1:length(m)) a[,,i] <- conf_to_prob(m[[i]], match(margin, c("user", "producer")))
      
      # fill state-level holes with national probabilities
      mn <- matrices_national[as.character(start_year:(start_year + 7))]
      an <- array(dim = c(dim(mn[[1]]), length(mn)))
      for(i in 1:length(mn)) an[,,i] <- conf_to_prob(mn[[i]], match(margin, c("user", "producer")))
      for(i in 1:(dim(a)[1])) for(j in 1:(dim(a)[3])) if(round(sum(a[i,,j]), 10) == 0) a[i,,j] <- an[i,,j] # fill
      
      # switching probability
      p <- switch_prob(d, a)
      
      # observed switch
      switch <- function(x) ! x[8] %in% c(x[1:7], 61L)
      s <- apply(d, 1, switch)
      
      # results
      i <- crp[crp[,3] == as.integer(state_code), 1:4]
      tibble(x = i[,1],
             y = i[,2],
             region = i[,4],
             state = state,
             year = ifelse(introduce, start_year + 7, start_year),
             mode = ifelse(introduce, "introduction", "discontinuation"),
             margin = margin,
             switch = s,
             prob = p)
}

if(from_scratch){
      d <- expand_grid(states, 
                       start_year = (2008:2015), 
                       margin = c("user", "producer"), 
                       introduce = c(T, F)) %>%
            sample_n(nrow(.)) %>%
            future_pmap_dfr(swp) %>%
            mutate(switch = ifelse(switch, "switch", "no switch"))
      write_csv(d, "big_data/confusion_results.csv")
}else{
      d <- read_csv("big_data/confusion_results.csv")
}

d <- d %>% mutate(p = ifelse(switch == "switch", prob, 1-prob)) %>%
      bind_rows(mutate(., switch = "overall"), .) %>%
      mutate(switch = factor(switch, levels = c("switch", "no switch", "overall")),
             margin = paste0(margin, "'s accuracy"),
             margin = factor(margin, levels = c("user's accuracy", "producer's accuracy")),
             mode = paste0("crop ", mode))

means <- d %>%
      group_by(switch, mode, margin) %>%
      summarize(mean = mean(p),
                median = median(p), 
                misclass = mean(p < .5)) %>% # prop pixels that were probably misclassified
      gather(stat, p, mean, median, misclass)

p <- ggplot(d, aes(p, fill = switch, color = switch)) +
      facet_grid(margin + mode ~ .) +
      geom_density(data = filter(d, switch == "overall"), alpha = .3) +
      geom_density(data = filter(d, switch != "overall"), alpha = .3) +
      geom_vline(data = means %>% filter(stat != "misclass"), 
                 aes(xintercept = p, color = switch, linetype = stat)) +
      scale_fill_manual(values = c("black", "red", "dodgerblue")) +
      scale_color_manual(values = c("black", "red", "dodgerblue")) +
      scale_x_continuous(breaks = seq(0, 1, .1), limits = 0:1, expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
      theme_bw() +
      theme(strip.text = element_text(color = "white"),
            strip.background = element_rect(fill = "black"),
            axis.text.y = element_blank(),
            axis.ticks.y = ) +
      labs(x = "accuracy",
           color = "category",
           fill = "category",
           linetype = "summary statistic",
           y = "relative frequency")
ggsave("figures/fig_s8.pdf", p, width = 8, height = 7, units = "in")


reg_means <- d %>%
      group_by(state, switch, mode, margin) %>%
      mutate(p = ifelse(switch == "switch", prob, 1-prob)) %>%
      summarize(mean = mean(p),
                median = median(p),
                n = n()) %>%
      mutate(level = "state")
regkey <- levels(region)[[1]]  
reg_means <- d %>%
      mutate(state = regkey$region[regkey$ID[region+1]+1]) %>%
      group_by(state, switch, mode, margin) %>%
      mutate(p = ifelse(switch == "switch", prob, 1-prob)) %>%
      summarize(mean = mean(p),
                median = median(p),
                n = n()) %>%
      mutate(level = "region") %>%
      bind_rows(reg_means) %>%
      gather(stat, p, mean, median) %>%
      group_by(state) %>%
      mutate(pm = mean(p)) %>% ungroup() %>%
      arrange(pm) %>%
      mutate(state = factor(state, levels = unique(state)))

p <- ggplot() +
      facet_grid(level ~ margin + mode, scales = "free_y", space = "free_y") +
      geom_point(data = reg_means, aes(p, state, color = switch, shape = stat)) +
      scale_x_continuous(breaks = seq(0, 1, .1)) +
      scale_color_manual(values = c("red", "dodgerblue", "black")) +
      scale_shape_manual(values = c(19, 1)) +
      theme_bw() +
      theme(strip.text = element_text(color = "white"),
            strip.background = element_rect(fill = "black")) +
      labs(x = "accuracy (probability that result was correct)",
           color = "result",
           y = NULL)
ggsave("figures/fig_s9.pdf", p, width = 8, height = 8, units = "in")


write_csv(means, "data/confusion_results.csv")
write_csv(reg_means, "data/confusion_results_regional.csv")


# summaries
d %>%
      group_by(switch, mode, margin) %>%
      summarize(mean = mean(p < .5)) %>%
      arrange(margin, switch, mode)
means %>%
      filter(stat == "mean") %>%
      group_by(switch, margin, stat) %>%
      summarize(p = mean(p))

