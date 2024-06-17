

library(terra)
library(tidyverse)
library(patchwork)
library(ggforce)
library(colormap)
library(GGally)
library(sf)
library(sfdep)
library(spdep)



# state boundaries
md <- map_data("state")
bg <- geom_polygon(data = md, aes(long, lat, group = group), fill = "gray80")
sb <- geom_path(data = md, aes(long, lat, group = group), color = "white", alpha = .25)



# load and format aggregated rasters =============

d <- list.files("data/aggregated", full.names = T) %>% 
      map(function(x) rast(x) %>% setNames(paste0(str_remove(basename(x), "\\.tif"), "__", 
                                                  1:nlyr(.), "_", names(.)))) %>%
      rast() %>%
      as.data.frame(xy = T) %>%
      as_tibble() %>%
      gather(var, value, -x, -y) %>%
      separate(var, c("var", "layer"), sep = "__")

d <- d %>%
      group_by(x, y, var) %>%
      mutate(n = length(layer)) %>%
      ungroup()

pw <- d %>% filter(str_detect(var, "pairwise"))
d <- d %>% filter(!str_detect(var, "pairwise|reclassified"))

ts <- d %>% filter(n == 8) %>%
      group_by(x, y, var) %>%
      mutate(year = case_when(str_detect(var, "abandon") ~ 2008:2015,
                              str_detect(var, "introduce") ~ 2015:2022)) %>%
      ungroup() %>%
      select(-n, -layer) %>%
      mutate(var = str_remove(var, "_8e_7r"))

conv <- ts %>%
      filter(str_detect(var, "convergence_")) %>%
      separate(var, c("level", "var", "phenom"), sep = "_") %>%
      select(-var) %>%
      rename(pc = value)

swch <- ts %>%
      filter(!str_detect(var, "convergence")) %>%
      separate(var, c("level", "phenom", "window", "outcome"), sep = "_")

swch <- swch %>% 
      unite(var, level, outcome) %>%
      spread(var, value) %>%
      group_by(phenom, window) %>%
      mutate(total = group_stasis + group_switch,
             t = group_switch - fallow_switch, # transformative
             i = crop_switch - group_switch, # incremental
             f = fallow_switch, # fallow
             pt = t/total,
             pi = i/total,
             pf = f/total) %>%
      filter(is.finite(total), total > 0) %>%
      left_join(conv)

y <- swch %>%
      group_by(phenom, window, year) %>%
      summarize(c = sum(pc*(t+i), na.rm = T),
                t = sum(t),
                i = sum(i),
                s = t + i,
                d = s - c,
                f = sum(f),
                total = sum(total),
                .groups = "drop") %>%
      mutate(pf = f/total,
             pt = t/total,
             pi = i/total,
             ps = s/total,
             pc = c/total,
             pd = d/total) %>%
      select(phenom, window, year, pt, pi, ps, pc, pd) %>%
      gather(stat, value, -year, -phenom, -window) %>%
      mutate(stat = factor(stat, levels = c("ps", "pi", "pt", "pc", "pd"),
                           labels = c("total switch", 
                                      "incremental switch", "transformative switch", 
                                      "convergent switch", "divergent switch")))

m <- swch %>%
      filter(window == "moving") %>%
      group_by(x, y) %>%
      summarize(a = sum(t[phenom == "abandon"]) + sum(i[phenom == "abandon"]),
                c = sum(pc*(t+i), na.rm = T),
                t = sum(t),
                i = sum(i),
                s = t + i,
                d = s - c,
                n = s - a,
                f = sum(f),
                total = sum(total),
                .groups = "drop") %>%
      mutate(pf = f/total,
             pt = t/total,
             pi = i/total,
             ps = s/total,
             pc = c/total,
             pd = d/total,
             pa = a/total,
             pn = n/total,
             
             ps_trans = pt/ps,
             ps_diverge = pd/ps,
             ps_intro = pn/ps
      ) %>%
      select(x, y, total, pc_switch = ps, 
             ps_trans, ps_diverge, ps_intro) %>%
      gather(stat, value, -x, -y, -total)

div <- d %>%
      filter(str_detect(var, "diversity"))


# time series plot (figure 2) ==========

p_ts <- y %>%
      mutate(phenom = ifelse(phenom == "abandon", "discontinuation", "introduction")) %>%
      filter(window == "moving") %>%
      ggplot(aes(year, value, 
                 color = stat,
                 group = paste(stat, phenom, window))) +
      scale_x_continuous(breaks = c(2008, 2015, 2022),
                         minor_breaks = unique(y$year)) +
      scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
      scale_color_manual(values = c("black", "firebrick2", "firebrick4", "dodgerblue", "dodgerblue4")) +
      scale_linetype_manual(values = c("twodash", "solid")) +
      geom_vline(xintercept = 2015, linewidth = .25) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      theme(axis.title.x = element_blank(),
            legend.title = element_blank()) +
      labs(y = "propotion of cropland nationwide")
ggsave("figures/fig_2.pdf", p_ts,
       width = 5, height = 3, units = "in")


# facet maps =========

### individual maps ==========

p <- m %>%
      ungroup() %>%
      mutate(stat = case_when(stat == "pc_switch" ~ "Prop cropland\nswitched",
                              stat == "ps_diverge" ~ "Prop switches\ndivergent",
                              stat == "ps_intro" ~ "Prop switches\nintroductions",
                              stat == "ps_trans" ~ "Prop switches\ntransformative")) %>%
      filter(is.finite(value)) %>%
      spread(stat, value) %>% 
      ggpairs(columns = 4:7,
              upper = list(continuous = wrap("cor", method = "spearman", size = 2.5, hjust = 0.7))) +
      theme_bw() +
      theme(strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white"))
ggsave("figures/fig_s5.png", p, width = 7, height = 7, units = "in")  


p <- ggplot() +
      facet_wrap(~stat) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "black", size = .25, alpha = .5) +
      geom_raster(data = m %>%
                        mutate(stat = factor(stat,
                                             levels = c("pc_switch", "ps_trans", "ps_intro", "ps_diverge"),
                                             labels = c("(a)  proportion cropland switched",
                                                        "(b)  proportion switches transformative",
                                                        "(c)  proportion switches introductions",
                                                        "(d)  proportion switches divergent"))), 
                  aes(x, y, fill = value)) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "white", size = .25, alpha = .25) +
      scale_alpha_continuous(trans = "log10",
                             range = c(0, 1)) +
      guides(fill = guide_colorbar(barwidth = .5, barheight = 25)) +
      theme_void(base_size = 16) +
      theme(legend.position = "right",
            strip.text = element_text(face = "bold")) +
      labs(fill = NULL) + 
      scale_fill_gradientn(colors = c("gray90", "#8bc9b3", "#3fa17d", "#136b4b", "#013d27", "black"), 
                           values = c(0, .2, .4, .55, .7, 1),
                           na.value = "gray90")
ggsave("figures/fig_3abcd.png", p, width = 10, height = 6, units = "in")



### rgb map =============

ppc <- m %>%
      group_by(stat) %>%
      mutate(value = value > median(value, na.rm = T)) %>%
      spread(stat, value) %>%
      mutate(combo = paste(ps_diverge, ps_trans, ps_intro),
             color = case_when(combo == "TRUE TRUE TRUE" ~ "black",
                               combo == "TRUE TRUE FALSE" ~ "#a49312", # purple
                               combo == "FALSE TRUE TRUE" ~ "#01589a", # orange
                               combo == "TRUE FALSE TRUE" ~ "#851f01", # lighblue
                               combo == "TRUE FALSE FALSE" ~ "#ebde71", # blue
                               combo == "FALSE TRUE FALSE" ~ "#6fcce9", # green
                               combo == "FALSE FALSE TRUE" ~ "#f29596", # red
                               combo == "FALSE FALSE FALSE" ~ "#cccccc"))
map <- ggplot() +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "black", linewidth = .25, alpha = .5) +
      geom_raster(data = ppc, aes(x, y, fill = color)) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "white", linewidth = .25, alpha = .25) +
      scale_fill_identity() +
      theme_void() +
      coord_fixed(ratio = 1.4)
ggsave("figures/fig_3e.png", map, width = 9, height = 6, units = "in")





### hot spot analyses ==================

# convert cell centers to polygons
s <- m %>%
      spread(stat, value) %>% 
      filter(is.finite(ps_diverge)) %>%
      select(x, y, pc_switch:ps_trans) %>%
      as.matrix() %>%
      rast(type = "xyz") %>%
      as.polygons(dissolve = F) %>%
      st_as_sf()

# compute hotspots via Getis Ord GI*
hs <- function(x, nn){
      x %>%
            gather(stat, value, pc_switch:ps_trans) %>%
            group_by(stat) %>% 
            mutate(nb = include_self(st_knn(geometry, k = nn)),
                   wt = st_weights(nb)) %>%
            mutate(Gi = local_g_perm(value, nb, wt, nsim = 999)) %>%
            unnest(Gi) %>%
            select(-nb, -wt) %>%
            mutate(k = paste(nn, "nn"))
}
hotspots <- bind_rows(hs(s, 4),
                      hs(s, 8),
                      hs(s, 36)) %>%
      mutate(k = factor(k, levels = c("4 nn", "8 nn", "36 nn"))) %>%
      mutate(classification = case_when(gi > 0 & p_folded_sim <= 0.01 ~ "Very hot",
                                        gi > 0 & p_folded_sim <= 0.05 ~ "Hot",
                                        gi > 0 & p_folded_sim <= 0.1 ~ "Somewhat hot",
                                        gi < 0 & p_folded_sim <= 0.01 ~ "Very cold",
                                        gi < 0 & p_folded_sim <= 0.05 ~ "Cold",
                                        gi < 0 & p_folded_sim <= 0.1 ~ "Somewhat cold",
                                        TRUE ~ "Not significant"),
             classification = factor(classification,
                                     levels = c("Very hot", "Hot", "Somewhat hot", "Not significant",
                                                "Somewhat cold", "Cold", "Very cold")),
             stat = factor(stat,
                           levels = c("pc_switch", "ps_diverge", "ps_intro", "ps_trans"),
                           labels = c("% cropland\nswitched", "% switches\ndivergent", 
                                      "% switches\nintroductions", "% switches\ntransformative"))) 

p <- ggplot() +
      facet_grid(stat ~ k, switch = "y") +
      geom_sf(data = hotspots, aes(fill = classification),
              color = NA, lwd = 0.1) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "black", size = .25, alpha = .25) +
      scale_fill_manual(values = c("darkred", "tomato", "#f4a582", "gray85", "#92c5de", "#0571b0", "darkblue")) +
      theme_void() +
      theme(strip.text.y = element_text(angle = 90)) +
      labs(fill = "Hot Spot\nClassification")
ggsave("figures/fig_s3.png", 
       p, width = 10, height = 6, units = "in")



# predictor regressions ==============

pd <- div %>%
      filter(var == "crop_stdiversity") %>%
      select(x, y, diversity = value) %>%
      left_join(m, .) %>%
      rename(cropland = total) %>%
      mutate(cropland = cropland / 1e6 / 15) %>%
      gather(pred, pred_value, cropland, diversity)

# load farm resource regions from 
# https://www.ers.usda.gov/data-products/arms-farm-financial-and-crop-production-practices/documentation.aspx
regions <- readxl::read_xls("data/reglink.xls", skip = 2) %>%
      janitor::clean_names() %>%
      mutate(fips = str_pad(fips, 5, "left", 0))
key <- regions %>%
      rename(key = key_to_ers_regions) %>%
      select(key) %>%
      filter(!is.na(key)) %>%
      separate(key, c("ers_resource_region", "region"), sep = "=") %>%
      mutate(ers_resource_region = as.numeric(ers_resource_region))
regions <- left_join(regions, key)

# join with county shapefile, limit to CONUS, fill a hole
regions <- st_read("data/Counties/cb_2018_us_county_500k.shp") %>%
      mutate(fips = paste0(STATEFP, COUNTYFP)) %>%
      left_join(regions) %>%
      filter(STATEFP <= 56, STATEFP != "02", STATEFP != "15") %>%
      mutate(region = ifelse(is.na(region), region[ers_resource_region == 3][1], region))
st_write(regions, "data/regions.shp")

# map 30 km cells to regions
xy <- swch %>% ungroup() %>% select(x, y) %>% distinct() 
xy$region_id <- xy %>%
      st_as_sf(coords = c("x", "y")) %>%
      st_set_crs(st_crs(regions))%>%
      st_within(regions) %>%
      map_int(function(x) x[1])
xy$region <- regions$region[xy$region_id]


rd <- pd %>%
      mutate(outcome = factor(stat, levels = c("pc_switch", "ps_intro", "ps_trans", "ps_diverge"),
                              labels = c("% cropland switched", 
                                         "% switches introductions", 
                                         "% switches transformative", 
                                         "% switches divergent")),
             predictor = factor(pred, levels = c("cropland", "diversity"),
                                labels = c("~ proportion landscape cropland", 
                                           "~ landscape crop diversity"))) %>%
      left_join(xy) %>%
      filter(!is.na(region)) %>%
      select(-predictor) %>%
      spread(pred, pred_value)


rds <- rd %>%
      group_by(outcome, region) %>%
      na.omit() %>%
      summarize(value = mean(value),
                cropland = mean(cropland),
                diversity = mean(diversity))
rdsw <- rd %>%
      group_by(outcome, region) %>%
      na.omit() %>%
      summarize(value = weighted.mean(value, cropland),
                diversity = weighted.mean(diversity, cropland),
                cropland = mean(cropland))

pal <- c('#332288','#88CCEE','#44AA99','#117733','#999933','#DDCC77','#CC6677','#882255','#AA4499')

reg <- regions %>% 
      group_by(region) %>% 
      summarise(m = sum(ALAND))
regmap <- ggplot(reg, aes(fill = region)) +
      geom_sf(color = "black") +
      scale_fill_manual(values = pal) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_void() +
      theme(legend.position = "top",
            legend.direction = "vertical") +
      labs(fill = "Farm resource region")


# predictor maps
pcult <- ggplot() +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "black", size = .25, alpha = .5) +
      geom_raster(data = rd, aes(x, y, fill = cropland)) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "white", size = .25, alpha = .25) +
      geom_sf(data = reg, color = "black", fill = NA) +
      scale_fill_gradientn(colors = c("black", "darkblue", "blue", "purple", "red", "orange", "gold"), 
                           values = c(0, .02, .1, .35, .6, .8, 1),
                           na.value = "black",
                           labels = scales::percent) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 5.2)) +
      theme_void() +
      theme(legend.direction = "vertical", 
            strip.text = element_blank()) +
      labs(fill = "% of landscape\ncultivated")
pdiv <- ggplot() +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "black", size = .25, alpha = .5) +
      geom_raster(data = rd, aes(x, y, fill = diversity)) +
      geom_path(data = md, aes(long, lat, group = group), 
                color = "white", size = .25, alpha = .25) +
      geom_sf(data = reg, color = "black", fill = NA) +
      scale_fill_gradientn(colors = c("black", "darkblue", "blue", "purple", "red", "orange", "gold"), 
                           #values = c(0, .02, .1, .35, .6, .8, 1),
                           na.value = "black") +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 5.2)) +
      theme_void() +
      labs(fill = "spatiotemporal\ncrop diversity")
p <- pcult / pdiv + 
      plot_annotation(tag_levels = 'a', 
                      tag_prefix = "(", tag_suffix = ")") 
ggsave("figures/fig_s2.png", 
       p, width = 7, height = 7, units = "in")



rdd <- rd %>% mutate(weight = cropland) %>% gather(pred, pred_value, cropland, diversity)
rdds <- rds %>% gather(pred, pred_value, cropland, diversity)
rddsw <- rdsw %>% mutate(weight = cropland) %>% gather(pred, pred_value, cropland, diversity)

strip <- function(prd = "cropland"){
      ggplot(mapping = aes(pred_value, value, group = region, color = region, fill = region)) +
            facet_wrap( ~ outcome, scales = "free", nrow = 1) +
            geom_smooth(data = rdd %>% filter(pred == prd),
                        linetype = "dotted",
                        aes(group = NULL), color = NA, fill = "black",
                        method = "glm", method.args = list(family = "binomial"),
                        alpha = .2, size = 1) +
            geom_smooth(data = rdd %>% filter(pred == prd),
                        method = "glm", method.args = list(family = "binomial"), se = F,
                        alpha = .2, size = .75) +
            geom_point(data = rdds %>% filter(pred == prd)) +
            geom_smooth(data = rdds %>% filter(pred == prd),
                        aes(group = NULL, color = NULL, fill = NULL),
                        color = "black",
                        method = "glm", method.args = list(family = "binomial"), se = F,
                        alpha = .2, size = 1.5) +
            scale_color_manual(values = pal) +
            scale_y_continuous(labels = scales::percent) +
            theme_bw(base_size = 11) +
            theme(strip.background = element_rect(fill = "black"),
                  strip.text = element_text(color = "white"),
                  strip.text.y = element_blank(),
                  legend.position = "none") +
            labs(x = ifelse(prd == "cropland",
                            "% landscape cropland",
                            "crop diversity"),
                 y = NULL)
}
p1 <- strip("diversity")
p2 <- strip("cropland") +
      scale_x_continuous(labels = scales::percent,
                         breaks = c(0, .5, 1))
p <- p1 / p2
ggsave("figures/fig_s6.pdf",
       p, width = 9, height = 6, units = "in")



# select relationships only
p1 <- ggplot(mapping = aes(pred_value, value, group = region, color = region, fill = region)) +
      geom_smooth(data = rdd %>% filter(pred == "diversity", outcome == "% cropland switched"), 
                  linetype = "dotted",
                  aes(group = NULL), color = NA, fill = "black",
                  method = "glm", method.args = list(family = "binomial"), 
                  alpha = .2, size = 1) +
      geom_smooth(data = rdd %>% filter(pred == "diversity", outcome == "% cropland switched"), 
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = .75) +
      geom_point(data = rdds %>% filter(pred == "diversity", outcome == "% cropland switched")) +
      geom_smooth(data = rdds %>% filter(pred == "diversity", outcome == "% cropland switched"), 
                  aes(group = NULL, color = NULL, fill = NULL),
                  color = "black",
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = 1.5) +
      scale_color_manual(values = pal) +
      scale_y_continuous(labels = scales::percent) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white"),
            strip.text.x = element_blank(),
            legend.position = "none") +
      labs(x = "crop diversity",
           y = "% cropland switched")
p2 <- ggplot(mapping = aes(pred_value, value, group = region, color = region, fill = region)) +
      geom_smooth(data = rdd %>% filter(pred == "cropland", outcome == "% cropland switched"), 
                  linetype = "dotted",
                  aes(group = NULL), color = NA, fill = "black",
                  method = "glm", method.args = list(family = "binomial"), 
                  alpha = .2, size = 1) +
      geom_smooth(data = rdd %>% filter(pred == "cropland", outcome == "% cropland switched"), 
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = .75) +
      geom_point(data = rdds %>% filter(pred == "cropland", outcome == "% cropland switched")) +
      geom_smooth(data = rdds %>% filter(pred == "cropland", outcome == "% cropland switched"), 
                  aes(group = NULL, color = NULL, fill = NULL),
                  color = "black",
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = 1.5) +
      scale_color_manual(values = pal) +
      scale_y_continuous(labels = scales::percent) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white"),
            strip.text.x = element_blank(),
            legend.position = "none") +
      labs(x = "% landscape cropland",
           y = "% cropland switched")
p3 <- ggplot(mapping = aes(pred_value, value, group = region, color = region, fill = region)) +
      geom_smooth(data = rdd %>% filter(pred == "cropland", outcome == "% switches divergent"), 
                  linetype = "dotted",
                  aes(group = NULL), color = NA, fill = "black",
                  method = "glm", method.args = list(family = "binomial"), 
                  alpha = .2, size = 1) +
      geom_smooth(data = rdd %>% filter(pred == "cropland", outcome == "% switches divergent"), 
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = .75) +
      geom_point(data = rdds %>% filter(pred == "cropland", outcome == "% switches divergent")) +
      geom_smooth(data = rdds %>% filter(pred == "cropland", outcome == "% switches divergent"), 
                  aes(group = NULL, color = NULL, fill = NULL),
                  color = "black",
                  method = "glm", method.args = list(family = "binomial"), se = F,
                  alpha = .2, size = 1.5) +
      scale_color_manual(values = pal) +
      scale_y_continuous(labels = scales::percent) +
      theme_bw(base_size = 11) +
      theme(strip.background = element_rect(fill = "black"),
            strip.text = element_text(color = "white"),
            strip.text.x = element_blank(),
            legend.position = "none") +
      labs(x = "% landscape cropland",
           y = "% switches divergent")
p <- p1 + p2 + p3 + 
      wrap_elements(plot = regmap) + 
      plot_layout(nrow = 1)
ggsave("figures/fig_4.pdf", 
       p, width = 9, height = 3.5, units = "in")







# spatial and temporal novelty and diversity =============

ai <- swch %>%
      ungroup() %>%
      filter(window == "moving") %>%
      group_by(x, y, phenom) %>%
      summarize(ps = mean(pt + pi), .groups = "drop") %>%
      spread(phenom, ps) %>%
      left_join(div %>% filter(str_detect(var, "crop_diversity")) %>% 
                      select(-n, -layer) %>% spread(var, value)) %>%
      left_join(select(m, x, y, total)) %>%
      mutate(ddiv = crop_diversity_2015_2022/total*8 - crop_diversity_2008_2015/total*8) %>% # delta temporal diversity (on all pixels, not just switched px)
      distinct()

spdiv <- div %>%
      filter(var == "crop_spdiversity") %>%
      separate(layer, c("layer", "year", "junk", "junk2"), sep = "_") %>%
      mutate(year = as.numeric(year)) %>%
      group_by(x, y) %>%
      summarize(t0 = mean(value[between(year, 2008, 2015)]),
                t1 = mean(value[between(year, 2015, 2022)]),
                diff = t1 - t0,
                ratio = t1 / t0)
spdiv <- m %>% 
      spread(stat, value) %>%
      left_join(spdiv)



pal <- c("#450078", "#9e67c7", "gray90", "#53ad5b", "#005407")
slim <- .75

sscat <- spdiv %>%
      mutate(ps_diverge = plyr::round_any(ps_diverge, .04),
             diff = plyr::round_any(diff, .06)) %>%
      group_by(ps_diverge, diff) %>%
      summarize(weight = sum(total * pc_switch)) %>%
      ggplot(aes(ps_diverge, diff, weight = weight, size = weight, color = diff)) +
      geom_point() +
      geom_smooth(method = lm, se = F, color = "black", size = .5) +
      scale_color_gradientn(colors = pal, limits = c(-slim, slim), values = c(0, .35, .5, .65, 1)) +
      scale_size_continuous(range = c(0, 6)) +
      coord_cartesian(ylim = c(-slim, slim)) +
      scale_x_continuous(labels = scales::percent) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(y = "spatial diversity trend",
           x = "% switches spatially novel (i.e. divergent)")
smap <- ggplot() +
      geom_raster(data = spdiv %>%
                        mutate(diff = case_when(diff > slim ~ slim,
                                                diff < -slim ~ -slim,
                                                TRUE ~ diff)), 
                  aes(x, y, fill = diff)) +
      geom_path(data = md, aes(long, lat, group = group), color = "black", alpha = .15) +
      scale_fill_gradientn(colors = pal, values = c(0, .35, .5, .65, 1)) +
      theme_void() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))

tlim <- .8
tscat <- ai %>% distinct() %>%
      mutate(ps_intro = introduce / (introduce + abandon),
             ps_intro = plyr::round_any(ps_intro, .04),
             diff = plyr::round_any(ddiv, tlim/10)) %>%
      group_by(ps_intro, diff) %>%
      summarize(weight = sum(total * (abandon + introduce))) %>%
      ggplot(aes(ps_intro, diff, weight = weight, size = weight, color = diff)) +
      geom_point() +
      geom_smooth(method = lm, se = F, color = "black", size = .5) +
      scale_color_gradientn(colors = pal, limits = c(-tlim, tlim), values = c(0, .35, .5, .65, 1)) +
      scale_size_continuous(range = c(0, 6)) +
      coord_cartesian(ylim = c(-tlim, tlim)) +
      scale_x_continuous(labels = scales::percent) +
      theme_minimal() +
      theme(legend.position = "none") +
      labs(y = "temporal diversity trend",
           x = "% switches temporally novel (i.e. introductions)")
tmap <- ggplot() +
      geom_raster(data = ai %>% 
                        mutate(ddiv = case_when(ddiv < -tlim ~ -tlim,
                                                ddiv > tlim ~ tlim,
                                                T ~ ddiv)), 
                  aes(x, y, fill = ddiv)) +
      geom_path(data = md, aes(long, lat, group = group), color = "black", alpha = .15) +
      scale_fill_gradientn(colors = pal, values = c(0, .35, .5, .65, 1)) +
      theme_void() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0))

p <- sscat + smap + tscat + tmap + 
      plot_layout(nrow = 2, widths = c(1, 1.7)) &
      theme(legend.position = "none")
ggsave("figures/fig_5.png", p,
       width = 8, height = 6, units = "in")




# pairwise =============

fct_crp <- function(x){
      factor(as.integer(x), 
             levels = c(301, 302, 303, 304, 305, 0),
             labels = c("field crops", "vegetables", "berries", "tree crops",
                        "fallow", "non-cultivated"))
}

x <- pw %>%
      split(.$var) %>%
      map_df(function(z){
            message(z$var[1])
            z %>%
                  separate(var, c("drop", "phenom", "year", "refyr"), sep = "_") %>%
                  separate(layer, c("drop2", "crops"), sep = "_") %>%
                  select(-n, -drop, -drop2) %>%
                  separate(crops, c("from", "to"), 3) %>%
                  filter(from != "000",
                         to != "000") %>%
                  filter(!(from == "305" & phenom == "abandon"),
                         !(to == "305" & phenom == "introduce"))
      })

n0 <- x %>% group_by(x, y, crop = from) %>% summarize(n0 = sum(value, na.rm = T), .groups = "drop")
n1 <- x %>% group_by(x, y, crop = to) %>% summarize(n1 = sum(value, na.rm = T), .groups = "drop")

f <- full_join(n0, n1) %>%
      group_by(x, y) %>%
      mutate(total = sum(n0)) %>%
      ungroup() %>%
      filter(total > 0)

notrans <- m %>% filter(stat == "ps_trans", value == 0)

p <- f %>%
      mutate(crop = fct_crp(crop)) %>%
      mutate(diff_tot = (n1-n0)/total) %>%
      ggplot(aes(x, y, fill = diff_tot)) +
      facet_wrap(~crop, nrow = 2) +
      geom_raster(data = notrans, aes(x, y), fill = "gray85") +
      geom_raster() +
      geom_polygon(data = md, aes(long, lat, group = group), 
                   fill = NA, color = "black", linewidth = .1) +
      scale_fill_gradientn(colors = c("darkred", "orange", "gray70", "dodgerblue", "darkblue"),
                           values = c(0, .4, .5, .6, 1),
                           limits = c(-1, 1)) +
      theme_void() +
      theme(
            legend.position = c(.825, .25),
            legend.direction = "horizontal",
            legend.title = element_text(size = 8),
            strip.text = element_text(size = 8,
                                      vjust = 1,
                                      face = "bold")) +
      guides(fill = guide_colorbar(barwidth = 10, barheight = .5,
                                   title.position="top", title.hjust = 0)) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(NA, max(f$y))) +
      labs(fill = "change in crop group cover as fraction of   \ntotal transformatively switched cropland")

sank <- function(ph){
      
      g <- x %>%
            filter(phenom == ph) %>%
            group_by(from, to) %>%
            summarize(n = sum(value), .groups = "drop") %>%
            mutate(from = fct_crp(from),
                   to = fct_crp(to)) %>%
            select(from, to, n) %>%
            gather_set_data(1:2)
      
      pal <- c("orange", "forestgreen", "red", "blue", "gray50")
      sankey <- ggplot(g %>% filter(from != to), 
                       aes(x, id = id, split = y, value = n)) +
            facet_wrap(~ base::switch(ph,
                                      "abandon" = "discontinued:                                                                            post-switch:  ",
                                      "introduce" = "  pre-switch:                                                                             introduced:")
            ) +
            geom_parallel_sets(aes(fill = from),
                               alpha = 0.6, axis.width = 0.03, sep = .1) +
            geom_parallel_sets_labels(angle = 0, lineheight = .7, 
                                      hjust = base::switch(ph,
                                                           "abandon" = c(rep(1, 4), rep(0, 5)),
                                                           "introduce" =c(rep(1, 5), rep(0, 4))),
                                      color = base::switch(ph,
                                                           "abandon" = c(rev(pal[1:4]), rev(pal)),
                                                           "introduce" = c(rev(pal), rev(pal[1:4]))),
                                      size = 2.5, sep = .1) +
            xlim(.8, 2.2) +
            theme_void() +
            theme(legend.position = "none",
                  plot.title = element_text(size = 8, face = "bold", hjust = .5),
                  strip.text = element_text(size = 8, face = "italic"),
                  plot.margin = unit(c(.05, 0, 0, 0), "npc")) +
            scale_fill_manual(values = pal)
      
      sankey
}

sa <- sank("abandon") + labs(title = "crop discontinuation")
si <- sank("introduce") + labs(title = "crop introduction")

ps <- (sa | si) / p +
      plot_layout(heights = c(1, 2.5)) +
      plot_annotation(tag_levels = 'a', 
                      tag_prefix = "(", tag_suffix = ")")
ggsave("figures/fig_6.png", ps,
       width = 8, height = 5, units = "in")




# trans-div odds ratios ==================

cmbi <- map(2015:2022, function(x) rast(paste0("data/aggregated/tidc_combo_introduce_", x, ".tif"))) %>%
      Reduce("+", .)
cmba <- map(2008:2015, function(x) rast(paste0("data/aggregated/tidc_combo_abandon_", x, ".tif"))) %>%
      Reduce("+", .)
cmb <- cmbi + cmba

pd <- as.data.frame(cmb, xy = T) %>%
      filter((incr_div + incr_conv + trans_div + trans_conv) > 0) %>%
      mutate_all(function(z) ifelse(z == 0, 1, z)) %>%
      mutate(odds_ratio = (trans_div / incr_div) / (trans_conv / incr_conv),
             odds_ratio = case_when(odds_ratio > 10 ~ 10,
                                    odds_ratio < .1 ~ .1,
                                    TRUE ~ odds_ratio))

p <- ggplot() +
      geom_raster(data = pd, aes(x, y, fill = odds_ratio)) +
      geom_path(data = md, aes(long, lat, group = group), color = "black", alpha = .25) +
      scale_fill_gradientn(trans = "log10",
                           colors = c("darkred", "tomato", "#f4a582", "gray85", "#92c5de", "#0571b0", "darkblue") %>% rev(),
                           breaks = c(.1, 1, 10),
                           labels = c("<=.1", 1, ">= 10")) +
      theme_void() +
      labs(fill = "odds ratio\n")
ggsave("figures/fig_s4.png", p, width = 8, height = 4, units = "in")



# reclassification ============

rd <- rast("data/aggregated/reclassified.tif") %>% as.data.frame(xy = T) %>% rename(p = sum)
p <- ggplot(rd, aes(x, y, fill = p)) + geom_raster() + 
      scale_fill_viridis_c(na.value = NA, trans = "sqrt") +
      guides(fill = guide_colorbar(barwidth = 12)) +
      theme_void() +
      theme(legend.position = "bottom") +
      labs(fill = "proportion of data reclassified  ")
ggsave("figures/fig_s7.png", p, width = 8, height = 5, units = "in")




# summary statistics =================

## crop switching prevalence ##

# total cropland area in analysis
swch %>%
      filter(window == "moving", phenom == "abandon", year == 2008) %>%
      summarize(total = sum(total) * 900 / 1e6)

# correlation between predictors
cor(rd$cropland, rd$diversity)

z <- y %>%
      group_by(phenom, window, stat) %>%
      summarize(value = mean(value))

# percent cropland switching
round(z$value[6], 3) # share of cropland with discontinuation in mean year
round(z$value[16], 3) # share of cropland with introduction in mean year


# switching rates by region
rd %>%
      filter(stat == "pc_switch") %>%
      group_by(stat, outcome, region) %>%
      summarize(value = weighted.mean(value, cropland)) %>%
      arrange(desc(value))

# pixels that switched at least once
bai1 <- swch %>% ungroup() %>%
      select(x, y, total) %>%
      distinct() %>%
      left_join(filter(d, str_detect(var, "switchONCE"))) %>%
      na.omit() %>%
      group_by(var) %>%
      summarize(p = sum(value) / sum(total))
bai1




## temporal novelty ##

# share of all crop switching that's introductions
z$value[16] / (z$value[6] + z$value[16])

# % landscapes with >2x diff between intro and discon
rd %>%
      filter(stat == "ps_intro",
             is.finite(value)) %>%
      summarize(p = 1 - mean(between(value, 1/3, 2/3)))
ai %>% mutate(r = 1 - between(introduce/abandon, .5, 2)) %>%
      summarize(r = mean(r))

# diversity trend vs intro-discon diff
cor(ai$ddiv, ai$introduce-ai$abandon, method = "spearman")


## categorical novelty ##

zm <- filter(z, window == "moving") %>%
      spread(stat, value) %>%
      janitor::clean_names() %>%
      mutate(incremental_switch = incremental_switch / total_switch,
             transformative_switch = transformative_switch / total_switch,
             convergent_switch = convergent_switch / total_switch,
             divergent_switch = divergent_switch / total_switch)

# percent switches transformative, for dicson an intro
select(zm, phenom, transformative_switch)




## spatial novelty ##

# share of all crop switching that's divergent
(z$value[10] + z$value[20]) / (z$value[6] + z$value[16])

# divergence rates by region
rd %>%
      filter(stat == "ps_diverge") %>%
      group_by(stat, outcome, region) %>%
      summarize(value = weighted.mean(value, cropland, na.rm = T)) %>%
      arrange(desc(value))

# correlate divergence vs diversity trends
cor(spdiv$ps_diverge, spdiv$diff, method = "spearman", use = "pairwise.complete.obs")




## relationships among novelty dimensions ##

# correlations, landscape scale
w <- m %>% spread(stat, value) %>% select(pc_switch:ps_trans) %>% na.omit()
cov.wt(select(w, pc_switch:ps_trans), wt = rep(1, nrow(w)), cor = TRUE)$cor
cor(w, method = "spearman")
cor(w, method = "spearman") ^ 2
cor(w) ^ 2

# odds ratio, spatial/temporal
log((zm$divergent_switch[2] / zm$convergent_switch[2]) / (zm$divergent_switch[1] / zm$convergent_switch[1]))

# odds ratio, categorical/temporal
log((zm$transformative_switch[2] / zm$incremental_switch[2]) / (zm$transformative_switch[1] / zm$incremental_switch[1]))

# odds ratio, categorical/spatial
as.data.frame(cmb) %>%
      summarize_all(sum) %>%
      summarize(lor = log((trans_div / incr_div) / (trans_conv / incr_conv)),
                lor2 = log((trans_div / trans_conv) / (incr_div / incr_conv)))

