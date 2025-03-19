
# Set Up -----------------------------------------------------------
# load packages
library(tidyverse)
library(ggplot2)
library(kohonen)

library(rnaturalearth)
library(ncdf4)
library(RColorBrewer)
library(reshape2)
library(here)

# set seed
set.seed(2025)
# land object for plotting
load(file = here("data/land.Rdata"))
life_history_info <- read_csv(here("data/species_life_history_info.csv"))
# Functions ---------------------------------------------------------------
# adapted plotting function
makeCovarianceMap <- function(input.season = input.season, variable, cmisst = cmisst, caption, ocean_variable) {
  # Covariance Map
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  myTitle <- switch(input.season,
                    "win" = paste("Winter", ocean_variable),
                    "sp" = paste("Spring", ocean_variable),
                    "sum" = paste("Summer", ocean_variable),
                    "aut" = paste("Autumn", ocean_variable))
  covMap<-cmisst
  limits <- c(-1, 1)
  extent <- c(25, 50, 225, 246) # min, max of lat, long
  
  
  gg <- ggplot() + ggtitle(myTitle) +
    geom_raster(data = melt(covMap), aes(x = Var1, y = Var2, fill=value)) +
    geom_sf(data=land, color="black", fill="grey", linewidth=0.25) +
    xlim(extent[3], extent[4]) + ylim(extent[1], extent[2]) +
    scale_fill_gradientn(colours = myPalette(100),limits=limits,name="Correlation", na.value = "white") +
    theme_classic() + theme(panel.border = element_rect(colour = "grey", fill=NA)) +
    labs(x = "Longitude", y = "Latitude", caption = paste0("In node:", "\n", caption))
  gg
}


# create lat/long df for cmisst data
makeLatLongDF <- function(grid.df, data.df, season){
  for(i in 1:length(data.df)){
    # pull out individual lists
    df1 <- data.df[[i]]
    
    # species/assessment specific 
    for(j in 1:length(df1)){
      spp.df <- df1[[j]][[season]]
      species <- names(df1)[j]
      
      spp.df <- as.data.frame(spp.df)
      spp.df <- spp.df %>% 
        rownames_to_column(var = "long") %>% 
        pivot_longer(!long, names_to = "lat", values_to = paste0(species)) %>% 
        mutate(lat_long = paste(lat, long, sep = "_")) %>% 
        select(lat_long, paste0(species))
      
      
      grid.df <- left_join(grid.df, spp.df)
    }
    
  }
  
  # remove empty cells
  grid.df <- grid.df %>% 
    filter(!if_any(-lat_long, is.na))
  
  # create matrix
  cor.matrix <- as.matrix(grid.df[2:ncol(grid.df)])
  rownames(cor.matrix) <- grid.df$lat_long
  cor.matrix <- t(cor.matrix)
  
  # return grid.df and cor.matrix
  return(list(grid.df, cor.matrix))
}


getNodeAvgs <- function(unit.classif, cor.df){
  
  # unit classification dataframe
  classif_df <- data.frame(classif = unit.classif, species = colnames(cor.df)[-1])
  
  # calculate node averages
  node_averages <- list()
  for(i in 1:max(classif_df$classif)){
    sp <- classif_df %>% 
      filter(classif == i)
    
    c <- which(colnames(cor.df) %in% sp$species)
    df <- cor.df[,c(1, c)]
    
    df$avg <- apply(df[,2:ncol(df)], 1, mean)
    
    df <- df %>% 
      select(lat_long, avg)
    
    node_averages[[i]] <- df
  }
  
  # format node averages to lat/long matrix to plot covariance maps
  node_matrix <- list()
  for(i in 1:max(classif_df$classif)){
    df <- node_averages[[i]]
    
    ll <- crossing(lats, longs) %>% 
      mutate(lat_long = paste(lats, longs, sep = "_"))
    
    df <- left_join(df, ll) %>% 
      select(lats, longs, avg)
    
    ll_matrix <- matrix(NA, nrow = length(longs), ncol = length(lats), dimnames = list(longs, lats))
    
    for(j in 1:nrow(df)){
      lat <- df$lats[j]
      long <- df$longs[j]
      
      x <- which(rownames(ll_matrix) == long)
      y <- which(colnames(ll_matrix) == lat)
      
      ll_matrix[x,y] <- df$avg[j]
    }
    
    node_matrix[[i]] <- ll_matrix
  }
  
  return(list(classif_df, node_averages, node_matrix))
  
}

SOMplots <- function(som, som_fp){
  # save plots
  png(paste0(som_fp, "/counts.png"))
  plot(som, type = "counts")
  dev.off()
  
  png(paste0(som_fp, "/codes.png"))
  plot(som, type = "codes")
  dev.off()
  
  png(paste0(som_fp, "/dist_neighbors.png"))
  plot(som, type = "dist.neighbours")
  dev.off()
  
}

NodeCovMaps <- function(som_fp, classif_df, node_matrix, season, variable){
  
  for(i in 1:max(classif_df$classif)){
    df <- node_matrix[[i]]
    
    spp_df <- classif_df %>% 
      filter(classif == i)
    
    cap <- paste(spp_df$species, collapse = "\n")
    
    # save as png
    makeCovarianceMap(input.season = season, cmisst = df, caption = cap, ocean_variable = variable)
    ggsave(paste0(som_fp, "/node_", i, ".png"), width = 4, height = 4)
  }
  
}


# ML  ---------------------------------------------------------------
## load data ---------------------------------------------------------------
recdev_cmisst_ml <- readRDS(here("data/processed/recdev_cmisst_ml.rds"))
yoy_cmisst_ml <- readRDS(here("data/processed/yoy_cmisst_ml.rds"))
yoy_index_cmisst_ml <- readRDS(here("data/processed/yoy_index_cmisst_ml.rds"))

# create grid for SOM
lats <- colnames(as.data.frame(yoy_cmisst_ml$bocaccio[[2]]))
longs <- rownames(as.data.frame(yoy_cmisst_ml$bocaccio[[2]]))

ml_grid <- crossing(lats, longs) %>% 
  mutate(lat_long = paste(lats, longs, sep = "_")) %>% 
  select(lat_long)

# create combined list for cmisst
# yoy ts series
names(yoy_cmisst_ml) <- c("bocaccio_core_index", "chilipepper_core_index", "halfbanded_core_index",
                          "shortbelly_core_index", "widow_core_index", "ytail_core_index")
# yoy index
names(yoy_index_cmisst_ml) <- c("canary_index", "chilipepper_index", "widow_index", "ytail_coast_index", "ytail_N_index")

# recruitment devs
names(recdev_cmisst_ml) <- c("black_WA_recdev", "black_CA_recdev", "blue_deacon_CA_recdev", "bocaccio_recdev", "cabezon_OR_recdev",
                             "cabezon_NCA_recdev", "cabezon_SCA_recdev", "california_scorpionfish_recdev", "chilipepper_recdev", "gopher_black_and_yellow_recdev",
                             "kelp_greenling_OR_recdev", "longspine_thornyhead_recdev", "pacific_ocean_perch_recdev", "rougheye_blackspotted_recdev",
                             "shortspine_thornyhead_recdev", "widow_recdev", "yelloweye_recdev", 
                             "ytail_N_recdev", "ytail_S_recdev", "black_OR_recdev")

ml_df <- list(yoy_cmisst_ml, yoy_index_cmisst_ml, recdev_cmisst_ml)


## ML Winter ---------------------------------------------------------------
# SOM grid and lat/long correlation df
ml_winter_cor_df <- makeLatLongDF(ml_grid, ml_df, 2)

### SOM ---------------------------------------------------------------------
#### 6 Node SOM --------------------------------------------------------------
som_grid_6_node <- somgrid(xdim = 3, ydim = 2, topo = "hexagonal")
som_6_node <- som(ml_winter_cor_df[[2]], grid = som_grid_6_node)

# plots
SOMplots(som_6_node, here("plots/soms/ml_winter/6_node_som"))

# Node averages
six_node_avg <- getNodeAvgs(som_6_node$unit.classif, ml_winter_cor_df[[1]])
NodeCovMaps(here("plots/soms/ml_winter/6_node_som"), six_node_avg[[1]], six_node_avg[[3]], season = "win", variable = "Mixed Layer Depth")


#### 4 Node SOM --------------------------------------------------------------
som_grid_4_node <- somgrid(xdim = 2, ydim = 2, topo = "hexagonal")
som_4_node <- som(ml_winter_cor_df[[2]], grid = som_grid_4_node)

# plots
SOMplots(som_4_node, here("plots/soms/ml_winter/4_node_som"))

# Node averages
four_node_avg <- getNodeAvgs(som_4_node$unit.classif, ml_winter_cor_df[[1]])
NodeCovMaps(here("plots/soms/ml_winter/4_node_som"), four_node_avg[[1]], four_node_avg[[3]], season = "win", variable = "Mixed Layer Depth")



### node qualities investigation --------------------------------------------
node6classif <- six_node_avg[[1]]
colnames(node6classif)[2] <- "index_stock"
node6classif <- left_join(node6classif, life_history_info)

node6classif %>%
  group_by(classif) %>% 
  summarise(avg_min_depth = mean(min_depth, na.rm = T),
            avg_max_depth = mean(max_depth, na.rm = T),
            avg_cog = mean(COG, na.rm = T),
            avg_peak_spawn = mean(peak_spawn, na.rm = T),
            avg_peak_spawn2 = mean(peak_spawn2, na.rm = T),
            avg_settle = mean(avg_settlement_age_mo, na.rm = T))


node4classif <- four_node_avg[[1]]
colnames(node4classif)[2] <- "index_stock"
node4classif <- left_join(node4classif, life_history_info)

node4classif %>%
  group_by(classif) %>% 
  summarise(avg_min_depth = mean(min_depth, na.rm = T),
            avg_max_depth = mean(max_depth, na.rm = T),
            avg_cog = mean(COG, na.rm = T),
            avg_peak_spawn = mean(peak_spawn, na.rm = T),
            avg_peak_spawn2 = mean(peak_spawn2, na.rm = T),
            avg_settle = mean(avg_settlement_age_mo, na.rm = T))

## ML Spring ---------------------------------------------------------------
# SOM grid and lat/long correlation df
ml_spring_cor_df <- makeLatLongDF(ml_grid, ml_df, 3)

### SOM ---------------------------------------------------------------------
#### 6 Node SOM --------------------------------------------------------------
som_grid_6_node <- somgrid(xdim = 3, ydim = 2, topo = "hexagonal")
som_6_node <- som(ml_spring_cor_df[[2]], grid = som_grid_6_node)

# plots
SOMplots(som_6_node, here("plots/soms/ml_spring/6_node_som"))

# Node averages
six_node_avg <- getNodeAvgs(som_6_node$unit.classif, ml_spring_cor_df[[1]])
NodeCovMaps(here("plots/soms/ml_spring/6_node_som"), six_node_avg[[1]], six_node_avg[[3]], season = "sp", variable = "Mixed Layer Depth")


#### 4 Node SOM --------------------------------------------------------------
som_grid_4_node <- somgrid(xdim = 2, ydim = 2, topo = "hexagonal")
som_4_node <- som(ml_spring_cor_df[[2]], grid = som_grid_4_node)

# plots
SOMplots(som_4_node, here("plots/soms/ml_spring/4_node_som"))

# Node averages
four_node_avg <- getNodeAvgs(som_4_node$unit.classif, ml_spring_cor_df[[1]])
NodeCovMaps(here("plots/soms/ml_spring/4_node_som"), four_node_avg[[1]], four_node_avg[[3]], season = "sp", variable = "Mixed Layer Depth")




