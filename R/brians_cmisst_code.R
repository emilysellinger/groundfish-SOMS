library(dplyr)
library(ggplot2)
library(readr)
library(rnaturalearth)
library(ncdf4)
library(RColorBrewer)
library(reshape2)
library(here)

# Index function ----------------------------------------------------------
get_CMISST_index <- function(response, oceanData=oceanData_ERSST,
                             years=NA, years.fit=years.fit,
                             months=1:12,
                             min.lon=225, max.lon=246,
                             min.lat=25, max.lat=50) {
  # Verify that the response is what we expect
  if (ncol(response)!=2) { print("incorrect data - requires a 2-column data frame with year and the response"); return(NA) }
  colnames(response)<-c("year","val")
  
  # 'years' will be considered 'all years'  If we need fit or pred, we can access them
  year_mo<-data.frame(year=rep(years, each=length(months)), month=rep(months, length(years)),
                      label=paste(rep(years, each=length(months)), rep(months, length(years)), sep = "_"))
  
  #***************************************************************
  # Extract and scale the spatial data
  #***************************************************************
  
  # Index the locations in the file
  lons <- as.numeric(dimnames(oceanData)[[1]])
  lats <- as.numeric(dimnames(oceanData)[[2]])
  yr_mo <- dimnames(oceanData)[[3]]
  lon.index<-which(lons >= min.lon & lons <= max.lon) 
  lat.index<-which(lats >= min.lat & lats <= max.lat)
  yr_mo.index<-which(yr_mo %in% year_mo$label)
  # Subset the ocean data with user-defined extent
  oceanData <- oceanData[lon.index, lat.index, yr_mo.index]
  
  # Create the function to calculate seasonal averages
  createSeasonalData<-function(oceanData,
                               years = years, months = months, year_mo=year_mo, season=1) {
    seasonal<-array(NA, c(dim(oceanData)[1], dim(oceanData)[2], length(years)), dimnames = list(dimnames(oceanData)[[1]], dimnames(oceanData)[[2]], years))
    for (yy in 1:length(years)) {
      if (season==1) seasonal[,,yy]<-(oceanData[,,year_mo$month == 1 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 2 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 3 & year_mo$year==years[yy]])/3
      if (season==2) seasonal[,,yy]<-(oceanData[,,year_mo$month == 4 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 5 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 6 & year_mo$year==years[yy]])/3
      if (season==3) seasonal[,,yy]<-(oceanData[,,year_mo$month == 7 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 8 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 9 & year_mo$year==years[yy]])/3
      if (season==4) seasonal[,,yy]<-(oceanData[,,year_mo$month == 10 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 11 & year_mo$year==years[yy]]+oceanData[,,year_mo$month == 12 & year_mo$year==years[yy]])/3
    } 
    # This scales (Z-score) the data cell-wise
    # The aperm is needed because for some reason the apply function returns the third dimension (time) as the first dimension
    oceanData.scl <- aperm(apply(seasonal, 1:2, scale), c(2,3,1))
    dimnames(oceanData.scl)[[3]]<-years
    return(oceanData.scl)
  }
  
  # Create the data by calling our function (returns scaled sst array and full dataset with fish)
  oceanData.s1.scl <- createSeasonalData(oceanData = oceanData, years = years, months = months, year_mo=year_mo, season = 1)
  oceanData.s2.scl <- createSeasonalData(oceanData = oceanData, years = years, months = months, year_mo=year_mo, season = 2)
  oceanData.s3.scl <- createSeasonalData(oceanData = oceanData, years = years, months = months, year_mo=year_mo, season = 3)
  oceanData.s4.scl <- createSeasonalData(oceanData = oceanData, years = years, months = months, year_mo=year_mo, season = 4)
  
  # Get covariance between each cell's temperature and survival (only for fit years!!)
  covs1<-apply(oceanData.s1.scl[,,as.character(years.fit)], 1:2, function(x) cor(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs2<-apply(oceanData.s2.scl[,,as.character(years.fit)], 1:2, function(x) cor(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs3<-apply(oceanData.s3.scl[,,as.character(years.fit)], 1:2, function(x) cor(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs4<-apply(oceanData.s4.scl[,,as.character(years.fit)], 1:2, function(x) cor(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  
  #********************************************************************
  # Create the index (how similar is each year to the covariance map)
  #********************************************************************
  coefs_cov<-NULL
  options(na.action="na.omit")
  for (tt in 1:dim(oceanData.s1.scl)[3])
    coefs_cov<-rbind(coefs_cov, c(lm(as.vector(oceanData.s1.scl[,,tt]) ~ -1 + as.vector(covs1))$coef,
                                  lm(as.vector(oceanData.s2.scl[,,tt]) ~ -1 + as.vector(covs2))$coef,
                                  lm(as.vector(oceanData.s3.scl[,,tt]) ~ -1 + as.vector(covs3))$coef,
                                  lm(as.vector(oceanData.s4.scl[,,tt]) ~ -1 + as.vector(covs4))$coef))
  coefs_cov<-data.frame(coefs_cov)
  coefs_cov$year<-years
  index_cov<-merge(coefs_cov, response[response$year %in% years.fit,], all.x=TRUE)
  colnames(index_cov)<-c("year","win.cov","spr.cov","sum.cov","aut.cov","val")
  
  # Returns index as a list
  #  cmisst[[1]] contains 6 columns (as one list item): year, 4 seasonal indices, response
  #  cmisst[[2]] winter spatial covariance values (for maps)
  #  cmisst[[3]] spring spatial covariance values (for maps)
  #  cmisst[[4]] summer spatial covariance values (for maps)
  #  cmisst[[5]] autumn spatial covariance values (for maps)
  #  cmisst[[6]] lat, long min and max, as 1 list item
  return(list(index_cov, covs1, covs2, covs3, covs4, c(min.lat, max.lat, min.lon, max.lon)))
}




# Plots -------------------------------------------------------------------
# Make a few plots from the results

# To get back to normal space
reverse_scale <- function(x, center = NULL, scale = NULL) {
  if (!is.null(attr(x, "scaled:scale"))) {
    x <- x * attr(x, "scaled:scale")
  } else { x <- x * scale }
  if (!is.null(attr(x, "scaled:center"))) {
    x <- x + attr(x, "scaled:center")
  } else { x <- x + center }
  x
}


makeCovarianceMap <- function(input.season = input.season, cmisst = cmisst) {
  # Covariance Map
  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
  season <- switch(input.season,
                   win = 2,
                   spr = 3,
                   sum = 4,
                   aut = 5)
  myTitle <- switch(input.season,
                    win = "Winter",
                    spr = "Spring",
                    sum = "Summer",
                    aut = "Autumn")
  covMap<-cmisst[[season]]
  lmt<-max(abs(covMap), na.rm=TRUE)
  # if(lmt > .5){
  #   limits<-c(-lmt, lmt)
  # }else{
  #   limits <- c(-.5, .5)
  # }
  limits <- c(-1, 1)
  
  extent <- cmisst[[6]] # min, max of lat, long
  
  gg <- ggplot() + ggtitle(myTitle) +
    geom_raster(data = melt(covMap), aes(x = Var1, y = Var2, fill=value)) +
    geom_sf(data=land, color="black", fill="grey", linewidth=0.25) +
    xlim(extent[3], extent[4]) + ylim(extent[1], extent[2]) +
    scale_fill_gradientn(colours = myPalette(100),limits=limits,name="Correlation", na.value = "white") +
    theme_classic() + theme(panel.border = element_rect(colour = "grey", fill=NA)) +
    labs(x = "Longitude", y = "Latitude")
  gg
}


makeBiplot <- function(input.season = input.season, cmisst = cmisst) {
  # Biplot with response
  index <- cmisst[[1]]
  season <- switch(input.season,
                   win = 2,
                   spr = 3,
                   sum = 4,
                   aut = 5)
  index$ind <- index[,season]
  myTitle <- switch(input.season,
                    win = "Winter",
                    spr = "Spring",
                    sum = "Summer",
                    aut = "Autumn")
  plot(index$ind, index$val, pch=20, cex=2, xlab=paste(myTitle, "CMISST Index"),
       ylab="Scaled (Z-score) Response", main=myTitle)
  lm1 <- lm(index$val~index$ind)
  abline(lm1)
  text(bquote(~ R^2 == .(round(summary(lm1)$adj.r.squared, 2))),
       x = par("usr")[1]*0.8, y=par("usr")[4]*0.80, cex=1.6, col="blue")
  if (input.loocv) {
    mae <- cmisst[[7]]
    colnames(mae)[colnames(mae)=="mae.mean"] <- "mae"
    text(paste("MAE =", round(mean(abs(mae[mae$season==input.season,"mae"])), 2)),
         x = par("usr")[1]*0.75, y=par("usr")[4]*0.60, cex=1.6, col="blue")
  }
}

makeTimeSeriesPlot <- function(input.season = input.season, cmisst = cmisst,
                               ylab="", yaxis_scaler=1) {
  # Time series plot in normal space
  response.tmp <- response
  response.tmp$year <- response.tmp$year - as.numeric(input.lag)
  response.tmp <- response.tmp[response.tmp$year %in% seq(input.years[1], input.years[2], 1), c('year', input.stock)]
  colnames(response.tmp) <- c('year','val')
  if(input.log) response.tmp$val <- log(response.tmp$val)
  response.tmp$val.scl <- scale(response.tmp$val)
  #reverse_scale(response.tmp$val.scl)
  
  index <- cmisst[[1]]
  season <- switch(input.season,
                   win = 2, spr = 3, sum = 4, aut = 5)
  index$ind <- index[,season]
  index$counts <- reverse_scale(index$val, attr(response.tmp$val.scl, "scaled:center"), attr(response.tmp$val.scl, "scaled:scale"))
  if (input.log) index$counts <- exp(index$counts)
  myTitle <- switch(input.season,
                    win = "Winter", spr = "Spring", sum = "Summer", aut = "Autumn")
  lm1 <- lm(index$val~index$ind)
  preds<-predict(lm1, newdata = index, interval = "confidence")
  preds<-reverse_scale(preds, attr(response.tmp$val.scl, "scaled:center"), attr(response.tmp$val.scl, "scaled:scale"))
  if (input.log) preds<-exp(preds)
  # Use prediction interval for predicted points
  preds_new<-predict(lm1, newdata = index, interval = "prediction")
  preds_new<-reverse_scale(preds_new, attr(response.tmp$val.scl, "scaled:center"), attr(response.tmp$val.scl, "scaled:scale"))
  if (input.log) preds_new<- exp(preds_new)
  # replace just the ones that were not used during fitting
  #preds[index$year %in% input.years.pred,]<-preds_new[index$year %in% input.years.pred,]
  preds[is.na(index$val),]<-preds_new[is.na(index$val),]
  
  preds<-data.frame(preds)
  # unlag the year to show the plot in return year
  index$year_return <- index$year + input.lag
  preds$year_return <- index$year_return
  # Plot for SOEM talk in 2024
  ggplot() +
    geom_line(data = index, aes(x=year_return, y=counts/yaxis_scaler)) +
    geom_point(data = index, aes(x=year_return, y=counts/yaxis_scaler)) +
    theme_classic() +
    ylab(label = ylab) + xlab("Response Year") +
    geom_line(data=preds, aes(x=year_return, y=fit/yaxis_scaler), color="deepskyblue2", linewidth=1.3) +
    geom_point(data=preds, aes(x=year_return, y=fit/yaxis_scaler), color="deepskyblue2") +
    geom_ribbon(data=preds, aes(x=year_return, ymin = lwr/yaxis_scaler, ymax = upr/yaxis_scaler), fill = "deepskyblue2", alpha = 0.2)
}

makeIndexPlot <- function(cmisst = cmisst) {
  # Output: Index time series
  index <- cmisst[[1]]
  plot(index$year, index$win.cov, type='b', pch=20, col="red4",
       xlab="", ylab="CMISST Index",
       ylim=c(min(index[,c("win.cov","spr.cov","sum.cov","aut.cov")], na.rm=TRUE),
              max(index[,c("win.cov","spr.cov","sum.cov","aut.cov")], na.rm=TRUE)))
  points(index$year, index$spr.cov, type='b', pch=20, col="blue")
  points(index$year, index$sum.cov, type='b', pch=20, col="green3")
  points(index$year, index$aut.cov, type='b', pch=20, col="purple")
  legend("topleft", legend = c("Win","Spr","Sum","Aut"), bty='n',
         col = c("red4","blue","green3","purple"), pch = 20, lty=1)
}


makeLOOplot <- function(cmisst = cmisst, season = "spr") {
  # Output: Observed and predicted time series from the LOO
  index <- cmisst[[1]] # This gets us the whole time series
  plot(index$year, index$val, type='b', pch=20, cex=2, col="black", xlab="", ylab="Scaled Response", main = input.stock)
  abline(0,0, lty=2)
  index <- cmisst[[7]] # this is just the loo results
  index2<-index[index$season==season & index$model=="cmisst",]
  lines(index2$year, index2$pred, lwd=3, col="deepskyblue2")
  text(labels = paste("LOO MAE CMISST =", round(mean(abs(index2$mae)),2)),
       x = par("usr")[1]+9, y=par("usr")[4]*0.80, cex=1.0, col="deepskyblue2")
}

makeTable <- function(cmisst = cmisst) {
  # Time series plot in normal space
  response.tmp <- response
  response.tmp$year <- response.tmp$year - as.numeric(input.lag)
  response.tmp <- response.tmp[response.tmp$year %in% seq(input.years[1], input.years[2], 1), c('year', input.stock)]
  colnames(response.tmp) <- c('year','val')
  if(input.log) response.tmp$val <- log(response.tmp$val)
  response.tmp$val.scl <- scale(response.tmp$val)
  index <- cmisst[[1]]
  index$response <- reverse_scale(index$val, attr(response.tmp$val.scl, "scaled:center"), attr(response.tmp$val.scl, "scaled:scale"))
  if (input.log) index$response <- exp(index$response)
  
  # Output: Table
  out<-cmisst[[1]]
  out$year <- as.integer(out$year)#out <- out[,c(5,6,1:4)]
  out$response <- index$response
  out
}



# Ocean data --------------------------------------------------------------
load(file = here("data/oceanSSHData.RData"))
load(file = here("data/oceanSSTData.RData"))
load(file = here("data/land.Rdata"))
load(file = here("data/oceanMLData.RData"))

# RREAS data --------------------------------------------------------------
yoy_rockfish_all <- read_csv(here("data/rreas/rreas_juv_rockfish_all.csv"))

## SSH ---------------------------------------------------------------------
yoy_rockfish_all_ssh <- get_CMISST_index(response = yoy_rockfish_all, oceanData = oceanSSHData, years = seq(2001, max(yoy_rockfish_all$year)),
                           years.fit = seq(2001, max(yoy_rockfish_all$year)))


### Plots -------------------------------------------------------------------
input.loocv <- FALSE
# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_ssh)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SSH/winter/yoy_all_2001-2022.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SSH/winter/yoy_all.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_ssh)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_ssh)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SSH/spring/yoy_all.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SSH/spring/yoy_all.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_ssh)
dev.off()


### Lag Plots ---------------------------------------------------------------
# try lagged responses
yoy_rockfish_all$year_lagged <- yoy_rockfish_all$year - 1

yoy_rockfish_all_ssh_lagged <- get_CMISST_index(response = yoy_rockfish_all[,c(3,2)], oceanData = oceanSSHData, years = yoy_rockfish_all$year_lagged,
                           years.fit = yoy_rockfish_all$year_lagged)

# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_ssh_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SSH/winter/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SSH/winter/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_ssh_lagged)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_ssh_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SSH/spring/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SSH/spring/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_ssh_lagged)
dev.off()


## SST ---------------------------------------------------------------------
yoy_rockfish_all_sst <- get_CMISST_index(response = yoy_rockfish_all[,1:2], oceanData = oceanSSTData, years = yoy_rockfish_all$year,
                                         years.fit = yoy_rockfish_all$year)


### Plots -------------------------------------------------------------------
# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_sst)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SST/winter/yoy_all.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SST/winter/yoy_all.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_sst)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_sst)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SST/spring/yoy_all.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SST/spring/yoy_all.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_sst)
dev.off()


### Lag Plots ---------------------------------------------------------------
# try lagged responses
yoy_rockfish_all_sst_lagged <- get_CMISST_index(response = yoy_rockfish_all[,c(3,2)], oceanData = oceanSSTData, years = yoy_rockfish_all$year_lagged,
                                                years.fit = yoy_rockfish_all$year_lagged)

# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_sst_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SST/winter/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SST/winter/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_sst_lagged)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_sst_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/SST/spring/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/SST/spring/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_sst_lagged)
dev.off()



## Mix Layer depth ---------------------------------------------------------
yoy_rockfish_all_ml <- get_CMISST_index(response = yoy_rockfish_all[c(9:38),c(1,2)], oceanData = oceanMLData, years = yoy_rockfish_all$year[9:38],
                                         years.fit = yoy_rockfish_all$year[9:38])

### Plots -------------------------------------------------------------------
# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_ml)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/ML/winter/yoy_all.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/ML/winter/yoy_all.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_ml)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_ml)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/ML/spring/yoy_all.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/ML/spring/yoy_all.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_ml)
dev.off()


### Lag Plots ---------------------------------------------------------------
# try lagged responses
#response$year_lagged <- response$year - 1

yoy_rockfish_all_ml_lagged <- get_CMISST_index(response = yoy_rockfish_all[c(10:39),c(3,2)], oceanData = oceanMLData, years = yoy_rockfish_all$year_lagged[10:39],
                                                years.fit = yoy_rockfish_all$year_lagged[10:39])

# winter 
makeCovarianceMap(input.season = 1, cmisst = yoy_rockfish_all_ml_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/ML/winter/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/ML/winter/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 1, cmisst = yoy_rockfish_all_ml_lagged)
dev.off()

# spring
makeCovarianceMap(input.season = 2, cmisst = yoy_rockfish_all_ml_lagged)
ggsave(here("plots/rreas_yoy_ts/covariance_maps/ML/spring/yoy_all_lag_1yr.png"), width = 4, height = 4)

png(here("plots/rreas_yoy_ts/cmisst_index/ML/spring/yoy_all_lag_1yr.png"))
makeBiplot(input.season = 2, cmisst = yoy_rockfish_all_ml_lagged)
dev.off()




# RREAS Individual species ------------------------------------------------
yoy_rockfish_spp <- read_csv(here("data/rreas/rreas_juv_rockfish_by_species.csv"))

yoy_species <- c("bocaccio", "chilipepper", "halfbanded", "shortbelly", "widow", "yellowtail")



## SSH ---------------------------------------------------------------------
yoy_cmisst_ssh <- list()
for(i in 1:length(yoy_species)){
  
  df <- yoy_rockfish_spp[,c(1, (i+2))]
  
  cmisst <- get_CMISST_index(response = df, oceanData = oceanSSHData,
                             years = df$year, years.fit = df$year)
  
  yoy_cmisst_ssh[[i]] <- cmisst
}
names(yoy_cmisst_ssh) <- yoy_species

saveRDS(yoy_cmisst_ssh, file = here("data/processed/yoy_cmisst_ssh.rds"))

## SST ---------------------------------------------------------------------
yoy_cmisst_sst <- list()
for(i in 1:length(yoy_species)){
  
  df <- yoy_rockfish_spp[,c(1, (i+2))]
  
  cmisst <- get_CMISST_index(response = df, oceanData = oceanSSTData,
                             years = df$year, years.fit = df$year)
  
  yoy_cmisst_sst[[i]] <- cmisst
}
names(yoy_cmisst_sst) <- yoy_species

saveRDS(yoy_cmisst_sst, file = here("data/processed/yoy_cmisst_sst.rds"))



# Mixed layer -------------------------------------------------------------
yoy_cmisst_ml <- list()
for(i in 1:length(yoy_species)){
    
    df <- yoy_rockfish_spp[c(9:38),c(1, (i+2))]
    
    cmisst <- get_CMISST_index(response = df, oceanData = oceanMLData,
                               years = df$year, years.fit = df$year)
    
    yoy_cmisst_ml[[i]] <- cmisst
}

names(yoy_cmisst_ml) <- yoy_species

saveRDS(yoy_cmisst_ml, file = here("data/processed/yoy_cmisst_ml.rds"))



## Plots -------------------------------------------------------------------
# create file path for plots
cov_map_ssh_fp <- here("plots/rreas_yoy_ts/covariance_maps/SSH")
cov_map_sst_fp <- here("plots/rreas_yoy_ts/covariance_maps/SST")
cov_map_ml_fp <- here("plots/rreas_yoy_ts/covariance_maps/ML")
index_ssh_fp <- here("plots/rreas_yoy_ts/cmisst_index/SSH")
index_sst_fp <- here("plots/rreas_yoy_ts/cmisst_index/SST")
index_ml_fp <- here("plots/rreas_yoy_ts/cmisst_index/ML")

# needed to remove error message
input.loocv <- F

# create plots
for(i in 1:length(yoy_species)){
  # stock name
  png_name <- yoy_species[i]
  
  # SSH plots
  df_ssh <- yoy_cmisst_ssh[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  
  png(paste0(index_ssh_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ssh)
  dev.off()
  
  png( paste0(index_ssh_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ssh)
  dev.off()
  
  
  
  # SST plots
  df_sst <- yoy_cmisst_sst[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  png(paste0(index_sst_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_sst)
  dev.off()
  
  png( paste0(index_sst_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_sst)
  dev.off()
  
  
  # Mixed layer plots
  # SST plots
  df_ml <- yoy_cmisst_ml[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  png(paste0(index_ml_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ml)
  dev.off()
  
  png( paste0(index_ml_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ml)
  dev.off()
}


# RREAS indices -----------------------------------------------------------
# instead of just CPUE from Mary's update, going to use indices Tanya created for some of the species from the rreas survey

## data --------------------------------------------------------------------
# load data
canary <- read_csv(here("data/rreas/canary_indices.csv"))
chilipepper <- read_csv(here("data/rreas/chilipepper_indices.csv"))
widow <- read_csv(here("data/rreas/widow_indices.csv"))
ytail_n <- read_csv(here("data/rreas/ytail_north_indices.csv"))
ytail_coast <- read_csv(here("data/rreas/ytail_coastwide_indices.csv"))

# combine to data frame
rreas_indices <- list(canary = canary, chilipepper = chilipepper, widow = widow, 
                      ytail_coast = ytail_coast, ytail_n = ytail_n)
yoy_species <- c("canary", "chilipepper", "widow", "ytail_coast", "ytail_n")
## SSH ---------------------------------------------------------------------
yoy_index_cmisst_ssh <- list()
for(i in 1:length(yoy_species)){
  
  df <- rreas_indices[[i]] %>% 
    filter(YEAR <= 2023) %>% 
    select(YEAR, est)
  
  cmisst <- get_CMISST_index(response = df, oceanData = oceanSSHData,
                             years = df$YEAR, years.fit = df$YEAR)
  
  yoy_index_cmisst_ssh[[i]] <- cmisst
}
names(yoy_index_cmisst_ssh) <- yoy_species

saveRDS(yoy_index_cmisst_ssh, file = here("data/processed/yoy_index_cmisst_ssh.rds"))

## SST ---------------------------------------------------------------------
yoy_index_cmisst_sst <- list()
oceanData_ERSST <- oceanSSTData
for(i in 1:length(yoy_species)){
  
  df <- rreas_indices[[i]]%>% 
    filter(YEAR <= 2023) %>% 
    select(YEAR, est)
  
  cmisst <- get_CMISST_index(response = df, oceanData = oceanData_ERSST,
                             years = df$YEAR, years.fit = df$YEAR)
  
  yoy_index_cmisst_sst[[i]] <- cmisst
}
names(yoy_index_cmisst_sst) <- yoy_species

saveRDS(yoy_index_cmisst_sst, file = here("data/processed/yoy_index_cmisst_sst.rds"))



## Mixed layer -------------------------------------------------------------
yoy_index_cmisst_ml <- list()
for(i in 1:length(yoy_species)){
  
  df <- rreas_indices[[i]] %>% 
    filter(YEAR <= 2019) %>% 
    select(YEAR, est)
  
  cmisst <- get_CMISST_index(response = df, oceanData = oceanMLData,
                             years = df$YEAR, years.fit = df$YEAR)
  
  yoy_index_cmisst_ml[[i]] <- cmisst
}

names(yoy_index_cmisst_ml) <- yoy_species

saveRDS(yoy_index_cmisst_ml, file = here("data/processed/yoy_index_cmisst_ml.rds"))



## Plots -------------------------------------------------------------------
# create file path for plots
cov_map_ssh_fp <- here("plots/rreas_yoy_indices/covariance_maps/SSH")
cov_map_sst_fp <- here("plots/rreas_yoy_indices/covariance_maps/SST")
cov_map_ml_fp <- here("plots/rreas_yoy_indices/covariance_maps/ML")
index_ssh_fp <- here("plots/rreas_yoy_indices/cmisst_index/SSH")
index_sst_fp <- here("plots/rreas_yoy_indices/cmisst_index/SST")
index_ml_fp <- here("plots/rreas_yoy_indices/cmisst_index/ML")

# create plots
for(i in 1:length(yoy_species)){
 
  # stock name
  png_name <- yoy_species[i]
  
  # SSH plots
  df_ssh <- yoy_index_cmisst_ssh[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  
  png(paste0(index_ssh_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ssh)
  dev.off()
  
  png( paste0(index_ssh_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ssh)
  dev.off()
  
  
  
  # SST plots
  df_sst <- yoy_index_cmisst_sst[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  png(paste0(index_sst_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_sst)
  dev.off()
  
  png(paste0(index_sst_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_sst)
  dev.off()
  
  
  # Mixed Layer depth plots
  df_ml <- yoy_index_cmisst_ml[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  
  png(paste0(index_ml_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ml)
  dev.off()
  
  png( paste0(index_ml_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ml)
  dev.off()
}





# Recruitment Deviations -------------------------------------
# (pulled from Ward et al. 2024)
species_df <- readRDS(here("data/rec_devs/species_assessment_sb_rec.rds"))
black_CA <- read_csv(here("data/rec_devs/black_rockfish_CA_2023.csv"))
black_CA <- black_CA %>% 
  filter(Yr >= 1980)

black_WA <- read_csv(here("data/rec_devs/black_rockfish_WA_2023.csv"))
black_WA <- black_WA %>% 
  filter(Yr >= 1980)


black_OR <- read_csv(here("data/rec_devs/black_rockfish_OR_2023.csv"))
black_OR <- black_OR %>% 
  filter(Yr >= 1980)


shortspine <- read_csv(here("data/rec_devs/shortspine_thornyhead_2023.csv"))
shortspine <- shortspine %>% 
  filter(Yr >= 1980)

# update rec devs to most recent assessment 
species_df$Black_rockfish_WA <- black_WA
species_df$Black_rockfish_CA <- black_CA
species_df$Black_rockfish_OR <- black_OR
species_df$Shortspine_thornyhead <- shortspine

species <- names(species_df)

## SSH ------------------------------------------------------
recdev_cmisst_ssh <- list()

for(i in 1:length(species)){
  df <- species_df[[i]]
  
  cmisst <- get_CMISST_index(response = df[,c("Yr", "dev")], oceanData = oceanSSHData,
                             years = df$Yr, years.fit = df$Yr)
  
  recdev_cmisst_ssh[[i]] <- cmisst
  
}

names(recdev_cmisst_ssh) <- species

# save
saveRDS(recdev_cmisst_ssh, file = here("data/processed/recdev_cmisst_ssh.rds"))

## SST ---------------------------------------------------------------------
# get cmisst sst
recdev_cmisst_sst <- list()
for(i in 1:length(species)){
  df <- species_df[[i]]
  
  cmisst <- get_CMISST_index(response = df[,c("Yr", "dev")], oceanData = oceanSSTData,
                             years = df$Yr, years.fit = df$Yr)
  
  recdev_cmisst_sst[[i]] <- cmisst
  
}

names(recdev_cmisst_sst) <- species

# save 
saveRDS(recdev_cmisst_sst, file = here("data/processed/recdev_cmisst_sst.rds"))

## Mixed Layer depth -----------------------------------------------------------
recdev_cmisst_ml <- list()
for(i in 1:length(species)){
  df <- species_df[[i]]
  df <- df %>% 
    filter(Yr >= 1991) %>% 
    filter(Yr <= 2020)
  
  cmisst <- get_CMISST_index(response = df[,c("Yr", "dev")], oceanData = oceanMLData,
                             years = df$Yr, years.fit = df$Yr)
  
  recdev_cmisst_ml[[i]] <- cmisst
  
}

names(recdev_cmisst_ml) <- species

# save 
saveRDS(recdev_cmisst_ml, file = here("data/processed/recdev_cmisst_ml.rds"))


## Plots -----------------------------------------------------------------
# create file path for plots
cov_map_ssh_fp <- here("plots/rec_dev_ts/covariance_maps/SSH")
cov_map_sst_fp <- here("plots/rec_dev_ts/covariance_maps/SST")
cov_map_ml_fp <- here("plots/rec_dev_ts/covariance_maps/ML")
index_ssh_fp <- here("plots/rec_dev_ts/cmisst_index/SSH")
index_sst_fp <- here("plots/rec_dev_ts/cmisst_index/SST")
index_ml_fp <- here("plots/rec_dev_ts/cmisst_index/ML")

# species image names (species list has back slashes, want to remove)
species_png_names <- c("black_rockfish_WA", "black_rockfish_CA", "blue_deacon_rockfish_CA", "bocaccio", "cabezon_OR",
                       "cabezon_NCA", "cabezon_SCA", "california_scorpionfish", "chilipepper", "gopher_black-and-yellow_rockfish",
                       "kelp_greenling_OR", "longspine_thornyhead", "pacific_ocean_perch", "rougheye_blackspotted_rockfish",
                       "shortspine_thornyhead", "widow_rockfish", "yelloweye_rockfish", 
                       "yellowtail_rockfish_N", "yellowtail_rockfish_S", "black_rockfish_OR")

# remove error message
input.loocv <- F

# create plots
for(i in 1:length(species)){
  # stock name
  png_name <- species_png_names[i]
  
  # SSH plots
  df_ssh <- recdev_cmisst_ssh[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ssh)
  ggsave(paste0(cov_map_ssh_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  
  png(paste0(index_ssh_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ssh)
  dev.off()
  
  png( paste0(index_ssh_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ssh)
  dev.off()
  
  
  
  # SST plots
  df_sst <- recdev_cmisst_sst[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_sst)
  ggsave(paste0(cov_map_sst_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  png(paste0(index_sst_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_sst)
  dev.off()
  
  png(paste0(index_sst_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_sst)
  dev.off()
  
  
  # Mixed Layer depth plots
  df_ml <- recdev_cmisst_ml[[i]]
  
  # save as png
  makeCovarianceMap(input.season = 1, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
  
  
  makeCovarianceMap(input.season = 2, cmisst = df_ml)
  ggsave(paste0(cov_map_ml_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
  
  
  png(paste0(index_ml_fp, "/winter/", png_name, ".png"))
  makeBiplot(input.season = 1, cmisst = df_ml)
  dev.off()
  
  png( paste0(index_ml_fp, "/spring/", png_name, ".png"))
  makeBiplot(input.season = 2, cmisst = df_ml)
  dev.off()
}



# look at chilipepper before/after
chili_before <- yoy_rockfish_spp %>% 
  select(year, rreas.yoy.chili) %>% 
  filter(year < 2014)
chili_after <- yoy_rockfish_spp %>% 
  select(year, rreas.yoy.chili) %>% 
  filter(year >= 2014)

# mix layer
chili_ml_before <- get_CMISST_index(response = chili_before[9:31,], oceanData = oceanMLData,
                           years = chili_before$year[9:31], years.fit = chili_before$year[9:31])


chili_ml_after <- get_CMISST_index(response = chili_after[1:7,], oceanData = oceanMLData,
                                    years = chili_after$year[1:7], years.fit = chili_after$year[1:7])
# plots
makeCovarianceMap(input.season = 1, cmisst = chili_ml_before)
makeCovarianceMap(input.season = 2, cmisst = chili_ml_before)

makeCovarianceMap(input.season = 1, cmisst = chili_ml_after)
makeCovarianceMap(input.season = 2, cmisst = chili_ml_after)


# ssh
chili_ssh_before <- get_CMISST_index(response = chili_before[9:31,], oceanData = oceanSSHData,
                                    years = chili_before$year[9:31], years.fit = chili_before$year[9:31])


chili_ssh_after <- get_CMISST_index(response = chili_after[1:7,], oceanData = oceanSSHData,
                                   years = chili_after$year[1:7], years.fit = chili_after$year[1:7])
# plots
# save as png
makeCovarianceMap(input.season = 1, cmisst = chili_ssh_before)
makeCovarianceMap(input.season = 2, cmisst = chili_ssh_before)

makeCovarianceMap(input.season = 1, cmisst = chili_ssh_after)
makeCovarianceMap(input.season = 2, cmisst = chili_ssh_after)


# Before/After Marine Heatwave --------------------------------------------
# going to look at covariance maps for before
for(i in 1:length(species)){
  df <- species_df[[i]]
  sp <- species[i]
  
  min_yr <- min(df$Yr)
  max_yr <- max(df$Yr)
  
  print(paste0(sp, ": ", min_yr, " - ", max_yr))
  }
# going to limit rec dev to stock assessments that have data through 2020
recdev_df <- list("black_CA_recdev" = species_df$Black_rockfish_CA, "blue_deacon_CA_recdev" = species_df$`Blue/Deacon_rockfish_CA`,
                    "cabezon_NCA_recdev" = species_df$Cabezon_NCA, "shortspine_thornyhead_recdev" = species_df$Shortspine_thornyhead, 
                    "widow_recdev" = species_df$Widow_rockfish)
# going to standardize columns between different lists
for(i in 1:length(names(recdev_df))){
  df <- recdev_df[[i]]
  
  df <- df %>% 
    ungroup %>% 
    mutate(year = Yr,
           response = dev) %>% 
    select(year, response)
  
  recdev_df[[i]] <- df
}


rreas_ts <- list()
for(i in 1:6){
  df <- yoy_rockfish_spp[,c(1, i+2)]
  
  colnames(df)[2] <- "response"
  
  rreas_ts[[i]] <- df
}
names(rreas_ts) <- c( "bocaccio_rreas", "chilipepper_rreas", "halfbanded_rreas", "shortbelly_rreas", "widow_rreas", "yellowtail_rreas")


names(rreas_indices) <- c("canary_index", "chilipepper_index", "widow_index", "ytail_coast_index", "ytail_n_index")
for(i in 1:length(names(rreas_indices))){
  df <- rreas_indices[[i]]
  df <- df %>%
    mutate(response = est,
           year = YEAR) %>% 
    select(year, response)
  
  rreas_indices[[i]] <- df
}



# combine into one list
heatwave_species_df <- c(rreas_ts, rreas_indices, recdev_df)


## ML ----------------------------------------------------------------------
heatwave_ml_before <- list()
heatwave_ml_after <- list()
for(i in 1:length(names(heatwave_species_df))){
  df <- heatwave_species_df[[i]]
  
  df_before <- df %>% 
    filter(year >= 1991) %>% 
    filter(year < 2014)
  
  df_after <- df %>% 
    filter(year >= 2014) %>% 
    filter(year <= 2020)
  
  heatwave_ml_before[[i]] <- get_CMISST_index(df_before, oceanData = oceanMLData,
                                              years = df_before$year, years.fit = df_before$year)
  
  heatwave_ml_after[[i]] <- get_CMISST_index(df_after, oceanData = oceanMLData,
                                              years = df_after$year, years.fit = df_after$year)
}

saveRDS(heatwave_ml_before, file = here("data/processed/heatwave_before_cmisst_ml.rds"))
saveRDS(heatwave_ml_after, file = here("data/processed/heatwave_after_cmisst_ml.rds"))

# # # Species spawning biomass ------------------------------------------------
# # ## Create data frame -------------------------------------------------------------
# clean_rec_devs <- readRDS(here("data/clean_rec_devs.rds"))
# 
# species <- unique(clean_rec_devs$short_name)
# 
# species_df <- list()
# for(i in 1:length(species)){
#   sp <- species[i]
# 
#   df <- clean_rec_devs %>%
#     filter(short_name == sp) %>%
#     mutate(SpawnBio.scaled = as.vector(scale(SpawnBio))) %>%
#     filter(Yr >= 1980)
# 
#   species_df[[i]] <- df
# 
# }
# names(species_df) <- species
# 
# # typo for cabezon NCA, end year 2030 instead of 2020
# species_df$Cabezon_NCA$Yr[40] <- 2020
# # going to trim some species time series because they jump around and that won't work with cmisst function
# 
# # cabezon OR to just 2015 since it skips from 2015 to 2020
# species_df$Cabezon_OR <- species_df$Cabezon_OR[-37,]
# species_df$`Gopher/black-and-yellow_rockfish` <- species_df$`Gopher/black-and-yellow_rockfish`[-40,]
# 
# 
# # save as separate rds file for later
# write_rds(species_df, file = here("data/species_assessment_sb_rec.rds"))
# # 
# # 
# ## SSH ---------------------------------------------------------------------
# # get cmisst ssh for each species
# species_cmisst_ssh <- list()
# 
# for(i in 1:length(species)){
#   df <- species_df[[i]]
#   
#   cmisst <- get_CMISST_index(response = df[,c(1,6)], oceanData = oceanData_SSH,
#                              years = df$Yr, years.fit = df$Yr)
#   
#   species_cmisst_ssh[[i]] <- cmisst
#   
# }
# 
# names(species_cmisst_ssh) <- species
# 
# # save
# save(species_cmisst_ssh, file = here("data/species_cmisst_ssh.rds"))
# 
# ## SST ---------------------------------------------------------------------
# # get cmisst sst
# species_cmisst_sst <- list()
# for(i in 1:length(species)){
#   df <- species_df[[i]]
#   
#   cmisst <- get_CMISST_index(response = df[,c(1,6)], oceanData = oceanData_ERSST,
#                              years = df$Yr, years.fit = df$Yr)
#   
#   species_cmisst_sst[[i]] <- cmisst
#   
# }
# 
# names(species_cmisst_sst) <- species
# 
# # save 
# save(species_cmisst_sst, file = here("data/species_cmisst_sst.rds"))
# 
# 
# ## Plots -------------------------------------------------------------------
# # create file path for plots
# cov_map_ssh_fp <- here("plots/spawning_biomass_ts/covariance_maps/SSH")
# cov_map_sst_fp <- here("plots/spawning_biomass_ts/covariance_maps/SST")
# index_ssh_fp <- here("plots/spawning_biomass_ts/cmisst_index/SSH")
# index_sst_fp <- here("plots/spawning_biomass_ts/cmisst_index/SST")
# 
# # species image names (species list has back slashes, want to remove)
# species_png_names <- c("black_rockfish_WA", "black_rockfish_CA", "blue_deacon_rockfish_CA", "bocaccio", "cabezon_OR",
#                       "cabezon_NCA", "cabezon_SCA", "california_scorpionfish", "chilipepper", "gopher_black-and-yellow_rockfish",
#                       "kelp_greenling_OR", "longspine_thornyhead", "pacific_ocean_perch", "rougheye_blackspotted_rockfish",
#                       "shortspine_thornyhead", "widow_rockfish", "yelloweye_rockfish", "yellowtail_rockfish_N", "yellowtail_rockfish_S")
# 
# # remove error message
# input.loocv <- F
# 
# # create plots
# for(i in 1:length(species)){
#   # stock name
#   png_name <- species_png_names[i]
#   
#   # SSH plots
#   df_ssh <- species_cmisst_ssh[[i]]
#   
#   # save as png
#   makeCovarianceMap(input.season = 1, cmisst = df_ssh)
#   ggsave(paste0(cov_map_ssh_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
#   
#   
#   makeCovarianceMap(input.season = 2, cmisst = df_ssh)
#   ggsave(paste0(cov_map_ssh_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
#   
#   
#   png(paste0(index_ssh_fp, "/winter/", png_name, ".png"))
#   makeBiplot(input.season = 1, cmisst = df_ssh)
#   dev.off()
#   
#   png( paste0(index_ssh_fp, "/spring/", png_name, ".png"))
#   makeBiplot(input.season = 2, cmisst = df_ssh)
#   dev.off()
#   
#   
#   
#   # SST plots
#   df_sst <- species_cmisst_sst[[i]]
#   
#   # save as png
#   makeCovarianceMap(input.season = 1, cmisst = df_sst)
#   ggsave(paste0(cov_map_sst_fp, "/winter/", png_name, ".png"), width = 4, height = 4)
#   
#   
#   makeCovarianceMap(input.season = 2, cmisst = df_sst)
#   ggsave(paste0(cov_map_sst_fp, "/spring/", png_name, ".png"), width = 4, height = 4)
#   
#   png(paste0(index_sst_fp, "/winter/", png_name, ".png"))
#   makeBiplot(input.season = 1, cmisst = df_sst)
#   dev.off()
#   
#   png( paste0(index_sst_fp, "/spring/", png_name, ".png"))
#   makeBiplot(input.season = 2, cmisst = df_sst)
#   dev.off()
# }
