library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(ncdf4)
library(RColorBrewer)
library(reshape2)

# Index function ----------------------------------------------------------
get_CMISST_index <- function(response, oceanData=oceanData_ERSST,
                             years=NA, years.fit=years.fit,
                             months=1:12,
                             min.lon=232, max.lon=244,
                             min.lat=30, max.lat=50) {
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
  covs1<-apply(oceanData.s1.scl[,,as.character(years.fit)], 1:2, function(x) cov(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs2<-apply(oceanData.s2.scl[,,as.character(years.fit)], 1:2, function(x) cov(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs3<-apply(oceanData.s3.scl[,,as.character(years.fit)], 1:2, function(x) cov(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  covs4<-apply(oceanData.s4.scl[,,as.character(years.fit)], 1:2, function(x) cov(x, response$val[response$year %in% years.fit], use="pairwise.complete.obs"))
  
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
  limits<-c(-lmt, lmt)
  extent <- cmisst[[6]] # min, max of lat, long
  
  gg <- ggplot() + ggtitle(myTitle) +
    geom_raster(data = melt(covMap), aes(x = Var1, y = Var2, fill=value)) +
    geom_sf(data=land, color="black", fill="grey", linewidth=0.25) +
    xlim(extent[3], extent[4]) + ylim(extent[1], extent[2]) +
    scale_fill_gradientn(colours = myPalette(100),limits=limits,name="Covariance", na.value = "white") +
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

# RREAS data --------------------------------------------------------------
yoy_rockfish_all <- read_csv(here("data/rreas_juv_rockfish_all.csv"))

response <- yoy_rockfish_all
response$response_scaled <- scale(response$response)
## SSH ---------------------------------------------------------------------
cmisst <- get_CMISST_index(response = response[,c(1,3)], oceanData = oceanData_SSH, years = response$year,
                           years.fit = response$year)
 
makeCovarianceMap(input.season = 1, cmisst = cmisst)
makeBiplot(input.season = 1, cmisst = cmisst)

makeCovarianceMap(input.season = 2, cmisst = cmisst)
makeBiplot(input.season = 2, cmisst = cmisst)

# try lagged responses
response$year_lagged <- response$year - 1

cmisst2 <- get_CMISST_index(response = response[,c(4,3)], oceanData = oceanData_SSH, years = response$year_lagged,
                           years.fit = response$year_lagged)

makeCovarianceMap(input.season = 1, cmisst = cmisst2)
makeBiplot(input.season = 1, cmisst = cmisst2)

makeCovarianceMap(input.season = 2, cmisst = cmisst2)
makeBiplot(input.season = 2, cmisst = cmisst2)
