#Loading libraries
library(terra)
library(spatialEco)
library(raster)

#Calculating bandwith for KDE - From Valavi et al 2021 - https://rvalavi.github.io/SDMwithRFs/
bandwith_calculation <- function(da){
  #Calculating latitudinal extent of data
  yext <- (max(da$latitude) - min(da$latitude))/5
  xext <- (max(da$longitude) - min(da$longitude))/5
  max_ext <- max(c(yext, xext))
  return(max_ext)
}


#This script alters the randomPoints function from the dismo package
#We are adding a "replace" argument to allow for sampling with replacement
#The function has also been adapted to work with terra.
randomPoints2 <- function (mask, n, p, ext = NULL, extf = 1.1, excludep = TRUE, replace = FALSE, 
                           prob = FALSE, cellnumbers = FALSE, tryf = 3, warn = 2, lonlatCorrection = TRUE) 
{
  if (nlayers(mask) > 1) {
    mask <- raster(mask, 1)
  }
  tryf <- max(round(tryf[1]), 1)
  if (missing(p)) {
    excludep <- FALSE
  }
  else {
    if (inherits(p, "SpatialPoints")) {
      p <- coordinates(p)
    }
  }
  if (inherits(ext, "character")) {
    if (!ext %in% c("points")) {
      stop("if ext is a character variable it should be 'points'")
    }
    else if (missing(p)) {
      warning("if p is missing, 'ext=points' is meaningless")
      ext <- extent(mask)
    }
    else {
      ext <- extent(min(p[, 1]), max(p[, 1]), min(p[, 
                                                    2]), max(p[, 2]))
    }
  }
  if (!is.null(ext)) {
    ext <- extent(ext)
    ext <- ext * extf
    ext <- intersect(ext, extent(mask))
    mask2 <- crop(raster(mask), ext)
  }
  else {
    mask2 <- raster(mask)
  }
  if (n > ncell(mask2)) {
    n <- ncell(mask2)
    if (warn > 0) {
      warning("changed n to ncell of the mask (extent)")
    }
  }
  nn <- n * tryf
  nn <- max(nn, 10)
  if (prob) {
    stopifnot(hasValues(mask))
    cells <- crop(mask, mask2)
    cells <- try(stats::na.omit(cbind(1:ncell(cells), getValues(cells))))
    if (inherits(cells, "try-error")) {
      stop("the raster is too large to be used with 'prob=TRUE'")
    }
    prob <- cells[, 2]
    cells <- cells[, 1]
    if (couldBeLonLat(mask) && lonlatCorrection) {
      rows <- rowFromCell(mask2, cells)
      y <- yFromRow(mask2, rows)
      dx <- pointDistance(cbind(0, y), cbind(xres(mask2), 
                                             y), longlat = TRUE)
      dx <- dx/max(dx)
      prob <- prob * dx
    }
    cells <- sample(cells, nn, prob = prob, replace = replace)
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
  }
  else if (canProcessInMemory(mask2)) {
    cells <- crop(mask, mask2)
    if (hasValues(cells)) {
      cells <- which(!is.na(getValues(cells)))
    }
    else {
      cells <- 1:ncell(cells)
    }
    nn <- min(length(cells), nn)
    if (lonlatCorrection & couldBeLonLat(mask)) {
      rows <- rowFromCell(mask2, cells)
      y <- yFromRow(mask2, rows)
      dx <- pointDistance(cbind(0, y), cbind(xres(mask2), 
                                             y), longlat = TRUE)
      cells <- sample(cells, nn, prob = dx)
    }
    else {
      cells <- sample(cells, nn)
    }
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
  }
  else {
    nn <- min(ncell(mask2), nn)
    if (couldBeLonLat(mask2)) {
      cells <- .randomCellsLonLat(mask2, nn)
    }
    else {
      if (nn >= ncell(mask2)) {
        cells <- 1:ncell(mask2)
      }
      else {
        cells <- sampleInt(ncell(mask2), nn)
      }
    }
    xy <- xyFromCell(mask2, cells)
    cells <- cellFromXY(mask, xy)
    if (hasValues(mask)) {
      vals <- cbind(cells, extract(mask, cells))
      cells <- stats::na.omit(vals)[, 1]
    }
  }
  if (excludep) {
    pcells <- cellFromXY(mask, p)
    cells <- cells[!(cells %in% pcells)]
  }
  if (length(cells) > n) {
    cells <- sample(cells, n)
  }
  else if (length(cells) < n) {
    frac <- length(cells)/n
    if (frac < 0.1) {
      stop("generated random points = ", frac, " times requested number; Use a higher value for tryf")
    }
    if (frac < 0.5 & warn == 1) {
      warning("generated random points = ", frac, " times requested number; Use a higher value for tryf")
    }
    else if (warn > 1) {
      warning("generated random points = ", frac, " times requested number")
    }
  }
  if (cellnumbers) {
    return(cells)
  }
  else {
    return(xyFromCell(mask, cells))
  }
}


#Creating background points
bg_pts <- function(da, ras, n, prob, try, replace){
  #Generating random samples from the KDE raster file
  set.seed(42)
  #Calculating bandwidth
  bw <- bandwith_calculation(da)
  #Calculating KDE
  tgb_kde <- sp.kde(x = da, bw = round(bw, 2), ref = ras, 
                    standardize = TRUE, scale.factor = 10000)
  #Getting 20 times the number of unique observations
  kde_samples <- randomPoints2(mask = raster(tgb_kde), n = nrow(da)*20, prob = T, 
                               tryf = 2, replace = T)
  #Getting background data points as data frame
  bg_samples <- kde_samples %>% 
    as.data.frame() %>% 
    rename("longitude" = "x", "latitude" = "y") 
  
  #Matching background points to observation dates  
  #Creating a vector with observation indices selected at random 20 times to match background
  ind <- sample(1:nrow(da), size = nrow(da)*20, replace = T)
  
  #Applying indices to data frame with unique obs
  bg_samples <- da[ind,] %>% 
    st_drop_geometry() %>% 
    dplyr::select(date, year:month, season_year:presence) %>% 
    mutate(presence = 0) %>% 
    cbind(bg_samples)
  
  return(bg_samples)
}