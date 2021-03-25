#*****************************************  Script description  ****************************************************************************
#
# This script is being modified to minimize the manual work required for the massive subsetting (for the second round of calibration)
# The only input required aside from the model domain files and properites is either a list of comIDs, or a list of gagesID, 
# If the gage does not exists in the Routlink file, it just throws a warning and continue with the subsetting the gauges which exists in the Routlink file
#
#*********************************************************************************************************************************************
library(data.table)
library(rwrfhydro)
library(ncdf4)
source("Utils_ReachFiles.R")
library(raster)
library(foreach)
#------------------------------------------ INPUTS : Please change the following parts as instructed ------------------------------------------

# Specify the path to your new subset domain files
myPath <- "/PATH/TO/DESTINATION/DIRECTORY"

# A charcter vector : list of the comIDs to be used for subsetting, 
# if gages is defined comIds should be NULL
# They should exists in the Routlink otherwise gives you a warning and proceed to the next one
comIds <- NULL

# A character vector : list of the gauges to be used for subsetting, 
# if comIds is specified, gages should be NULL
# if specified then should exists on the Routlink, otherwise only gives you a warning and proceed to the next one
gageIds <- c("07106500")

#******  Specify the path to the ORIGINAL (full extent) domain files: *****************
# If you do not want to subset any of the following files, assign it as NULL.
# The fullGeoFile (Geogrid domain) is required at all times.

# Path to the Routing domain file
fullHydFile <- "/PATH/TO/Fulldom_hires.nc or NULL"

# Path to the Geogrid file
fullGeoFile <- "/PATH/TO/geo_em.d0x.nc"

# Path to the Wrfinput file
fullWrfFile <- "/PATH/TO/wrfinput.nc or wrfinput_d0x or NULL"

# Path to the Routelink file
fullRtlinkFile <- "/PATH/TO/Route_Link.nc or NULL"

# Path to the Spatial Weights file
fullSpwtFile <- "/PATH/TO/spatialweights.nc or NULL"

# Path to the Ground Water Bucket parameter file
fullGwbuckFile <- "/PATH/TO/GWBUCKPARM.nc or NULL"

# Path to the Soil Properties parameter file
fullSoilparmFile <- "/PATH/TO/soil_properties.nc or NULL"

# Path to the Lake parameter file
fullLakeparmFile <- "/PATH/TO/LAKEPARM.nc or NULL"

# Path to the Lake Tyep file
fullLakeTypesFile <- "/PATH/TO/LAKETYPE.nc"

# Path to the hydro2D file , set this to NULL if you do not have a hydro 2D file
fullHydro2dFile <- "/PATH/TO/hydro2dtbl.nc or NULL"

# Path to the geo spatial file required for the new outputting option
geoSpatialFile <- "/PATH/TO/GEOGRID_LDASOUT_Spatial_Metadata.nc or NULL"

# Path to the nudging parameter file
fullNudgeParamFile <- "/PATH/TO/nudgingParams.nc or NULL"

# Path to the re expression files
downstreamReExpFile <- "/PATH/TO/RouteLink.reExpTo.Rdb"
upstreamReExpFile   <- "/PATH/TO/RouteLink.reExpFrom.Rdb"
reIndFile           <- "/PATH/TO/RouteLink.reInd.Rdb"

# Projection for boundaing coordinates. This needs to be a PROJ4 string
# (e.g., "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs").
coordProj <- "+proj=lcc +lat_1=30 +lat_2=60 +lat_0=40.0000076293945 +lon_0=-97 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs"

# Multiplier between routing grid and LSM grid
# (e.g., 1-km LSM and 250-m routing means a value of 4)
dxy <- 4

# **** perform some checks on the matching length of gageMeta and comIds/gages in case gageMeta is not NULL

if (!is.null(gageIds) & !is.null(comIds)) print("Only one of the two variables (gauges and comIds) is required. The gages would be ignored")

#------------------------------ Extract those ComIds and relvalnt information form the Routlink file to perform the subsetting for them ------------------------------

# load the reexpression files
print(load(downstreamReExpFile))
print(load(upstreamReExpFile))
print(load(reIndFile))

routlinkInfo <- as.data.table(GetNcdfFile(fullRtlinkFile,c('gages','link','order'), q=TRUE))
routlinkInfo$ind <- 1:nrow(routlinkInfo)
#setnames(routlinkInfo, "index", "ind") # If index was one of the variables in the Routlink, we do not need the line above and this line should be commented out ...

if (!is.null(comIds)) {
  routlinkInfoSub <- routlinkInfo[link %in% comIds]
  if (length(setdiff(comIds, routlinkInfoSub$link)) > 0 ) {
    print("Warning : The following ComIds are not in the Routlink, therefore no subsetting will be performed for them : ")
    print(setdiff(comIds, routlinkInfo$link))
  }
} else if (!is.null(gageIds)) {
  routlinkInfoSub <- routlinkInfo[trimws(gages) %in% gageIds]
  if (length(setdiff(gageIds, trimws(routlinkInfoSub$gages))) > 0 ) {
    print("Warning : The following gages are not in the Routlink, therefore no subsetting will be performed for them : ")
    print(setdiff(gageIds, trimws(routlinkInfoSub$gages)))
  }
}

#------------------------------ Find all the COMIDs above the specified comIds ---------------------------------------

gageInds <- routlinkInfoSub$ind
names(gageInds) <- routlinkInfoSub$link

## make this a list of pairs: gageId, gageInd
dumList <- 1:length(gageInds)
names(dumList) <- trimws(names(gageInds))

gageIndsList <- plyr::llply(dumList, function(ii) list(gageId=names(gageInds)[ii],
                                                       gageInd=gageInds[[ii]]))
rm(gageInds, dumList)

gather <- function(gageIdInd) {
  upBranches <- rwrfhydro:::GatherStreamInds(from, start = gageIdInd$gageInd, linkLength=reInd$length)
  indsAll <- c(upBranches$ind, upBranches$startInd)
  upComIds <- routlinkInfo[ind %in% indsAll, link]
}

upComIdsAll <- plyr::llply(gageIndsList, gather, .parallel = FALSE)
rm(gageIndsList, gather, routlinkInfo)

#------------------  To find the min, max i and j values from the spatial weight files -------------------------------------------------------------- 

# build a datatable using the spatial weight file created by Kevin for NHDPlus
IDmask <- ncdump(fullSpwtFile, "IDmask", quiet = TRUE) # ID mask is the name of the polygon which is the same as the link in Routlink file
i_index  <- ncdump(fullSpwtFile, "i_index", quiet = TRUE) # index in the x dimension of the raster grid (starting from 1,1 in the LL corner)
j_index <- ncdump(fullSpwtFile, "j_index", quiet = TRUE) # index in the y dimension of the raster grid (starting from 1,1 in the LL corner)

swDT <- data.table(i_index = i_index, j_index = j_index, IDmask = IDmask)

infoDT <- as.data.table(plyr::ldply(upComIdsAll, function(x) {
  # print(x)
  sub <- swDT[IDmask %in% x]
  data.frame(hyd_w = min(sub$i_index), hyd_e = max(sub$i_index), hyd_s = min(sub$j_index), hyd_n = max(sub$j_index))
}))
rm(IDmask, swDT)

infoDT$hyd_w <- floor((infoDT$hyd_w - 0.00001*dxy)/dxy)*dxy + 1
infoDT$hyd_s <- floor((infoDT$hyd_s - 0.00001*dxy)/dxy)*dxy + 1

infoDT$hyd_e <- ceiling((infoDT$hyd_e)/dxy)*dxy
infoDT$hyd_n <- ceiling((infoDT$hyd_n)/dxy)*dxy

#******** Then we need to add the geo_min, geo_max (the min and max row number when the 1,1 is upperRight corner) and the hydro information
# grab the geo_domain information from hydro boundaries
infoDT[, `:=` (geo_w = floor((hyd_w-0.00001*dxy)/dxy) + 1 ,
               geo_e = floor((hyd_e-0.00001*dxy)/dxy) + 1,
               geo_s = floor((hyd_s-0.00001*dxy)/dxy) + 1,
               geo_n = floor((hyd_n-0.00001*dxy)/dxy) + 1)]


# Create temp geogrid tif
tmpfile <- tempfile(fileext=".tif")
ExportGeogrid(fullGeoFile, "HGT_M", tmpfile)
geohgt <- raster::raster(tmpfile)
file.remove(tmpfile)
rm(tmpfile) 
infoDT[, `:=` (geo_min = dim(geohgt)[1] - geo_n + 1,
               geo_max = dim(geohgt)[1] - geo_s + 1 )]

infoDT[, `:=` (hyd_min = (geo_min-1)*dxy+1,
               hyd_max = geo_max*dxy)]

save(infoDT , file = paste0(myPath, "/infoDT.Rdata"))

#----------------------- Now is time to start the subsetting -------------------------------------------------------------------

# Here I (Arezoo) have not changes the infromation from the script Aubrey had, eventhough we do not need buffer anymore, 
# I kept it to be similar to what Aubrey had 

# wrapping all the cutout parts into a function so we can parallize it easily
CutOut <- function(gage) {
  
  source("Utils_ReachFiles.R")
  if (!is.null(gageIds)) {
    folderName <- trimws(routlinkInfoSub[link == gage]$gages)
  }else{
    folderName <- gage
  }
  
  system(paste0("mkdir -p ", myPath, "/", folderName))
  
  # =====================================================================
  # Specify the NEW (subset extent) domain files:
  
  # Routing domain file
  subHydFile <- paste0(myPath, "/", folderName, "/Fulldom_hires.nc")
  
  # Geogrid domain file
  subGeoFile <- paste0(myPath, "/", folderName, "/geo_em.d0x.nc")
  
  # Wrfinput domain file
  subWrfFile <- paste0(myPath, "/", folderName, "/wrfinput_d0x.nc")
  
  # Route link file
  subRtlinkFile <- paste0(myPath, "/", folderName, "/Route_Link.nc")
  
  # Spatial weights file
  subSpwtFile <- paste0(myPath, "/", folderName, "/spatialweights.nc")
  
  # GW bucket parameter file
  subGwbuckFile <- paste0(myPath, "/", folderName, "/GWBUCKPARM.nc")
  
  # Soil parameter file
  subSoilparmFile <- paste0(myPath, "/", folderName, "/soil_properties.nc")
  
  # Lake parameter file
  subLakeparmFile <- paste0(myPath, "/", folderName,  "/LAKEPARM.nc")
  
  # Hydro 2D file
  subHydro2dFile <- paste0(myPath, "/", folderName, "/hydro2dtbl.nc")
  
  # geo Spatial file 
  subGeoSpatialFile <- paste0(myPath, "/", folderName, "/GEOGRID_LDASOUT_Spatial_Metadata.nc")
  
  # Coordinate parameter text file
  subCoordParamFile <- paste0(myPath, "/", folderName, "/params.txt")
  
  # Nudging parameter 
  subNudgeParamFile <- paste0(myPath, "/", folderName, "/nudgingParams.nc")
  
  # Forcing clip script file
  subScriptFile <- paste0(myPath, "/", folderName, "/script_forcing_subset.txt")
  
  #=================================================================================
  
  # Get subsetting dimensions
  geo_w <- infoDT[.id == gage]$geo_w
  geo_e <- infoDT[.id == gage]$geo_e
  geo_s <- infoDT[.id == gage]$geo_s
  geo_n <- infoDT[.id == gage]$geo_n
  hyd_w <- infoDT[.id == gage]$hyd_w
  hyd_e <- infoDT[.id == gage]$hyd_e
  hyd_s <- infoDT[.id == gage]$hyd_s
  hyd_n <- infoDT[.id == gage]$hyd_n
  hyd_min <- infoDT[.id == gage]$hyd_min
  hyd_max <- infoDT[.id == gage]$hyd_max
  geo_min <- infoDT[.id == gage]$geo_min
  geo_max <- infoDT[.id == gage]$geo_max
  
  # Get relevant real coords for new bounds
  geo_min_col <- infoDT[.id == gage]$geo_w
  geo_max_col <- infoDT[.id == gage]$geo_e
  geo_min_row <- infoDT[.id == gage]$geo_min
  geo_max_row <- infoDT[.id == gage]$geo_max
  
  rowcol_new <- data.frame(id=c(1,2,3,4), row=c(geo_max_row, geo_min_row, geo_min_row, geo_max_row),
                           col=c(geo_min_col, geo_min_col, geo_max_col, geo_max_col))
  sp_new_proj <- xyFromCell(geohgt, raster::cellFromRowCol(geohgt, rowcol_new$row, rowcol_new$col), spatial=TRUE)
  sp_new_wrf <- sp::coordinates(sp::spTransform(sp_new_proj, "+proj=longlat +a=6370000 +b=6370000 +no_defs"))
  
  #---------------------------- SUBSET DOMAINS  
  
  # NCO starts with 0 index for dimensions, so we have to subtract 1
  
  # Geo Spatial File
  if (!is.null(geoSpatialFile)) {
    
    if (!file.exists(geoSpatialFile)) stop(paste0("The geoSpatialFile :", geoSpatialFile, " does not exits"))
    cmd <- paste0("ncks -O -d x,", geo_w-1, ",", geo_e-1, " -d y,", geo_min-1, ",", geo_max-1, " ", geoSpatialFile, " ", subGeoSpatialFile)
    
    print(cmd)
    system(cmd)
    
    # read in the x and y values of the cutout domain
    xy_data <- rwrfhydro::GetNcdfFile(subGeoSpatialFile, q = TRUE)
    
    # let s modify the GeoTransform attribute of variable crs in the geoSpatialFile
    cmd <- paste0("ncatted -O -a GeoTransform,crs,o,c,'", xy_data$x[1] - attributes(xy_data)$x$resolution/2, " ",
                  attributes(xy_data)$x$resolution, " 0 ", xy_data$y[1] + attributes(xy_data)$y$resolution/2, " 0 ",
                  attributes(xy_data)$y$resolution,"' ", subGeoSpatialFile)
    print(cmd)
    system(cmd)
  }
  
  # ROUTING GRID
  if (!is.null(fullHydFile)) {
    
    if  (!file.exists(fullHydFile)) stop(paste0("The fullHydFile : ", fullHydFile, " does not exits"))
    cmd <- paste0("ncks -O -d x,", hyd_w-1, ",", hyd_e-1, " -d y,", hyd_min-1, ",", hyd_max-1, " ", fullHydFile, " ", subHydFile)
    print(cmd)
    system(cmd)
    
    # read in the x and y values of the cutout domain
    xy_data2 <- rwrfhydro::GetNcdfFile(subHydFile, q = TRUE)
    
    # let s modify the GeoTransform attribute of variable crs in the geoSpatialFile
    cmd <- paste0("ncatted -O -a GeoTransform,crs,o,c,'", xy_data2$x[1] - attributes(xy_data)$x$resolution/2/dxy, " ", # Move from Center to the upper left corner
                  attributes(xy_data)$x$resolution/dxy, " 0 ", xy_data2$y[1] + attributes(xy_data)$y$resolution/2/dxy, " 0 ",
                  attributes(xy_data)$y$resolution/dxy, "' ", subHydFile)
    print(cmd)
    system(cmd)
  }
  
  # GEO GRID
  # Dimension subsetting
  cmd <- paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1,
                " -d west_east_stag,", geo_w-1, ",", geo_e, " -d south_north_stag,",geo_s-1, ",", geo_n, " ",
                fullGeoFile, " ", subGeoFile)
  print(cmd)
  system(cmd)
  
  # Read the 2D corner coordinates
  corner_lats <- c()
  for (ncVarName in c('XLAT_M', 'XLAT_U', 'XLAT_V', 'XLAT_C')) {
    if (ncVarName %in% names(rwrfhydro::ncdump(subGeoFile, quiet = TRUE)$var)) {
      a = rwrfhydro::ncdump(subGeoFile, variable = ncVarName, quiet = TRUE)
      corners = c(a[1,1], a[1, ncol(a)], a[nrow(a), ncol(a)], a[nrow(a), 1])
      rm(a)
    }else{
      corners = c(0,0,0,0)
    }
    corner_lats = c(corner_lats, corners)                           # Populate corner_lats lis
  }
  
  corner_lons <- c()
  for (ncVarName in c('XLONG_M', 'XLONG_U', 'XLONG_V', 'XLONG_C')) {
    if (ncVarName %in% names(rwrfhydro::ncdump(subGeoFile, quiet = TRUE)$var)) {
      a = rwrfhydro::ncdump(subGeoFile, variable = ncVarName, quiet = TRUE)
      corners = c(a[1,1], a[1, ncol(a)], a[nrow(a), ncol(a)], a[nrow(a), 1])
      rm(a)
    }else{
      corners = c(0,0,0,0)
    }
    corner_lons = c(corner_lons, corners)                           # Populate corner_lats lis
  }
  
  # Attribute updates
  cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -O -a corner_lons,global,o,f,", paste(corner_lons, collapse  = ","), " ", subGeoFile)
  system(cmd)
  cmd <- paste0("ncatted -O -a corner_lats,global,o,f,", paste(corner_lats, collapse  = ","), " ", subGeoFile)
  system(cmd)
  
  #HYDRO_TBL_2D GRID
  if (!is.null(fullHydro2dFile)) {
    
    if (!file.exists(fullHydro2dFile)) stop(paste0("The fullHydro2dFile : ", fullHydro2dFile, " does not exits"))
    
    # Dimension subsetting
    #DIMENSION IS CURRENTLY north_south.. may change this after talking with Wei 04/28/2017
    cmd <- paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullHydro2dFile, " ", subHydro2dFile)
    print(cmd)
    system(cmd)
    
    # Attribute updates
    cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lons,global,o,f,", paste(corner_lons, collapse  = ","), " ", subHydro2dFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lats,global,o,f,", paste(corner_lats, collapse  = ","), " ", subHydro2dFile)
    system(cmd)
  }
  
  # WRFINPUT GRID
  
  if (!is.null(fullWrfFile)) {
    
    if (!file.exists(fullWrfFile)) stop(paste0("The fullWrfFile : ", fullWrfFile, " does not exits"))
    
    cmd <- paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullWrfFile, " ", subWrfFile)
    print(cmd)
    system(cmd)
    
    # Attribute updates
    cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lons,global,o,f,", paste(corner_lons, collapse  = ","), " ", subWrfFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lats,global,o,f,", paste(corner_lats, collapse  = ","), " ", subWrfFile)
    system(cmd)
  }
  
  
  #---------------------------- SUBSET PARAMS -----------------------------------
  
  if (!is.null(fullSpwtFile) & !is.null(fullRtlinkFile)) {
    
    if (!file.exists(fullSpwtFile)) stop(paste0("The fullSpwtFile : ", fullSpwtFile, " does not exits"))
    if (!file.exists(fullRtlinkFile)) stop(paste0("The fullRtlinkFile : ", fullRtlinkFile, " does not exits"))
    
    # Identify catchments to keep
    fullWts <- ReadWtFile(fullSpwtFile)
    keepIdsPoly <- subset(fullWts[[1]], fullWts[[1]]$i_index >= hyd_w & fullWts[[1]]$i_index <= hyd_e &
                            fullWts[[1]]$j_index >= hyd_s & fullWts[[1]]$j_index <= hyd_n)
    
    # keep only those basins that are fully within the cutout doamin
    FullBasins <- as.data.table(keepIdsPoly)
    FullBasins <- FullBasins[, .(sumBas = sum(weight)), by = IDmask] # find How much of the basins are in the cutout domain
    keepIdsPoly <- FullBasins[sumBas > .999]$IDmask
    
    fullRtlink <- ReadLinkFile(fullRtlinkFile)
    keepIdsLink <- upComIdsAll[[as.character(gage)]]
    keepIds <- unique(c(keepIdsPoly, keepIdsLink))
    
    # SPATIAL WEIGHT
    subWts <- SubsetWts(fullWts, keepIdsPoly, hyd_w, hyd_e, hyd_s, hyd_n)
    file.copy(fullSpwtFile, subSpwtFile, overwrite = TRUE)
    UpdateWtFile(subSpwtFile, subWts[[1]], subWts[[2]], subDim=TRUE)
    
    # ROUTE LINK
    subRtlink <- subset(fullRtlink, fullRtlink$link %in% keepIds)
    subRtlink$to <- ifelse(subRtlink$to %in% unique(subRtlink$link), subRtlink$to, 0)
    
    # reorder the ascendingIndex if ascendingIndex exists in the variables
    if ("ascendingIndex" %in% names(subRtlink)) subRtlink$ascendingIndex <- (rank(subRtlink$ascendingIndex) - 1)
    file.copy(fullRtlinkFile, subRtlinkFile, overwrite = TRUE)
    UpdateLinkFile(subRtlinkFile, subRtlink, subDim=TRUE)
  }
  
  # GWBUCK PARAMETER
  
  if (!is.null(fullGwbuckFile)) {
    
    if (is.null(fullSpwtFile) | is.null(fullRtlinkFile)) {
      stop("To subset the fullSpwtFile, you need fullSpwtFile and fullRtlinkFile")
    }
    if (!file.exists(fullGwbuckFile)) stop(paste0("the fullGwbuckFile : ", fullGwbuckFile, " does not exits"))
    fullGwbuck <- GetNcdfFile(fullGwbuckFile, quiet=TRUE)
    subGwbuck <- subset(fullGwbuck, fullGwbuck$ComID %in% keepIdsPoly)
    subGwbuck$Basin <- seq(1, nrow(subGwbuck), 1)
    file.copy(fullGwbuckFile, subGwbuckFile, overwrite = TRUE)
    UpdateGwbuckFile(subGwbuckFile, subGwbuck, subDim=TRUE)
  }
  
  # SOIL PARAMETER
  
  if (!is.null(fullSoilparmFile)) {
    
    if (!file.exists(fullSoilparmFile)) stop(paste0("the fullSoilparmFile : ", fullSoilparmFile, " does not exits"))
    cmd <- paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ", fullSoilparmFile, " ", subSoilparmFile)
    system(cmd)
    
    # Attribute updates
    cmd <- paste0("ncatted -h -a WEST-EAST_GRID_DIMENSION,global,o,l,", geo_e-geo_w+2, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_GRID_DIMENSION,global,o,l,", geo_n-geo_s+2, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_UNSTAG,global,o,l,", geo_e-geo_w+1, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_UNSTAG,global,o,l,", geo_n-geo_s+1, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_START_STAG,global,d,,, ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_START_STAG,global,d,,, ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a WEST-EAST_PATCH_END_STAG,global,d,,, ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a SOUTH-NORTH_PATCH_END_STAG,global,d,,, ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a i_parent_end,global,o,l,", geo_e-geo_w+2, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -h -a j_parent_end,global,o,l,", geo_n-geo_s+2, " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lons,global,o,f,", paste(corner_lons, collapse  = ","), " ", subSoilparmFile)
    system(cmd)
    cmd <- paste0("ncatted -O -a corner_lats,global,o,f,", paste(corner_lats, collapse  = ","), " ", subSoilparmFile)
    system(cmd)
  }
  
  # LAKE PARAMETER
  
  # Make a copy of the LAKEPARM file, Read the subsetted Routlink file
  # and finds out which lakes fall into the domain and keep only those in the LAKEPARM.nc file
  
  if (!is.null(fullLakeparmFile)) {
    
    if (is.null(fullSpwtFile) | is.null(fullRtlinkFile)) {
      stop("To subset the fullLakeparmFile, you need fullSpwtFile and fullRtlinkFile")
    }
    if (!file.exists(fullLakeparmFile)) stop(paste0("the fullLakeparmFile : ", fullLakeparmFile, " does not exits"))
    rl <- ReadRouteLink(subRtlinkFile)
    rll <- subset(rl, rl$NHDWaterbodyComID>=0)
    lk <- GetNcdfFile(fullLakeparmFile)
    lkl <- subset(lk, lk$lake_id %in% rll$NHDWaterbodyComID)
    
    # Identify lake outlets for full domain
    rl.full <- ReadRouteLink(fullRtlinkFile)
    rl.full$ind <- 1:nrow(rl.full)
    lakeTypes <- GetNcdfFile(fullLakeTypesFile, quiet=TRUE)
    names(lakeTypes)[which(names(lakeTypes)=="LINKID")] <- "link"
    rl.full <- plyr::join(rl.full, lakeTypes, by="link")
    lakeOutlets <- subset(rl.full, rl.full$NHDWaterbodyComID > 0 & rl.full$TYPEL==1)
    rm(rl.full)
    
    routlinkInfo <- as.data.table(GetNcdfFile(fullRtlinkFile,c('gages','link','order'), q=TRUE))
    routlinkInfo$ind <- 1:nrow(routlinkInfo)
    gather <- function(gageIdInd) {
      upBranches <- rwrfhydro:::GatherStreamInds(from, start = gageIdInd$gageInd, linkLength=reInd$length)
      indsAll <- c(upBranches$ind, upBranches$startInd)
      upComIds <- routlinkInfo[ind %in% indsAll, link]
    }
    
    lakeIds <- unique(rll$NHDWaterbodyComID)
    lakesToRemove <- c()
    for (lake in lakeIds) {
      lakecom <- subset(lakeOutlets$link, lakeOutlets$NHDWaterbodyComID == lake)
      if (length(lakecom)>0) {
        message(paste0("Processing lake ", lake, " outlet ", lakecom))
        
        # Taken from Arezoo's tracing code. Overkill here, but it works!
        # Extract those ComIds and relvalnt information form the Routlink file to perform the subsetting for them
        routlinkInfoSub <- routlinkInfo[link == lakecom]
        # Find all the COMIDs above the specified comIds
        gageInds <- routlinkInfoSub$ind
        names(gageInds) <- routlinkInfoSub$link
        # make this a list of pairs: gageId, gageInd
        dumList <- 1:length(gageInds)
        names(dumList) <- trimws(names(gageInds))
        gageIndsList <- plyr::llply(dumList, function(ii) list(gageId=names(gageInds)[ii],
                                                               gageInd=gageInds[[ii]]))
        # Trace up
        upComIdsAll <- plyr::llply(gageIndsList, gather, .parallel = FALSE)
        rm(lakecom, routlinkInfoSub, gageInds, dumList, gageIndsList)
        
        # Check if all upstream reaches are in route link
        if (!all(unlist(upComIdsAll) %in% rl$link)) lakesToRemove <- c(lakesToRemove, lake)
      } else {
        message(paste0("No outlet found for ", lake, "so skipping"))
      }
    }
    
    # Adjust lakeparm and routelink
    lkl <- subset(lkl, !(lkl$lake_id %in% lakesToRemove))
    rl[rl$NHDWaterbodyComID %in% lakesToRemove, "NHDWaterbodyComID"] <- -9999
    UpdateLinkFile(subRtlinkFile, rl, subDim=FALSE)
    
    if (nrow(lkl) != 0) {
      if ("ascendingIndex" %in% names(lkl)) lkl$ascendingIndex <- (rank(lkl$ascendingIndex) - 1)
      file.copy(fullLakeparmFile, subLakeparmFile, overwrite = TRUE)
      UpdateLakeFile(subLakeparmFile, lkl, subDim=TRUE)
    }
  }
  
  # Nudging PARAMETER FILE, not tested yet
  
  if (!is.null(fullNudgeParamFile)) {
    
    if (is.null(fullSpwtFile) | is.null(fullRtlinkFile)) {
      stop("To subset the fullNudgeParamFile, you need fullSpwtFile and fullRtlinkFile")
    }
    if (!file.exists(fullNudgeParamFile)) stop(paste0("The fullNudgeParamFile : ", fullNudgeParamFile, " does not exist"))
    
    file.copy(fullNudgeParamFile, subNudgeParamFile, overwrite = TRUE)
    # read the list of gages from the subsetted Routlink
    RtGages <- unique(rwrfhydro::ncdump(subRtlinkFile, "gages", quiet = TRUE))
    # read the list of all gages from the full nudging param file
    nudgeGages <- rwrfhydro::ncdump(fullNudgeParamFile, "stationId", quiet = TRUE)
    # Find the indices for which will be used from subsetting
    subRtlinkInds <- which(nudgeGages %in% RtGages)
    # read the nudging param file in, and subset them
    nudgeInfo <- rwrfhydro::GetNcdfFile(fullNudgeParamFile, quiet=TRUE)
    for (varName in c("stationId", "R", "G", "tau")){ nudgeInfo[[varName]] <- nudgeInfo[[varName]][subRtlinkInds]}
    nudgeInfo$qThresh <- nudgeInfo$qThresh[, subRtlinkInds]
    nudgeInfo$expCoeff <- nudgeInfo$expCoeff[,, subRtlinkInds]
    
    # reduce the size of the file
    cmdtxt <- paste0("ncks -O -d stationIdInd,1,",  length(subRtlinkInds), " ", subNudgeParamFile, " ", subNudgeParamFile)
    system(cmdtxt)
    
    # read all the variables with link as their dimension (only dimension) and subset the link dimension
    ncin <- ncdf4::nc_open(subNudgeParamFile, write=TRUE)
    for (i in names(nudgeInfo)) {
      #      print(i)
      ncdf4::ncvar_put(ncin, i, nudgeInfo[[i]])
    }
    ncdf4::nc_close(ncin)
  }
  
  #--------------------------------- CREATE SCRIPT FILES 
  
  # Save the coordinate parameter file
  coordsExport <- data.frame(grid=c("hyd_sn", "hyd_ns", "geo_sn", "geo_ns"),
                             imin=c(hyd_w, hyd_w, geo_w, geo_w),
                             imax=c(hyd_e, hyd_e, geo_e, geo_e),
                             jmin=c(hyd_s, hyd_min, geo_s, geo_min),
                             jmax=c(hyd_n, hyd_max, geo_n, geo_max),
                             index_start=c(1,1,1,1))
  write.table(coordsExport, file=subCoordParamFile, row.names=FALSE, sep="\t")
  
  # Save the forcing subset script file
  #ncksCmd <- paste0("ncks -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ${OLDFORCPATH}/${i} ${NEWFORCPATH}/${i}")
  ncksCmd <- paste0("ncks -O -d west_east,", geo_w-1, ",", geo_e-1, " -d south_north,", geo_s-1, ",", geo_n-1, " ${OLDFORCPATH}/${i} ${NEWFORCPATH}/${i##*/}")
  fileConn <- file(subScriptFile)
  writeLines(c("#!/bin/bash",
               "OLDFORCPATH='PATH_TO_OLD_FORCING_DATA_FOLDER'",
               "NEWFORCPATH='PATH_TO_NEW_FORCING_DATA_FOLDER'",
               "for i in `ls $OLDFORCPATH`; do",
               "echo ${i##*/}",
               ncksCmd,
               "done"),
             fileConn)
  close(fileConn)
  
  
  ncksCmd <- paste0("ncks -O -d x,", geo_w-1, ",", geo_e-1, " -d y,", geo_s-1, ",", geo_n-1, " ${OLDFORCPATH}/${i} ${NEWFORCPATH}/${i##*/}")
  fileConn <- file(paste0(subScriptFile, "_NWM_REALTIME"))
  writeLines(c("#!/bin/bash",
               "OLDFORCPATH='PATH_TO_OLD_FORCING_DATA_FOLDER'",
               "NEWFORCPATH='PATH_TO_NEW_FORCING_DATA_FOLDER'",
               "for i in `ls $OLDFORCPATH`; do",
               "echo ${i##*/}",
               ncksCmd,
               "done"),
             fileConn)
  close(fileConn)
  
  # copy this script to the myPath dir, so the user haev from what files this cutout has been generated, and what are the options.
  
  file.copy(from = paste0(getwd(), "/subset_domain.R"), to = paste0(myPath, "/subset_domain.R"), overwrite = TRUE)
}

foreach (gage = unique(infoDT$.id), .export=ls(.GlobalEnv), .packages = c("data.table", "rwrfhydro", "ncdf4", "raster")) %do% {CutOut(gage)}

