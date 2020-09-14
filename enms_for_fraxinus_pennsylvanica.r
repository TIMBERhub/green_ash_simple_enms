### modeling green ash
### Adam B. Smith | adam.smith@mobot.org | Missouri Botanical Garden | Spring 2019
###
### INTRODUCTION
### The goal of this project is to generate a "bank" of ~10 ENMs for Green Ash (*Fraxinus pennsylvanica*) for use in TIMBER. The models will be projected back in time to ~21 Kypb.
### 
### We will define models based on combinations of:
### 
### * Occurrence set (broad, narrow)  
### * SDM algorithm (Maxent, GLM)  
### * Global circulation model (ECBilt, CCSM)
### 
### for a total of 2 * 2 * 2 = 8 models. This script draws on several data sources that are not downloaded or generated in the script:
### * Paleo and current climate data from Lorenz et al 2016 Scientific Data
### * Little's range map of Fraxinus pennsylvanica
### 
### CONTENTS
### setup ###
### constants ###
### obtain and process basemap data ###
### obtain and clean occurrences ###
### define calibration regions ###
### select climate variables ###
### collate calibration/evaluation data ###
### create plot demonstrating local SAC and LOO folds ###
### calibrate models using leave-one-out cross-validation ###
### calibrate final models ###
### project final models ###
### convert prediction rasters to matrices for TIMBER ###
### figure of predictions ###
### analyze similarity between model projections ###
### calculate biotic velocity ###
### plot potential total population size ###

### NOTES
### Variable naming conventions:
### * objects that are spatial (ie polygons/points) will be suffixed with "Sp" as in "speciesRecordsSp"
### * spatial objects that are in an equal-area projection will be further suffixed with "Alb" as in "speciesRecordsSpAlb"

#############
### setup ###
#############

	# source('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/code/enms_for_fraxinus_pennsylvanica.r')

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	options(keep.source=FALSE) # manage memory
	gc()
	print('')
	print(date())

	options(stringsAsFactors=FALSE)
	raster::rasterOptions(format='GTiff', overwrite=TRUE)

	library(BIEN)
	library(sp)
	library(rgdal)
	library(raster)
	library(RColorBrewer)
	library(geosphere)
	library(rgeos)
	library(dismo)
	library(blockCV)
	library(scales)
	library(omnibus) # Adam's grab-bag library (https://github.com/adamlilith/omnibus)
	library(enmSdm) # Adam's SDM library (https://github.com/adamlilith/enmSdm)
	library(statisfactory) # Adam's statistics library (https://github.com/adamlilith/statisfactory)
	library(legendary) # Adam's plotting library (https://github.com/adamlilith/legenday)
	
	setwd('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms')

	dirCreate('./figures_and_tables')
	dirCreate('./species_records')
	dirCreate('./regions')
	dirCreate('./environmental_rasters')
	dirCreate('./environmental_rasters/lorenz_et_al_2016')
	dirCreate('./models')
	dirCreate('./predictions')

	### return stack of climate rasters for North America from Lorenz et al. 2016 Scientific Reports data set
	#########################################################################################################
	# note that these are assumed to live in a folder outside the folder for green ash
	
	getClimRasts <- function(gcm, year, variables, rescale = TRUE, fillCoasts = FALSE) {

		# gcm		'ccsm' or 'ecbilt'
		# year		year BP (from 0 to 21000 for ccsm or 22000 for ecbilt)
		# variables names of variables
		# rescale	TRUE ==> rescale rasters to [0, 1] using present-day values for min/max
		# fillCoasts FALSE ==> use rasters as-is; TRUE ==> extrapolate to NA cells immediately adjacent to non-NA cells (typically coastal cells)
	
		gcmFolder <- if (gcm == 'ccsm') {
			'ccsm3_22-0k_all_tifs'
		} else if (gcm == 'ecbilt') {
			'ecbilt_21-0k_all_tifs'
		}
	
		# get current version of each variable
		for (variable in variables) {
		
			rast <- stack(paste0('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/V2/', gcmFolder, '/', year, 'BP/', variable, '.tif'))
			names(rast) <- variable

			rasts <- if (exists('rasts', inherits=FALSE)) {
				stack(rasts, rast)
			} else {
				rast
			}
		
		}
		
		# rescale
		if (rescale) {

			# not present... rescale by present-day values
			if (year != 0 | file.exists ('./environmental_rasters/lorenz_et_al_2016/variable_statistics_0bp_across_north_america.csv')) {
			
				variableStats <- read.csv('./environmental_rasters/lorenz_et_al_2016/variable_statistics_0bp_across_north_america.csv')
				
				for (variable in variables) {
					minVal <- variableStats$min[variableStats$gcm == gcm & variableStats$variable == variable]
					maxVal <- variableStats$max[variableStats$gcm == gcm & variableStats$variable == variable]
					rasts[[variable]] <- (rasts[[variable]] - minVal) / (maxVal - minVal)
				}
		
			# present-day... just rescale to [0, 1]
			} else {
		
				for (i in 1:nlayers(rasts)) rasts[[i]] <- stretch(rasts[[i]], 0, 1)
				
			}
			
		}

		if (fillCoasts) {
			
			name <- names(rasts)
			
			# fill NA cells near coasts to account for fact that some records may not fall in a cell near a coast
			for (i in 1:nlayers(rasts)) {
				rasts[[i]] <- focal(rasts[[i]], w=matrix(1, nrow=3, ncol=3), fun=mean, na.rm=TRUE, NAonly=TRUE)
			}
			
			names(rasts) <- name
		}
			
		names(rasts) <- paste0(gcm, '_', variables)
		rasts
		
	}

#################
### constants ###
#################

	# species name
	species <- 'Fraxinus pennsylvanica'

	# buffer around occurrences ("accessible" area) used to define calibration area
	exts <- c(80, 160, 320) # in km

	# time periods represented by rasters
	rastTimes <- seq(0, 21000, by=500)
	
	# generation time
	genTime_yr <- 30 # generation time in year
	genTimes <- seq(0, 21000, by=genTime_yr)

	# maximum pairwise correlation allowable between variables
	maxCor <- 0.7
	
	# longitude/latitude field names
	ll <- c('longitude', 'latitude')
	
	# minimum number of occurrences in a k-fold
	minFoldSize <- 200
	
	# buffer around each training/test fold to be excluded (in meters)
	foldBuffer_m <- 8 * 40 * 1000
	
	# GCMs
	gcms <- c('ccsm', 'ecbilt')
	
	# decided in "### process current/paleoclimate layers ###"
	if (file.exists('./figures_and_tables/predictors.rda')) load('./figures_and_tables/predictors.rda')

	# study region polygon
	studyRegionSpAlb <- shapefile('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/study_region_spatial_polygons/study_region_mask_without_glaciers')
	
	# mask raster for present-day land (1 = land, NA = not)
	if (file.exists('./regions/maskRaster.tif')) {
	
		maskRast <- raster('./regions/maskRaster.tif')
		
	} else {
	
		# mask raster for entire projection region
		maskRast <- raster('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/V2/ccsm3_22-0k_all_tifs/0BP/an_avg_ETR.tif')
		maskRast <- 0 * maskRast + 1
		names(maskRast) <- 'mask'
		
		writeRaster(maskRast, './regions/maskRaster', datatype='INT2U')
		
	}

	maskRastAlb <- projectRaster(maskRast, crs=getCRS('albersNA'))
	
	# states that are undersampled for Fraxinus pennsylvanica species in BIEN 4.1
	undersampled <- c('Louisiana', 'Kentucky')

	# calibration regions
	if (file.exists('./regions/calibration_regions.rda')) load('./regions/calibration_regions.rda')
	
	# North America spatial polygons
	if (file.exists('./regions/greatLakesSpAlb.rda')) load('./regions/greatLakesSpAlb.rda')

	if (file.exists('./regions/gadm36_northAmerica_sp_level_0.rda')) load('./regions/gadm36_northAmerica_sp_level_0.rda')
	if (file.exists('./regions/gadm36_northAmerica_spAlb_level_0.rda')) load('./regions/gadm36_northAmerica_spAlb_level_0.rda')
	
	if (file.exists('./regions/gadm36_northAmerica_sp_level_1.rda')) load('./regions/gadm36_northAmerica_sp_level_1.rda')
	if (file.exists('./regions/gadm36_northAmerica_spAlb_level_1.rda')) load('./regions/gadm36_northAmerica_spAlb_level_1.rda')
	
	# range maps
	if (file.exists('./regions/bien_range_map/Fraxinus_pennsylvanica.shp')) {
		bienRange <- shapefile('./regions/bien_range_map/Fraxinus_pennsylvanica.shp')
		bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA', TRUE))
	}
	
	if (file.exists('./regions/little_range_map/littleRangeSpAlb.rda')) load('./regions/little_range_map/littleRangeSpAlb.rda')
	
say('#######################################', pre=1)
say('### obtain and process basemap data ###')
say('#######################################', post=2)

	say('Basemap data includes shapefiles of North American countries and the range maps of the species from BIEN and Little. We will create equal-area versions of each in an Albers projection (object names will be suffixed with "Alb")', breaks=80)

	### North America Level 2
		
		mex <- getData('GADM', country='MEX', level=2, path='D:/ecology/!Scratch')
		usa <- getData('GADM', country='USA', level=2, path='D:/ecology/!Scratch')
		can <- getData('GADM', country='CAN', level=2, path='D:/ecology/!Scratch')

		nam2Sp <- rbind(can, usa, mex)

		# get lakes for removal from other geographies... add small buffer because lake borders don't exactly align
		greatLakesSp2 <- nam2Sp[nam2Sp$ENGTYPE_2 == 'Water body', ]
		greatLakesSpAlb2 <- sp::spTransform(greatLakesSp2, getCRS('albersNA', TRUE))
		greatLakesSpAlb2 <- gBuffer(greatLakesSpAlb2, width=10)
		greatLakesSp2 <- sp::spTransform(greatLakesSpAlb2, getCRS('wgs84', TRUE))
		greatLakesSp <- gUnaryUnion(greatLakesSp2)
		
		greatLakesSpAlb <- sp::spTransform(greatLakesSp2, getCRS('albersNA', TRUE))
		
		# names of level 1 areas with lakes
		lakeNames0 <- unique(nam2Sp@data$NAME_0[nam2Sp@data$ENGTYPE_2 == 'Water body'])
		lakeNames1 <- unique(nam2Sp@data$NAME_1[nam2Sp@data$ENGTYPE_2 == 'Water body'])
		
		# remove lakes
		nam2Sp <- nam2Sp[nam2Sp@data$ENGTYPE_2 != 'Water body', ]
		
		# not saving because only need for identifying lake areas
		
		save(greatLakesSpAlb, file='./regions/greatLakesSpAlb.rda')

	### North America Level 1
	
		mex <- getData('GADM', country='MEX', level=1, path='D:/ecology/!Scratch')
		usa <- getData('GADM', country='USA', level=1, path='D:/ecology/!Scratch')
		can <- getData('GADM', country='CAN', level=1, path='D:/ecology/!Scratch')

		nam1Sp <- rbind(can, usa, mex)
		
		for (thisLevel in lakeNames1) {
		
			thisLargerSp <- nam1Sp[nam1Sp@data$NAME_1 == thisLevel, ]
			df <- thisLargerSp@data
			nam1Sp <- nam1Sp[-which(nam1Sp@data$NAME_1 == thisLevel), ]
			thisLargerSansLakesSp <- gDifference(thisLargerSp, greatLakesSp)
			thisLargerSansLakesSp <- as(thisLargerSansLakesSp, 'SpatialPolygonsDataFrame')
			projection(thisLargerSansLakesSp) <- projection(nam1Sp)
			thisLargerSansLakesSp@data <- df
			nam1Sp <- rbind(nam1Sp, thisLargerSansLakesSp)
			
		}

		nam1SpAlb <- sp::spTransform(nam1Sp, getCRS('albersNA', TRUE))
		save(nam1Sp, file='./regions/gadm36_northAmerica_sp_level_1.rda')
		save(nam1SpAlb, file='./regions/gadm36_northAmerica_spAlb_level_1.rda')
	
	### North America level 0

		mex <- getData('GADM', country='MEX', level=0, path='D:/ecology/!Scratch')
		usa <- getData('GADM', country='USA', level=0, path='D:/ecology/!Scratch')
		can <- getData('GADM', country='CAN', level=0, path='D:/ecology/!Scratch')

		nam0Sp <- rbind(can, usa, mex)

		for (thisLevel in lakeNames0) {
		
			thisLargerSp <- nam0Sp[nam0Sp@data$NAME_0 == thisLevel, ]
			df <- thisLargerSp@data
			nam0Sp <- nam0Sp[-which(nam0Sp@data$NAME_0 == thisLevel), ]
			thisLargerSansLakesSp <- gDifference(thisLargerSp, greatLakesSp)
			thisLargerSansLakesSp <- as(thisLargerSansLakesSp, 'SpatialPolygonsDataFrame')
			projection(thisLargerSansLakesSp) <- projection(nam1Sp)
			thisLargerSansLakesSp@data <- df
			nam0Sp <- rbind(nam0Sp, thisLargerSansLakesSp)
			
		}

		nam0SpAlb <- sp::spTransform(nam0Sp, getCRS('albersNA', TRUE))
		save(nam0Sp, file='./regions/gadm36_northAmerica_sp_level_0.rda')
		save(nam0SpAlb, file='./regions/gadm36_northAmerica_spAlb_level_0.rda')
			

	### range maps
	dirCreate('./regions/bien_range_map')
	BIEN_ranges_species(species, directory='./regions/bien_range_map')
	bienRange <- BIEN_ranges_load_species(species)
	bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA', TRUE))

	littleRangeFileName <- 'C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/range_maps/little/fraxpenn.shp'
	littleRange <- shapefile(littleRangeFileName)
	projection(littleRange) <- getCRS('wgs84')
	littleRange <- littleRange[littleRange$CODE == 1, ] # remove holes
	littleRangeSpAlb <- sp::spTransform(littleRange, getCRS('albersNA', TRUE))
	dirCreate('./regions/little_range_map')
	save(littleRangeSpAlb, file='./regions/little_range_map/littleRangeSpAlb.rda')
			
say('####################################', pre=1)
say('### obtain and clean occurrences ###')
say('####################################', post=2)

	### obtain records
	##################
		
		speciesFileName <- paste0(
			'./species_records/00_',
			gsub(tolower(species), pattern=' ', replacement='_'),
			'_bien_all_occurrences.rda'
		)

		if (file.exists(speciesFileName)) {
			load(speciesFileName)
		} else {
			occsRaw <- BIEN_occurrence_species(
				species=species,
				cultivated = FALSE,
				only.new.world = TRUE,
				all.taxonomy = FALSE,
				native.status = FALSE,
				natives.only = TRUE,
				observation.type = TRUE,
				political.boundaries = TRUE,
				collection.info = TRUE
			)

			sink('./species_records/bien_download.txt')
				say('Data representing species records were downloaded from BIEN on', date(), '.')
			sink()
			save(occsRaw, file=speciesFileName)
		}

		### remove occurrences with missing coordinates and dates:
		occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
		dim(occs)

		### remove records <1950 (period covered by Lorenz et al. 2016 climate data is 1950-2005):
		occs$date_collected <- as.Date(occs$date_collected)
		occs <- occs[!is.na(occs$date_collected), ]
		occs <- occs[which(occs$date_collected >= as.Date(paste0(1950 + genTime_yr, '-01-01'))), ]

	### remove records in states where green ash has been naturalized (according to BONAP)
	######################################################################################
	
		occs <- occs[!(occs$state_province %in% c('Washington', 'Oregon', 'Idaho', 'Utah', 'Colorado', 'Arizona', 'New Mexico', 'Chihuahua', 'British Columbia')), ]
		occs <- occs[-which(occs$state_province == 'Texas' & occs$county == 'El Paso'), ]
		
	### inspect occurrences
	#######################
		
		### convert occurrences to spatial object
		occsSp <- SpatialPointsDataFrame(
			occs[ , ll],
			data=occs,
			proj4=getCRS('wgs84', TRUE)
		)

		occsSpAlb <- sp::spTransform(occsSp, getCRS('albersNA', TRUE))

		png('./figures_and_tables/occurrences_raw.png', width=1200, height=1000)
		
			par(mfrow=c(1, 2))
			
			# plot occurrences with BIEN range map
			plot(occsSpAlb, col=NA, main='BIEN range map')
			plot(bienRangeAlb, col=alpha('green', 0.4), border='green', add=TRUE)
			points(occsSpAlb, pch=16, cex=0.9, col='black')
			plot(nam1SpAlb, add=TRUE, border='gray')

			# plot occurrences with Little's range map
			plot(occsSpAlb, col=NA, main='Little range map')
			plot(littleRangeSpAlb, col=alpha('red', 0.4), border='red', add=TRUE)
			points(occsSpAlb, pch=16, cex=0.9, col='black')
			plot(nam1SpAlb, add=TRUE, border='gray')

			title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)
			
		dev.off()

		say('Things to note:')
		say('* There is at least one occurrence in South America or its vicinity.')
		say('* There are many occurrences outside Little\'s range polygon.')
		say('* ...Especially in the northern part of peninsular Florida', post=2)
		say('* ...and the western US', post=2)
		
		### examine areas in vicinity of Kentucky and Louisiana that appear to lack data

		png('./figures_and_tables/occurrences_near_louisiana_kentucky.png', width=1200, height=1000)
		
			par(mfrow=c(1, 2))

			plot(nam1SpAlb[nam1SpAlb@data$NAME_1 == 'Kentucky', ])
			points(occsSpAlb, pch=16, cex=0.9, col='black')

			plot(nam1SpAlb[nam1SpAlb@data$NAME_1 == 'Louisiana', ])
			points(occsSpAlb, pch=16, cex=0.9, col='black')
			
			title(main=paste0('Missing records in ', paste(undersampled, collapse=' and ')), sub=date(), outer=TRUE, cex.sub=0.7, line=-2)

		dev.off()
			
		say('So, yes, there is a bias in the central and Appalachian region of Kentucky... Note that the species genuinely does not appear in the "high" Appalachians. This is disturbing since the usual methods of correcting for bias like inverse-p weighting or spatial thinning would remove signal from its absence in these areas *plus* in areas that are undersampled like Kentucky and Louisiana.', post=2, breaks=80)

		say('So, we are going to remove Louisiana and Kentucky from the training set. This means occurrences from these areas will be be use to train the models, nor will background sites drawn from these areas.', post=2, breaks=80)

	### remove occurrences in undersampled areas and in South America
	#################################################################
	
		nrowStart <- nrow(occsSpAlb)
		occsSpAlb <- occsSpAlb[(occsSpAlb@data$country %in% c('Canada', 'United States', 'Mexico')) & !(occsSpAlb@data$state_province %in% undersampled), ]
		
		sink('./species_records/cleaning_notes_outliers_undersampled.txt', split=TRUE)
			say('By removing occurrences outside Mexico/US/Canada we reduced the number of occurrences from ', nrowStart, ' to ', nrow(occsSpAlb), '.')
		sink()

	### correct for hypothesized spatial sampling bias
	##################################################
		
		say('To remove sampling bias we will spatially thin occurrences so that there is but one occurrence per cell.')

		# get occurrences in WGS84 CRS
		priority <- rep(NA, nrow(occsSpAlb))
		priority[occsSpAlb$observation_type == 'plot'] <- 1
		priority[occsSpAlb$observation_type == 'specimen'] <- 2
		priority[occsSpAlb$observation_type == 'trait occurrence'] <- 3
		priority[occsSpAlb$observation_type == 'human observation'] <- 4
		priority[occsSpAlb$observation_type == 'unknown'] <- 5
		
		nrowStart <- nrow(occsSpAlb)
		occsSpAlb <- elimCellDups(occsSpAlb, r=maskRastAlb, priority=priority)

		sink('./species_records/cleaning_notes_thinning.txt', split=TRUE)
			say('Through spatial thinning we reduced the number of occurrences from ', nrowStart, ' to ', nrow(occsSpAlb), '.')
		sink()
		
		# inspect
		png('./figures_and_tables/occurrences_thinned_one_per_cell.png', width=1200, height=1000)
		
			plot(occsSpAlb, pch=16, cex=1.8, main='All cleaned and thinnned occurrences')
			plot(nam1SpAlb, add=TRUE, border='gray')
			
			title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)
			
		dev.off()
		
		# save
		save(occsSpAlb, file=paste0('./species_records/01_', (gsub(tolower(species), pattern=' ', replacement='_')), '_bien_cleaned_occurrences_thinned.rda'))
	
say('##################################', pre=1)
say('### define calibration regions ###')
say('##################################', post=2)	
	
	say('Currently, best practices assert that one should model the niche using only data from an area considered accessible to the species over long time. However, in practice this is very difficult to estimate, so we will use buffers of different size (defined above in the "### constants ###" section above).', post=2, breaks=80)

	### load data
	
	load(paste0('./species_records/01_', gsub(tolower(species), pattern=' ', replacement='_'), '_bien_cleaned_occurrences_thinned.rda'))
	
	### create buffers
	
	# buffer around occurrences
	calibRegionsSpAlb <- list()
	for (i in seq_along(exts)) {
		calibRegionsSpAlb[[i]] <- gBuffer(occsSpAlb, width=exts[i] * 1000)
	}

	# crop area outside North America, remove undersampled states, remove Great Lakes
	undersampledSpAlb <- nam1SpAlb[nam1SpAlb$NAME_1 %in% undersampled, ]
	
	for (i in seq_along(exts)) {
		# calibRegionsSpAlb[[i]] <- gIntersection(calibRegionsSpAlb[[i]], sr) # not doing this because some points outside study region
		calibRegionsSpAlb[[i]] <- rgeos::gIntersection(calibRegionsSpAlb[[i]], nam0SpAlb)
		calibRegionsSpAlb[[i]] <- rgeos::gDifference(calibRegionsSpAlb[[i]], undersampledSpAlb, checkValidity=TRUE)
		calibRegionsSpAlb[[i]] <- rgeos::gDifference(calibRegionsSpAlb[[i]], greatLakesSpAlb, checkValidity=TRUE)
	}
	
	studyRegionCropSpAlb <- crop(studyRegionSpAlb, nam0SpAlb)
	
	# inspect
	png('./figures_and_tables/regions.png', width=1200, height=1000)

		par(cex.main=1.4)
		plot(calibRegionsSpAlb[[length(calibRegionsSpAlb)]], border=NA, col=NA, main='Calibration Regions')
		plot(studyRegionCropSpAlb, add=TRUE, lwd=0.8, col='khaki1')
		
		for (i in rev(seq_along(exts))) {
			plot(calibRegionsSpAlb[[i]], col=paste0('aquamarine', 4 + 1 - i), lwd=0.8, add=TRUE)
		}
		
		plot(nam1SpAlb, add=TRUE)
		points(occsSpAlb, pch=16, cex=1.4)

		legend <- c('Records', paste('Calibration region:', rev(exts), 'km'), 'Study region (modern)')
		fill <- c(NA, paste0('aquamarine', 4 + 1 - rev(seq_along(exts))), 'khaki1')
		pch <- c(16, rep(NA, length(exts)), NA)
		pt.cex <- c(1.4, rep(NA, length(exts)), NA)
		border <- c(NA, rep(0.8, length(exts)), 0.8)
		col <- c('black', rep(NA, length(exts)), NA)
		
		legend('bottomright',
			   legend=legend,
			   fill=fill,
			   pch=pch,
			   pt.cex=pt.cex,
			   border=border,
			   col=col,
			   bty='n', cex=1.6
		)
		
	dev.off()

	names(calibRegionsSpAlb) <- paste0('bufferAroundOccs_', exts, 'km')
	save(calibRegionsSpAlb, file='./regions/calibration_regions.rda')

say('################################', pre=1)
say('### select climate variables ###')
say('################################', post=2)

	say('We will be using the paleo climate layers prepared by Lorenz, D.J., Nieto-Lugilde, D., Blois, J.L., Fitzpatrick, M.C., and Williams, J.W.  2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100 AD.  Scientific Data 3:160048. The DOI of the original data is 10.5061/dryad.1597g. This data set contains future, current, and historical GCM projections. We will be using the current and historical model output from the ECBilt and CCSM GCMs (other models have only current/future projections associated with them). The spatial resolution is 0.5deg for all of North America.', breaks=80, post=2)

	say('Available variables include min/max temperature, precipitation, growing degree days, AET, PET, ETR (= AET/PET), and WDI (PET - precipitation), summarized across an annual, quarterly, or monthly basis. Variables are projected back to 22 Kybp (CCSM) or 21 Kybp (ECBilt) in 500-yr timeslices. The current period covers averages across 1950-2005. Past climate encompasses averages across 200-yr intervals centered on every 500 yr (e.g., 500 ypb, 1000 ypb, 1500 ypb, etc.).', breaks=80, post=2)
	
	say('We will *not* use variable "an_cv_ETR" because it has NAs in Canada during the present in large portions of Canada (and possibly through time).', breaks=80, post=2)
	say('We will *not* use variable "an_cv_WDI" because it has some very extreme values in just a few cells in Mexico at 0 ybp (and possibly through time).', breaks=80, post=2)
	
	say('Tasks to complete include:')
	say('* Selecting variables.')
	say('* Cropping/masking rasters to study region *extent* (versus the areas within the study region).')
	say('* Interpolating rasters to 30-yr intervals (approximate generation time for green ash).')

	say('### select variables', level=2)
	####################################
		
		calibRegionBroadestSp <- calibRegionsSpAlb[[paste0('bufferAroundOccs_', max(exts), 'km')]]
		
		### get set of randomly located points for assessing collinearity between variables
		
		climCcsm_0ypb <- raster::stack(listFiles('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/V2/ccsm3_22-0k_all_tifs/0BP', pattern='an_'))
		climEcbilt_0ypb <- raster::stack(listFiles('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/V2/ecbilt_21-0k_all_tifs/0BP', pattern='an_'))
		
		# remove an_cv_ETR and an_cv_WDI
		keeps <- names(climCcsm_0ypb)[-which(names(climCcsm_0ypb) %in% c('an_cv_ETR', 'an_cv_WDI'))]
		climCcsm_0ypb <- subset(climCcsm_0ypb, keeps)
		climEcbilt_0ypb <- subset(climEcbilt_0ypb, keeps)
		
		maskBroad0ypb <- rasterize(calibRegionBroadestSp, climCcsm_0ypb[[1]])
		
		randBg <- randomPoints(maskBroad0ypb, 11000)
		randBg <- as.data.frame(randBg)
		names(randBg) <- ll
		
		envCcsm_0ybp <- extract(climCcsm_0ypb, randBg)
		envEcbilt_0ybp <- extract(climEcbilt_0ypb, randBg)

		completeCases <- complete.cases(envCcsm_0ybp) & complete.cases(envEcbilt_0ybp)
		
		randBg <- randBg[completeCases, ]
		randBg <- randBg[1:max(nrow(randBg), 10000), ]
		envCcsm_0ybp <- envCcsm_0ybp[completeCases, ]
		envEcbilt_0ybp <- envEcbilt_0ybp[completeCases, ]
		
		env <- rbind(envCcsm_0ybp, envEcbilt_0ybp)
		
		env <- as.data.frame(env)
		envScaled <- scale(env)
		variableSel <- usdm::vifcor(envScaled, maxCor, maxobservations=Inf)

		sink('./figures_and_tables/variable_selection.txt', split=TRUE)
		
			say('Based on a threshold maximum correlation of ', maxCor, ' and VIF, we will use the following variables:')
			print(variableSel)
			
		sink()
		
		predictors <- variableSel@results$Variables
		save(predictors, file='./figures_and_tables/predictors.rda')
		
	say('calculate statistics for rescaling rasters from present-day rasters', level=2)
	###################################################################################
		
		# saves statistics on each variable (min/max/etc)
		variableStats <- data.frame()
		
		for (gcm in gcms) {
		
			rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale = FALSE)
			mins <- minValue(rasts)
			maxs <- maxValue(rasts)
			
			# remember
			variableStats <- rbind(
				variableStats,
				data.frame(
					gcm = gcm,
					variable = predictors,
					min = mins,
					max = maxs
				)
			)
	
		} # next GCM
				
		write.csv(variableStats, './environmental_rasters/lorenz_et_al_2016/variable_statistics_0bp_across_north_america.csv', row.names=FALSE)

say('###########################################', pre=1)
say('### collate calibration/evaluation data ###')
say('###########################################', post=2)

	say('In this step we will:')
	say('* Extract environmental data.')
	say('* Locate background sites.')
	say('* Define geographically distinct cross-validation folds for model evaluation.')

	calibRegionsSp <- calibRegionsSpAlb
	for (i in seq_along(calibRegionsSp)) {
		calibRegionsSp[[i]] <- sp::spTransform(calibRegionsSp[[i]], getCRS('wgs84', TRUE))
	}
	
	load(paste0('./species_records/01_', gsub(tolower(species), pattern=' ', replacement='_'), '_bien_cleaned_occurrences_thinned.rda'))
	
		occsSp <- sp::spTransform(occsSpAlb, getCRS('wgs84', TRUE))
	
	say('extract environmental data', level=2)
	##########################################

		# by GCM
		for (gcm in gcms) {

			if (exists('rasts', inherits=FALSE)) rm(rasts)
		
			# get current version of each variable raster
			rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
			areas <- area(rasts[[1]])

			# extract
			env <- raster::extract(rasts, occsSp)
			env <- as.data.frame(env)
			
			thisArea <- raster::extract(areas, occsSp)
			thisArea <- as.data.frame(thisArea)
			names(thisArea) <- 'cellArea_km2'
				
			thisEnv <- cbind(thisArea, env)
			
			# remember
			occsSpAlb@data <- cbind(occsSpAlb@data, thisEnv)
			
		} # next GCM
		
		# remove occurrences with incomplete environmental data
		vars <- character()
		for (gcm in gcms) vars <- c(vars, paste0(gcm, '_', predictors))
		noNas <- complete.cases(occsSpAlb@data[ , vars])
		occsSpAlb <- occsSpAlb[noNas, ]

	say('generate test background sites', level=2)
	##############################################
	
		calibRegionSpAlb <- calibRegionsSpAlb[[paste0('bufferAroundOccs_', max(exts), 'km')]]
		
		# get random background sites
		bgTestSpAlb <- spsample(calibRegionSpAlb, n=11000, type='random', iter=10)
		bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
		bgTest <- as.data.frame(coordinates(bgTestSp))
		names(bgTest) <- ll
		
		if (exists('env')) rm(env)

		# by GCM
		for (gcm in gcms) {

			# get current version of each variable raster
			rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
			areas <- area(rasts[[1]])

			# extract
			thisEnv <- raster::extract(rasts, bgTest)
			thisEnv <- as.data.frame(thisEnv)
			
			thisArea <- raster::extract(areas, bgTest)
			thisArea <- as.data.frame(thisArea)
			names(thisArea) <- 'cellArea_km2'
			
			thisEnv <- cbind(thisArea, thisEnv)
			
			# remember
			env <- if (exists('env', inherits=FALSE)) {
				cbind(env, thisEnv)
			} else {
				thisEnv
			}
			
		} # next GCM

		# cull NA cases
		bgTest <- cbind(bgTest, env)
		bgTest <- bgTest[complete.cases(bgTest), ]
		bgTest <- bgTest[1:min(nrow(bgTest), 10000), ]
			
	say('generate training background sites', level=2)
	##################################################

		for (countExt in seq_along(exts)) {
				
			ext <- exts[countExt]
			calibRegionSpAlb <- calibRegionsSpAlb[[countExt]]
			
			# get random background sites
			bgCalibSpAlb <- spsample(calibRegionSpAlb, n=11000, type='random', iter=10)
			bgCalibSp <- sp::spTransform(bgCalibSpAlb, getCRS('wgs84', TRUE))
			bgCalib <- as.data.frame(coordinates(bgCalibSp))
			names(bgCalib) <- ll
			
			if (exists('env')) rm(env)

			# by GCM
			for (gcm in gcms) {

				# get current version of each variable raster
				rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
				areas <- area(rasts[[1]])

				# extract
				thisEnv <- raster::extract(rasts, bgCalib)
				thisEnv <- as.data.frame(thisEnv)
				
				thisArea <- raster::extract(areas, bgCalib)
				thisArea <- as.data.frame(thisArea)
				names(thisArea) <- 'cellArea_km2'
				
				thisEnv <- cbind(thisArea, thisEnv)
				
				# remember
				env <- if (exists('env', inherits=FALSE)) {
					cbind(env, thisEnv)
				} else {
					thisEnv
				}
				
			} # next GCM

			# cull NA cases
			bgCalib <- cbind(bgCalib, env)
			bgCalib <- bgCalib[complete.cases(bgCalib), ]
			bgCalib <- bgCalib[1:min(nrow(bgCalib), 10000), ]
			
			# remember
			assign(paste0('bgCalib_extent', exts[countExt], 'km'), bgCalib)
			
		} # next calibration region

	say('assign spatially exclusive training/calibration/evaluation cross-validation folds', level=2)
	#################################################################################################
		
		bgTestSp <- SpatialPoints(bgTest[ , ll], getCRS('wgs84', TRUE))
		bgTestSpAlb <- sp::spTransform(bgTestSp, getCRS('albersNA', TRUE))
		
		### !!!!! NEXT LINE CAN TAKE HOURS !!!!! ###
		studyRegionCropSpAlb <- crop(studyRegionSpAlb, nam0SpAlb)
		
		blockRast <- rastWithSquareCells(studyRegionCropSpAlb, res=foldBuffer_m)

		# create folds based on occurrences such that there is at least a minimum number of occurrences per fold
		set.seed(123)
		minNumInFold <- -Inf
		numFolds <- 50
		maxTries <- 20
		tries <- 1
		
		while (minNumInFold < minFoldSize) {
		
			values(blockRast) <- sample(1:numFolds, ncell(blockRast), replace=TRUE)
			folds <- extract(blockRast, occsSpAlb)
			minNumInFold <- min(table(folds))
			
			tries <- tries + 1
			if (tries == maxTries) {
				tries <- 1
				numFolds <- numFolds - 1
			}
			
		}
	
		say('Final number of folds: ', numFolds)
		
		occs <- occsSpAlb@data
		folds <- data.frame(fold=folds)
		occs <- cbind(occs, folds)
		
		# assign folds to test background sites
		bgTestSp <- SpatialPointsDataFrame(bgTest[ , ll], data=bgTest, proj4string=getCRS('wgs84', TRUE))
		bgTestSpAlb <- sp::spTransform(bgTestSp, getCRS('albersNA', TRUE))
		bgTestFolds <- extract(blockRast, bgTestSpAlb)
		bgTestSpAlb$fold <- bgTestFolds
		bgTest <- bgTestSpAlb@data
		
		# assign folds to calibration background sites
		for (ext in exts) {
		
			x <- get(paste0('bgCalib_extent', ext, 'km'))
			xSp <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			xSpAlb <- sp::spTransform(xSp, getCRS('albersNA', TRUE))
			bgCalibFolds <- extract(blockRast, xSpAlb)
			x$fold <- bgCalibFolds
			x <- as.data.frame(x)
			
			assign(paste0('bgCalib_extent', ext, 'km'), x)
			
		}

	say('collate calibration and evaluation	occurrences and background sites', level=2)
	##################################################################################
	
		occsBg <- list()
		occsBg$note <- 'This list object contains spatially thinned occurrence data and associated background sites. $testBg is a set of background sites drawn from the largest training region extent. $calibEvalOccsBg[[i]] is a list of lists, one per training region extent. Each sublist has two data frames. $calibEvalOccsBg[[i]]occsBg is a data frame representing calibration occurrences and background sites. Weights are 1 for occurrences and a fractional value for background sites such that their combined weight is equal to the total weight of occurrences. Folds are geographically exclusive and comprise a set of intermingled spatial blocks. $calibEvalOccsBg[[i]]$folds is a data frame with one row per row in $calibEvalOccsBg[[i]]$occsBg and one column per fold. Values in each column are 1 (training data), 2 (calibration data), or NA (masked... these represent evaluation data censored from the training and calibration sets).'
		
		occsBg$meta <- data.frame(
			item = c(
				'numOccs',
				'numFolds',
				'minFoldSize_numOccs'
			),
			value = c(
				nrow(occs),
				numFolds,
				minFoldSize
			)
		)
		
		columnsJustBg <- c('fold', ll, 'cellArea_km2', paste0(rep(gcms, each=length(predictors)), '_', predictors))
		columnsOccsAndBg <- c('presBg', 'fold', 'weight', ll, 'cellArea_km2', paste0(rep(gcms, each=length(predictors)), '_', predictors))
		
		occsBg$allOccsThinned <- occs
		occsBg$testBg <- bgTest[ , columnsJustBg]
		occsBg$calibEvalOccsBg <- list()
		
		for (countExents in seq_along(exts)) {
		
			ext <- exts[countExents]
		
			# collate occurrences and background sites
			thisOccs <- occs
			thisOccs$presBg <- 1
			thisOccs$weight <- 1
			thisOccs <- thisOccs[ , columnsOccsAndBg]
			
			thisBg <- get(paste0('bgCalib_extent', ext, 'km'))
			thisBg$presBg <- 0
			thisBg$weight <- nrow(thisOccs) / nrow(thisBg)
			thisBg <- thisBg[ , columnsOccsAndBg]
			
			thisOccsBg <- rbind(thisOccs, thisBg)
			
			# make matrix indicating which occurrences/background sites are in which fold set
			if (exists('folds')) rm(folds)
			
			for (calibFold in 1:numFolds) {
			
				for (evalFold in (1:numFolds)[-calibFold]) {
			
					designation <- rep(1, nrow(thisOccsBg)) # start with all training
					designation[thisOccsBg$fold == calibFold] <- 2 # calibration folds
					designation[thisOccsBg$fold == evalFold] <- NA # mask evaluation folds
					
					designation <- matrix(designation, ncol=1)
					colnames(designation) <- paste0('calib', calibFold, '_vs_eval', evalFold)
					
					folds <- if (exists('folds')) {
						cbind(folds, designation)
					} else {
						designation
					}
			
				}
			
			}
			
			occsBg$calibEvalOccsBg[[countExents]] <- list()
			occsBg$calibEvalOccsBg[[countExents]]$occsBg <- thisOccsBg
			occsBg$calibEvalOccsBg[[countExents]]$folds <- folds
		
		}
		
		names(occsBg$calibEvalOccsBg) <- paste0('occsBg_extent', exts, 'km')
	
		### make a map of calibration/evaluation sites

		calibFold <- 1
		evalFold <- 2
		cex <- 0.6
		pch <- 3

		png('./figures_and_tables/example_of_training_calibration_evaluation_folds.png', width=1800, height=600)
		
			par(mfrow=c(1, 3), oma=rep(0, 4), cex.main=2.2)
			
			plot(studyRegionCropSpAlb, main=('Training, calibration,\nand evaluation occurrences'))
			x <- occsBg$calibEvalOccsBg[[1]]$occsBg
			x <- x[x$presBg == 1, ]
			x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			points(x, pch=pch, cex=cex)
			points(x[x$fold==calibFold, ], pch=pch, cex=cex, col='cyan')
			points(x[x$fold==evalFold, ], pch=pch, cex=cex, col='magenta')
			
			legend('bottomright', inset=0.1, legend=c('Training', 'Calibration', 'Evaluation'), pch=pch, cex=2.2, col=c('black', 'cyan', 'magenta'), bty='n')
			
			plot(studyRegionCropSpAlb, main=('Training and calibration\nbackground sites (narrowest extent)'))
			x <- occsBg$calibEvalOccsBg[[1]]$occsBg
			x <- x[x$presBg == 0, ]
			x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			points(x, pch=pch, cex=cex)
			points(x[x$fold==calibFold, ], pch=pch, cex=cex, col='cyan')
			
			plot(studyRegionCropSpAlb, main=('Evaluation background sites'))
			x <- occsBg$testBg
			x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			points(x, pch=pch, cex=cex)
			points(x[x$fold==evalFold, ], pch=pch, cex=cex, col='magenta')
			
		dev.off()
		
		save(occsBg, file='./species_records/02_fraxinus_pennsylvanica', gsub(tolower(species), pattern=' ', replacement='_'), '_collated_occurrence_and_background_sites.rda')
		
# say('########################')
# say('### calibrate models ###')
# say('########################')

	# say('Note that the BRTs and GLMs use weighted calibration occurrences sites while MaxEnt does not. LOO calibration statistics are calculated for all models using weighted occurrences.  Background sites have weights equal to one another with the sum of weights such that it equals the weight of all occurrences.', breaks=80)

	# load(paste0('./species_records/02_', gsub(species, pattern=' ', replacement='_'), '_occurrences_background_env_kFolds_broad.rda'))
	# load(paste0('./species_records/02_', gsub(species, pattern=' ', replacement='_'), '_occurrences_background_env_kFolds_narrow.rda'))

	# metrics = c('logLossEqualWeight')
	
	# # for (ext in exts) {
	# # for (ext in exts[1]) { # broad
	# for (ext in exts[2]) { # narrow

		# # collate data
		# data <- get(paste0('occsBg', ext))
		# rownames(data) <- 1:nrow(data)
		# vars1and2 <- c(paste0('ccsm_', predictors), paste0('ecbilt_', predictors))
		# data <- data[complete.cases(data[ , vars1and2]), ]
		# data$presBg <- as.numeric(data$presBg)
		# folds <- data[ , grepl(names(data), pattern='fold')]
		# weight <- data$weight

		# # subsample folds to make modeling feasible
		# set.seed(123)
		# useFolds <- sort(sample(1:ncol(folds), round(0.2 * ncol(folds))))
		# folds <- folds[ , useFolds]
		
		# for (gcm in gcms) {
		# # for (gcm in gcms[1]) { # ccsm
		# # for (gcm in gcms[2]) { # ecbilt

			# vars <- paste0(gcm, '_', predictors)
			
# source('C:/Ecology/Drive/R/enmSdm/R/trainByCrossValid.r')
# # vars <- vars[1:2]

			# # # GLM
			# # algo <- 'glm'
			# # say(gcm, ' ', ext, ' tuning with ', algo, level=2)
			# # calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainGlm, metrics=metrics, w=weight, na.rm=TRUE, out='tuning', verbose=Inf)
			# # save(calib, file=paste0('./models/loo_calibration_for_', tolower(ext), '_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # # MaxEnt
			# # algo <- 'maxent'
			# # say(gcm, ' ', ext, ' tuning with ', algo, level=2)
			# # calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainMaxEnt, metrics=metrics, w=weight, na.rm=TRUE, out='tuning', jackknife=FALSE, scratchDir='D:/ecology/!Scratch/_temp', cores=4)
			# # save(calib, file=paste0('./models/loo_calibration_for_', tolower(ext), '_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # BRTs
			# algo <- 'brt'
			# say(gcm, ' ', ext, ' tuning with ', algo, level=2)
			# calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainBrt, metrics=metrics, w=weight, na.rm=TRUE, cores=4, out='tuning', maxTrees=8000)
			# save(calib, file=paste0('./models/loo_calibration_for_', tolower(ext), '_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
		# } # next GCM
			
	# } # next extent

# say('##############################')
# say('### calibrate final models ###')
# say('##############################')

	# # proportion of models for which a variable/feature must be included for the variable/feature to be included in the final model
	# bestThreshGlm <- 0.2
	# bestThreshMx <- 0.4

	# load(paste0('./species_records/02_', gsub(species, pattern=' ', replacement='_'), '_occurrences_background_env_kFolds_broad.rda'))
	# load(paste0('./species_records/02_', gsub(species, pattern=' ', replacement='_'), '_occurrences_background_env_kFolds_narrow.rda'))

	# sink('./models/final_model_calibration.txt', split=TRUE)
	
	# for (ext in exts) {

		# # collate data
		# data <- get(paste0('occsBg', ext))
		# folds <- data[ , grepl(names(data), pattern='fold')]
		# weight <- data$weight

		# for (gcm in gcms) {

			# vars <- paste0(gcm, '_', predictors)
			
			# ### GLM
			# say(tolower(ext), ' ', gcm, ' glm', level=1)

			# load(paste0('./models/loo_calibration_for_', tolower(ext), '_extent_with_glm_', gcm, '_gcm.rda'))
			# calibSummary <- summaryByCrossValid(calib, 'trainGlm', metric='logLossTestEqualWeight', descending=FALSE)
			
			# say('Summary of LOO calibration:')
			# print(calibSummary)

			# bestPreds <- calibSummary$term[calibSummary$proportionOfModels >= bestThreshGlm]
			# say('Best predictors:')
			# print(bestPreds)
			
			# form <- paste0('presBg ~ ', paste(bestPreds, collapse=' + '))
			# model <- glm(form, data=data, weights=weight, family='binomial')
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), '_extent_with_glm_', gcm, '_gcm.rda'))

			# ### MaxEnt
			# say(tolower(ext), ' ', gcm, ' maxent', level=1)

			# load(paste0('./models/loo_calibration_for_', tolower(ext), '_extent_with_maxent_', gcm, '_gcm.rda'))
			# calibSummary <- summaryByCrossValid(calib, 'trainMaxEnt', metric='logLossTestEqualWeight', descending=FALSE)

			# say('Summary of LOO calibration:')
			# print(calibSummary)
			
			# bestPreds <- calibSummary$term[calibSummary$proportionOfModels >= bestThreshMx]
			# linear <- ('linear' %in% bestPreds)
			# quadratic <- ('quadratic' %in% bestPreds)
			# hinge <- ('hinge' %in% bestPreds)
			# product <- ('product' %in% bestPreds)
			# threshold <- ('threshold' %in% bestPreds)
			
			# regMultMin <- calibSummary$value[calibSummary$term == 'regMult_25PercQuant']
			# regMultMax <- calibSummary$value[calibSummary$term == 'regMult_75PercQuant']
			# regMult <- seq(regMultMin, regMultMax, length.out=10)

			# say('Features included:', pre=1)
			# say('linear ........ ', linear)
			# say('quadratic ..... ', quadratic)
			# say('hinge ..........', hinge)
			# say('product ........', product)
			# say('threshold ......', threshold)
			# say('regMult ........', paste(round(regMult, 3), collapse=', '), pre=1)

			# model <- trainMaxEnt(data=data, resp='presBg', preds=vars, linear=linear, quadratic=quadratic, hinge=hinge, product=product, threshold=threshold, regMult=regMult, out=c('model', 'tuning'), jackknife=FALSE, scratchDir='C:/ecology/!Scratch/_temp')
			# model <- model$model
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), '_extent_with_maxent_', gcm, '_gcm.rda'))
			
		# } # next GCM
			
	# } # next extent
	
	# sink()
	
# say('############################')
# say('### project final models ###')
# say('############################')

	# # mask from Great Lakes and Canadian large lakes
	# load('D:/Ecology/Water Bodies/Great Lakes - USGS/greatLakes.rda')
	# load('D:/Ecology/Water Bodies/Major Lakes of Canada - USGS/canadaMajorLakes.rda')
	# canadaMajorLakes <- sp::spTransform(canadaMajorLakes, CRS(projection(greatLakes)))

	# !!!
	# mask <- raster('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/!ENSEMBLE 1950-2005 across 12 ESMs/an_avg_ETR.tif')
	# mask1 <- rasterize(greatLakes, mask)
	# mask2 <- rasterize(canadaMajorLakes, mask)
	
	# mask <- stack(mask1, mask2)
	# mask <- max(mask, na.rm=TRUE)
	# mask <- calc(mask, fun=function(x) ifelse(is.na(x), 1, NA))
	
	# for (ext in exts) {
# # for (ext in 'Broad') {

		# for (gcm in gcms) {
# # for (gcm in 'ccsm') {

			# for (algo in c('glm', 'maxent')) {
			
				# say(tolower(ext), ' ', gcm, ' ', algo)
				# load(paste0('./models/final_model_for_', tolower(ext), '_extent_with_', algo, '_', gcm, '_gcm.rda'))
							
				# if (exists('timeStack')) rm(timeStack)
				
				# for (year in rastTimes) {
					
					# clim <- getClimRasts(gcm=gcm, year=year, variables=predictors, rescale=TRUE)
					
					# pred <- if (algo == 'glm') {
						# predict(clim, model, type='response')
					# } else if (algo == 'maxent') {
						# raster::predict(clim, model, fun=enmSdm::predictMaxEnt, type='cloglog')
					# }
					
					# pred <- pred * mask
					# names(pred) <- paste0('year', year, 'BP')
					
					# timeStack <- if (exists('timeStack')) {
						# raster::stack(timeStack, pred)
					# } else {
						# pred
					# }
					
				# }
				
				# say('   interpolating...')
				# predByGen <- enmSdm::interpolateRasters(timeStack, interpFrom=rastTimes, interpTo=genTimes)
				
				# predByGen <- round(100 * predByGen)
				# names(predByGen) <- paste0('year', genTimes, 'BP')
				
				# writeRaster(predByGen, paste0('./predictions/predictions_', gcm, '_', algo, '_for_', tolower(ext), '_extent'), datatype='INT1U')
				
			# } # next algorithm
			
		# } # next GCM
			
	# } # next extent

# say('#########################################################')
# say('### convert prediction rasters to matrices for TIMBER ###')
# say('#########################################################')
	
	# # load extent for TIMBER... manually selected based on current occurrences and ENM and GLM ENM
	# load('./regions/timber_extent_manually_selected.rda')

	# # for (ext in exts) {
# for (ext in 'Broad') {

		# # for (gcm in gcms) {
# for (gcm in 'ccsm') {

			# for (algo in c('glm', 'maxent')) {

				# say(tolower(ext), ' ', gcm, ' ', algo)
		
				# # predictions
				# pred <- brick(paste0('./predictions/predictions_', gcm, '_', algo, '_for_', tolower(ext), '_extent.tif'), datatype='INT1U')
				
				# # crop to TIMBER extent
				# pred <- aggregate(pred, fact=2)
				# pred <- crop(pred, timberExtent)
				
				# pred <- projectRaster(pred, crs=getCRS('albersNA'))
				# pred <- round(pred)
				# names(pred) <- paste0('year', genTimes, 'BP')

				# predArray <- raster::as.array(pred)
				# saveRDS(predArray, paste0('./predictions/prediction_array_for_', algo, '_', gcm, '_gcm_', tolower(ext), '_extent.rds'))
				
			# } # next algorithm
					
		# } # next GCM
			
	# } # next extent

# say('#############################')
# say('### figure of predictions ###')
# say('#############################')
	
	# # americas <- c('Canada', 'Mexico', 'United States', 'Guatemala', 'Belize', 'Honduras', 'El Salvador', 'Nicaragua', 'Costa Rica', 'Panama', 'Colombia', 'Venezuela', 'Guyana', 'Ecuador', 'Brazil', 'Peru', 'Suriname', 'French Guiana', 'Cuba', 'Haiti', 'Dominican Republic', 'Jamaica', 'Chile', 'Paraguay', 'Uruguay', 'Trinidad and Tobago', 'Bahamas', 'Barbados', 'Saint Lucia', ' Saint Vincent and the Grenadines', 'Grenada', 'Antigua and Barbuda', 'Dominica', 'Saint Kitts and Nevis')
	# # americas <- c('Canada', 'Mexico', 'United States', 'Guatemala', 'Belize', 'Honduras', 'El Salvador', 'Nicaragua', 'Costa Rica', 'Panama', 'Colombia', 'Ecuador', 'Cuba', 'Haiti', 'Dominican Republic', 'Jamaica', 'Trinidad and Tobago', 'Bahamas', 'Barbados', 'Saint Lucia', ' Saint Vincent and the Grenadines', 'Grenada', 'Antigua and Barbuda', 'Dominica', 'Saint Kitts and Nevis')
	
	# # load('D:/Ecology/Political Geography/GADM/ver3pt6/gadm36_0.RData')
	# # gadm0Am <- gadm[gadm$NAME_0 %in% americas, ]
	# # gadm0AmAlb <- sp::spTransform(gadm0Am, getCRS('albersNA', TRUE))
	
	# load('./species_records/01_Fraxinus_pennsylvanica_bien_cleaned_occurrences_thinned.rda')
	# figExt <- extent(occsSpAlb)
	# figExt <- as(figExt, 'SpatialPolygons')
	# projection(figExt) <- getCRS('albersNA')
	# figExt <- gBuffer(figExt, width=300000)

	# # colors
	# breaks <- seq(0, 100, by=10)
	# cols <- c('gray', 'gray', '#e5f5e0', '#c7e9c0', '#a1d99b', '#74c476', '#41ab5d', '#238b45', '#006d2c', '#00441b')
	
	# for (ext in exts) {
	# # for (ext in exts[1]) {
	# # for (ext in exts[2]) {

		# # for (gcm in gcms[1]) {
		# for (gcm in gcms) {
		# # for (gcm in gcms[2]) {

			# algo <- 'glm'
			# predGlm <- brick(paste0('./predictions/predictions_', gcm, '_', algo, '_for_', tolower(ext), '_extent.tif'), datatype='INT1U')
			# predGlm <- projectRaster(predGlm, crs=getCRS('albersNA'))
			# names(predGlm) <- paste0('year', genTimes, 'BP')
			
			# algo <- 'maxent'
			# predMx <- brick(paste0('./predictions/predictions_', gcm, '_', algo, '_for_', tolower(ext), '_extent.tif'), datatype='INT1U')
			# predMx <- projectRaster(predMx, crs=getCRS('albersNA'))
			# names(predMx) <- paste0('year', genTimes, 'BP')
			
			# dirCreate(paste0('./figures_and_tables/prediction_figures_for_', gcm, '_gcm_', tolower(ext), '_ext'))
			
			# for (period in seq_along(genTimes)) {
				
				# genTime_yr <- genTimes[period]
				# say(tolower(ext), ' ', gcm, ' ', genTime_yr)
				
				# png(paste0('./figures_and_tables/prediction_figures_for_', gcm, '_gcm_', tolower(ext), '_ext/', prefix(genTime_yr, 5), 'YPB.png'), width=1200, height=800)

					# par(mfrow=c(1, 2), oma=c(0, 0, 0, 0), mar=rep(0.25, 4), cex.main=2)

					# plot(figExt, border=NA, ann=FALSE, xpd=NA)
					# plot(predGlm[[period]], col=cols, breaks=breaks, legend=FALSE, add=TRUE)
					# # plot(gadm0AmAlb, border='gray30', add=TRUE)
					# # title(main=paste0('GLM with ', ext, ' Background and ', toupper(gcm), ' AOGCM\n', genTime_yr, ' YBP'), line=-7)
					# text(-2474842, -2936141, labels=paste0(genTime_yr, ' YBP\n', ext, ' Extent | ', toupper(gcm), ' GCM | GLM'), cex=3, pos=4)
				
					# plot(figExt, border=NA, ann=FALSE, xpd=NA)
					# plot(predMx[[period]], col=cols, breaks=breaks, legend=FALSE, add=TRUE)
					# # plot(gadm0AmAlb, border='gray30', add=TRUE)
					# # title(main=paste0('Maxent with ', ext, ' Background and ', toupper(gcm), ' AOGCM\n', genTime_yr, ' YBP'), line=-7)
					# text(-2474842, -2936141, labels=paste0(genTime_yr, ' YBP\n', ext, ' Extent | ', toupper(gcm), ' GCM | MAXENT'), cex=3, pos=4)
					
					# # title(sub=date(), outer=TRUE, cex.sub=0.8, line=-2)
					
				# dev.off()
				
			# }
					
		# } # next GCM
			
	# } # next extent
	
# say('####################################################')
# say('### analyze similarity between model projections ###')
# say('####################################################')

	# say('To reduce the number of ABC iterations required, we will try to reduce the number of ENM model outputs that are included in the analysis. To do this, we want to cluster ENM output and select a subset of models based on similarity.')
	
	# say('Metrics to consider for clustering:')
	# say('Similarity between present-day predictions.')
	# say('Similarity between LGM predictions.')

	# algos <- c('glm', 'maxent')
	
	# ### get prediction rasters for present and LGM

	# # raster stacks for holding predictions
	# if (exists('sq')) rm(sq)
	# if (exists('lgm')) rm(lgm)
	
	# for (ext in exts) {
		# for (gcm in gcms) {
			# for (algo in algos) {

				# idString <- paste0(gcm, '_', algo, '_for_', tolower(ext))
				# thisPreds <- stack(paste0('./predictions/predictions_', idString, '_extent.tif'))

				# thisSq <- thisPreds[[1]]
				# thisLgm <- thisPreds[[nlayers(thisPreds)]]
				
				# maxPred <- max(c(cellStats(thisSq, 'max'), cellStats(thisLgm, 'max')))
				
				# thisSq <- thisSq / maxPred
				# thisLgm <- thisLgm / maxPred
				
				# names(thisSq) <- names(thisLgm) <- idString

				# sq <- if (exists('sq')) {
					# stack(sq, thisSq)
				# } else {
					# thisSq
				# }
				
				# lgm <- if (exists('lgm')) {
					# stack(lgm, thisLgm)
				# } else {
					# thisLgm
				# }

			# } # algo
		# } # GCM
	# } # extent
	
	# ### calculate similarity between predictions
	
	# sqArea <- area(sq, na.rm=TRUE, weights=TRUE)
	# lgmArea <- area(lgm, na.rm=TRUE, weights=TRUE)
				
	# sqArea <- values(sqArea)
	# lgmArea <- values(lgmArea)
				
	# sqAreaWeights <- sum(sqArea, na.rm=TRUE)
	# lgmAreaWeights <- sum(lgmArea, na.rm=TRUE)

	# # status quo
	# sqDiffs <- matrix(NA, nlayers(sq), nlayers(sq))
	# colnames(sqDiffs) <- rownames(sqDiffs) <- names(sq)
	
	# for (i in 1L:nlayers(sq)) {
		# for (j in 1L:nlayers(sq)) {
		
			# diff <- (sq[[i]] - sq[[j]])
			# diff <- abs(diff)
			# diff <- values(diff)
			# diff <- sum(diff * sqArea, na.rm=TRUE) / sqAreaWeights
			
			# sqDiffs[i, j] <- diff
		
		# }
	# }
	
	# # LGM
	# lgmDiffs <- matrix(NA, nlayers(lgm), nlayers(lgm))
	# colnames(lgmDiffs) <- rownames(lgmDiffs) <- names(lgm)
	
	# for (i in 1L:nlayers(lgm)) {
		# for (j in 1L:nlayers(lgm)) {
		
			# diff <- (lgm[[i]] - lgm[[j]])
			# diff <- abs(diff)
			# diff <- values(diff)
			# diff <- sum(diff * lgmArea, na.rm=TRUE) / lgmAreaWeights
			
			# lgmDiffs[i, j] <- diff
		
		# }
	# }
	
	# sqDiffs <- as.dist(sqDiffs)
	# lgmDiffs <- as.dist(lgmDiffs)

	# sqClust <- hclust(sqDiffs)
	# lgmClust <- hclust(lgmDiffs)
	
	# png('./figures_and_tables/clustering_of_prediction_rasters.png', width=1600, height=800)
		# par(mfrow=c(1, 2), cex=1.6)
		# plot(sqClust, main='Present-day differences')
		# plot(lgmClust, main='LGM differences')
	# dev.off()
	
	# png('./figures_and_tables/prediction_rasters_present.png', width=1600, height=800)
		# par(cex.main=3, mfrow=c(2, 4))
		# for (i in 1:nlayers(sq)) plot(sq[[i]], main=names(sq)[i])
	# dev.off()
	
	# png('./figures_and_tables/prediction_rasters_lgm.png', width=1600, height=800)
		# par(cex.main=3, mfrow=c(2, 4))
		# for (i in 1:nlayers(lgm)) plot(lgm[[i]], main=names(lgm)[i])
	# dev.off()
	
# say('#################################')
# say('### calculate biotic velocity ###')
# say('#################################')

	# algos <- c('glm', 'maxent')
	
	# # # NB times are 0 at LGM, 21000 at present
	# times <- c(0, 21000)
	# atTimes <- c(0, 21000)
	
	# # times <- seq(0, 21000, by=30)
	# # atTimes <- seq(0, 21000, by=30)
	
	# metrics <- c('centroid', 'nsCentroid', 'prevalence')

	# velocities <- data.frame()
	
	# for (ext in exts) {
		# for (gcm in gcms) {
			# for (algo in algos) {
			
				# say(paste(ext, gcm, algo))

				# # get predictions
				# idString <- paste0(gcm, '_', algo, '_for_', tolower(ext))
				# preds <- stack(paste0('./predictions/predictions_', idString, '_extent.tif'))
				# preds <- subset(preds, c(1, nlayers(preds)))
				
				# # project to equal-area
				# preds <- projectRaster(preds, crs=getCRS('albersNA'))
				# preds <- round(preds)

				# # biotic velocity
				# thisVelocity <- bioticVelocity(preds, times=times, atTimes=atTimes, metrics=metrics)
				
				# # remember
				# velocities <- rbind(
					# velocities,
					# cbind(
						# data.frame(
							# ext = tolower(ext),
							# gcm = gcm,
							# algo = algo
						# ),
						# thisVelocity
					# )
				# )

			# }
		# }
	# }
	
	# write.csv(velocities, './figures_and_tables/biotic_velocities_0_to_21Kybp.csv', row.names=FALSE)
	
# say('############################################')
# say('### plot potential total population size ###')
# say('############################################')

	# (velocities <-  read.csv('./figures_and_tables/biotic_velocities_0_to_21Kybp_by_30yr.csv')

	# years <- seq(21000 - 30, 0, by=-30)
	
	# plot(1, 1, col='white', ylim=c(0, 1), xlim=c(21000, 0), ylab='Mean suitability', xlab='Ybp')
	# lines(years, rev(velocities$prevalence[velocities$ext == 'broad' & velocities$gcm == 'ccsm' & velocities$algo == 'glm']), col='black', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'narrow' & velocities$gcm == 'ccsm' & velocities$algo == 'glm']), col='red', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'broad' & velocities$gcm == 'ccsm' & velocities$algo == 'maxent']), col='purple', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'narrow' & velocities$gcm == 'ccsm' & velocities$algo == 'maxent']), col='orange', lwd=2)
	
	# lines(years, rev(velocities$prevalence[velocities$ext == 'broad' & velocities$gcm == 'ecbilt' & velocities$algo == 'glm']), col='black', lty='dotted', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'narrow' & velocities$gcm == 'ecbilt' & velocities$algo == 'glm']), col='red', lty='dotted', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'broad' & velocities$gcm == 'ecbilt' & velocities$algo == 'maxent']), col='purple', lty='dotted', lwd=2)
	# lines(years, rev(velocities$prevalence[velocities$ext == 'narrow' & velocities$gcm == 'ecbilt' & velocities$algo == 'maxent']), col='orange', lty='dotted', lwd=2)
	
	# legend('bottomright',
		# legend=c(
			# 'CCSM GLM Broad',
			# 'CCSM GLM Narrow',
			# 'CCSM Maxent Broad',
			# 'CCSM Maxent Narrow',
			# 'ECBilt GLM Broad',
			# 'ECBilt GLM Narrow',
			# 'ECBilt Maxent Broad',
			# 'ECBilt Maxent Narrow'
		# ),
		# col=c(
			# 'black',
			# 'red',
			# 'purple',
			# 'orange'
		# ),
		# lty=c(
			# 'solid',
			# 'solid',
			# 'solid',
			# 'solid',
			# 'dotted',
			# 'dotted',
			# 'dotted',
			# 'dotted'
		# ),
		# pch=NA,
		# bty='n',
		# lwd=3
	# )	
	
#############################################
say('DONE', deco='~', pre=2, post=2, level=1)
