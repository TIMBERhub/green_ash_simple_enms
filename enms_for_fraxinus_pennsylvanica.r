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
### calibrate models ###
### calibrate and evaluate final models ###
### plot model performance against present-day data ###
### assess differences between model output ###
### project models back in time ###
### make maps of unthresholded predictions ###
### make maps of thresholded predictions ###
### determine thresholds that best recreate Little range ###
### calculate biotic velocity ###
### plot biotic velocity for periods < 21 Ka ###
### plot biotic velocity for periods of 21 Ka ###

### NOTES
### Variable naming conventions:
### * objects that are spatial (ie polygons/points) will be suffixed with "Sp" as in "speciesRecordsSp"
### * spatial objects that are in an equal-area projection will be further suffixed with "Alb" as in "speciesRecordsSpAlb"

#############
### setup ###
#############

	memory.limit(memory.limit() * 2^30)
	rm(list=ls()) # reproducibility!
	options(keep.source=FALSE) # manage memory
	gc()
	print('')
	print(date())

	# ### setup to run on Adam's computer
		
		# ### source('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/code/enms_for_fraxinus_pennsylvanica.r')

		# setwd('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms')
		# lorenzPath <- 'D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/'
		# studyRegionRastsFileName <- 'C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
		# elevationRastFileName <- 'D:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_elevationInMeters_fromEtop.tif'
		# tempDir <- 'D:/ecology/!Scratch/_temp'
	
	### setup to run on SHIGOTO computer
		
		### source('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/code/enms_for_fraxinus_pennsylvanica.r')

		setwd('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms')
		lorenzPath <- 'E:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/'
		studyRegionRastsFileName <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
		demoGeneticRasterTemplate <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_resampled_to_genetic_demographic_simulation_resolution.tif'
		elevationRastFileName <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_elevationInMeters_fromEtop.tif'
		tempDir <- 'E:/ecology/!Scratch/_temp'
	
	# ### setup to run on POWERBANK
	
		# ### source('H:/Global Change Program/Research/ABC for Biogeographic History of Trees/code/enms_for_fraxinus_pennsylvanica.r')
		
		# setwd('H:/Global Change Program/Research/ABC for Biogeographic History of Trees')
		# lorenzPath <- './!lorenz_et_al/Version 2017-06-16/'
		# studyRegionRastsFileName <- './!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif' # for POWERBANK computers
		# elevationRastFileName <- 'H:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_elevationInMeters_fromEtop.tif'
		# tempDir <- 'E:/ecology/!Scratch/_temp'

	options(stringsAsFactors=FALSE)
	raster::rasterOptions(format='GTiff', overwrite=TRUE)

	library(BIEN)
	library(brglm2)
	library(cluster)
	library(dismo)
	library(geosphere)
	library(rgdal)
	library(raster)
	library(RColorBrewer)
	library(rgeos)
	library(rJava)
	library(scales)
	library(sp)

	library(omnibus) # Adam's grab-bag library (https://github.com/adamlilith/omnibus)
	library(enmSdm) # Adam's SDM library (https://github.com/adamlilith/enmSdm)
	library(statisfactory) # Adam's statistics library (https://github.com/adamlilith/statisfactory)
	library(legendary) # Adam's plotting library (https://github.com/adamlilith/legenday)
	
	source('./code/assign_refugia_from_abundance_raster.r')
	
	dirCreate(tempDir)
	
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
		
			rast <- stack(paste0(lorenzPath, '/', gcmFolder, '/', year, 'BP/', variable, '.tif'))
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

	set.seed(pi)

	# species name
	species <- 'Fraxinus pennsylvanica'

	# buffer around occurrences ("accessible" area) used to define calibration area
	exts <- c(80, 160, 320) # in km

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
	
	# SDM algorithms
	algos <- c('brt', 'glm', 'maxent', 'ns')
	
	# desired number of final model clusters... used to simplify output
	numModelClusts <- 8

	# colors for clusters... only the first "numModelClusts" will be used
	clustCols <- c('#a6cee3','#1f78b4','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f')

	# decided in "### process current/paleoclimate layers ###"
	if (file.exists('./figures_and_tables/predictors.rda')) load('./figures_and_tables/predictors.rda')

	# study region polygon
	if (file.exists('./regions/studyRegion.rda')) load('./regions/studyRegion.rda')
	
	# mask raster for present-day land (1 = land, NA = not)
	if (file.exists('./regions/maskRaster.tif')) {
	
		maskRast <- raster('./regions/maskRaster.tif')
		
	} else {
	
		# mask raster for entire projection region
		maskRast <- raster('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/ccsm3_22-0k_all_tifs/0BP/an_avg_ETR.tif')
		maskRast <- 0 * maskRast + 1
		names(maskRast) <- 'mask'
		
		writeRaster(maskRast, './regions/maskRaster', datatype='INT2U')
		
	}

	maskRastAlb <- projectRaster(maskRast, crs=getCRS('albersNA'))
	
	# states that are undersampled for Fraxinus pennsylvanica species in BIEN 4.1
	undersampled <- c('Louisiana', 'Kentucky')

	# calibration regions
	if (file.exists('./regions/calibration_regions.rda')) load('./regions/calibration_regions.rda')
	if (file.exists('./regions/studyRegion_croppedToPresentLand.rda')) load('./regions/studyRegion_croppedToPresentLand.rda')
	
	# North America spatial polygons
	if (file.exists('./regions/greatLakesSpAlb.rda')) load('./regions/greatLakesSpAlb.rda')

	if (file.exists('./regions/gadm36_northAmerica_sp_level_0.rda')) load('./regions/gadm36_northAmerica_sp_level_0.rda')
	if (file.exists('./regions/gadm36_northAmerica_spAlb_level_0.rda')) load('./regions/gadm36_northAmerica_spAlb_level_0.rda')
	
	if (file.exists('./regions/gadm36_northAmerica_sp_level_1.rda')) load('./regions/gadm36_northAmerica_sp_level_1.rda')
	if (file.exists('./regions/gadm36_northAmerica_spAlb_level_1.rda')) load('./regions/gadm36_northAmerica_spAlb_level_1.rda')
	
	if (file.exists('./regions/gadm36_northAmerica_spAlb_continental.rda')) load('./regions/gadm36_northAmerica_spAlb_continental.rda')
	
	if (file.exists('./regions/gadm36_northAmerica_spAlb_continental_cropToStudyRegion.rda')) load('./regions/gadm36_northAmerica_spAlb_continental_cropToStudyRegion.rda')
	
	# range maps
	if (file.exists('./regions/bien_range_map/Fraxinus_pennsylvanica.shp')) {
		bienRange <- shapefile('./regions/bien_range_map/Fraxinus_pennsylvanica.shp')
		bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA', TRUE))
	}
	
	if (file.exists('./regions/little_range_map/littleRangeSpAlb.rda')) load('./regions/little_range_map/littleRangeSpAlb.rda')
	
# say('#######################################', pre=1)
# say('### obtain and process basemap data ###')
# say('#######################################', post=2)

	# say('Basemap data includes shapefiles of North American countries and the range maps of the species from BIEN and Little. We will create equal-area versions of each in an Albers projection (object names will be suffixed with "Alb")', breaks=80)

	# # study region
	# if (!file.exists('./regions/studyRegion.rda')) {
	
		# studyRegionSpAlb <- shapefile('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/study_region_spatial_polygons/study_region_mask_without_glaciers')
		# save(studyRegionSpAlb, file='./regions/studyRegion.rda')
		
	# }

	# ### North America Level 2
		
		# mex <- getData('GADM', country='MEX', level=2, path='D:/ecology/!Scratch')
		# usa <- getData('GADM', country='USA', level=2, path='D:/ecology/!Scratch')
		# can <- getData('GADM', country='CAN', level=2, path='D:/ecology/!Scratch')

		# nam2Sp <- rbind(can, usa, mex)

		# # get lakes for removal from other geographies... add small buffer because lake borders don't exactly align
		# greatLakesSp2 <- nam2Sp[nam2Sp$ENGTYPE_2 == 'Water body', ]
		# greatLakesSpAlb2 <- sp::spTransform(greatLakesSp2, getCRS('albersNA', TRUE))
		# greatLakesSpAlb2 <- gBuffer(greatLakesSpAlb2, width=10)
		# greatLakesSp2 <- sp::spTransform(greatLakesSpAlb2, getCRS('wgs84', TRUE))
		# greatLakesSp <- gUnaryUnion(greatLakesSp2)
		
		# greatLakesSpAlb <- sp::spTransform(greatLakesSp2, getCRS('albersNA', TRUE))
		
		# # names of level 1 areas with lakes
		# lakeNames0 <- unique(nam2Sp@data$NAME_0[nam2Sp@data$ENGTYPE_2 == 'Water body'])
		# lakeNames1 <- unique(nam2Sp@data$NAME_1[nam2Sp@data$ENGTYPE_2 == 'Water body'])
		
		# # remove lakes
		# nam2Sp <- nam2Sp[nam2Sp@data$ENGTYPE_2 != 'Water body', ]
		
		# # not saving because only need for identifying lake areas
		
		# save(greatLakesSpAlb, file='./regions/greatLakesSpAlb.rda')

	# ### North America Level 1
	
		# mex <- getData('GADM', country='MEX', level=1, path='D:/ecology/!Scratch')
		# usa <- getData('GADM', country='USA', level=1, path='D:/ecology/!Scratch')
		# can <- getData('GADM', country='CAN', level=1, path='D:/ecology/!Scratch')

		# nam1Sp <- rbind(can, usa, mex)
		
		# for (thisLevel in lakeNames1) {
		
			# thisLargerSp <- nam1Sp[nam1Sp@data$NAME_1 == thisLevel, ]
			# df <- thisLargerSp@data
			# nam1Sp <- nam1Sp[-which(nam1Sp@data$NAME_1 == thisLevel), ]
			# thisLargerSansLakesSp <- gDifference(thisLargerSp, greatLakesSp)
			# thisLargerSansLakesSp <- as(thisLargerSansLakesSp, 'SpatialPolygonsDataFrame')
			# projection(thisLargerSansLakesSp) <- projection(nam1Sp)
			# thisLargerSansLakesSp@data <- df
			# nam1Sp <- rbind(nam1Sp, thisLargerSansLakesSp)
			
		# }

		# nam1SpAlb <- sp::spTransform(nam1Sp, getCRS('albersNA', TRUE))
		# save(nam1Sp, file='./regions/gadm36_northAmerica_sp_level_1.rda')
		# save(nam1SpAlb, file='./regions/gadm36_northAmerica_spAlb_level_1.rda')
	
	# ### North America level 0

		# mex <- getData('GADM', country='MEX', level=0, path='D:/ecology/!Scratch')
		# usa <- getData('GADM', country='USA', level=0, path='D:/ecology/!Scratch')
		# can <- getData('GADM', country='CAN', level=0, path='D:/ecology/!Scratch')

		# nam0Sp <- rbind(can, usa, mex)

		# for (thisLevel in lakeNames0) {
		
			# thisLargerSp <- nam0Sp[nam0Sp@data$NAME_0 == thisLevel, ]
			# df <- thisLargerSp@data
			# nam0Sp <- nam0Sp[-which(nam0Sp@data$NAME_0 == thisLevel), ]
			# thisLargerSansLakesSp <- gDifference(thisLargerSp, greatLakesSp)
			# thisLargerSansLakesSp <- as(thisLargerSansLakesSp, 'SpatialPolygonsDataFrame')
			# projection(thisLargerSansLakesSp) <- projection(nam1Sp)
			# thisLargerSansLakesSp@data <- df
			# nam0Sp <- rbind(nam0Sp, thisLargerSansLakesSp)
			
		# }

		# nam0SpAlb <- sp::spTransform(nam0Sp, getCRS('albersNA', TRUE))
		# save(nam0Sp, file='./regions/gadm36_northAmerica_sp_level_0.rda')
		# save(nam0SpAlb, file='./regions/gadm36_northAmerica_spAlb_level_0.rda')
		
		# # continent
		# namSpAlb <- gUnaryUnion(nam0SpAlb)
		# save(namSpAlb, file='./regions/gadm36_northAmerica_spAlb_continental.rda')
		
		# # # continent cropped to study region
		# studyRegionRasts <- brick('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif')
					
		# namSpAlbStudyRegion <- crop(namSpAlb, studyRegionRasts[[1]])
		# save(namSpAlbStudyRegion, file='./regions/gadm36_northAmerica_spAlb_continental_cropToStudyRegion.rda')
			
	# ### study region cropped to present-day land
	# if (!exists('studyRegionCroppedToPresentSpAlb')) {
		# studyRegionCroppedToPresentSpAlb <- crop(studyRegionSpAlb, nam0SpAlb)
		# save(studyRegionCroppedToPresentSpAlb, file='./regions/studyRegion_croppedToPresentLand.rda')
	# }
		

	# ### range maps
	# dirCreate('./regions/bien_range_map')
	# BIEN_ranges_species(species, directory='./regions/bien_range_map')
	# bienRange <- BIEN_ranges_load_species(species)
	# bienRangeAlb <- sp::spTransform(bienRange, getCRS('albersNA', TRUE))

	# littleRangeFileName <- 'C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/range_maps/little/fraxpenn.shp'
	# littleRange <- shapefile(littleRangeFileName)
	# projection(littleRange) <- getCRS('wgs84')
	# littleRange <- littleRange[littleRange$CODE == 1, ] # remove holes
	# littleRangeSpAlb <- sp::spTransform(littleRange, getCRS('albersNA', TRUE))
	# dirCreate('./regions/little_range_map')
	# save(littleRangeSpAlb, file='./regions/little_range_map/littleRangeSpAlb.rda')
			
# say('####################################', pre=1)
# say('### obtain and clean occurrences ###')
# say('####################################', post=2)

	# ### obtain records
	# ##################
		
		# speciesFileName <- paste0(
			# './species_records/00_',
			# gsub(tolower(species), pattern=' ', replacement='_'),
			# '_bien_all_occurrences.rda'
		# )

		# if (file.exists(speciesFileName)) {
			# load(speciesFileName)
		# } else {
			# occsRaw <- BIEN_occurrence_species(
				# species=species,
				# cultivated = FALSE,
				# only.new.world = TRUE,
				# all.taxonomy = FALSE,
				# native.status = FALSE,
				# natives.only = TRUE,
				# observation.type = TRUE,
				# political.boundaries = TRUE,
				# collection.info = TRUE
			# )

			# sink('./species_records/bien_download.txt')
				# say('Data representing species records were downloaded from BIEN on', date(), '.')
			# sink()
			# save(occsRaw, file=speciesFileName)
		# }

		# ### remove occurrences with missing coordinates and dates:
		# occs <- occsRaw[!(is.na(occsRaw$longitude) | is.na(occsRaw$latitude) | is.na(occsRaw$date_collected)), ]
		# dim(occs)

		# ### remove records <1950 (period covered by Lorenz et al. 2016 climate data is 1950-2005):
		# occs$date_collected <- as.Date(occs$date_collected)
		# occs <- occs[!is.na(occs$date_collected), ]
		# occs <- occs[which(occs$date_collected >= as.Date(paste0(1950 + genTime_yr, '-01-01'))), ]

	# ### remove records in states where green ash has been naturalized (according to BONAP)
	# ######################################################################################
	
		# occs <- occs[!(occs$state_province %in% c('Washington', 'Oregon', 'Idaho', 'Utah', 'Colorado', 'Arizona', 'New Mexico', 'Chihuahua', 'British Columbia')), ]
		# occs <- occs[-which(occs$state_province == 'Texas' & occs$county == 'El Paso'), ]
		
	# ### inspect occurrences
	# #######################
		
		# ### convert occurrences to spatial object
		# occsSp <- SpatialPointsDataFrame(
			# occs[ , ll],
			# data=occs,
			# proj4=getCRS('wgs84', TRUE)
		# )

		# occsSpAlb <- sp::spTransform(occsSp, getCRS('albersNA', TRUE))

		# png('./figures_and_tables/occurrences_raw.png', width=1200, height=1000)
		
			# par(mfrow=c(1, 2))
			
			# # plot occurrences with BIEN range map
			# plot(occsSpAlb, col=NA, main='BIEN range map')
			# plot(bienRangeAlb, col=alpha('green', 0.4), border='green', add=TRUE)
			# points(occsSpAlb, pch=16, cex=0.9, col='black')
			# plot(nam1SpAlb, add=TRUE, border='gray')

			# # plot occurrences with Little's range map
			# plot(occsSpAlb, col=NA, main='Little range map')
			# plot(littleRangeSpAlb, col=alpha('red', 0.4), border='red', add=TRUE)
			# points(occsSpAlb, pch=16, cex=0.9, col='black')
			# plot(nam1SpAlb, add=TRUE, border='gray')

			# title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)
			
		# dev.off()

		# say('Things to note:')
		# say('* There is at least one occurrence in South America or its vicinity.')
		# say('* There are many occurrences outside Little\'s range polygon.')
		# say('* ...Especially in the northern part of peninsular Florida', post=2)
		# say('* ...and the western US', post=2)
		
		# ### examine areas in vicinity of Kentucky and Louisiana that appear to lack data

		# png('./figures_and_tables/occurrences_near_louisiana_kentucky.png', width=1200, height=1000)
		
			# par(mfrow=c(1, 2))

			# plot(nam1SpAlb[nam1SpAlb@data$NAME_1 == 'Kentucky', ])
			# points(occsSpAlb, pch=16, cex=0.9, col='black')

			# plot(nam1SpAlb[nam1SpAlb@data$NAME_1 == 'Louisiana', ])
			# points(occsSpAlb, pch=16, cex=0.9, col='black')
			
			# title(main=paste0('Missing records in ', paste(undersampled, collapse=' and ')), sub=date(), outer=TRUE, cex.sub=0.7, line=-2)

		# dev.off()
			
		# say('So, yes, there is a bias in the central and Appalachian region of Kentucky... Note that the species genuinely does not appear in the "high" Appalachians. This is disturbing since the usual methods of correcting for bias like inverse-p weighting or spatial thinning would remove signal from its absence in these areas *plus* in areas that are undersampled like Kentucky and Louisiana.', post=2, breaks=80)

		# say('So, we are going to remove Louisiana and Kentucky from the training set. This means occurrences from these areas will be be use to train the models, nor will background sites drawn from these areas.', post=2, breaks=80)

	# ### remove occurrences in undersampled areas and in South America
	# #################################################################
	
		# nrowStart <- nrow(occsSpAlb)
		# occsSpAlb <- occsSpAlb[(occsSpAlb@data$country %in% c('Canada', 'United States', 'Mexico')) & !(occsSpAlb@data$state_province %in% undersampled), ]
		
		# sink('./species_records/cleaning_notes_outliers_undersampled.txt', split=TRUE)
			# say('By removing occurrences outside Mexico/US/Canada we reduced the number of occurrences from ', nrowStart, ' to ', nrow(occsSpAlb), '.')
		# sink()

	# ### correct for hypothesized spatial sampling bias
	# ##################################################
		
		# say('To remove sampling bias we will spatially thin occurrences so that there is but one occurrence per cell.')

		# # get occurrences in WGS84 CRS
		# priority <- rep(NA, nrow(occsSpAlb))
		# priority[occsSpAlb$observation_type == 'plot'] <- 1
		# priority[occsSpAlb$observation_type == 'specimen'] <- 2
		# priority[occsSpAlb$observation_type == 'trait occurrence'] <- 3
		# priority[occsSpAlb$observation_type == 'human observation'] <- 4
		# priority[occsSpAlb$observation_type == 'unknown'] <- 5
		
		# nrowStart <- nrow(occsSpAlb)
		# occsSpAlb <- elimCellDups(occsSpAlb, r=maskRastAlb, priority=priority)

		# sink('./species_records/cleaning_notes_thinning.txt', split=TRUE)
			# say('Through spatial thinning we reduced the number of occurrences from ', nrowStart, ' to ', nrow(occsSpAlb), '.')
		# sink()
		
		# # inspect
		# png('./figures_and_tables/occurrences_thinned_one_per_cell.png', width=1200, height=1000)
		
			# plot(occsSpAlb, pch=16, cex=1.8, main='All cleaned and thinnned occurrences')
			# plot(nam1SpAlb, add=TRUE, border='gray')
			
			# title(sub=date(), outer=TRUE, cex.sub=0.7, line=-2)
			
		# dev.off()
		
		# # save
		# save(occsSpAlb, file=paste0('./species_records/01_', (gsub(tolower(species), pattern=' ', replacement='_')), '_bien_cleaned_occurrences_thinned.rda'))
	
# say('##################################', pre=1)
# say('### define calibration regions ###')
# say('##################################', post=2)	
	
	# say('Currently, best practices assert that one should model the niche using only data from an area considered accessible to the species over long time. However, in practice this is very difficult to estimate, so we will use buffers of different size (defined above in the "### constants ###" section above).', post=2, breaks=80)

	# ### load data
	
	# load(paste0('./species_records/01_', gsub(tolower(species), pattern=' ', replacement='_'), '_bien_cleaned_occurrences_thinned.rda'))
	
	# ### create buffers
	
	# if (!exists('calibRegionsSpAlb')) {
		
		# # buffer around occurrences
		# calibRegionsSpAlb <- list()
		# for (i in seq_along(exts)) {
			# calibRegionsSpAlb[[i]] <- gBuffer(occsSpAlb, width=exts[i] * 1000)
		# }

		# # crop area outside North America, remove undersampled states, remove Great Lakes
		# undersampledSpAlb <- nam1SpAlb[nam1SpAlb$NAME_1 %in% undersampled, ]
		
		# for (i in seq_along(exts)) {
			# # calibRegionsSpAlb[[i]] <- gIntersection(calibRegionsSpAlb[[i]], sr) # not doing this because some points outside study region
			# calibRegionsSpAlb[[i]] <- rgeos::gIntersection(calibRegionsSpAlb[[i]], nam0SpAlb)
			# calibRegionsSpAlb[[i]] <- rgeos::gDifference(calibRegionsSpAlb[[i]], undersampledSpAlb, checkValidity=TRUE)
			# calibRegionsSpAlb[[i]] <- rgeos::gDifference(calibRegionsSpAlb[[i]], greatLakesSpAlb, checkValidity=TRUE)
		# }
		
		# studyRegionCroppedToPresentSpAlb <- crop(studyRegionSpAlb, nam0SpAlb)
		
		# # inspect
		# png('./figures_and_tables/regions.png', width=1200, height=1000)

			# par(cex.main=1.4)
			# plot(calibRegionsSpAlb[[length(calibRegionsSpAlb)]], border=NA, col=NA, main='Calibration Regions')
			# plot(studyRegionCroppedToPresentSpAlb, add=TRUE, lwd=0.8, col='khaki1')
			
			# for (i in rev(seq_along(exts))) {
				# plot(calibRegionsSpAlb[[i]], col=paste0('aquamarine', 4 + 1 - i), lwd=0.8, add=TRUE)
			# }
			
			# plot(nam1SpAlb, add=TRUE)
			# points(occsSpAlb, pch=16, cex=1.4)

			# legend <- c('Records', paste('Calibration region:', rev(exts), 'km'), 'Study region (modern)')
			# fill <- c(NA, paste0('aquamarine', 4 + 1 - rev(seq_along(exts))), 'khaki1')
			# pch <- c(16, rep(NA, length(exts)), NA)
			# pt.cex <- c(1.4, rep(NA, length(exts)), NA)
			# border <- c(NA, rep(0.8, length(exts)), 0.8)
			# col <- c('black', rep(NA, length(exts)), NA)
			
			# legend('bottomright',
				   # legend=legend,
				   # fill=fill,
				   # pch=pch,
				   # pt.cex=pt.cex,
				   # border=border,
				   # col=col,
				   # bty='n', cex=1.6
			# )
			
		# dev.off()

		# names(calibRegionsSpAlb) <- paste0('bufferAroundOccs_', exts, 'km')
		# save(calibRegionsSpAlb, file='./regions/calibration_regions.rda')
		
	# }

# say('################################', pre=1)
# say('### select climate variables ###')
# say('################################', post=2)

	# say('We will be using the paleo climate layers prepared by Lorenz, D.J., Nieto-Lugilde, D., Blois, J.L., Fitzpatrick, M.C., and Williams, J.W.  2016.  Downscaled and debiased climate simulations for North America from 21,000 years ago to 2100 AD.  Scientific Data 3:160048. The DOI of the original data is 10.5061/dryad.1597g. This data set contains future, current, and historical GCM projections. We will be using the current and historical model output from the ECBilt and CCSM GCMs (other models have only current/future projections associated with them). The spatial resolution is 0.5deg for all of North America.', breaks=80, post=2)

	# say('Available variables include min/max temperature, precipitation, growing degree days, AET, PET, ETR (= AET/PET), and WDI (PET - precipitation), summarized across an annual, quarterly, or monthly basis. Variables are projected back to 22 Kybp (CCSM) or 21 Kybp (ECBilt) in 500-yr timeslices. The current period covers averages across 1950-2005. Past climate encompasses averages across 200-yr intervals centered on every 500 yr (e.g., 500 ypb, 1000 ypb, 1500 ypb, etc.).', breaks=80, post=2)
	
	# say('We will *not* use variable "an_cv_ETR" because it has NAs in Canada during the present in large portions of Canada (and possibly through time).', breaks=80, post=2)
	# say('We will *not* use variable "an_cv_WDI" because it has some very extreme values in just a few cells in Mexico at 0 ybp (and possibly through time).', breaks=80, post=2)
	
	# say('Tasks to complete include:')
	# say('* Selecting variables.')
	# say('* Cropping/masking rasters to study region *extent* (versus the areas within the study region).')
	# say('* Interpolating rasters to 30-yr intervals (approximate generation time for green ash).')

	# say('### select variables', level=2)
	# ####################################
		
		# calibRegionBroadestSp <- calibRegionsSpAlb[[paste0('bufferAroundOccs_', max(exts), 'km')]]
		
		# ### get set of randomly located points for assessing collinearity between variables
		
		# climCcsm_0ypb <- raster::stack(listFiles('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/ccsm3_22-0k_all_tifs/0BP', pattern='an_'))
		# climEcbilt_0ypb <- raster::stack(listFiles('D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/ecbilt_21-0k_all_tifs/0BP', pattern='an_'))
		
		# # remove an_cv_ETR and an_cv_WDI
		# keeps <- names(climCcsm_0ypb)[-which(names(climCcsm_0ypb) %in% c('an_cv_ETR', 'an_cv_WDI'))]
		# climCcsm_0ypb <- subset(climCcsm_0ypb, keeps)
		# climEcbilt_0ypb <- subset(climEcbilt_0ypb, keeps)
		
		# maskBroad0ypb <- rasterize(calibRegionBroadestSp, climCcsm_0ypb[[1]])
		
		# randBg <- randomPoints(maskBroad0ypb, 11000)
		# randBg <- as.data.frame(randBg)
		# names(randBg) <- ll
		
		# envCcsm_0ybp <- extract(climCcsm_0ypb, randBg)
		# envEcbilt_0ybp <- extract(climEcbilt_0ypb, randBg)

		# completeCases <- complete.cases(envCcsm_0ybp) & complete.cases(envEcbilt_0ybp)
		
		# randBg <- randBg[completeCases, ]
		# randBg <- randBg[1:max(nrow(randBg), 10000), ]
		# envCcsm_0ybp <- envCcsm_0ybp[completeCases, ]
		# envEcbilt_0ybp <- envEcbilt_0ybp[completeCases, ]
		
		# env <- rbind(envCcsm_0ybp, envEcbilt_0ybp)
		
		# env <- as.data.frame(env)
		# envScaled <- scale(env)
		# variableSel <- usdm::vifcor(envScaled, maxCor, maxobservations=Inf)

		# sink('./figures_and_tables/variable_selection.txt', split=TRUE)
		
			# say('Based on a threshold maximum correlation of ', maxCor, ' and VIF, we will use the following variables:')
			# print(variableSel)
			
		# sink()
		
		# predictors <- variableSel@results$Variables
		# save(predictors, file='./figures_and_tables/predictors.rda')
		
	# say('calculate statistics for rescaling rasters from present-day rasters', level=2)
	# ###################################################################################
		
		# # saves statistics on each variable (min/max/etc)
		# variableStats <- data.frame()
		
		# for (gcm in gcms) {
		
			# rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale = FALSE)
			# mins <- minValue(rasts)
			# maxs <- maxValue(rasts)
			
			# # remember
			# variableStats <- rbind(
				# variableStats,
				# data.frame(
					# gcm = gcm,
					# variable = predictors,
					# min = mins,
					# max = maxs
				# )
			# )
	
		# } # next GCM
				
		# write.csv(variableStats, './environmental_rasters/lorenz_et_al_2016/variable_statistics_0bp_across_north_america.csv', row.names=FALSE)

# say('###########################################', pre=1)
# say('### collate calibration/evaluation data ###')
# say('###########################################', post=2)

	# say('In this step we will:')
	# say('* Extract environmental data.')
	# say('* Locate background sites.')
	# say('* Define geographically distinct cross-validation folds for model evaluation.')

	# calibRegionsSp <- calibRegionsSpAlb
	# for (i in seq_along(calibRegionsSp)) {
		# calibRegionsSp[[i]] <- sp::spTransform(calibRegionsSp[[i]], getCRS('wgs84', TRUE))
	# }
	
	# load(paste0('./species_records/01_', gsub(tolower(species), pattern=' ', replacement='_'), '_bien_cleaned_occurrences_thinned.rda'))
	
		# occsSp <- sp::spTransform(occsSpAlb, getCRS('wgs84', TRUE))
	
	# say('extract environmental data', level=2)
	# ##########################################

		# # by GCM
		# for (gcm in gcms) {

			# if (exists('rasts', inherits=FALSE)) rm(rasts)
		
			# # get current version of each variable raster
			# rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
			# areas <- area(rasts[[1]])

			# # extract
			# env <- raster::extract(rasts, occsSp)
			# env <- as.data.frame(env)
			
			# thisArea <- raster::extract(areas, occsSp)
			# thisArea <- as.data.frame(thisArea)
			# names(thisArea) <- 'cellArea_km2'
				
			# thisEnv <- cbind(thisArea, env)
			
			# # remember
			# occsSpAlb@data <- cbind(occsSpAlb@data, thisEnv)
			
		# } # next GCM
		
		# # remove occurrences with incomplete environmental data
		# vars <- character()
		# for (gcm in gcms) vars <- c(vars, paste0(gcm, '_', predictors))
		# noNas <- complete.cases(occsSpAlb@data[ , vars])
		# occsSpAlb <- occsSpAlb[noNas, ]

	# say('generate test background sites', level=2)
	# ##############################################
	
		# calibRegionSpAlb <- calibRegionsSpAlb[[paste0('bufferAroundOccs_', max(exts), 'km')]]
		
		# # get random background sites
		# bgTestSpAlb <- spsample(calibRegionSpAlb, n=11000, type='random', iter=10)
		# bgTestSp <- sp::spTransform(bgTestSpAlb, getCRS('wgs84', TRUE))
		# bgTest <- as.data.frame(coordinates(bgTestSp))
		# names(bgTest) <- ll
		
		# if (exists('env')) rm(env)

		# # by GCM
		# for (gcm in gcms) {

			# # get current version of each variable raster
			# rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
			# areas <- area(rasts[[1]])

			# # extract
			# thisEnv <- raster::extract(rasts, bgTest)
			# thisEnv <- as.data.frame(thisEnv)
			
			# thisArea <- raster::extract(areas, bgTest)
			# thisArea <- as.data.frame(thisArea)
			# names(thisArea) <- 'cellArea_km2'
			
			# thisEnv <- cbind(thisArea, thisEnv)
			
			# # remember
			# env <- if (exists('env', inherits=FALSE)) {
				# cbind(env, thisEnv)
			# } else {
				# thisEnv
			# }
			
		# } # next GCM

		# # cull NA cases
		# bgTest <- cbind(bgTest, env)
		# bgTest <- bgTest[complete.cases(bgTest), ]
		# bgTest <- bgTest[1:min(nrow(bgTest), 10000), ]
			
	# say('generate training background sites', level=2)
	# ##################################################

		# for (countExt in seq_along(exts)) {
				
			# ext <- exts[countExt]
			# calibRegionSpAlb <- calibRegionsSpAlb[[countExt]]
			
			# # get random background sites
			# bgCalibSpAlb <- spsample(calibRegionSpAlb, n=11000, type='random', iter=10)
			# bgCalibSp <- sp::spTransform(bgCalibSpAlb, getCRS('wgs84', TRUE))
			# bgCalib <- as.data.frame(coordinates(bgCalibSp))
			# names(bgCalib) <- ll
			
			# if (exists('env')) rm(env)

			# # by GCM
			# for (gcm in gcms) {

				# # get current version of each variable raster
				# rasts <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=TRUE)
				# areas <- area(rasts[[1]])

				# # extract
				# thisEnv <- raster::extract(rasts, bgCalib)
				# thisEnv <- as.data.frame(thisEnv)
				
				# thisArea <- raster::extract(areas, bgCalib)
				# thisArea <- as.data.frame(thisArea)
				# names(thisArea) <- 'cellArea_km2'
				
				# thisEnv <- cbind(thisArea, thisEnv)
				
				# # remember
				# env <- if (exists('env', inherits=FALSE)) {
					# cbind(env, thisEnv)
				# } else {
					# thisEnv
				# }
				
			# } # next GCM

			# # cull NA cases
			# bgCalib <- cbind(bgCalib, env)
			# bgCalib <- bgCalib[complete.cases(bgCalib), ]
			# bgCalib <- bgCalib[1:min(nrow(bgCalib), 10000), ]
			
			# # remember
			# assign(paste0('bgCalib_extent', exts[countExt], 'km'), bgCalib)
			
		# } # next calibration region

	# say('assign spatially exclusive training/calibration/evaluation cross-validation folds', level=2)
	# #################################################################################################
		
		# bgTestSp <- SpatialPoints(bgTest[ , ll], getCRS('wgs84', TRUE))
		# bgTestSpAlb <- sp::spTransform(bgTestSp, getCRS('albersNA', TRUE))
		
		# blockRast <- rastWithSquareCells(studyRegionCroppedToPresentSpAlb, res=foldBuffer_m)

		# # create folds based on occurrences such that there is at least a minimum number of occurrences per fold
		# minNumInFold <- -Inf
		# numFolds <- 50
		# maxTries <- 20
		# tries <- 1
		
		# while (minNumInFold < minFoldSize) {
		
			# values(blockRast) <- sample(1:numFolds, ncell(blockRast), replace=TRUE)
			# folds <- extract(blockRast, occsSpAlb)
			# minNumInFold <- min(table(folds))
			
			# tries <- tries + 1
			# if (tries == maxTries) {
				# tries <- 1
				# numFolds <- numFolds - 1
			# }
			
		# }
	
		# say('Final number of folds: ', numFolds)
		
		# occs <- occsSpAlb@data
		# folds <- data.frame(fold=folds)
		# occs <- cbind(occs, folds)
		
		# # assign folds to test background sites
		# bgTestSp <- SpatialPointsDataFrame(bgTest[ , ll], data=bgTest, proj4string=getCRS('wgs84', TRUE))
		# bgTestSpAlb <- sp::spTransform(bgTestSp, getCRS('albersNA', TRUE))
		# bgTestFolds <- extract(blockRast, bgTestSpAlb)
		# bgTestSpAlb$fold <- bgTestFolds
		# bgTest <- bgTestSpAlb@data
		
		# # assign folds to calibration background sites
		# for (ext in exts) {
		
			# x <- get(paste0('bgCalib_extent', ext, 'km'))
			# xSp <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			# xSpAlb <- sp::spTransform(xSp, getCRS('albersNA', TRUE))
			# bgCalibFolds <- extract(blockRast, xSpAlb)
			# x$fold <- bgCalibFolds
			# x <- as.data.frame(x)
			
			# assign(paste0('bgCalib_extent', ext, 'km'), x)
			
		# }

	# say('collate calibration and evaluation	occurrences and background sites', level=2)
	# ##################################################################################
	
		# occsBg <- list()
		# occsBg$note <- 'This list object contains spatially thinned occurrence data and associated background sites. $testBg is a set of background sites drawn from the largest training region extent. $calibEvalOccsBg[[i]] is a list of lists, one per training region extent. Each sublist has two data frames. $calibEvalOccsBg[[i]]occsBg is a data frame representing calibration occurrences and background sites. Weights are 1 for occurrences and a fractional value for background sites such that their combined weight is equal to the total weight of occurrences. Folds are geographically exclusive and comprise a set of intermingled spatial blocks. $calibEvalOccsBg[[i]]$folds is a data frame with one row per row in $calibEvalOccsBg[[i]]$occsBg and one column per fold. Values in each column are 1 (training data), 2 (calibration data), or NA (masked... these represent evaluation data censored from the training and calibration sets).'
		
		# occsBg$meta <- data.frame(
			# item = c(
				# 'numOccs',
				# 'numFolds',
				# 'minFoldSize_numOccs'
			# ),
			# value = c(
				# nrow(occs),
				# numFolds,
				# minFoldSize
			# )
		# )
		
		# columnsJustBg <- c('fold', ll, 'cellArea_km2', paste0(rep(gcms, each=length(predictors)), '_', predictors))
		# columnsOccsAndBg <- c('presBg', 'fold', 'weight', ll, 'cellArea_km2', paste0(rep(gcms, each=length(predictors)), '_', predictors))
		
		# occsBg$allOccsThinned <- occs
		# occsBg$testBg <- bgTest[ , columnsJustBg]
		# occsBg$calibEvalOccsBg <- list()
		
		# for (countExents in seq_along(exts)) {
		
			# ext <- exts[countExents]
		
			# # collate occurrences and background sites
			# thisOccs <- occs
			# thisOccs$presBg <- 1
			# thisOccs$weight <- 1
			# thisOccs <- thisOccs[ , columnsOccsAndBg]
			
			# thisBg <- get(paste0('bgCalib_extent', ext, 'km'))
			# thisBg$presBg <- 0
			# thisBg$weight <- nrow(thisOccs) / nrow(thisBg)
			# thisBg <- thisBg[ , columnsOccsAndBg]
			
			# thisOccsBg <- rbind(thisOccs, thisBg)
			
			# # make matrix indicating which occurrences/background sites are in which fold set
			# if (exists('folds')) rm(folds)
			
			# for (calibFold in 1:numFolds) {
			
				# for (evalFold in (1:numFolds)[-calibFold]) {
			
					# designation <- rep(1, nrow(thisOccsBg)) # start with all training
					# designation[thisOccsBg$fold == calibFold] <- 2 # calibration folds
					# designation[thisOccsBg$fold == evalFold] <- NA # mask evaluation folds
					
					# designation <- matrix(designation, ncol=1)
					# colnames(designation) <- paste0('calib', calibFold, '_vs_eval', evalFold)
					
					# folds <- if (exists('folds')) {
						# cbind(folds, designation)
					# } else {
						# designation
					# }
			
				# }
			
			# }
			
			# occsBg$calibEvalOccsBg[[countExents]] <- list()
			# occsBg$calibEvalOccsBg[[countExents]]$occsBg <- thisOccsBg
			# occsBg$calibEvalOccsBg[[countExents]]$folds <- folds
		
		# }
		
		# names(occsBg$calibEvalOccsBg) <- paste0('occsBg_extent', exts, 'km')
	
		# ### make a map of calibration/evaluation sites

		# calibFold <- 1
		# evalFold <- 2
		# cex <- 0.6
		# pch <- 3

		# png('./figures_and_tables/example_of_training_calibration_evaluation_folds.png', width=1800, height=600)
		
			# par(mfrow=c(1, 3), oma=rep(0, 4), cex.main=2.2)
			
			# plot(studyRegionCroppedToPresentSpAlb, main=('Training, calibration,\nand evaluation occurrences'))
			# x <- occsBg$calibEvalOccsBg[[1]]$occsBg
			# x <- x[x$presBg == 1, ]
			# x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			# x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			# points(x, pch=pch, cex=cex)
			# points(x[x$fold==calibFold, ], pch=pch, cex=cex, col='cyan')
			# points(x[x$fold==evalFold, ], pch=pch, cex=cex, col='magenta')
			
			# legend('bottomright', inset=0.1, legend=c('Training', 'Calibration', 'Evaluation'), pch=pch, cex=2.2, col=c('black', 'cyan', 'magenta'), bty='n')
			
			# plot(studyRegionCroppedToPresentSpAlb, main=('Training and calibration\nbackground sites (narrowest extent)'))
			# x <- occsBg$calibEvalOccsBg[[1]]$occsBg
			# x <- x[x$presBg == 0, ]
			# x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			# x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			# points(x, pch=pch, cex=cex)
			# points(x[x$fold==calibFold, ], pch=pch, cex=cex, col='cyan')
			
			# plot(studyRegionCroppedToPresentSpAlb, main=('Evaluation background sites'))
			# x <- occsBg$testBg
			# x <- SpatialPointsDataFrame(x[ , ll], data=x, proj4string=getCRS('wgs84', TRUE))
			# x <- sp::spTransform(x, getCRS('albersNA', TRUE))
			# points(x, pch=pch, cex=cex)
			# points(x[x$fold==evalFold, ], pch=pch, cex=cex, col='magenta')
			
		# dev.off()
		
		# save(occsBg, file=paste0('./species_records/02_', gsub(tolower(species), pattern=' ', replacement='_'), '_collated_occurrence_and_background_sites.rda'))
		
# say('########################')
# say('### calibrate models ###')
# say('########################')

	# say('This portion of the script calibrates ENMs using withheld calibration data.  Multiple models are trained using the same set of occurrence and background sites, then evaluated against a set of withheld occurrence and background sites.  The best model is chosen among these.  The best model is then evaluated using a second set of withheld sites.', breaks=80)

	# load(paste0('./species_records/02_', gsub(tolower(species), pattern=' ', replacement='_'), '_collated_occurrence_and_background_sites.rda'))

	# metrics <- c('logLoss', 'cbi', 'trainSe95')
	
	# set.seed(123)
	
	# for (countExt in seq_along(exts)) {

		# ext <- exts[countExtent]
	
		# # collate data
		# data <- occsBg$calibEvalOccsBg[[countExtent]]$occsBg
		# rownames(data) <- 1:nrow(data)
		# vars1and2 <- c(paste0('ccsm_', predictors), paste0('ecbilt_', predictors))
		# data$presBg <- as.numeric(data$presBg)
		# folds <- occsBg$calibEvalOccsBg[[countExtent]]$folds
		# weight <- data$weight
	
		# for (gcm in gcms) {
		# # for (gcm in gcms[1]) { # ccsm

			# vars <- paste0(gcm, '_', predictors)

			# # GLM
			# algo <- 'glm'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainGlm, metrics=metrics, w=weight, na.rm=TRUE, out='tuning')
			# save(calib, file=paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # MaxEnt
			# algo <- 'maxent'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainMaxEnt, metrics=metrics, w=weight, na.rm=TRUE, out='tuning', jackknife=FALSE, scratchDir=tempDir, cores=6)
			# save(calib, file=paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # BRTs
			# algo <- 'brt'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainBrt, metrics=metrics, w=weight, na.rm=TRUE, cores=6, out='tuning', maxTrees=8000, anyway=TRUE)
			# save(calib, file=paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
			# # NS
			# algo <- 'ns'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# calib <- trainByCrossValid(data=data, resp='presBg', preds=vars, folds=folds, trainFx=trainNs, metrics=metrics, w=weight, na.rm=TRUE, out='tuning')
			# save(calib, file=paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
		# } # next GCM
			
	# } # next extent

# say('###########################################')
# say('### calibrate and evaluate final models ###')
# say('###########################################')

	# say('Make one model per extent, GCM, and algorithm using all occurrence data. Call this the "final" model. Make k-fold models with all data except that from each fold (combine training and calibration data from previous step and use just evaluation data as the withheld portion.', breaks=80)

	# # proportion of models for which a variable/feature must be included for the variable/feature to be included in the final model
	# bestThreshGlm <- 1 / 3

	# set.seed(123)
	
	# load(paste0('./species_records/02_', gsub(tolower(species), pattern=' ', replacement='_'), '_collated_occurrence_and_background_sites.rda'))

	# evals <- data.frame() # k-fold evaluation statistics
	
	# for (countExtent in seq_along(exts)) {

		# ext <- exts[countExtent]
	
		# # collate data
		# data <- occsBg$calibEvalOccsBg[[countExtent]]$occsBg
		# rownames(data) <- 1:nrow(data)
		# vars1and2 <- c(paste0('ccsm_', predictors), paste0('ecbilt_', predictors))
		# data$presBg <- as.numeric(data$presBg)
		# weight <- data$weight

		# # get test folds... training data will consists of (prior) training data plus calibration fold
		# # test data is the same... complicated because "folds" matrix was intended to cycle across all possible
		# # combinations of calibration/evaluation folds
		# folds <- occsBg$calibEvalOccsBg[[countExtent]]$folds
		# foldsNames <- colnames(folds)
		# foldsNames <- strsplit(foldsNames, '_vs_eval')
		# evalFolds <- integer()
		# for (i in seq_along(foldsNames)) evalFolds <- c(evalFolds, as.integer(foldsNames[[i]][2]))
		
		# foldColumnsToModel <- integer()
		# numFolds <- max(evalFolds)
		# for (i in 1:numFolds) {
			# foldColumnsToModel <- c(foldColumnsToModel, min(which(evalFolds == i)))
		# }
		# foldInteger <- evalFolds[foldColumnsToModel]

		# for (gcm in gcms) {

			# vars <- paste0(gcm, '_', predictors)

			# ### GLM
			# algo <- 'glm'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			
			# load(paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			# synopsis <- summaryByCrossValid(calib, trainFxName='trainGlm', metric='cbiTest')
			# terms <- synopsis$term[synopsis$proportionOfModels >= bestThreshGlm]
			# form <- paste0('presBg ~ 1 + ', paste(terms, collapse=' + '))
			
			# # full model
			# model <- glm(form, data=data, family='binomial', weights=weight)
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
			# # k-fold models
			# for (subK in 1:numFolds) {

				# # get column of folds matrix to use and the actual number identifying this fold
				# thisFoldCol <- foldColumnsToModel[subK]
				# thisFoldInteger <- foldInteger[subK]
			
				# # training and test data
				# thisFolds <- folds[ , thisFoldCol]
				# trains <- which(!is.na(thisFolds))
				# tests <- which(is.na(thisFolds))
				# thisTrainData <- data[trains, ]
				# thisTestPres <- data[tests, ]
				# thisTestPres <- thisTestPres[thisTestPres$presBg == 1, ]
				# thisTestBg <- occsBg$testBg[occsBg$testBg$fold == foldInteger[subK], ]
				
				# # model and evaluation
				# kModel <- glm(form, data=thisTrainData, family='binomial', weights=thisTrainData$weight)
				# predsAtPres <- predict(kModel, thisTestPres, type='response')
				# predsAtBg <- predict(kModel, thisTestBg, type='response')
				
				# cbi <- contBoyce(predsAtPres, predsAtBg, na.rm=TRUE)
				
				# evals <- rbind(
					# evals,
					# data.frame(
						# algo = algo,
						# gcm = gcm,
						# extent_km = ext,
						# fold = thisFoldInteger,
						# cbi = cbi
					# )
				# )
				
			# }
			
			# ### Maxent
			# algo <- 'maxent'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			
			# load(paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			# synopsis <- summaryByCrossValid(calib, trainFxName='trainMaxEnt', metric='cbiTest')
			
			# feats <- synopsis$featureSet[which.max(synopsis$frequencyInBestModels)]
			# regMult <- synopsis$meanRegMult[which.max(synopsis$frequencyInBestModels)]
			
			# args <- c(
				# paste0('betamultiplier=', regMult),
				# paste0('linear=', ifelse(grepl('l', feats), 'true', 'false')),
				# paste0('product=', ifelse(grepl('p', feats), 'true', 'false')),
				# paste0('quadratic=', ifelse(grepl('q', feats), 'true', 'false')),
				# paste0('hinge=', ifelse(grepl('h', feats), 'true', 'false')),
				# paste0('threshold=', ifelse(grepl('t', feats), 'true', 'false')),
				# paste0('jackknife=true')
			# )
			
			# model <- maxent(data[ , vars], p=data$presBg, args=args)
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # k-fold models
			# for (subK in 1:numFolds) {

				# # get column of folds matrix to use and the actual number identifying this fold
				# thisFoldCol <- foldColumnsToModel[subK]
				# thisFoldInteger <- foldInteger[subK]
			
				# # training and test data
				# thisFolds <- folds[ , thisFoldCol]
				# trains <- which(!is.na(thisFolds))
				# tests <- which(is.na(thisFolds))
				# thisTrainData <- data[trains, ]
				# thisTestPres <- data[tests, ]
				# thisTestPres <- thisTestPres[thisTestPres$presBg == 1, ]
				# thisTestBg <- occsBg$testBg[occsBg$testBg$fold == foldInteger[subK], ]
				
				# # model and evaluation
				# kModel <- maxent(thisTrainData[ , vars], p=thisTrainData$presBg, args=args)
				# predsAtPres <- predictMaxEnt(kModel, thisTestPres, type='cloglog')
				# predsAtBg <- predictMaxEnt(kModel, thisTestBg, type='cloglog')
				
				# cbi <- contBoyce(predsAtPres, predsAtBg, na.rm=TRUE)
				
				# evals <- rbind(
					# evals,
					# data.frame(
						# algo = algo,
						# gcm = gcm,
						# extent_km = ext,
						# fold = thisFoldInteger,
						# cbi = cbi
					# )
				# )
				
			# }

			# ### BRTs
			# algo <- 'brt'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# load(paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			# synopsis <- summaryByCrossValid(calib, trainFxName='trainBrt', metric='cbiTest')
			
			# learningRate <- synopsis$learningRate[synopsis$value == 'mean']
			# treeComplexity <- synopsis$treeComplexity[synopsis$value == 'mean']
			# bagFraction <- synopsis$bagFraction[synopsis$value == 'mean']
			
			# model <- trainBrt(data=data, resp='presBg', preds=vars, learningRate=learningRate, treeComplexity=treeComplexity, bagFraction=bagFraction, nTrees=8000, w=weight, anyway=TRUE)
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))

			# # k-fold models
			# for (subK in 1:numFolds) {

				# # get column of folds matrix to use and the actual number identifying this fold
				# thisFoldCol <- foldColumnsToModel[subK]
				# thisFoldInteger <- foldInteger[subK]
			
				# # training and test data
				# thisFolds <- folds[ , thisFoldCol]
				# trains <- which(!is.na(thisFolds))
				# tests <- which(is.na(thisFolds))
				# thisTrainData <- data[trains, ]
				# thisTestPres <- data[tests, ]
				# thisTestPres <- thisTestPres[thisTestPres$presBg == 1, ]
				# thisTestBg <- occsBg$testBg[occsBg$testBg$fold == foldInteger[subK], ]
				
				# # model and evaluation
				# kModel <- trainBrt(data=thisTrainData, resp='presBg', preds=vars, learningRate=learningRate, treeComplexity=treeComplexity, bagFraction=bagFraction, nTrees=8000, w=thisTrainData$weight, anyway=TRUE)

				# if (!is.null(kModel)) {

					# predsAtPres <- predictEnmSdm(kModel, thisTestPres)
					# predsAtBg <- predictEnmSdm(kModel, thisTestBg)
					
					# cbi <- contBoyce(predsAtPres, predsAtBg, na.rm=TRUE)

					# evals <- rbind(
						# evals,
						# data.frame(
							# algo = algo,
							# gcm = gcm,
							# extent_km = ext,
							# fold = thisFoldInteger,
							# cbi = cbi
						# )
					# )
					
				# # model did not converge or was insufficient
				# } else {
				
					# evals <- rbind(
						# evals,
						# data.frame(
							# algo = algo,
							# gcm = gcm,
							# extent_km = ext,
							# fold = thisFoldInteger,
							# cbi = NA
						# )
					# )
					
				# }
					
			# }
			
			# ### NS
			# algo <- 'ns'
			# say(gcm, ' ', ext, '-km extent tuning with ', algo, level=2)
			# load(paste0('./models/external_calibration_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			# synopsis <- summaryByCrossValid(calib, trainFxName='trainNs', metric='cbiTest')

			# nsPreds <- synopsis$term
			# nsDfs <- pmax(1, round(synopsis$meanDf))
			
			# form <- paste0('presBg ~ 1 + ', paste(paste('splines::ns(', nsPreds, ', df=', nsDfs, ')', sep=''), collapse=' + '))
			# model <- glm(form, data=data, family='binomial', weights=weight)
			# save(model, file=paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
			# # k-fold models
			# for (subK in 1:numFolds) {

				# # get column of folds matrix to use and the actual number identifying this fold
				# thisFoldCol <- foldColumnsToModel[subK]
				# thisFoldInteger <- foldInteger[subK]
			
				# # training and test data
				# thisFolds <- folds[ , thisFoldCol]
				# trains <- which(!is.na(thisFolds))
				# tests <- which(is.na(thisFolds))
				# thisTrainData <- data[trains, ]
				# thisTestPres <- data[tests, ]
				# thisTestPres <- thisTestPres[thisTestPres$presBg == 1, ]
				# thisTestBg <- occsBg$testBg[occsBg$testBg$fold == foldInteger[subK], ]
				
				# # model and evaluation
				# kModel <- glm(form, data=thisTrainData, family='binomial', weights=thisTrainData$weight)
				# predsAtPres <- predict(kModel, thisTestPres, type='response')
				# predsAtBg <- predict(kModel, thisTestBg, type='response')
				
				# cbi <- contBoyce(predsAtPres, predsAtBg, na.rm=TRUE)
				
				# evals <- rbind(
					# evals,
					# data.frame(
						# algo = algo,
						# gcm = gcm,
						# extent_km = ext,
						# fold = thisFoldInteger,
						# cbi = cbi
					# )
				# )
				
			# }

		# } # next GCM
			
	# } # next extent
	
	# write.csv(evals, './figures_and_tables/final_model_evaluations.csv', row.names=FALSE)
	
# say('#######################################################')
# say('### plot model performance against present-day data ###')
# say('#######################################################')

	# evals <- read.csv('./figures_and_tables/final_model_evaluations.csv')

	# evals$combos <- paste0(evals$gcm, ' ', evals$ext, '-km extent ', evals$algo)
	
	# n <- length(unique(evals$combos))

	# png('./figures_and_tables/final_model_evaluations.png', width=1200, height=800)

		# ylim <- c(roundTo(min(0, evals$cbi), 0.1, floor), 1)

		# par(oma=c(15, 2, 1, 1), mar=c(5, 5, 4, 2), cex.lab=2, cex.axis=1.8)
		# boxplot(cbi ~ combos, data=evals, ylim=ylim, ylab='CBI', xlab='', col='cornflowerblue', names=NA, xpd=NA)

		# labels <- sort(unique(evals$combos))
		# ylen <- diff(ylim)
		# text(1:n, rep(ylim[1] - 0.07 * ylen, n), labels=labels, srt=90, xpd=NA, cex=1.8, adj=c(1, 0.5))
		
	# dev.off()
	
# say('###############################################')
# say('### assess differences between model output ###')
# say('###############################################')
	
	# say('To help select the final models, we will do a cluster analysis on the rasters for the predictions for suitability in the present-day and 21 Kybp. Distance between rasters is calculated as RMSD.', breaks=80)
	
	# if (exists('sqPred')) rm(sqPred)
	# if (exists('lgmPred')) rm(lgmPred)
	
	# # get study region rasters for present and 21 Kybp, rescale so fully-available land cells are 1 and fully-covered glacial cells are NA
	# studyRegionRasts <- brick(studyRegionRastsFileName)
	# sqMask <- studyRegionRasts[[nlayers(studyRegionRasts)]]
	# lgmMask <- studyRegionRasts[[1]]

	# # rescaled so cell with partial glacier range from 0 to 1 (1 = fully open, <1 = some glacier)
	# sqGradMask <- calc(sqMask, fun=function(x) ifelse(x == 1, NA, 1 - x))
	# lgmGradMask <- calc(lgmMask, fun=function(x) ifelse(x == 1, NA, 1 - x))

	# # remove all cells with full glacier
	# sqMask <- calc(sqMask, fun=function(x) ifelse(x == 1, NA, 1))
	# lgmMask <- calc(lgmMask, fun=function(x) ifelse(x == 1, NA, 1))
	
	# for (gcm in gcms) {
	
		# sqClim <- getClimRasts(gcm=gcm, year=0, variables=predictors, rescale=TRUE, fillCoasts=FALSE)
		# lgmClim <- getClimRasts(gcm=gcm, year=21000, variables=predictors, rescale=TRUE, fillCoasts=FALSE)
	
		# sqClim <- projectRaster(sqClim, sqMask)
		# lgmClim <- projectRaster(lgmClim, lgmMask)

		# sqNames <- names(sqClim)
		# sqClim <- sqClim * sqMask
		# names(sqClim) <- sqNames
	
		# lgmNames <- names(lgmClim)
		# lgmClim <- lgmClim * lgmMask
		# names(lgmClim) <- lgmNames
	
		# for (ext in exts) {
	
			# for (algo in c('brt', 'glm', 'maxent', 'ns')) {
			
				# say(gcm, ' ', ext, '-km extent ', algo)
			
				# load(paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
			
				# thisSqPred <- raster::predict(sqClim, model, fun=enmSdm::predictEnmSdm)
				# thisLgmPred <- raster::predict(lgmClim, model, fun=enmSdm::predictEnmSdm)
				
				# thisSqPred <- thisSqPred * sqGradMask
				# thisLgmPred <- thisLgmPred * lgmGradMask
				
				# names(thisSqPred) <- names(thisLgmPred) <- paste0(gcm, '_', ext, 'kmExtent_', algo)
				
				# sqPred <- if (exists('sqPred')) {
					# stack(sqPred, thisSqPred)
				# } else {
					# thisSqPred
				# }
				
				# lgmPred <- if (exists('lgmPred')) {
					# stack(lgmPred, thisLgmPred)
				# } else {
					# thisLgmPred
				# }
			
			# } # next algorithm
	
		# } # extent
	
	# } # next GCM

	# ### cluster analysis
	# # take sum of square differences as measure of difference between rasters
	# sqDiffs <- lgmDiffs <- matrix(NA, nrow=nlayers(sqPred), ncol=nlayers(sqPred))
	# colnames(sqDiffs) <- colnames(lgmDiffs) <- names(sqPred)

	# for (i in 1:nlayers(sqPred)) {
	
		# for (j in 1:nlayers(sqPred)) {
	
			# # present-day
			# diffs <- sqPred[[i]] - sqPred[[j]]
			# n <- cellStats(diffs * 0 + 1, 'sum')
			# diffs <- diffs^2
			# diffs <- cellStats(diffs, 'sum')
			# diffs <- diffs / n
			# diffs <- sqrt(diffs)
			# sqDiffs[j, i] <- diffs
			
			# # LGM
			# diffs <- lgmPred[[i]] - lgmPred[[j]]
			# n <- cellStats(diffs * 0 + 1, 'sum')
			# diffs <- diffs^2
			# diffs <- cellStats(diffs, 'sum')
			# diffs <- diffs / n
			# diffs <- sqrt(diffs)
			# lgmDiffs[j, i] <- diffs
			
		# }
	
	# }
		
	# sqDiffs <- lgmDiffs <- matrix(NA, nrow=nlayers(sqPred), ncol=nlayers(sqPred))
	# colnames (sqDiffs) <- rownames(sqDiffs) <- colnames(lgmDiffs) <- rownames(lgmDiffs) <- names(sqPred)
		
	# for (i in 1:nlayers(sqPred)) {
		# for (j in 1:nlayers(sqPred)) {
			
			# diffs <- sqPred[[i]] - sqPred[[j]]
			# diffs <- diffs^2
			# n <- cellStats(diffs * 0 + 1, 'sum')
			# sqDiffs[i, j] <- cellStats(diffs, 'sum') / n
			
			# diffs <- lgmPred[[i]] - lgmPred[[j]]
			# diffs <- diffs^2
			# n <- cellStats(diffs * 0 + 1, 'sum')
			# lgmDiffs[i, j] <- cellStats(diffs, 'sum') / n
			
		# }
	# }
				
	
	# clust0ybp <- agnes(as.dist(sqDiffs), diss=TRUE, method='average')
	# clust21000ybp <- agnes(as.dist(lgmDiffs), diss=TRUE, method='average')
	
	# save(clust0ybp, file='./figures_and_tables/agnes_cluster_of_models_based_on_rasters_0ybp.rda')
	# save(clust21000ybp, file='./figures_and_tables/agnes_cluster_of_models_based_on_rasters_21000ybp.rda')
	
	# n <- ncol(sqDiffs)

	# png('./figures_and_tables/clustering_of_prediction_rasters.png', width=2200, height=800)

		# par(mfrow=c(1, 2), cex.main=2.2, cex.lab=1.6, cex.axis=1.4)

		# # present-day
		# clustMembership <- cutree(clust21000ybp, k=numModelClusts)
		# names(clustMembership) <- colnames(as.matrix(clust21000ybp$diss))

		# group <- clustMembership[match(clust0ybp$order.lab, names(clustMembership))]
		# cols <- clustCols[group]

		# plot(clust0ybp, main='Present-day differences', which.plot=2, ylab='RMSD', xlab=NA, cex=1.8)
		# usr <- par('usr')
		# y <- usr[3] - 0.03 * (usr[4] - usr[3])
		# points(1:n, rep(y, n), pch=16, col=cols, xpd=NA, cex=4)
		# y <- usr[3] - 0.07 * (usr[4] - usr[3])
		# text(n / 2, y, labels='LGM Group', xpd=NA, cex=2.2)

		# # LGM
		# clustMembership <- cutree(clust21000ybp, k=numModelClusts)
		# names(clustMembership) <- colnames(as.matrix(clust21000ybp$diss))

		# group <- clustMembership[match(clust21000ybp$order.lab, names(clustMembership))]
		# cols <- clustCols[group]

		# plot(clust21000ybp, main='LGM differences', which.plot=2, ylab='RMSD', xlab=NA, cex=1.8)
		# usr <- par('usr')
		# y <- usr[3] - 0.03 * (usr[4] - usr[3])
		# points(1:n, rep(y, n), pch=16, col=cols, xpd=NA, cex=4)
		# y <- usr[3] - 0.07 * (usr[4] - usr[3])
		# text(n / 2, y, labels='LGM Group', xpd=NA, cex=2.2)
		
	# dev.off()

# say('###################################')
# say('### project models back in time ###')
# say('###################################')

	# say('Write prediction rasters. Climate layers are interpolated linearly to 30-yr time periods. Predictions are made to these layers then re-projected to an equal-area projection and masked by the study region rasters.', breaks=80)
	
	# # get study region rasters for present and 21 Kybp, rescale so fully-available land cells are 1 and fully-covered glacial cells are NA
	# studyRegionRasts <- brick(studyRegionRastsFileName)
	# names(studyRegionRasts) <- paste0('year', seq(21000, 0, by=-30), 'ybp')
	
	# # time periods represented by rasters
	# climYears <- seq(21000, 0, by=-500)
	
	# for (gcm in gcms) {

		# for (ext in exts) {

			# for (algo in algos) {
			
				# say(paste(gcm, ext, algo))
			
				# load(paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
		
				# if (exists('preds')) rm(preds)
				# for (climYear in climYears) {

					# # get climate data
					# clim <- getClimRasts(gcm=gcm, year=climYear, variables=predictors, rescale=TRUE, fillCoasts=FALSE)
					
					# thisPred <- raster::predict(clim, model, fun=enmSdm::predictEnmSdm)
					
					# preds <- if (exists('preds')) {
						# stack(preds, thisPred)
					# } else {
						# thisPred
					# }
		
				# } # next year

				# interpFrom <- -1 * climYears
				# interpTo <- seq(-21000, 0, by=30)
				# preds <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')

				# # project
				# preds <- projectRaster(preds, studyRegionRasts)
				# preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
				# preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))
				
				# # mask by study region and force values to be within [0, 1] (can get pushed outside this during re-projection)
				# for (i in 1:nlayers(preds)) {
					
					# landMask <- (1 - studyRegionRasts[[i]])
					# preds[[i]] <- preds[[i]] * landMask
					
				# }
				
				# names(preds) <- paste0('ybp', seq(21000, 0, by=-30))
				# writeRaster(preds, paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo))
					
			# } # next algorithm
			
		# } # next extent
		
	# } # next GCM

# say('############################################################')
# say('### determine thresholds that best recreate Little range ###')
# say('############################################################')

	# # get threshold for sensitivity at this value
	# seThold <- 0.9

	# thresholds <- data.frame()

	# for (gcm in gcms) {

		# for (ext in exts) {

			# for (algo in algos) {
			
				# say(paste(gcm, ext, algo))
			
				# rasts <- stack(paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))
				# present <- subset(rasts, 701)
				
				# # extract predictions from cells inside or overlapping Little's range
				# insideLittleList <- raster::extract(present, littleRangeSpAlb, weights=TRUE, normalizeWeights=FALSE, cellnumbers=TRUE)
				
				# inside <- data.frame()
				# for (i in seq_along(insideLittleList)) inside <- rbind(inside, insideLittleList[[i]])
				# aggSum <- aggregate(inside, by=list(inside$cell), 'sum')
				# aggMean <- aggregate(inside, by=list(inside$cell), 'mean', na.rm=TRUE)
				
				# inside <- data.frame(cell=aggSum$cell, pred=aggMean$value, weight=aggSum$weight)
				# if (any(inside$weight > 1)) inside$weight[inside$weight > 1] <- 1
				
				# # extract predictions to cells outside or overlapping Little's range
				# outside <- as.vector(present)
				# outside <- data.frame(
					# cell=1:ncell(present),
					# pred=outside,
					# weight=1
				# )
				
				# # adjust weights of predictions outside Little range by proportion of cell inside the range
				# for (i in 1:nrow(inside)) {
					# outside$weight[outside$cell == i] = 1 - inside$weight[i]
				# }
				
				# # AUC
				# auc <- aucWeighted(inside$pred, outside$pred, presWeight=inside$weight, contrastWeight=outside$weight, na.rm=TRUE)
				
				# # CBI
				# cbi <- contBoyce(inside$pred, outside$pred, presWeight=inside$weight, na.rm=TRUE)
				
				# # determine threshold that maximizes TSS
				# best <- c(NA, -Inf)
				# names(best) <- c('threshold', 'tss')
				# for (thold in seq(0, 1, by=0.0001)) {
				
					# tss <- tssWeighted(inside$pred, outside$pred, presWeight=inside$weight, contrastWeight=outside$weight, na.rm=TRUE, thresholds=thold)
					
					# if (tss > best[['tss']]) {
						# best[['threshold']] <- thold
						# best[['tss']] <- tss
					# }
				
				# }
				
				# stats <- thresholdStats(best[['threshold']], inside$pred, outside$pred, presWeight=inside$weight, contrastWeight=outside$weight, na.rm=TRUE)
				
				# thresholds <- rbind(
					# thresholds,
					# data.frame(
						# gcm = gcm,
						# extent = ext,
						# algorithm = algo,
						# threshold = best[['threshold']],
						# auc = auc,
						# cbi = cbi,
						# tss = best[['tss']],
						# se = stats[['sensitivity']],
						# sp = stats[['specificity']]
					# )
				# )
				
				
			# } # next algorithm
			
		# } # next extent
		
	# } # next GCM
	
	# write.csv(thresholds, './figures_and_tables/thresholds.csv', row.names=FALSE)

# say('##############################################')
# say('### make maps of unthresholded predictions ###')
# say('##############################################')

	# # study region
	# studyRegionRasts <- brick(studyRegionRastsFileName)

	# # map extent
	# plotExtent <- extent(namSpAlbStudyRegion)
	# plotExtent <- as(plotExtent, 'SpatialPolygons')
	# projection(plotExtent) <- projection(namSpAlbStudyRegion)
	
	# load('./figures_and_tables/biotic_velocities.rda')
	
	# ### load prediction stacks
	# preds <- list()
	# for (algo in algos) {
		# for (gcm in gcms) {
			# for (ext in exts) {
				# thesePreds <- brick(paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))
				# names(thesePreds) <- paste0('year', seq(21000, 0, by=-30), 'ybp')
				# preds[[length(preds) + 1]] <- thesePreds
				# names(preds)[[length(preds)]] <- paste0(gcm, '_', ext, 'kmExtent_', algo)
			# }
		# }
	# }

	# # clustering
	# load('./figures_and_tables/agnes_cluster_of_models_based_on_rasters_21000ybp.rda')
	# clustMembership <- cutree(clust21000ybp, k=numModelClusts)
	# names(clustMembership) <- colnames(as.matrix(clust21000ybp$diss))

	# # suitability raster colors
	
	# # cols <- c('gray83', '#f7fcfd','#e5f5f9','#ccece6','#99d8c9','#66c2a4','#41ae76','#238b45','#006d2c','#00441b')
	# cols <- c('gray83', '#ccece6', '#99d8c9', '#66c2a4', '#41ae76', '#238b45', '#006d2c', '#00441b')
	# breaks <- seq(0, 1, length.out=length(cols) + 1)
	
	# # colors for centroids and quantile markers
	# colAllCells <- 'black'
	# colConstantCells <- 'red'
	
	# currentLty <- 'dotted'
	# centroidPathLwd <- 2
	
	# ### plot
	# dirCreate('./figures_and_tables/series')
	
	# conts <- list() # contours of LGM rasters
	# interval <- 30
	# years <- seq(21000, 0, by=-30)
	# for (countYear in seq_along(years)) {
		
		# year <- years[countYear]
		# if (year %% 210 == 0) {
			
			# say(year)
			
			# png(paste0('./figures_and_tables/series/predicted_suitable_area_all_models_', prefix(21000 - year, 5), 'yr_after_21000ybp_', prefix(year, 5), 'ybp.png'), width=2400, height=1600)
				
				# par(mfrow=c(4, 7), oma=c(2, 2, 4, 2), mar=c(0, 0, 6, 0))
			
				# count <- 1
				# plot(0, 0, col='white', fg='white', xaxt='n', yaxt='n', main='', xlim=c(0, 1), ylim=c(0, 1), ann=FALSE)

				# # land
				# year500 <- 500 * ceiling(year / 500)
				# land <- getClimRasts('ccsm', year=year500, variables=predictors[1], rescale=FALSE)
				# land <- land * 0
				# land <- projectRaster(land, crs=projection(plotExtent))
				# land <- crop(land, plotExtent)

				# for (algo in algos) {
					# for (gcm in gcms) {
						# for (ext in exts) {

							# # empty plots for time slider
							# if (count %% 7 == 0) {
								# plot(NA, fg='white', bg=NA, xaxt='n', yaxt='n', main='', xlim=c(0, 1), ylim=c(0, 1), ann=FALSE)
								# count <- count + 1
							# }
							
							# # time slider
							# if (count == 22) {
							
								# mult <- 1.15 # relative height
								# x <- 0.7
								# y <- (1 - (year / 21000)) * 3.95 * mult * 1
								# lines(c(x, x), c(0, 3.95 * mult * 1), lwd=50, col='gray70', xpd=NA)
								# lines(c(x, x), c(0, y), lwd=50, col='gray20', xpd=NA)
								# text(x, 0.17 + y, labels=paste(year, 'ybp'), cex=5, xpd=NA)
							
								# # legend
								# y <- 0
								# x <- -0.02
								# points(x, y + 0, pch=1, col=colAllCells, cex=6, xpd=NA)
								# lines(c(x, x), y + c(0, 0.20), col=colAllCells, xpd=NA, lwd=centroidPathLwd)
								# text(x, y + 0.25, labels='Starting centroid and centroid path of all cells', srt=90, cex=4, adj=c(0, 0.5), xpd=NA, col=colAllCells)
							
								# x <- 0.1
								# points(x, y + 0, pch=1, col=colConstantCells, cex=6, xpd=NA)
								# lines(c(x, x), y + c(0, 0.20), col=colConstantCells, xpd=NA, lwd=centroidPathLwd)
								# text(x, y + 0.25, labels='Starting centroid and centroid path of constantly-exposed, never-ice cells', srt=90, cex=4, adj=c(0, 0.5), xpd=NA, col=colConstantCells)
							
								# x <- 0.22
								# lines(c(x, x), y + c(0, 0.06), col=colAllCells, xpd=NA, lwd=3)
								# lines(c(x, x), y + c(0.14, 0.20), col=colAllCells, xpd=NA, lwd=3, lty=currentLty)
								# text(x, y + 0.25, labels='Starting/current latitude of extreme 5th quantiles using all cells', srt=90, cex=4, adj=c(0, 0.5), xpd=NA, col=colAllCells)
							
								# x <- 0.34
								# lines(c(x, x), y + c(0, 0.06), col=colConstantCells, xpd=NA, lwd=3)
								# lines(c(x, x), y + c(0.14, 0.20), col=colConstantCells, xpd=NA, lwd=3, lty=currentLty)
								# text(x, y + 0.25, labels='Starting/current latitude of extreme 5th quantiles using constantly-exposed, never-ice cells', srt=90, cex=4, adj=c(0, 0.5), xpd=NA, col=colConstantCells)
							
							# }
							
							# # plot extent
							# thisOne <- paste0(gcm, '_', ext, 'kmExtent_', algo)
							# boxCol <- clustCols[clustMembership[[thisOne]]]
							# plot(plotExtent, border=NA, bg=NA, fg=NA, col=NA, ann=FALSE, main='')

							# # land
							# plot(land, col='gray90', legend=FALSE, add=TRUE)
							
							# # scenario class indicator (border)
							# plot(plotExtent, border=boxCol, lwd=10, add=TRUE, fg=NA, bg=NA)
							
							# # predictions
							# x <- preds[[thisOne]][[countYear]]
							# plot(x, legend=FALSE, add=TRUE, col=cols, breaks=breaks)
							
							# ice <- studyRegionRasts[[countYear]]
							# ice <- calc(ice, fun=function(x) ifelse(x == 1, 1, NA))
							# plot(ice, col='steelblue1', legend=FALSE, add=TRUE)
							
							# plot(namSpAlbStudyRegion, add=TRUE, lwd=0.2, border='gray20')
							
							# # contours for LGM densities
							# if (year == 21000) {
								# quant <- quantile(x, 0.95)
								# conts[[length(conts) + 1]] <- rasterToContour(x, level=quant)
								# names(conts)[length(conts)] <- thisOne
							# }
							
							# plot(conts[[thisOne]], add=TRUE)
							
							# #### quantile/centroid indicators

								# ### centroids

								# # starting centroids
								# startCentroid_anyCells <- velocities[velocities$timeFrom == -21000 & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('centroidLong', 'centroidLat')]
								# startCentroid_constantCells <- velocities[velocities$timeFrom == -21000 & velocities$timeTo == -21000 + interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, c('centroidLong', 'centroidLat')]

								# cex <- 3
								# points(startCentroid_constantCells, pch=1, cex=cex, col=colConstantCells)
								# points(startCentroid_anyCells, pch=1, cex=cex, col=colAllCells)
								
								# ### centroid paths

								# startCentroid_constantCells_path <- velocities[velocities$timeFrom <= -1 * year & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInContinuouslyExposedLand, c('centroidLong', 'centroidLat')]
								# startCentroid_allCells_path <- velocities[velocities$timeFrom <= -1 * year & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('centroidLong', 'centroidLat')]

								# lines(startCentroid_constantCells_path$centroidLong, startCentroid_constantCells_path$centroidLat, col=colConstantCells, lwd=2)
								# lines(startCentroid_allCells_path$centroidLong, startCentroid_constantCells_path$centroidLat, col=colAllCells, lwd=2)
								
								# ### northern/southern latitude markers
								
								# # quantile locations
								# right <- xmax(plotExtent)
								# left <- xmin(plotExtent) + 0.02 * (xmax(plotExtent) - xmin(plotExtent))
								# width <- right - left
								# extension <- 0.06
								
								# # starting latitude of 5th quantile
								# latConstant <- velocities[velocities$timeFrom == -21000 & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p05')]

								# latAll <- velocities[velocities$timeFrom == -21000 & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p05')]
								
								# lines(c(left, left + width * extension), c(latAll, latAll), col=colConstantCells, lwd=cex)
								# lines(c(left, left + width * extension), c(latAll, latAll), col=colAllCells, lwd=cex)

								# # current latitude of 5th quantile
								# if (year < 21000) {
									
									# latConstant <- velocities[velocities$timeFrom == -1 * year - interval & velocities$timeTo == -1 * year & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p05')]
								
									# latAll <- velocities[velocities$timeFrom == -1 * year - interval & velocities$timeTo == -1 * year & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p05')]
									
									# lines(c(left, left + width * extension), c(latConstant, latConstant), col=colConstantCells, lwd=cex, lty=currentLty)
									# lines(c(left, left + width * extension), c(latAll, latAll), col=colAllCells, lwd=cex, lty=currentLty)

								# }

								# # starting latitude of 95th quantile
								# latConstant <- velocities[velocities$timeFrom == -21000 & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p95')]

								# latAll <- velocities[velocities$timeFrom == -21000 & velocities$timeSpan == interval & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p95')]
								
								# lines(c(left, left + width * extension), c(latConstant, latConstant), col=colConstantCells, lwd=cex)
								
								# lines(c(left, left + width * extension), c(latAll, latAll), col=colAllCells, lwd=cex)

								# # current latitude of 95th quantile
								# if (year < 21000) {
								
									# latConstant <- velocities[velocities$timeFrom == -1 * year - interval & velocities$timeTo == -1 * year & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p95')]
								
									# latAll <- velocities[velocities$timeFrom == -1 * year - interval & velocities$timeTo == -1 * year & velocities$algo == algo & velocities$gcm == gcm & velocities$ext == ext & !velocities$onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, c('nsQuantLat_quant0p95')]
								
									# lines(c(left, left + width * extension), c(latConstant, latConstant), col=colConstantCells, lwd=cex, lty=currentLty)
									# lines(c(left, left + width * extension), c(latAll, latAll), col=colAllCells, lwd=cex, lty=currentLty)

								# }
								
							# main <- paste(toupper(gcm), '\n', ext, '-km extent ', toupper(algo), collapse='', sep='')
							# title(main, line=0, cex.main=3.3, xpd=NA)
							
							# count <- count + 1

						# }
					# }
				# }
				
			# # mtext(paste(year, 'YBP'), outer=TRUE, xpd=NA, cex=4, font=2)
			# mtext(date(), side=1, cex=1, outer=TRUE)

			# dev.off()
			
		# } # if plotting this year

	# } # next year
		
say('############################################')
say('### make maps of thresholded predictions ###')
say('############################################')

	# study region
	studyRegionRasts <- brick(studyRegionRastsFileName)

	# map extent
	plotExtent <- extent(namSpAlbStudyRegion)
	plotExtent <- as(plotExtent, 'SpatialPolygons')
	projection(plotExtent) <- projection(namSpAlbStudyRegion)
	
	### thresholds
	thresholds <- read.csv('./figures_and_tables/thresholds.csv')
	
	### load prediction stacks
	preds <- list()
	for (algo in algos) {
		for (gcm in gcms) {
			for (ext in exts) {
				thesePreds <- brick(paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))
				names(thesePreds) <- paste0('year', seq(21000, 0, by=-30), 'ybp')
				preds[[length(preds) + 1]] <- thesePreds
				names(preds)[[length(preds)]] <- paste0(gcm, '_', ext, 'kmExtent_', algo)
			}
		}
	}

	### plot
	dirCreate('./figures_and_tables/refugia_and_current_range')
	
	# land
	lgmLand <- getClimRasts('ccsm', year=21000, variables=predictors[1], rescale=FALSE)
	lgmLand <- lgmLand * 0
	lgmLand <- projectRaster(lgmLand, crs=projection(plotExtent))
	lgmLand <- crop(lgmLand, plotExtent)

	sqLand <- getClimRasts('ccsm', year=0, variables=predictors[1], rescale=FALSE)
	sqLand <- sqLand * 0
	sqLand <- projectRaster(sqLand, crs=projection(plotExtent))
	sqLand <- crop(sqLand, plotExtent)

	# simulation raster (same resolution/extent as raster used for genetic/demographic simulations)
	simRast <- raster(demoGeneticRasterTemplate)

	for (gcm in gcms) {

		for (ext in exts) {

			say(gcm, ' ', ext)

			png(paste0('./figures_and_tables/refugia_and_current_range/', gcm, '_', ext, '_extent.png'), width=2400, height=1350)
				
			par(mfrow=c(2, 4), oma=c(2, 2, 4, 2), mar=c(0, 0, 6, 0))

			for (algo in algos) {
		
				pred <- preds[[paste0(gcm, '_', ext, 'kmExtent_', algo)]]
		
				### LGM
				x <- pred[[1]]
				land <- lgmLand
				
				# thresholded predictions
				thold <- thresholds$threshold[thresholds$gcm == gcm & thresholds$extent == ext & thresholds$algorithm == algo]

				x <- assign_refugia_from_abundance_raster(x, simRast, threshold=thold)
				numRefugia <- cellStats(x[['id']], 'max')
				cols <- rainbow(numRefugia)
				cols <- alpha(cols, 0.8)

				# plot extent
				plot(plotExtent, border=NA, bg=NA, fg=NA, col=NA, ann=FALSE, main='')
				plot(land, col='gray80', legend=FALSE, add=TRUE)
				plot(x[['id']], legend=FALSE, add=TRUE, col=cols)
				
				ice <- studyRegionRasts[[1]]
				ice <- calc(ice, fun=function(x) ifelse(x == 1, 1, NA))
				# plot(ice, col=alpha('steelblue1', 0.5), legend=FALSE, add=TRUE)
				plot(ice, col='darkgray', legend=FALSE, add=TRUE)

				labelFig(paste0(toupper(algo), ' using ', toupper(gcm), ' with ', ext, '-km extent'), adj=c(0.65, 0.01), cex=4.2)
				
				labelFig('LGM', adj=c(0.08, -0.08), cex=3.8)
				text <- if (numRefugia == 1) { '1 refugium' } else { paste(numRefugia, 'refugia') }
				labelFig(text, adj=c(0.08, -0.14), cex=3.8)
				labelFig('ice', adj=c(0.5, -0.3), cex=3)
		
				### SQ
				x <- pred[[701]]
				land <- sqLand
				
				# thresholded predictions
				thold <- thresholds$threshold[thresholds$gcm == gcm & thresholds$extent == ext & thresholds$algorithm == algo]

				x <- x >= thold
				x <- calc(x, fun=function(x) ifelse(x == 0, NA, x))
				cols <- alpha('forestgreen', 0.7)

				# plot extent
				plot(plotExtent, border=NA, bg=NA, fg=NA, col=NA, ann=FALSE, main='')
				plot(land, col='gray80', legend=FALSE, add=TRUE)
				plot(littleRangeSpAlb, col='orange', border='darkorange4', lwd=2, add=TRUE)
				plot(x, legend=FALSE, add=TRUE, col=cols)
				
				tss <- thresholds$tss[thresholds$gcm == gcm & thresholds$extent == ext & thresholds$algorithm == algo]
				se <- thresholds$se[thresholds$gcm == gcm & thresholds$extent == ext & thresholds$algorithm == algo]
				sp <- thresholds$sp[thresholds$gcm == gcm & thresholds$extent == ext & thresholds$algorithm == algo]
				
				legend('bottomright', inset=c(0.08, 0.17), legend=c('Little', 'Predicted'), fill=c('orange', alpha('forestgreen', 0.7)), border=c('darkorange4', NA), bty='n', cex=2.5, xpd=NA)
				
				labelFig('Present', adj=c(0.08, -0.08), cex=3.8)
				labelFig(paste0('TSS: ', sprintf('%.2f', tss), ' | Se = ', sprintf('%.2f', se), ' | Sp = ', sprintf('%.2f', sp)), adj=c(0.08, -0.14), cex=3.8)
		
			} # next algorithm
			
			mtext('LGM and Present-Day Predicted Distributions', side=3, cex=3.8, outer=TRUE, line=-0.5)
			mtext(date(), side=1, cex=1, outer=TRUE)

			dev.off()
			
		} # next extent

	} # next GCM
		
# say('#################################')
# say('### calculate biotic velocity ###')
# say('#################################')

	# say('Cycle through: time intervals (30, ..., 21000 yr); whether or not to consider velocity in only shared cells or all cells; and whether or not to examine velocity in only cells that never had ice and were always land across all time periods.', breaks=80, post=2)

	# # times represented by suitability rasters
	# times <- seq(-21000, 0, by=30)

	# # raster with same resolution/extent as that used for demographic/genetic simulations
	# demoGeneticTemplate <- raster(demoGeneticRasterTemplate)
	
	# # time intervals at which to calculate velocities
	# # intervals <- c(30, 60, 120, 240, 480, 990, 21000)
	# intervals <- c(30, 990, 21000)

	# # to store it all
	# velocities <- data.frame()
	
	# # allowing land and glaciers to shift
	# for (interval in intervals) {
	
		# atTimes <- seq(-21000, 0, by=interval)
			
		# for (ext in exts) {
			# for (gcm in gcms) {
				# for (algo in algos) {

					# # get predictions
					# preds <- brick(paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))
					
					# # resample to same resolution as demographic/genetic simulations then ensure values are in [0, 1]
					# preds <- resample(preds, demoGeneticTemplate)
					# preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
					# preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))
					
					# for (onlyInSharedCells in c(TRUE, FALSE)) {
			
						# say('ext ', ext, ' | gcm ', gcm, ' | algo ', algo, ' | interval ', interval, ' | shared cells ', onlyInSharedCells, ' | dynamic land | ', date())

						# # biotic velocity
						# thisVelocity <- bioticVelocity(preds, times=times, atTimes=atTimes, onlyInSharedCells=onlyInSharedCells, cores=4)

						# gc()

						# # remember
						# velocities <- rbind(
							# velocities,
							# cbind(
								# data.frame(
									# ext = tolower(ext),
									# gcm = gcm,
									# algo = algo,
									# onlyInSharedCells = onlyInSharedCells,
									# onlyInContinuouslyExposedLand = FALSE
								# ),
								# thisVelocity
							# )
						# )
								
					# } # next in shared cells
		
				# } # next algo
			# } # next GCM
		# } # next extent
			
	# } # next interval

	# studyRegion <- brick(studyRegionRastsFileName)
	
	# studyRegionExposedLandMask <- sum(studyRegion)
	# studyRegionExposedLandMask <- calc(studyRegionExposedLandMask, fun=function(x) ifelse(x %==na% 0, 1, NA))
	
	# # using only cells that were never covered by glaciers and always land
	# for (interval in intervals) {
	
		# atTimes <- seq(-21000, 0, by=interval)
			
		# for (ext in exts) {
			# for (gcm in gcms) {
				# for (algo in algos) {

					# say('ext ', ext, ' | gcm ', gcm, ' | algo ', algo, ' | interval ', interval, ' | shared cells ', onlyInSharedCells, ' | static land | ', date())

					# # get predictions
					# preds <- brick(paste0('./predictions/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))
					# preds <- studyRegionExposedLandMask * preds

					# preds <- resample(preds, demoGeneticTemplate)
				
					# preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
					# preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))

					# # biotic velocity
					# thisVelocity <- bioticVelocity(preds, times=times, atTimes=atTimes, onlyInSharedCells=onlyInSharedCells, cores=4)
					
					# gc()

					# # remember
					# velocities <- rbind(
						# velocities,
						# cbind(
							# data.frame(
								# ext = tolower(ext),
								# gcm = gcm,
								# algo = algo,
								# onlyInSharedCells = TRUE,
								# onlyInContinuouslyExposedLand = TRUE
							# ),
							# thisVelocity
						# )
					# )
							
				# } # next algorithm
			# } # next GCM
		# } # next extent
	# } # next interval

	# save(velocities, file='./figures_and_tables/biotic_velocities.rda')
	
# say('################################################')
# say('### plot biotic velocity for periods < 21 Ka ###')
# say('################################################')

	# # clustering
	# load('./figures_and_tables/agnes_cluster_of_models_based_on_rasters_21000ybp.rda')
	# clustMembership <- cutree(clust21000ybp, k=numModelClusts)
	# names(clustMembership) <- colnames(as.matrix(clust21000ybp$diss))

	# load('./figures_and_tables/biotic_velocities.rda')
	
	# # metrics <- c('centroidVelocity', 'nsCentroidVelocity', 'ewCentroidVelocity', 'nsQuantVelocity_quant0p05', 'nsQuantVelocity_quant0p95')
	# metrics <- c('centroidVelocity', 'nsCentroidVelocity', 'nsQuantVelocity_quant0p05', 'nsQuantVelocity_quant0p95')
	
	# for (metric in metrics) {
			
		# metricNice <- if (metric == 'centroidVelocity') {
			# 'Centroid Velocity'
		# } else if (metric == 'nsQuantVelocity_quant0p95') {
			# 'Northern Range Edge Velocity (95th quantile)'
		# } else if (metric == 'nsQuantVelocity_quant0p05') {
			# 'Southern Range Edge Velocity (5th quantile)'
		# } else if (metric == 'nsCentroidVelocity') {
			# 'North- vs Southward Centroid Movement'
		# } else if (metric == 'ewCentroidVelocity') {
			# 'East- vs Westward Centroid Movement'
		# }

		# png(paste0('./figures_and_tables/bioticVelocity_', metric, '.png'), width=1.5 * 1280, height=1.5 * 720)
			
			# par(mfrow=c(2, 3), oma=c(1, 1, 3, 1), cex.main=2.2, cex.lab=2, cex.axis=1.8)

			# for (interval in c(30, 990)) {
				
				# # get maximum velocity across models
				# minVel <- min(velocities[velocities$timeSpan == interval, metric])
				# maxVel <- max(velocities[velocities$timeSpan == interval, metric])

				# ### for cases where cells are not necessarily continuously-exposed and available
				# for (onlyInSharedCells in c(FALSE, TRUE)) {

					# # plot
					# main <- paste0(
						# # metricNice,
						# if (onlyInSharedCells) { 'Only Shared non-NA Cells' } else { 'All Cells' },
						# ' | ', interval, '-yr Intervals'
					# )

					# plot(0, 0, col='white', xlim=c(-21000, 0), ylim=c(minVel, maxVel), xlab='YBP', ylab='Velocity (m / y)', main=main)
					
					# for (i in seq(-21000, 0, by=500)) {
						# lines(c(i, i), c(minVel, maxVel), lwd=0.2, col='gray')
					# }
				
					# for (gcm in gcms) {
						# for (ext in exts) {
							# for (algo in algos) {
								
								# thisVel <- velocities[velocities$timeSpan == interval & velocities$gcm == gcm & velocities$ext == ext & velocities$algo == algo & velocities$onlyInSharedCells == onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, ]
								# modelGroup <- paste0(gcm, '_', ext, 'kmExtent_', algo)
								# group <- clustMembership[[modelGroup]]
								
								# col <- clustCols[group]
								
								# years <- rowMeans(thisVel[ , c('timeFrom', 'timeTo')])
								# lines(years, thisVel[ , metric], col=col, lwd=2)
								
							# }
						# }
					# }
					
					# # legend('topleft', inset=-0.01, legend=paste0('Group', 1:numModelClusts), col=clustCols[1:numModelClusts], lwd=1, cex=0.5, box.col='white')
					
				# } # next only in shared cells

				# # plot
				# main <- paste0('Continuously Terrestrial, Ice-free Cells', ' | ', interval, '-yr Intervals')

				# plot(0, 0, col='white', xlim=c(-21000, 0), ylim=c(minVel, maxVel), xlab='YBP', ylab='Velocity (m / y)', main=main)
				
				# for (i in seq(-21000, 0, by=500)) {
					# lines(c(i, i), c(minVel, maxVel), lwd=0.2, col='gray')
				# }
			
				# for (gcm in gcms) {
					# for (ext in exts) {
						# for (algo in algos) {
							
							# thisVel <- velocities[velocities$timeSpan == interval & velocities$gcm == gcm & velocities$ext == ext & velocities$algo == algo & velocities$onlyInSharedCells == onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, ]
							# modelGroup <- paste0(gcm, '_', ext, 'kmExtent_', algo)
							# group <- clustMembership[[modelGroup]]
							
							# col <- clustCols[group]
							
							# years <- rowMeans(thisVel[ , c('timeFrom', 'timeTo')])
							# lines(years, thisVel[ , metric], col=col, lwd=2)
							
						# }
					# }
				# }
				
					# # legend('topleft', inset=-0.01, legend=paste0('Group', 1:numModelClusts), col=clustCols[1:numModelClusts], lwd=1, cex=0.5, box.col='white')
					
			# } # next interval
		
			# title(sub=date(), outer=TRUE, line=-1)
			# title(main=metricNice, outer=TRUE, line=0.8, cex.main=2.9)
		
		# dev.off()
			
	# } # next metric
			
# say('#################################################')
# say('### plot biotic velocity for periods of 21 Ka ###')
# say('#################################################')

	# # clustering
	# load('./figures_and_tables/agnes_cluster_of_models_based_on_rasters_21000ybp.rda')
	# clustMembership <- cutree(clust21000ybp, k=numModelClusts)
	# names(clustMembership) <- colnames(as.matrix(clust21000ybp$diss))

	# load('./figures_and_tables/biotic_velocities.rda')
	
	# # metrics <- c('centroidVelocity', 'nsCentroidVelocity', 'ewCentroidVelocity', 'nsQuantVelocity_quant0p05', 'nsQuantVelocity_quant0p95')
	# metrics <- c('centroidVelocity', 'nsCentroidVelocity', 'nsQuantVelocity_quant0p05', 'nsQuantVelocity_quant0p95')
	
	# interval <- 21000
	
	# for (metric in metrics) {
			
		# if (metric == 'centroidVelocity') {
			# metricNice <- 'Centroid Velocity'
			# zeroLine <- FALSE
		# } else if (metric == 'nsQuantVelocity_quant0p95') {
			# metricNice <- 'Northern Range Edge Velocity (95th quantile)'
			# zeroLine <- TRUE
		# } else if (metric == 'nsQuantVelocity_quant0p05') {
			# metricNice <- 'Southern Range Edge Velocity (5th quantile)'
			# zeroLine <- TRUE
		# } else if (metric == 'nsCentroidVelocity') {
			# metricNice <- 'North- vs Southward Centroid Movement'
			# zeroLine <- TRUE
		# } else if (metric == 'ewCentroidVelocity') {
			# metricNice <- 'East- vs Westward Centroid Movement'
			# zeroLine <- TRUE
		# }

		# png(paste0('./figures_and_tables/bioticVelocity_21Ka_', metric, '.png'), width=1900, height=800, res=300)
		
			# par(mfrow=c(1, 3), oma=c(6.5, 0, 1.6, 1), mar=c(0, 4, 1.2, 0), mgp=c(1.2, 0.23, 0), cex.main=0.8, cex.lab=0.7, cex.axis=0.6, tck=-0.02, lwd=0.8)

			# # get extreme velocities across models
			# minVel <- min(velocities[velocities$timeSpan == interval, metric])
			# maxVel <- max(velocities[velocities$timeSpan == interval, metric])
			
			# ylim <- c(min(0, minVel), maxVel)

			# ### for cases where cells are not necessarily continuously-exposed and available
			# for (onlyInSharedCells in c(FALSE, TRUE)) {

				# # plot
				# main <- ifelse(onlyInSharedCells, 'Only Shared non-NA Cells', 'All Cells')

				# thisVel <- velocities[velocities$timeSpan == interval & velocities$onlyInSharedCells == onlyInSharedCells & !velocities$onlyInContinuouslyExposedLand, ]

				# # color and sort by group
				# thisVel$group <- clustMembership[match(paste0(thisVel$gcm, '_', thisVel$ext, 'kmExtent_', thisVel$algo), names(clustMembership))]
				# thisVel <- thisVel[order(thisVel$group), ]
				# thisVel$col <- clustCols[thisVel$group]

				# # plot
				# n <- nrow(thisVel)
				# labels <- paste0(thisVel$gcm, ' ', thisVel$ext, '-km extent ', thisVel$algo)
				# ylab <- 'Biotic velocity (m / y)'
				# plot(1:n, thisVel[ , metric], pch=21, cex=1, bg=thisVel$col, xaxt='n', ylab=ylab, xlab='', main=main, ylim=ylim, bty='n', xpd=NA)
				# if (zeroLine) lines(c(1, n), c(0, 0), col=1, lwd=0.6)
				# points(1:n, thisVel[ , metric], pch=21, cex=1, bg=thisVel$col)
				# text(1:n, thisVel[ , metric], labels=labels, adj=c(1.1, 0.2), xpd=NA, srt=90, col=thisVel$col, cex=0.6)
				
			# }
			
			# # plot
			# main <- 'Continuously-exposed, non-ice cells'

			# thisVel <- velocities[velocities$timeSpan == interval & velocities$onlyInSharedCells & velocities$onlyInContinuouslyExposedLand, ]

			# # color and sort by group
			# thisVel$group <- clustMembership[match(paste0(thisVel$gcm, '_', thisVel$ext, 'kmExtent_', thisVel$algo), names(clustMembership))]
			# thisVel <- thisVel[order(thisVel$group), ]
			# thisVel$col <- clustCols[thisVel$group]

			# # plot
			# n <- nrow(thisVel)
			# labels <- paste0(thisVel$gcm, ' ', thisVel$ext, '-km extent ', thisVel$algo)
			# ylab <- 'Biotic velocity (m / y)'
			# plot(1:n, thisVel[ , metric], pch=21, cex=1, bg=thisVel$col, xaxt='n', ylab=ylab, xlab='', main=main, ylim=ylim, bty='n', xpd=NA)
			# if (zeroLine) lines(c(1, n), c(0, 0), col=1, lwd=0.6)
			# points(1:n, thisVel[ , metric], pch=21, cex=1, bg=thisVel$col)
			# text(1:n, thisVel[ , metric], labels=labels, adj=c(1.1, 0.2), xpd=NA, srt=90, col=thisVel$col, cex=0.6)

			# title(main=paste('Biotic velocity for', metricNice, 'across a 21-Ka interval'), outer=TRUE)
			# title(sub=date(), outer=TRUE, cex.sub=0.3, line=5.5)

		# dev.off()
			
	# } # next metric
			
#############################################
say('DONE', deco='~', pre=2, post=2, level=1)
