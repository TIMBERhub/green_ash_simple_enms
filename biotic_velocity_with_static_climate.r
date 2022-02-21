### BIOTIC VELOCITY WITH STATIC CLIMATE
### Adam B. Smith | adam.smith@mobot.org | Missouri Botanical Garden | Spring 2021
###
### CONTENTS
### setup ###
### constants ###
### project models back in time using period-static (non-interpolated) climate ###
### calculate biotic velocity using period-static (non-interpolated) climate ###



#############
### setup ###
#############

	memory.limit(memory.limit() * 2^30)
	rm(list=ls()) # reproducibility!
	options(keep.source=FALSE) # manage memory
	gc()
	print('')
	print(date())

	### source('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/code/biotic_velocity_with_static_climate.r')

	setwd('E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms')
	lorenzPath <- 'E:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/Version 2017-06-16/'
	studyRegionRastsFileName <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'
	demoGeneticRasterTemplate <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_resampled_to_genetic_demographic_simulation_resolution.tif'
	elevationRastFileName <- 'E:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_elevationInMeters_fromEtop.tif'
	tempDir <- 'E:/ecology/!Scratch/_temp'
	
	options(stringsAsFactors=FALSE)
	raster::rasterOptions(format='GTiff', overwrite=TRUE)

	library(cowplot)
	library(ggplot2)
	library(raster)
	library(scales)
	library(sp)

	library(omnibus) # Adam's grab-bag library (https://github.com/adamlilith/omnibus)
	library(enmSdm) # Adam's SDM library (https://github.com/adamlilith/enmSdm)
	library(statisfactory) # Adam's statistics library (https://github.com/adamlilith/statisfactory)
	library(legendary) # Adam's plotting library (https://github.com/adamlilith/legenday)
	
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
	# exts <- c(80, 160, 320) # in km
	exts <- c(160) # in km

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
	# algos <- c('brt', 'glm', 'maxent', 'ns')
	algos <- c('maxent')
	
# say('##################################################################################')
# say('### project models back in time using period-static (non-interpolated) climate ###')
# say('##################################################################################')

	# say('Write prediction rasters. Climate layers are NOT interpolated but held constant. Each 30-yr time step is assigned the climate layers that are closest in time to it. Predictions are made to these layers then re-projected to an equal-area projection and masked by the study region rasters.', breaks=80)
	
	# # get study region rasters for present and 21 Kybp, rescale so fully-available land cells are 1 and fully-covered glacial cells are NA
	# studyRegionRasts <- brick(studyRegionRastsFileName)
	# names(studyRegionRasts) <- paste0('year', seq(21000, 0, by=-30), 'ybp')
	
	# # time periods represented by rasters
	# climYears <- seq(21000, 0, by=-500)
	
	# for (gcm in gcms) {
	# # for (gcm in 'ccsm') {

		# for (ext in exts) {
		# # for (ext in 160) {

			# for (algo in algos) {
			# # for (algo in 'maxent') {
			
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

				# # # interpFrom <- -1 * climYears
				# # # interpTo <- seq(-21000, 0, by=30)
				# # # preds <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')

				# # project
				# preds <- projectRaster(preds, studyRegionRasts)
				# preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
				# preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))
				
				# # mask by study region and force values to be within [0, 1] (can get pushed outside this during re-projection)
				# for (i in 1:nlayers(preds)) {
					
					# landMask <- (1 - studyRegionRasts[[i]])
					# preds[[i]] <- preds[[i]] * landMask
					
				# }
				
				# names(preds) <- paste0('ybp', climYears)
				# dirCreate('./predictions/climate_static_within_period')
				# writeRaster(preds, paste0('./predictions/climate_static_within_period/', gcm, '_', ext, 'kmExtent_', algo))
					
			# } # next algorithm
			
		# } # next extent
		
	# } # next GCM

# say('################################################################################')
# say('### calculate biotic velocity using period-static (non-interpolated) climate ###')
# say('################################################################################')

	# say('Cycle through: time period over which to calculate velocity; whether or not to consider velocity in only shared cells or all cells; and whether or not to examine velocity in only cells that never had ice and were always land across all time periods.', breaks=80, post=2)

	# # times represented by suitability rasters
	# genTimes <- seq(-21000, 0, by=30)
	# climTimes <- seq(-21000, 0, by=500)

	# # raster with same resolution/extent as that used for demographic/genetic simulations
	# demoGeneticTemplate <- raster(demoGeneticRasterTemplate)
	
	# # time intervals at which to calculate velocities
	# intervals <- c(30, 990)

	# # to store it all
	# velocities <- data.frame()
	
	# # allowing land and glaciers to shift
	# for (ext in exts) {
		# for (gcm in gcms) {
			# for (algo in algos) {

				# # get predictions
				# preds500yr <- brick(paste0('./predictions/climate_static_within_period/', gcm, '_', ext, 'kmExtent_', algo, '.tif'))

				# # resample to same resolution as demographic/genetic simulations then ensure values are in [0, 1]
				# preds500yr <- resample(preds500yr, demoGeneticTemplate)
				# preds500yr <- calc(preds500yr, fun=function(x) ifelse(x < 0, 0, x))
				# preds500yr <- calc(preds500yr, fun=function(x) ifelse(x > 1, 1, x))
				# names(preds500yr) <- paste0('ybp', abs(climTimes))
				
				# preds <- preds500yr[[1]]
				# names(preds) <- 'ybp21000'
				# for (genTime in genTimes[2:length(genTimes)]) {
					# closestIndex <- which.min(abs(climTimes - genTime))
					# preds <- stack(preds, preds500yr[[closestIndex]])
					# names(preds)[nlayers(preds)] <- paste0('ybp', climTimes[closestIndex], '_gen', abs(genTime))
				# }
				
				# for (onlyInSharedCells in c(TRUE, FALSE)) {
		
					# for (interval in intervals) {
			
						# say('ext ', ext, ' | gcm ', gcm, ' | algo ', algo, ' | shared cells ', onlyInSharedCells, ' | dynamic land | interval ', interval, ' | ', date())

						# intervalTimes <- seq(-21000, 0, by=interval)

						# # biotic velocity
						# thisVelocity <- bioticVelocity(preds, times=genTimes, atTimes=intervalTimes, onlyInSharedCells=onlyInSharedCells, cores=1)

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
						
					# } # next interval
							
				# } # next in shared cells
	
			# } # next algo
		# } # next GCM
	# } # next extent

	# save(velocities, file='./figures_and_tables/biotic_velocities_climate_static_within_period.rda')

say('####################################################################################################')
say('### plot biotic velocity using interpolated and non-interpolated climate (static period climate) ###')
say('####################################################################################################')

	load('./figures_and_tables/biotic_velocities_climate_static_within_period.rda')
	static <- velocities
	load('./figures_and_tables/biotic_velocities.rda')
	dynamic <- velocities

	staticExts <- unique(static$ext)
	staticGcms <- unique(static$gcm)
	staticAlgos <- unique(static$algo)
	staticOnlyInSharedCells <- unique(static$onlyInSharedCells)
	staticOnlyInContinuouslyExposedLand <- unique(static$onlyInContinuouslyExposedLand)
	staticTimeSpan <- unique(static$timeSpan)
	
	dynamic <- dynamic[dynamic$ext %in% staticExts & dynamic$gcm %in% staticGcms & dynamic$algo %in% staticAlgos & dynamic$onlyInSharedCells %in% staticOnlyInSharedCells & dynamic$onlyInContinuouslyExposedLand %in% staticOnlyInContinuouslyExposedLand & dynamic$timeSpan %in% staticTimeSpan, ]
	
	static$type <- 'static'
	dynamic$type <- 'dynamic'
	vels <- rbind(static, dynamic)
	vels$timeAt <- rowMeans(vels[ , c('timeFrom', 'timeTo')])
	
	data <- vels[vels$timeSpan == 30 & vels$onlyInSharedCells, ]
	sharedInterval30 <- ggplot(data, aes(x=timeAt, y=centroidVelocity, color=type)) +
		geom_line() +
		xlab('YBP') + ylab('Centroid velocity (m/y)') + ggtitle('Centroid velocity: Shared cells only, 30-yr periods') +
		facet_wrap(~gcm)
		
	sharedInterval30

	data <- vels[vels$timeSpan == 30 & !vels$onlyInSharedCells, ]
	allCellsInterval30 <- ggplot(data, aes(x=timeAt, y=centroidVelocity, color=type)) +
		geom_line() +
		xlab('YBP') + ylab('Centroid velocity (m/y)') + ggtitle('Centroid velocity: All cells, 30-yr periods') +
		facet_wrap(~gcm)
		
	allCellsInterval30

	data <- vels[vels$timeSpan == 990 & vels$onlyInSharedCells, ]
	sharedInterval990 <- ggplot(data, aes(x=timeAt, y=centroidVelocity, color=type)) +
		geom_line() +
		xlab('YBP') + ylab('Centroid velocity (m/y)') + ggtitle('Centroid velocity: Shared cells only, 990-yr periods') +
		facet_wrap(~gcm)
		
	sharedInterval990

	data <- vels[vels$timeSpan == 990 & !vels$onlyInSharedCells, ]
	allCellsInterval990 <- ggplot(data, aes(x=timeAt, y=centroidVelocity, color=type)) +
		geom_line() +
		xlab('YBP') + ylab('Centroid velocity (m/y)') + ggtitle('Centroid velocity: All cells, 990-yr periods') +
		facet_wrap(~gcm)
		
	allCellsInterval990

	main <- plot_grid(sharedInterval30, allCellsInterval30, sharedInterval990, allCellsInterval990, labels='auto', label_size=14, ncol=1, rel_widths=1)
	
	main
	
	ggsave('./figures_and_tables/exploring_biotic_velocity/velocities_with_period-static_and_dynamic_climate.pdf', width=8.5, height=11, units='in')
	


			
#############################################
say('DONE', deco='~', pre=2, post=2, level=1)
