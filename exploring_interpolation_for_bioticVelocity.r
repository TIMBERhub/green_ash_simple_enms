# source('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms/code/exploring_interpolation_for_bioticVelocity.r')

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	
	library(raster)
	library(omnibus)
	library(statisfactory)
	library(enmSdm)

	exts <- c(80, 160, 320) # in km
	gcms <- c('ccsm', 'ecbilt')
	# algos <- c('brt', 'glm', 'maxent', 'ns')
	algos <- c('glm')

	setwd('C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/enms')
	
	if (file.exists('./figures_and_tables/predictors.rda')) load('./figures_and_tables/predictors.rda')

	lorenzPath <- 'D:/Ecology/Climate/Lorenz et al 2016 North America 21Kybp to 2100 CE/'
	
		studyRegionRastsFileName <- 'C:/Ecology/Drive/Research/ABC vs Biogeography/NSF_ABI_2018_2021/data_and_analyses/green_ash/study_region/!study_region_raster_masks/study_region_daltonIceMask_lakesMasked_linearIceSheetInterpolation.tif'

	
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
		
			rast <- stack(paste0(lorenzPath, '/V2/', gcmFolder, '/', year, 'BP/', variable, '.tif'))
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
	
	
###################################
###################################
###################################	
	
	# get study region rasters for present and 21 Kybp, rescale so fully-available land cells are 1 and fully-covered glacial cells are NA
	studyRegionRasts <- brick(studyRegionRastsFileName)
	names(studyRegionRasts) <- paste0('year', seq(21000, 0, by=-30), 'ybp')
	
	# time periods represented by rasters
	climYears <- seq(21000, 0, by=-500)
	
	for (gcm in gcms) {

		for (ext in exts) {

			for (algo in algos) {
			
				say(paste(gcm, ext, algo))
			
				load(paste0('./models/final_model_for_', tolower(ext), 'km_extent_with_', algo, '_', gcm, '_gcm.rda'))
		
				if (exists('preds')) rm(preds)
				for (climYear in climYears) {

					# get climate data
					clim <- getClimRasts(gcm=gcm, year=climYear, variables=predictors, rescale=TRUE, fillCoasts=FALSE)
					
					thisPred <- raster::predict(clim, model, fun=enmSdm::predictEnmSdm)
					
					preds <- if (exists('preds')) {
						stack(preds, thisPred)
					} else {
						thisPred
					}
		
				} # next year

				interpFrom <- -1 * climYears
				interpTo <- seq(-21000, 0, by=30)
				
				predsLogit <- logitAdj(preds, epsilon=0)

				source('C:/Ecology/Drive/R/enmSdm/R/interpolateRasters.r')

				say('linear')
				predsLinear <- interpolateRasters(preds, interpFrom=interpFrom, interpTo=interpTo, type='linear')
				say('gam')
				predsGamLogit <- interpolateRasters(predsLogit, interpFrom=interpFrom, interpTo=interpTo, type='gam', family='gaussian', onFail='linear')
				say('poly')
				predsPolyLogit <- interpolateRasters(predsLogit, interpFrom=interpFrom, interpTo=interpTo, type='poly', family='gaussian', onFail='linear')
				say('bs')
				predsBsLogit <- interpolateRasters(predsLogit, interpFrom=interpFrom, interpTo=interpTo, type='bs', family='gaussian', onFail='linear')
				say('smooth')
				predsSmoothLogit <- interpolateRasters(predsLogit, interpFrom=interpFrom, interpTo=interpTo, type='smooth.spline', onFail='linear')

				predsGam <- probitAdj(predsGamLogit, epsilon=0)
				predsPoly <- probitAdj(predsPolyLogit, epsilon=0)
				predsBs <- probitAdj(predsBsLogit, epsilon=0)
				predsSmooth <- probitAdj(predsSmoothLogit, epsilon=0)

				xy <- rbind(c(-87.10108, 32.27589))
				
				linear <- c(extract(predsLinear, xy))
				gam <- c(extract(predsGam, xy))
				poly <- c(extract(predsPoly, xy))
				bs <- c(extract(predsBs, xy))
				ss <- c(extract(predsSmooth, xy))
				
				plot(interpTo, linear, type='l', ylim=c(0, 1), lwd=2)
				lines(interpTo, gam, col='red', lwd=2)
				lines(interpTo, poly, col='blue')
				lines(interpTo, bs, col='orange')
				lines(interpTo, ss, col='cyan')
				
				
				legend('topleft', inset=0.01, legend=c('linear', 'gam', 'poly', 'bs', 'ss'), col=c('black', 'red', 'blue', 'orange', 'cyan'), lwd=1)
				
				if (exists('bv')) rm(bv)
				bv <- data.frame()
				
				# biotic velocity for each type of interpolation
				for (this in c('predsLinear', 'predsGam', 'predsPoly', 'predsBs', 'predsSmooth')) {
				
					preds <- get(this)
				
					# project
					preds <- projectRaster(preds, studyRegionRasts)
					preds <- calc(preds, fun=function(x) ifelse(x < 0, 0, x))
					preds <- calc(preds, fun=function(x) ifelse(x > 1, 1, x))
					
					# mask by study region and force values to be within [0, 1] (can get pushed outside this during re-projection)
					for (i in 1:nlayers(preds)) {
						
						landMask <- (1 - studyRegionRasts[[i]])
						preds[[i]] <- preds[[i]] * landMask
						
					}
					
					names(preds) <- paste0('ybp', seq(21000, 0, by=-30))
						
					source('C:/Ecology/Drive/R/enmSdm/R/bioticVelocity.r')
					thisBv <- bioticVelocity(preds, time=interpTo, onlyInSharedCells=TRUE)
					
					thisBv$interp <- this
					
					bv <- rbind(bv, thisBv)
					
				}
						
				print(NON)	
					
			} # next algorithm
			
		} # next extent
		
	} # next GCM

	# maxs <- max(bv$centroidVelocity)
	maxs <- max(1000)
	
	plot(1, xlim=c(-21000, 0), ylim=c(0, maxs), col='white', main='Different methods for interpolating suitability layers', ylab='Centroid velocity (m / y)', xlab='Year')
	
	for (i in seq(-21000, 0, by=500)) lines(c(i, i), c(0, maxs), col='gray70')
	
	interps <- c('predsLinear', 'predsGam', 'predsPoly', 'predsBs', 'predsSmooth')
	for (this in seq_along(interps)) {
	
		interp <- interps[this]
		these <- bv[bv$interp == interp, ]
		
		lines(these$timeTo, these$centroidVelocity, col=this, lwd=2)
	
	}
	
	legend('topright', legend=interps, lwd=1, col=seq_along(interps), bg='white')
	
	
	