	metrics <- c('centroidVelocity', 'nsQuantVelocity_quant0p05', 'nsQuantVelocity_quant0p95')
	
	for (metric in metrics) {
			
		metricNice <- if (metric == 'centroidVelocity') {
			'Centroid Velocity'
		} else if (metric == 'nsQuantVelocity_quant0p95') {
			'Velocity of Northern Range Edge (95th quantile)'
		} else if (metric == 'nsQuantVelocity_quant0p05') {
			'Velocity of Southern Range Edge (5th quantile)'
		}
	
		png(paste0('./figures_and_tables/bioticVelocity_', metric, '.png'), width=1.5 * 1280, height=1.5 * 720)
			
			par(mfrow=c(2, 2), oma=c(1, 1, 3, 1), cex.main=2.2, cex.lab=2, cex.axis=1.8)

			for (interval in c(30, 990)) {
				
				# get maximum velocity across models
				minVel <- Inf
				maxVel <- -Inf
				for (onlyInSharedCells in c(TRUE, FALSE)) {
			
					cells <- if (onlyInSharedCells) { 'allCells'} else { 'sharedCells' }
					cells <- capIt(cells)
					x <- get(paste0('bv', interval, 'Yr', cells))
					minVel <- min(c(minVel, x[ , metric]))
					maxVel <- max(c(maxVel, x[ , metric]))
					
					assign(paste0('bv', interval, 'Yr', cells), x)
					
				}
				
				minVel <- min(0, minVel)

				for (onlyInSharedCells in c(FALSE, TRUE)) {
				
					cells <- if (onlyInSharedCells) { 'allCells'} else { 'sharedCells' }
					cells <- capIt(cells)
					x <- get(paste0('bv', interval, 'Yr', cells))

					# plot
					main <- paste0(
						metricNice,
						if (onlyInSharedCells) { ' | Shared Cells' } else { ' | All Cells' },
						' | ', interval, '-yr Time Steps'
					)

					plot(0, 0, col='white', xlim=c(-21000, 0), ylim=c(minVel, maxVel), xlab='YBP', ylab='Velocity (m / y)', main=main)
					
					for (i in seq(-21000, 0, by=500)) {
						lines(c(i, i), c(minVel, maxVel), lwd=0.2, col='gray')
					}
				
					for (gcm in gcms) {
						for (ext in exts) {
							for (algo in algos) {
								
								thisVel <- x[x$gcm == gcm & x$ext == ext & x$algo == algo & x$onlyInSharedCells == onlyInSharedCells, ]
								modelGroup <- paste0(gcm, '_', ext, 'kmExtent_', algo)
								group <- clustMembership[[modelGroup]]
								
								col <- clustCols[group]
								
								years <- rowMeans(thisVel[ , c('timeFrom', 'timeTo')])
								lines(years, thisVel[ , metric], col=col, lwd=2)
								
							}
						}
					}
					
					# legend('topleft', inset=-0.01, legend=paste0('Group', 1:numModelClusts), col=clustCols[1:numModelClusts], lwd=1, cex=0.5, box.col='white')
					
				} # next only in shared cells
				
			} # next interval
			
			title(sub=date(), outer=TRUE, line=-1)
			title(main=metricNice, outer=TRUE, line=0.8, cex.main=2.9)
			
		dev.off()
		
	} # next metric
