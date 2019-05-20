require(data.table)

PSANOVA.plots <- function(x, which = c("components","perspective","global"), main,...) {
	
	persp.function <- function(z, zlim, x1, x2, ...) {
		nrz <- nrow(z)
		ncz <- ncol(z)
		jet.colors <- colorRampPalette( c("blue", "green") )
		# Generate the desired number of colors from this palette
		nbcol <- 100
		#color <- jet.colors(nbcol)
		color <- topo.colors(nbcol)
		# Compute the z-value at the facet centres
		zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
		# Recode facet z-values into color indices
		facetcol <- cut(zfacet, nbcol)
		persp(x1, x2, z, expand = 0.5, theta = 45, phi = 30, ltheta = 120, shade = 0.75, col = color[facetcol], ...)
	}
	
	which <- match.arg(which)
	
	terms.formula <- x$terms$spatial$terms.formula
	
	xlab <- terms.formula$x.coord
	ylab <- terms.formula$y.coord
	
	x.coord <- x$data[,xlab]
	y.coord <- x$data[,ylab]
	response <- x$data[,x$model$response]
	
	columns <- sort(unique(x.coord))
	rows <- sort(unique(y.coord))
	
	xy.coord <- data.table(expand.grid(columns = columns, rows = rows))
	setkeyv(xy.coord, c("rows", "columns"))
	ONE <- rep(1, length(x.coord))
	ONE[x$data$weights == 0] <- NA
	df <- data.table(columns = x.coord, rows = y.coord, response = response, ONE = ONE)
	setkeyv(df, c("rows", "columns"))
	df <- df[xy.coord]
	df <- df[order(df$columns, df$rows),]
	
	# Grid for prediction
	p1 <- if(length(columns) > 100) 1 else 100%/%length(columns) + 1
	p2 <- if(length(rows) > 100) 1 else 100%/%length(rows) + 1
	
	grid.new <- c(length(columns)*p1, length(rows)*p2)
	
	col.p <- seq(min(x.coord), max(x.coord), l = grid.new[1])
	row.p <- seq(min(y.coord), max(y.coord), l = grid.new[2])

	B1p <- SpATS:::spline.bbase(x$terms$spatial$MM$MM1$knots, col.p, terms.formula$degree[1])
	B2p <- SpATS:::spline.bbase(x$terms$spatial$MM$MM2$knots, row.p, terms.formula$degree[2])

	X1p <- B1p%*%x$terms$spatial$MM$MM1$U.X
	X2p <- B2p%*%x$terms$spatial$MM$MM2$U.X

	Z1p <- B1p%*%x$terms$spatial$MM$MM1$U.Z
	Z2p <- B2p%*%x$terms$spatial$MM$MM2$U.Z

	Xp = X2p%x%X1p
	Xp <- Xp[,-1,drop = FALSE]
	
	if(terms.formula$type != "PSANOVA") {
		error("This function can only be used for PSANOVA decomposition")
	} else {
		smooth.comp <- names(x$dim)[attr(x$dim, "spatial") & attr(x$dim, "random")]
		
		B1pn <- SpATS:::spline.bbase(x$terms$spatial$MMn$MM1$knots, col.p, terms.formula$degree[1])
		B2pn <- SpATS:::spline.bbase(x$terms$spatial$MMn$MM2$knots, row.p, terms.formula$degree[2])
		Z1pn <- B1pn%*%x$terms$spatial$MMn$MM1$U.Z
		Z2pn <- B2pn%*%x$terms$spatial$MMn$MM2$U.Z
		
		# Coefficients associated to the spatial component
		fixed.spat.coef <- x$coeff[x$terms$spatial$fixed$pos]
		random.spat.coef <- x$coeff[x$terms$spatial$random$pos]
		
		Zp1 <- X2p[,1, drop = FALSE]%x%Z1p
		Zp2 <- Z2p%x%X1p[,1, drop = FALSE]
		Zp3 <- X2p[,-1, drop = FALSE]%x%Z1p
		Zp4 <- Z2p%x%X1p[,-1, drop = FALSE]
		Zp5 <- Z2pn%x%Z1pn
		
		Zp <- cbind(Zp1, Zp2, Zp3, Zp4, Zp5)
		
		# Parametric part
		eta0 <- matrix(Xp%*%fixed.spat.coef, nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		# Smooth function
		n.r.i <- 1
		n.r.f <- ncol(Zp1)
		eta1 <- matrix(Zp1%*%random.spat.coef[n.r.i:n.r.f], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		n.r.i <- n.r.f + 1
		n.r.f <- n.r.f + ncol(Zp2)
		eta2 <- matrix(Zp2%*%random.spat.coef[n.r.i:n.r.f], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		n.r.i <- n.r.f + 1
		n.r.f <- n.r.f + ncol(Zp3)
		eta3 <- matrix(Zp3%*%random.spat.coef[n.r.i:n.r.f], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		n.r.i <- n.r.f + 1
		n.r.f <- n.r.f + ncol(Zp4)
		eta4 <- matrix(Zp4%*%random.spat.coef[n.r.i:n.r.f], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		n.r.i <- n.r.f + 1
		n.r.f <- n.r.f + ncol(Zp5)
		eta5 <- matrix(Zp5%*%random.spat.coef[n.r.i:n.r.f], nrow = length(row.p), ncol = length(col.p), byrow = TRUE)
		
		eta <- cbind(Xp,Zp)%*%c(fixed.spat.coef, random.spat.coef)
		
		zlim <- range(eta)
		colors = topo.colors(100)		
		if (which == 'perspective') {
			# Perspective plots
			#dev.new(width = 24, height = 16)
			#op <- par(mfrow = c(2,3), oma = c(2, 1, 3, 2), mar = c(2.5, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))
			op <- par(mfrow = c(2,3))
			persp.function(z = t(eta0), zlim = zlim, x1 = col.p, x2 = row.p, main = 'Parametric', xlab = xlab, ylab = ylab, zlab = "", ticktype = "detailed", ...)
			plot(col.p, eta1[1,], main = smooth.comp[1], type = "l", xlab = xlab, ylab = "", ylim = zlim, ...)
			plot(row.p, eta2[,1], main = smooth.comp[2], type = "l", xlab = ylab, ylab = "", ylim = zlim, ...)
			persp.function(z = t(eta3), zlim = zlim, x1 = col.p, x2 = row.p, main = smooth.comp[3], xlab = xlab, ylab = ylab, zlab = "", ticktype = "detailed", ...)
			persp.function(z = t(eta4), zlim = zlim, x1 = col.p, x2 = row.p, main = smooth.comp[4], xlab = xlab, ylab = ylab, zlab = "", ticktype = "detailed", ...)
			persp.function(z = t(eta5), zlim = zlim, x1 = col.p, x2 = row.p, main = smooth.comp[5], xlab = xlab, ylab = ylab, zlab = "", ticktype = "detailed", ...)		
			par(op)
			#dev.off()
		}
		# Plot the components (image plots)
		Mf = kronecker(matrix(df$ONE, ncol = length(columns), nrow = length(rows)), matrix(1, p2, p1))
		if (which == 'components') {
			#dev.new(width = 24, height = 16)
			op <- par(mfrow = c(2,3), oma = c(2, 1, 3, 2), mar = c(2.5, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))
			fields::image.plot(col.p, row.p, t(eta0*Mf), main = 'Parametric', xlab = xlab, ylab = ylab, graphics.reset = TRUE, col = colors, zlim = zlim, ...)
			plot(col.p, eta1[1,], main = smooth.comp[1], type = "l", xlab = xlab, ylab = "", ylim = zlim, ...)
			plot(row.p, eta2[,1], main = smooth.comp[2], type = "l", xlab = ylab, ylab = "", ylim = zlim, ...)
			fields::image.plot(col.p, row.p, t(eta3*Mf), main = smooth.comp[3], xlab = xlab, ylab = ylab, graphics.reset = TRUE, col = colors, zlim = zlim, ...)
			fields::image.plot(col.p, row.p, t(eta4*Mf), main = smooth.comp[4], xlab = xlab, ylab = ylab, graphics.reset = TRUE, col = colors, zlim = zlim, ...)
			fields::image.plot(col.p, row.p, t(eta5*Mf), main = smooth.comp[5], xlab = xlab, ylab = ylab, graphics.reset = TRUE, col = colors, zlim = zlim, ...)
			par(op)
			#dev.off()
		}
		if (which == 'global') {
			# Raw data
			#dev.new(width = 16, height = 8)
			#op <- par(mfrow = c(1,2), oma = c(2, 1, 3, 2), mar = c(2.5, 4, 2.5, 2.5), mgp = c(1.7, 0.5, 0))
			#image.plot(columns, rows, t(matrix(df$response, ncol = length(columns), nrow = length(rows))), main = "Raw data", col = colors, xlab = xlab, ylab = ylab, zlim = range(df$response, na.rm = TRUE), graphics.reset = TRUE)
			# Smooth surface
		  fields::image.plot(col.p, row.p, t(matrix(eta, nrow = length(row.p), ncol = length(col.p), byrow = TRUE)*Mf), 
                 main = main, xlab = xlab, ylab = ylab, graphics.reset = TRUE, col = colors, zlim = zlim)
			#par(op)
			#dev.off()
		}
	}
}