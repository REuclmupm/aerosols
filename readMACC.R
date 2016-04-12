readMACC <- function(file, loc,
                     vars = c('tcwv','gtco3',
                         'aod469','aod550','aod670',
                         'aod865','aod1240')){
    ## Location definition
    lon <- loc$lon
    lat <- loc$lat
    elev <- loc$elev
    ## lon ranges from 0 to 360: correction for convention of lon<0
    if (lon < 0) lon <- lon + 360
    myPoint <- SpatialPoints(cbind(lon, lat),
                             proj4string = CRS('+proj=longlat +datum=WGS84 '))
    ## Read NetCDF file and extract values at the location
    readAndExtract <- function(varname){
        b <- brick(file, varname = varname)
        tt <- getZ(b)
        val <- extract(b, myPoint)
        z <- zoo(t(val), order = as.POSIXct(tt))
        names(z) <- varname
        z
    }
    ## Use this function with each variable
    vals <- lapply(vars, readAndExtract)
    ## Combine results in a multivariate zoo
    vals <- do.call(cbind, vals)
    
    ## Precipitable water vapor in kg/m2
    ## Conversion to atm-cm, divided by water density 1 g/cc
    vals$pwc <- vals$tcwv/10
    ## Ozone in Kg/m2
    ## Conversion to atm-cm (1 atm-cm=1000 Dobson Units)
    vals$uo <- vals$gtco3/2.1414e-2
    
    ## Data for linear fit
    aodLogs <- log(vals[, c("aod469", "aod550", "aod670", "aod865")])
    wlLogs <- log(c(0.469, 0.550, 0.670, 0.865))
    
    ## Linear fit to Angstrom law
    lmCoefs <- function(x){
        if (any(is.na(x))) coefs <- rep(NA, 2)
        else {
            lmLogs <- lm(x ~ wlLogs)
            coefs <- coefficients(lmLogs)
        }
    }
    
    coefs <- apply(aodLogs, 1, lmCoefs)
    alpha <- -coefs[2,]
    beta <- exp(coefs[1,])
    alpha2 <- with(vals, -(log(aod1240/aod865))/0.36013)
    
    ## aod 700
    aod700 <- beta*0.7^(-alpha)
    ## aod 380
    aod380 <- beta*0.38^(-alpha)
    ## aod 500
    aod500 <- beta*0.5^(-alpha)
    
    ## Parameters for TL
    p0 <- exp(-elev/8434.5)
    p0in <- 1/p0
    resto <- 2.0 + 0.54*p0in - 0.5*p0in^2 + 0.16*p0in^3
    
    ## Linke Turbidity
    TL <- with(vals, 3.91 * exp(0.689*p0in) * aod550 + 0.376*log(pwc) + resto)
    
    ## definition of zoo objects with all variables
    z <- cbind(vals[, c('pwc', 'uo', 'aod550')],
               aod380, aod500, aod700,
               beta, alpha, alpha2, TL)
}

