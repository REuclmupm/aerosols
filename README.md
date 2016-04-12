# Aerosols

Utilities for working with aerosols data in the context of solar radiation

## MACC

    lat <- 37
    lon <- -2
    elev <- 500
    loc <- list(lon = lon, 
	            lat = lat, 
				elev = elev)
				
    vals <- readMACC('MACC_2003_jan.nc', loc)

