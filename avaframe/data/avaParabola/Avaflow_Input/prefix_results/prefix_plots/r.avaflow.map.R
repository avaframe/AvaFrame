#!/usr/bin/R

##############################################################################
#
# MODULE:       r.avaflow.map.R
# AUTHOR:       Martin Mergili
#
# PURPOSE:      The simulation model for avalanche and debris flows
#               Script for the creation of map plots of simulation results
#
# COPYRIGHT:    (c) 2013 - 2022 by the author
#               (c) 2020 - 2022 by the University of Graz
#               (c) 2013 - 2021 by the BOKU University, Vienna
#               (c) 2015 - 2020 by the University of Vienna
#               (c) 1993 - 2022 by the R Development Core Team
#
#               This program is free software under the GNU General Public
#               License (>=v2).
#
##############################################################################

# Loading libraries
options('rgdal_show_exportToProj4_warnings'='none')
options(warn = -1)
library(sp, quietly = T)
library(maptools, quietly = T)
library(rgdal, quietly = T)
library(raster, quietly = T)

# Defining arguments
wkdir <- './data/avaParabola/Avaflow_Input/'  #working directory
prefix <- 'prefix_'  #prefix for file names
aflag <- 1  #control for additional output
tflag <- 0  #control for tsunami
model <- 7  #model
ymax0 <- -2997.500000  #coordinates defining map boundaries
ymin0 <- -4997.500000
xmin0 <- 997.500000
xmax0 <- 5997.500000
cellsize <- 20.000000
tint <- 10.000000
tstop <- 300.000000
impdef <- 0
depdef <- 0
thresholdsh <- 0.100000
thresholdst <- 10000.000000
thresholdsp <- 10000.000000
maxhflow <- 3.188508
maxtflow <- 1186548.500000
maxpflow <- 701102.375000
maxvflow <- 22.788925
maxtreach <- 292.928619
maxhtsun <- 0.000000
maxhentr <- -0.000000
maxhdep <- 0.000000
mapprof <- 0
hydrograph <- 0
releasemass <- 1
ninhyd <- 0
nouthyd <- 0
ntimemax <- 30
ctrlpts <- 0
ortho1 <- 'None'
ortho2 <- 'None'
ortho3 <- 'None'
phase1 <- 1  #first phase
phase2 <- 2  #second phase
phase3 <- 3  #third phase
slomo <- 1.000000  #factor for slow motion
mult <- 0  #control for multiple simulations

# Preparing vectors for looping over different sets of maps
if ( mult == 0 ) {
  jrange <- c(0)
  mstringlist <- c('hflow', 'None', 'None', 'None', 'None', 'None', 'None')
  if ( aflag == 1 ) {
      jrange <- append(jrange, 1)
      jrange <- append(jrange, 2)
      mstringlist[2] <- 'tflow'
      mstringlist[3] <- 'pflow'
  }

  if ( maxhentr > 0 || maxhdep > 0 ) {
    jrange <- append(jrange, 3)
    mstringlist[4] <- 'basechange'
  }
  if ( model == 7 && tflag == 1 ) {
    jrange <- append(jrange, 5)
    mstringlist[6] <- 'htsun'
  }
  jrange <- append(jrange, 6)
  mstringlist[7] <- 'treach'

} else {
  jrange <- c(0, 1, 2, 3)
  mstringlist <- c('iii_hflow', 'iii_tflow', 'iii_pflow', 'dii')
}
  for ( j in jrange ) {  # loop over all sets of maps to be created

    mstring <- mstringlist[j+1]
    if ( mult != 0 ) {
      d_hdisp <- 1  # absolute maximum raster value
      d_hdispn <- 0  # absolute minimum raster value
      thresholds <- 0.0
    } else if ( j == 6 ) {
      d_hdisp <- maxtreach  # absolute maximum raster value
      d_hdispn <- 0  # absolute minimum raster value
      thresholds <- 0.0
    } else if ( j == 0 ) {
      d_hdisp <- maxhflow  # absolute maximum raster value
      d_hdispn <- 0  # absolute minimum raster value
      thresholds <- thresholdsh
    } else if ( j == 1 ) {
      d_hdisp <- maxtflow  # absolute maximum raster value
      d_hdispn <- 0  # absolute minimum raster value
      thresholds <- thresholdst
    } else if ( j == 2 ) {
      d_hdisp <- maxpflow  # absolute maximum raster value
      d_hdispn <- 0  # absolute minimum raster value                
      thresholds <- thresholdsp
    } else if ( j == 6 ) {
      d_hdisp <- maxbasechange
      d_hdispn <- minbasechange
      thresholds <- thresholdsh
    } else {
      d_hdisp <- maxhtsun
      d_hdispn <- maxhtsun * -1
      thresholds <- thresholdsh
    }

    # Unit conversion and number of digits for legend according to maximum value:
    if ( j == 6 ) {
      mconv <- 1
    } else if ( max(d_hdisp, -d_hdispn) != 0 ) {
      mconv <- log10(max(d_hdisp, -d_hdispn))
    } else {
      mconv <- 3
    }

    if ( mconv < 3 ) {
        mconv2 <- '1'
        if ( mconv < 0 ) {
            mdig <- 3
        } else if ( mconv < 1 ) {
            mdig <- 2
        } else if ( mconv < 2 ) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    } else if ( mconv < 6 ) {
        mconv2 <- '0.001'
        if ( mconv < 4 ) {
            mdig <- 2
        } else if ( mconv < 5 ) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    } else if ( mconv < 9 ) {
        mconv2 <- '0.000001'
        if (mconv < 7) {
            mdig <- 2
        } else if (mconv < 8) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    } else if ( mconv < 12 ) {
        mconv2 <- '0.000000001'
        if ( mconv < 10 ) {
            mdig <- 2
        } else if ( mconv < 11 ) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    } else if ( mconv < 15 ) {
        mconv2 <- '0.000000000001'
        if ( mconv < 13 ) {
            mdig <- 2
        } else if ( mconv < 14 ) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    } else {
        mconv2 <- '0.000000000000001'
        if ( mconv < 16 ) {
            mdig <- 2
        } else if ( mconv < 17 ) {
            mdig <- 1
        } else {
            mdig <- 0
        }
    }
    d_hdisp <- d_hdisp * as.numeric(mconv2)
  
    if ( mult != 0 ) {
        thrs5 <- 1.00
        thrs4 <- 0.80
        thrs3 <- 0.60
        thrs2 <- 0.40
        thrs1 <- 0.20
        thrsm <- 0.01

    } else if ( model <= 3 && j == 3 ) {  # for entrainment and deposition map

        d_hdispn <- d_hdispn * as.numeric(mconv2)
        jimp <- 0  # index of impact parameter
        thrs5 <- 0
        thrs4 <- round( d_hdisp, digits <- mdig) # breaks for raster maps
        thrs3 <- round( d_hdisp * 10 / 20, digits <- mdig )
        thrs2 <- round( d_hdisp * 6 / 20, digits <- mdig )
        thrs1 <- round( d_hdisp * 3 / 20, digits <- mdig )
        if (d_hdisp == 0) {
            thrsm <- 0
        } else {
            thrsm <- min( as.numeric(thresholds) * float(mconv2), as.numeric(thrs1) * 0.75) # minimum value displayed
        }

        thrs5n <- 0
        thrs4n <- round( d_hdispn, digits <- mdig )
        thrs3n <- round( d_hdispn * 10 / 20, digits <- mdig )
        thrs2n <- round( d_hdispn * 6 / 20, digits <- mdig )
        thrs1n <- round( d_hdispn * 3 / 20, digits <- mdig )
        if ( d_hdispn == 0 ) {
            thrsmn <- 0
        } else {
            thrsmn <- max(-1 * as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1n) * 0.75) # minimum value displayed
        }

    } else if ( model == 7 && j == 5 ) {  # for tsunami map
        d_hdispn <- d_hdispn * as.numeric(mconv2)
        thrs5 <- 0
        thrs4 <- round( d_hdisp, digits <- mdig ) # breaks for raster maps
        thrs3 <- round( d_hdisp * 15 / 20, digits <- mdig )
        thrs2 <- round( d_hdisp * 10 / 20, digits <- mdig )
        thrs1 <- round( d_hdisp * 5 / 20, digits <- mdig )
        thrsm <- round( d_hdisp * 1 / 25, digits <- mdig )

        thrs5n <- 0
        thrs4n <- round( d_hdispn, digits <- mdig )
        thrs3n <- round( d_hdispn * 15 / 20, digits <- mdig )
        thrs2n <- round( d_hdispn * 10 / 20, digits <- mdig )
        thrs1n <- round( d_hdispn * 5 / 20, digits <- mdig )
        thrsmn <- round( d_hdispn * 1 / 25, digits <- mdig )

    } else if ( j == 6 ) {  # for time of reach map
        jimp <- j
        thrs5 <- round( d_hdisp, digits <- mdig ) # breaks for raster maps
        thrs4 <- round( d_hdisp * 10 / 20, digits <- mdig )
        thrs3 <- round( d_hdisp * 7 / 20, digits <- mdig )
        thrs2 <- round( d_hdisp * 4 / 20, digits <- mdig )
        thrs1 <- round( d_hdisp * 2 / 20, digits <- mdig )
        thrsm <- 0

        thrs5n <- 0
        thrs4n <- 0
        thrs3n <- 0
        thrs2n <- 0
        thrs1n <- 0
        thrsmn <- 0

    } else {  # for all other maps
        jimp <- j
        thrs5 <- round( d_hdisp, digits <- mdig )  # breaks
        thrs4 <- round( d_hdisp * 12 / 20, digits <- mdig )
        thrs3 <- round( d_hdisp * 8 / 20, digits <- mdig )
        thrs2 <- round( d_hdisp * 5 / 20, digits <- mdig )
        thrs1 <- round( d_hdisp * 2 / 20, digits <- mdig )
        if ( d_hdisp == 0) {
            thrsm <- 0
        } else {
            thrsm <- min(as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1) * 0.75) # minimum value displayed
        }

        if (model == 7 && j == 3) {

            d_hdispn <- d_hdispn * as.numeric(mconv2)
            thrs5n <- round( d_hdispn, mdig )  # breaks (negative values)
            thrs4n <- round( d_hdispn * 12 / 20, digits <- mdig )
            thrs3n <- round( d_hdispn * 8 / 20, digits <- mdig )
            thrs2n <- round( d_hdispn * 5 / 20, digits <- mdig )
            thrs1n <- round( d_hdispn * 2 / 20, digits <- mdig )
            if (d_hdisp == 0) {
                thrsmn <- 0
            } else {
                thrsmn <- max(-as.numeric(thresholds) * as.numeric(mconv2), as.numeric(thrs1n) * 0.75) # minimum value displayed
            }

        } else {

            thrs5n <- 0
            thrs4n <- 0
            thrs3n <- 0
            thrs2n <- 0
            thrs1n <- 0
            thrsmn <- 0
        }
    }

    if (j == 0 && impdef == 1 ) {

        impactareaname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'impactarea.asc', sep = '')
        impactarea0 <- raster(impactareaname)
        impactarea <- t(as.matrix(impactarea0))
        impactarea <- impactarea[, ncol(impactarea):1]
    }

    if ( j == 0 && depdef == 1 ) {

        hdepositname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hdeposit.asc', sep = '')
        hdeposit0 <- raster(hdepositname)
        hdeposit0 <- reclassify(hdeposit0,c(-Inf,as.numeric(thresholds),0,as.numeric(thresholds),Inf,1))
        hdeposit <- t(as.matrix(hdeposit0))
        hdeposit <- hdeposit[, ncol(hdeposit):1]
    }

  if ( j == 6) {
      ntimemin <- ntimemax+1
  } else {
      ntimemin <- 0
  }  
  if ( mult != 0 ) {
      ntimemax <- 0
  }

  for ( ntimesteps in ntimemin:( ntimemax+1 )) { # loop over all time steps plus one for maps of maximum values

    if ( mult != 0 ) {
        fill <- 'iscore'
    } else if ( ntimesteps < 10 ) {
        fill <- paste('000', as.character(ntimesteps), sep='')  # formatting time step string
    } else if ( ntimesteps < 100 ) {
        fill <- paste('00', as.character(ntimesteps), sep='')
    } else if ( ntimesteps < 1000 ) {
        fill <- paste('0', as.character(ntimesteps), sep='')
    } else {
        fill <- as.character(ntimesteps)
    }  

    if ( mult != 0 ) {  # strings for names of maps:
        mstringt <- paste(prefix, mstring, sep='')
    } else if ( ntimesteps <= ntimemax && model == 7 && j == 5 ) {
        mstringt <- paste(prefix, 'htsun', fill, sep='')
    } else if ( ntimesteps <= ntimemax && j != 6 ) {
        mstringt <- paste( prefix, mstring, fill, sep='')
        mstrings <- paste( prefix, mstring, '1', fill, sep='')
        if ( model == 7 ) {
            mstringf <- paste( prefix, mstring, '2', fill, sep='')
            mstringw <- paste( prefix, mstring, '3', fill, sep='')
        }
    } else if ( j < 3 ) {
        mstringt <- paste( prefix, mstring, '_max', sep='')
        mstrings <- paste( prefix, mstring, '1_max', sep='')
        if ( model == 7 ) {
            mstringf <- paste( prefix, mstring, '2_max', sep='')
            mstringw <- paste( prefix, mstring, '3_max', sep='')
        }
    } else if ( j == 3 ) {
        mstringt <- paste( prefix, mstring, '_fin', sep='')
        mstrings <- paste( prefix, mstring, '1_fin', sep='')
        if ( model == 7 ) {
            mstringf <- paste( prefix, mstring, '2_fin', sep='')
            mstringw <- paste( prefix, mstring, '3_fin', sep='')
        }
    } else if ( model == 7 && j == 5 ) {
        mstringt <- paste( prefix, 'htsun_max', sep='')
    } else if ( ntimesteps == ntimemax + 1 && j == 6 ) {
        mstringt <- paste( prefix, 'treach', sep='')
    }

    # Reading flow velocity and direction data from file defining control variable:
    if (ntimesteps > 0 && ntimesteps <= ntimemax) {
        velocity <- 1
    } else {
        velocity <- 0
    }
    if (velocity == 1) {

        dirline <- ntimesteps - 1

        intable = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions1.txt', 
            sep = '')
        dirx <- scan(intable, nlines = 1, quiet = TRUE)  #(PHASE 1) x coordinate
        diry <- scan(intable, skip = 1, nlines = 1, quiet = TRUE)  #(PHASE 1) y coordinate

        if (model <= 3) {

            dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow height
            dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow velocity in x direction
            dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #flow velocity in y direction

        } else if (model == 7) {

            intable2 = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions2.txt', 
                sep = '')
            dirxf <- scan(intable2, nlines = 1, quiet = TRUE)  #PHASE 2 x coordinate
            diryf <- scan(intable2, skip = 1, nlines = 1, quiet = TRUE)  #PHASE 2 y coordinate

            dirh <- scan(intable, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow height
            dirvx <- scan(intable, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow velocity in x direction
            dirvy <- scan(intable, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 1 flow velocity in y direction
            dirhf <- scan(intable2, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow height
            dirvxf <- scan(intable2, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow velocity in x direction
            dirvyf <- scan(intable2, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 2 flow velocity in y direction

            dirhf[dirhf > 0] <- 1

            intable3 = paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'directions3.txt', 
                sep = '')
            dirxw <- scan(intable3, nlines = 1, quiet = TRUE)  #PHASE 3 x coordinate
            diryw <- scan(intable3, skip = 1, nlines = 1, quiet = TRUE)  #PHASE 3 y coordinate

            dirhw <- scan(intable3, skip = 2 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow height
            dirvxw <- scan(intable3, skip = 3 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow velocity in x direction
            dirvyw <- scan(intable3, skip = 4 + 3 * dirline, nlines = 1, quiet = TRUE)  #PHASE 3 flow velocity in y direction

            dirhw[dirhw > 0] <- 1
        }

        dirh[dirh > 0] <- 1
    }

    # Computing map extent
    xmin <- xmin0
    xmax <- xmax0
    ymin <- ymin0
    ymax <- ymax0
    xdiff <- xmax - xmin  #extent in x direction
    ydiff <- ymax - ymin  #extent in y direction

    # if ( ydiff > xdiff ) { #ensuring minimum extent in x direction

    # xmin<-xmin-(ydiff-xdiff)/2 xmax<-xmax+(ydiff-xdiff)/2 xdiff<-xmax-xmin

    # } else

    if (xdiff > 1.6 * ydiff) {
        # ensuring minimum extent in y direction

        ymin <- ymin - (xdiff/1.6 - ydiff)/2
        ymax <- ymax + (xdiff/1.6 - ydiff)/2
        ydiff <- ymax - ymin
    }

    vunit <- 0.05 * min(xdiff, ydiff)/maxvflow

    # Building map geometry
    asprat <- xdiff/(1.4 * ydiff)  #x/y ratio
    margx <- 6  #margin in x direction
    margy <- 2.5  #margin in y direction
    dispx <- 15  #total width of output image
    dispy <- (dispx - margx)/asprat + margy  #total height of output image

    # Adapting boundaries for vector layers
    xsmin <- xmin + 0.035 * xdiff
    xsmax <- xmax - 0.035 * xdiff
    ysmin <- ymin + 0.035 * ydiff
    ysmax <- ymax - 0.035 * ydiff

    # Defining graphic parameters for raster plot
    if (model <= 3 || fill == 'iscore' || mstring == 'htsun' || mstring == 'treach') {

        if (mstring == 'hflow' || mstring == 'treach' || mstring == 'iii_hflow' || 
            mstring == 'dii') {

            if (fill == 'iscore') {
                munit <- ''
                if (mstring == 'iii_hflow') 
                    mlabel <- 'Impact indicator index'
                if (mstring == 'dii') 
                    mlabel <- 'Deposition indicator index'
                rcol1 <- rgb(0, 0, 0.8, 0.4)
                rcol7 <- rgb(0.15, 0.1, 0.6, 0.5)
                rcol13 <- rgb(0.3, 0.2, 0.4, 0.6)
                rcol19 <- rgb(0.45, 0.3, 0.2, 0.7)
                rcol25 <- rgb(0.6, 0.4, 0, 0.8)
            } else if (mstring == 'hflow') {
                if (mconv2 == '1') 
                    munit <- 'm'
                if (mconv2 == '0.001') 
                    munit <- 'km'
                mlabel <- 'Flow height'
                rcol1 <- rgb(0.5, 0.5, 0, 0.2)
                rcol7 <- rgb(0.75, 0.25, 0, 0.4)
                rcol13 <- rgb(1, 0, 0, 0.6)
                rcol19 <- rgb(0.8, 0, 0.2, 0.8)
                rcol25 <- rgb(0.6, 0, 0.4, 1)
            } else if (mstring == 'treach') {
                munit <- 's'
                mlabel <- 'Time of reach'
                rcol1 <- rgb(0.7, 0.3, 0.1, 1)
                rcol7 <- rgb(0.55, 0.275, 0.2, 0.8)
                rcol13 <- rgb(0.4, 0.25, 0.3, 0.6)
                rcol19 <- rgb(0.25, 0.225, 0.4, 0.4)
                rcol25 <- rgb(0.1, 0.2, 0.5, 0.2)
            }
        } else if (mstring == 'tflow' || mstring == 'iii_tflow') {

            if (fill != 'iscore') {
                if (mconv2 == '1') 
                    munit <- 'J'
                if (mconv2 == '0.001') 
                    munit <- 'kJ'
                if (mconv2 == '0.000001') 
                    munit <- 'MJ'
                if (mconv2 == '0.000000001') 
                    munit <- 'GJ'
                if (mconv2 == '0.000000000001') 
                    munit <- 'TJ'
                mlabel <- 'Flow kinetic energy'
                rcol1 <- rgb(0, 0.5, 0.5, 0.2)
                rcol7 <- rgb(0.15, 0.45, 0.4, 0.4)
                rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)
                rcol19 <- rgb(0.45, 0.35, 0.2, 0.8)
                rcol25 <- rgb(0.6, 0.3, 0.1, 1)
            } else {
                munit <- ''
                mlabel <- 'Impact indicator index'
                rcol1 <- rgb(0, 0.5, 0.5, 0.4)
                rcol7 <- rgb(0.15, 0.45, 0.4, 0.5)
                rcol13 <- rgb(0.3, 0.4, 0.3, 0.6)
                rcol19 <- rgb(0.45, 0.35, 0.2, 0.7)
                rcol25 <- rgb(0.6, 0.3, 0.1, 0.8)
            }
        } else if (mstring == 'pflow' || mstring == 'iii_pflow') {

            if (fill != 'iscore') {
                if (mconv2 == '1') 
                    munit <- 'Pa'
                if (mconv2 == '0.001') 
                    munit <- 'kPa'
                if (mconv2 == '0.000001') 
                    munit <- 'MPa'
                if (mconv2 == '0.000000001') 
                    munit <- 'GPa'
                if (mconv2 == '0.000000000001') 
                    munit <- 'TPa'
                mlabel <- 'Flow pressure'
                rcol1 <- rgb(0.2, 0, 0.6, 0.2)
                rcol7 <- rgb(0.4, 0.3, 0.35, 0.4)
                rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)
                rcol19 <- rgb(0.5, 0.45, 0.15, 0.8)
                rcol25 <- rgb(0.6, 0.6, 0, 1)
            } else {
                munit <- ''
                mlabel <- 'Impact indicator index'
                rcol1 <- rgb(0.2, 0, 0.6, 0.4)
                rcol7 <- rgb(0.4, 0.3, 0.35, 0.5)
                rcol13 <- rgb(0.4, 0.3, 0.3, 0.6)
                rcol19 <- rgb(0.5, 0.45, 0.15, 0.7)
                rcol25 <- rgb(0.6, 0.6, 0, 0.8)
            }
        } else if (mstring == 'basechange') {

            if (mconv2 == '1') 
                munit <- 'm'
            if (mconv2 == '0.001') 
                munit <- 'km'
            mlabel <- 'Change of basal topography'

            rcol3 <- rgb(0.1, 0.2, 0.6, 1)
            rcol7 <- rgb(0.1, 0.2, 0.6, 0.75)
            rcol11 <- rgb(0.1, 0.2, 0.6, 0.5)
            rcol15 <- rgb(0.1, 0.2, 0.6, 0.25)
            rcol17 <- rgb(1, 1, 1, 0)
            rcol20 <- rgb(0.2, 0.6, 0.1, 0.25)
            rcol24 <- rgb(0.2, 0.6, 0.1, 0.5)
            rcol28 <- rgb(0.2, 0.6, 0.1, 0.75)
            rcol32 <- rgb(0.2, 0.6, 0.1, 1)

        } else if (mstring == 'htsun') {

            if (mconv2 == '1') 
                munit <- 'm'
            if (mconv2 == '0.001') 
                munit <- 'km'
            mlabel <- 'Tsunami height'

            rcol3 <- rgb(0, 0.6, 0.1, 1)
            rcol7 <- rgb(0, 0.6, 0.1, 0.75)
            rcol11 <- rgb(0, 0.6, 0.1, 0.5)
            rcol15 <- rgb(0, 0.6, 0.1, 0.25)
            rcol17 <- rgb(1, 1, 1, 0)
            rcol20 <- rgb(0, 0.1, 0.6, 0.25)
            rcol24 <- rgb(0, 0.1, 0.6, 0.5)
            rcol28 <- rgb(0, 0.1, 0.6, 0.75)
            rcol32 <- rgb(0, 0.1, 0.6, 1)
        }

        if (mstring == 'basechange' || mstring == 'htsun') {

            ctext = vector('expression', 10)  #vector for colour bar labels

            if (mstring == 'basechange') {
                ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), 
                    mdig), nsmall = mdig), fmunit = format(munit)))[2]
            } else {
                ctext[1] <- ''
            }
            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsmn)), 
                fmunit = format(munit)))[2]
            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
                fmunit = format(munit)))[2]
            ctext[7] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[8] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[9] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[10] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]

            rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5)
            rlabels <- ctext
            rcolours <- c(rcol3, rcol7, rcol11, rcol15, rcol17, rcol20, rcol24, rcol28, 
                rcol32)

        } else {

            ctext = vector('expression', 6)  #vector for colour bar labels

            ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
                fmunit = format(munit)))[2]
            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]

            rthresholds <- c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5)
            rlabels <- ctext
            rcolours <- c(rcol1, rcol7, rcol13, rcol19, rcol25)
        }

    } else if (model == 7) {

        if (mstring == 'hflow') {

            if (mconv2 == '1') 
                munit <- 'm'
            if (mconv2 == '0.001') 
                munit <- 'km'
            mlabel <- 'Flow height'

            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
            rthresholds <- rep(0:625) + 0.5

        } else if (mstring == 'tflow') {

            if (mconv2 == '1') 
                munit <- 'J'
            if (mconv2 == '0.001') 
                munit <- 'kJ'
            if (mconv2 == '0.000001') 
                munit <- 'MJ'
            if (mconv2 == '0.000000001') 
                munit <- 'GJ'
            if (mconv2 == '0.000000000001') 
                munit <- 'TJ'
            mlabel <- 'Flow kinetic
energy'

            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
            rthresholds <- rep(0:625) + 0.5

        } else if (mstring == 'pflow') {

            if (mconv2 == '1') 
                munit <- 'Pa'
            if (mconv2 == '0.001') 
                munit <- 'kPa'
            if (mconv2 == '0.000001') 
                munit <- 'MPa'
            if (mconv2 == '0.000000001') 
                munit <- 'GPa'
            if (mconv2 == '0.000000000001') 
                munit <- 'TPa'
            mlabel <- 'Flow pressure'

            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
            rthresholds <- rep(0:625) + 0.5

        } else if (mstring == 'basechange') {

            if (mconv2 == '1') 
                munit <- 'm'
            if (mconv2 == '0.001') 
                munit <- 'km'
            mlabel <- 'Change of
topography'

            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
            rthresholds <- rep(0:625) + 0.5

        } else if (mstring == 'htsun') {

            if (mconv2 == '1') 
                munit <- 'm'
            if (mconv2 == '0.001') 
                munit <- 'km'
            mlabel <- 'Height of tsunami'

            rcolours <- rgb(rep(0:4, 25, each = 25) * 0.25, rep(0:4, 125, each = 5) * 
                0.25, rep(0:4, 625, each = 1) * 0.25, rep(1:5, 5, each = 125) * 0.2)
            rthresholds <- rep(0:625) + 0.5
        }

        if (mstring == 'basechange') {

            ctext = vector('expression', 6)  #vector for contour line legend labels

            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
                fmunit = format(munit)))[2]
            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]

            ctext_neg = vector('expression', 6)  #vector for negative contour line legend labels

            ctext_neg[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrs5n)), 
                fmunit = format(munit)))[2]
            ctext_neg[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext_neg[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext_neg[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext_neg[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1n), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext_neg[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrsmn), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]

        } else {

            ctext = vector('expression', 6)  #vector for contour line legend labels

            ctext[6] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(as.numeric(thrsm)), 
                fmunit = format(munit)))[2]
            ctext[5] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs1), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[4] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs2), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[3] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs3), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[2] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs4), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
            ctext[1] <- substitute(expression(fthrs ~ fmunit), list(fthrs = format(round(as.numeric(thrs5), 
                mdig), nsmall = mdig), fmunit = format(munit)))[2]
        }
    }

    # Importing raster layers
    cat(paste('Plotting ', mstringt, ' ...', sep=''))
    mstringtname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringt, '.asc', sep = '')
    mstringtr <- raster(mstringtname)  #raster of total value

    if (ortho1 != 'None') {
        orthophoto <- stack(ortho1, ortho2, ortho3)  #orthophoto
    } else if ( mult != 0 ) {
        hillshadename <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hillshade0000.asc', sep = '')
        hillshade <- raster(hillshadename)  #raster of hillshade
    } else if ( ntimesteps <= ntimemax ) {
        hillshadename <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hillshade', fill, '.asc', sep = '')
        hillshade <- raster(hillshadename)  #raster of hillshade
    }

    # Computing vectors
    if ( mult == 0 ) { # for single model run
        mstringrname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', prefix, 'hflow0000.asc', sep = '')
        startshape0 <- raster(mstringrname)  #raster of release area
        startshape <- t(as.matrix(startshape0)) # raster map for contour creation
        startshape <- startshape[, ncol(startshape):1]
    }

    # Computing composite raster and defining contour line parameters
    if ( mult != 0 ) { # for multiple model runs

      a2 <- reclassify(mstringtr, c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs5,5))

      a2c <- t(as.matrix(a2)) # raster map for contour creation
      a2c <- a2c[, ncol(a2c):1]

      contstep = 1  # contour interval
      contmin = 1  # minimum contour
      contmax = 5  # maximum contour  

    } else if ( model == 7 && j == 5 ) { # for tsunami map

      a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrs3n,1,thrs3n,thrs2n,2,thrs2n,thrs1n,3,thrs1n,thrsmn,4,thrsmn,thrsm,5,
          thrsm,thrs1,6,thrs1,thrs2,7,thrs2,thrs3,8,thrs3,Inf,9))

      a2c <- t(as.matrix(a2)) # raster map for contour creation
      a2c <- a2c[, ncol(a2c):1]

      contstep = 1  # contour interval
      contmin = 0  # minimum contour
      contmax = 9  # maximum contour  

    } else if ( model <= 3 || (ntimesteps == ntimemax && j == 6 )) {

      if ( j == 3 ) { # for entrainment and deposition map:

        a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrs3n,1,thrs3n,thrs2n,2,thrs2n,thrs1n,3,thrs1n,thrsmn,4,thrsmn,thrsm,5,
            thrsm,thrs1,6,thrs1,thrs2,7,thrs2,thrs3,8,thrs3,Inf,9))
    
        a2c <- t(as.matrix(a2)) # raster map for contour creation
        a2c <- a2c[, ncol(a2c):1]

        contstep = 1  # contour interval
        contmin = 0  # minimum contour
        contmax = 9  # maximum contour
  
      } else { # for all other maps:

        a2 <- reclassify(mstringtr * as.numeric(mconv2), c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,Inf,5))
    
        a2c <- t(as.matrix(a2)) # raster map for contour creation
        a2c <- a2c[, ncol(a2c):1]

        contstep = 1  # contour interval
        contmin = 1  # minimum contour
        contmax = 5  # maximum contour
      }
    } else if ( j != 6 ) {

      mstringsname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstrings, '.asc', sep = '')
      mstringsr <- raster(mstringsname)  #raster of phase 1 value

      mstringfname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringf, '.asc', sep = '')
      mstringfr <- raster(mstringfname)  #raster of phase 2 value

      mstringwname <- paste(wkdir, prefix, 'results/', prefix, 'ascii/', mstringw, '.asc', sep = '')
      mstringwr <- raster(mstringwname)  #raster of phase 3 value

      a21 <- mstringtr * as.numeric(mconv2)
  
      a22r <- abs(mstringsr/mstringtr)
      a23r <- abs(mstringfr/mstringtr)
      a24r <- abs(mstringwr/mstringtr)
      a22 <- mask(a22r, mstringsr, maskvalue=0, updatevalue=0) # ratio of PHASE 1 component
      a23 <- mask(a23r, mstringfr, maskvalue=0, updatevalue=0) # ratio of PHASE 2 component
      a24 <- mask(a24r, mstringwr, maskvalue=0, updatevalue=0) # ratio of PHASE 3 component

      a2 <- reclassify(a21, c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,Inf,5)) # reclass for magnitude

      if ( j == 3 ) {
  
        a2n <- reclassify(a21, c(-Inf,thrsm,0,thrsm,thrs1,1,thrs1,thrs2,2,thrs2,thrs3,3,thrs3,thrs4,4,thrs4,Inf,5)) # reclass for magnitude

      } else {
  
        a2n <- 0

      }

      a2x <- max( a2, a2n)

      a2sr <- reclassify(a22, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))
      a2fr <- reclassify(a23, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))
      a2wr <- reclassify(a24, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5))

      a2s <- mask(a2sr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 1 ratio
      a2f <- mask(a2fr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 2 ratio
      a2w <- mask(a2wr, a2, maskvalue=0, updatevalue=0) # reclass for PHASE 3 ratio
      a2 <- (a2-1)*125+(a2s-1)*25+(a2f-1)*5+a2w #combined ratio
      a2c <- t(as.matrix(a2+124)) # raster map for contour creation
      a2c <- a2c[, ncol(a2c):1]

      if ( j == 3 ) {

          a2snr <- reclassify(a22, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5)) # reclass for PHASE 1 ratio
          a2fnr <- reclassify(a23, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5)) # reclass for PHASE 2 ratio
          a2wnr <- reclassify(a24, c(-Inf,0.2,1,0.2,0.4,2,0.4,0.6,3,0.6,0.8,4,0.8,Inf,5)) # reclass for PHASE 3 ratio

          a2sn <- mask(a2snr, a2n, maskvalue=0, updatevalue=0) # reclass for PHASE 1 ratio
          a2fn <- mask(a2fnr, a2n, maskvalue=0, updatevalue=0) # reclass for PHASE 2 ratio
          a2wn <- mask(a2wnr, a2n, maskvalue=0, updatevalue=0) # reclass for PHASE 3 ratio
          a2n <- (a2n-1)*125+(a2sn-1)*25+(a2fn-1)*5+a2wn # combined ratio
          a2cn <- t(as.matrix(a2n+124)) # raster map for contour creation
          a2cn <- a2cn[, ncol(a2cn):1]
          a2 <- (a2x-1)*125+(a2s+a2sn-1)*25+(a2f+a2fn-1)*5+a2w+a2wn # combined ratio
      }

      contstep = 125  # contour interval
      contmin = 125  # minimum contour
      contmax = 625  # maximum contour
    }

    rastplot <- a2

    # Creating plot file
    if (fill == 'iscore') {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '.png', 
            sep = '')
    } else if (ntimesteps <= ntimemax) {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, 'maps_timesteps/', 
            prefix, mstring, fill, '.png', sep = '')
    } else if (mstring == 'basechange') {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_fin.png', 
            sep = '')
    } else if (mstring == 'htsun') {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_max.png', 
            sep = '')
    } else if (mstring == 'treach') {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '.png', 
            sep = '')
    } else {
        mapplot = paste(wkdir, prefix, 'results/', prefix, 'plots/', prefix, mstring, '_max.png', 
            sep = '')
    }
    png(filename = mapplot, width = dispx, height = dispy, units = 'cm', res = 300)

    # Defining margins
    par(mar = c(2, 1.3, 0.5, 2.5))
    par(oma = c(0, 0, 0, 1.5))
    clip(xmin, xmax, ymin, ymax)  #constraining drawing area

    if (ydiff > xdiff * 1.2) {
        # shrink factor and scaling of labels for colour bar legend
        lshrink <- 0.7
        lcex <- 1
    } else if (ydiff > xdiff/1.2) {
        lshrink <- 0.85
        lcex <- 1
    } else {
        lshrink <- 1
        lcex <- 0.8
    }

    # Plotting raster layers
    if (ortho1 != 'None') {
        plot(rastplot, legend.width = 1, legend = FALSE, col = rgb(0.75, 0.75, 0.75, 
            1), axes = FALSE, box = FALSE, xlab = NA, ylab = NA, xlim = c(xmin, xmax), 
            ylim = c(ymin, ymax), useRaster = TRUE)  #hillshade
        plotRGB(orthophoto, add = T, r = 3, g = 2, b = 1, stretch = 'hist', alpha = 150)  #orthophoto
    } else {
        plot(hillshade, legend.width = 1, col = gray((max(0, cellStats(hillshade, 'min')):min(255, 
            cellStats(hillshade, 'max')))/255), legend = FALSE, axes = FALSE, box = FALSE, 
            xlab = NA, ylab = NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), useRaster = TRUE)  #hillshade
    }

    par(new = TRUE)
    par(cex = lcex)

    if (model < 7 || fill == 'iscore' || mstring == 'htsun' || mstring == 'treach') {
        plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, 
            xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, 
            box = FALSE, legend.args = list(text = mlabel, side = 4, line = -2.1, cex = lcex), 
            axis.args = list(labels = rlabels), legend.shrink = lshrink)  #flow parameter
    } else {
        plot(rastplot, legend.width = 1, useRaster = TRUE, col = rcolours, breaks = rthresholds, 
            xlim = c(xmin, xmax), ylim = c(ymin, ymax), axes = FALSE, xlab = NA, ylab = NA, 
            box = FALSE, legend = FALSE)  #flow parameter
    }

    # Plotting vector layers
    if (model <= 3 && releasemass > 0 && mult == 0 ) {
        contour(x = seq(xmin0, xmax0, length.out = nrow(startshape)), y = seq(ymin0, ymax0, length.out = ncol(startshape)), z = startshape, 
          levels=c(thrsm), drawlabels=FALSE, col = 'red', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), 
          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # release area
    } else if (model == 7 && mult == 0 ) {
        if (releasemass > 3) {
            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape3)), y = seq(ymin0, ymax0, length.out = ncol(startshape3)), z = startshape3, 
              levels=c(thrsm), drawlabels=FALSE, col = 'blue', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 3 release area
        }
        if (releasemass == 2 || releasemass == 3 || releasemass == 7) {
            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape2)), y = seq(ymin0, ymax0, length.out = ncol(startshape2)), z = startshape2, 
              levels=c(thrsm), drawlabels=FALSE, col = 'green', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 2 release area
        }
        if (releasemass == 1 || releasemass == 3 || releasemass == 5 || releasemass == 7) {
            contour(x = seq(xmin0, xmax0, length.out = nrow(startshape)), y = seq(ymin0, ymax0, length.out = ncol(startshape)), z = startshape, 
              levels=c(thrsm), drawlabels=FALSE, col = 'red', lty = 3, lwd = 1.5, xlim = c(xmin, xmax), 
              ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) # PHASE 1 release area
        }
    }

    if (j == 0 && impdef == 1 && (fill != 'iscore' || mlabel == 'Impact indicator index')) {

        contour(x = seq(xmin0, xmax0, length.out = nrow(impactarea)), y = seq(ymin0, ymax0, length.out = ncol(impactarea)), z = impactarea, 
          levels=c(1), drawlabels=FALSE, col = 'red', lty = 1, lwd = 1.5, xlim = c(xmin, xmax), 
          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #observed impact area
    }

    if ( j == 0 && depdef == 1 && (fill != 'iscore' || mlabel == 'Deposition indicator index')) {

        contour(x = seq(xmin0, xmax0, length.out = nrow(hdeposit)), y = seq(ymin0, ymax0, length.out = ncol(hdeposit)), z = hdeposit, 
          levels=c(1), drawlabels=FALSE, col = 'orange', lty = 1, lwd = 1.5, xlim = c(xmin, xmax), 
          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #height of observed deposit
    }

    if (mapprof != 0) {

        par(new = TRUE)
        intable <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'mapprof.txt', sep='')
        profilex <- scan(intable, nlines = 1, quiet = TRUE)  #profile x coordinate
        profiley <- scan(intable, skip = 1, nlines = 1, quiet = TRUE)  #profile y coordinate
        lines(x = profilex, y = profiley, col = 'yellow', lwd = 1.5, lty = 2)  #profile
    }

    if (max(replace(a2c, is.na(a2c), 0)) != min(replace(a2c, is.na(a2c), 0))) {

        contour(x = seq(xmin0, xmax0, length.out = nrow(a2c)), y = seq(ymin0, ymax0, length.out = ncol(a2c)), z = a2c, 
          levels=seq(contmin, contmax, by=contstep), drawlabels=FALSE, col = 'black', lty = 1, lwd = 0.5, xlim = c(xmin, xmax), 
          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #contour lines
    }

    if (model == 7 && mstring == 'basechange' )  {
      if ( max(replace(a2cn, is.na(a2cn), 0)) != min(replace(a2cn, is.na(a2cn), 0))) {

        contour(x = seq(xmin0, xmax0, length.out = nrow(a2cn)), y = seq(ymin0, ymax0, length.out = ncol(a2cn)), z = a2cn, 
          levels=seq(contmin, contmax, by=contstep), drawlabels=FALSE, col = 'gainsboro', lty = 1, lwd = 0.5, xlim = c(xmin, xmax), 
          ylim = c(ymin, ymax), axes = F, xlab = NA, ylab = NA, add = TRUE) #negative contour lines
      }
    }

    if (velocity == 1) {

        if (model <= 3) {
            par(new = TRUE)
            arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * 
                dirvx, length = 0.025, code = 2, angle = 45, col = 'red', lwd = 0.5 * 
                dirh)  #directions
        } else if (model == 7) {
            par(new = TRUE)
            arrows(x0 = dirxw, y0 = diryw, x1 = dirxw + vunit * dirvyw, y1 = diryw - 
                vunit * dirvxw, length = 0.025, code = 2, angle = 45, col = rgb(0, 0, 
                1), lwd = 0.5 * dirhw)  #PHASE 3 directions
            par(new = TRUE)
            arrows(x0 = dirxf, y0 = diryf, x1 = dirxf + vunit * dirvyf, y1 = diryf - 
                vunit * dirvxf, length = 0.025, code = 2, angle = 45, col = rgb(0, 1, 
                0), lwd = 0.5 * dirhf)  #PHASE 2 directions
            par(new = TRUE)
            arrows(x0 = dirx, y0 = diry, x1 = dirx + vunit * dirvy, y1 = diry - vunit * 
                dirvx, length = 0.025, code = 2, angle = 45, col = rgb(1, 0, 0), lwd = 0.5 * 
                dirh)  #PHASE 1 directions
        }
    }

    # Plotting control points
    if (ctrlpts > 0) {

        par(new = TRUE)
        inctrl <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'ctrlpoints.txt', 
            sep = '')

        jz <- 1  #initializing counter for loop over all control points
        repeat {
            # starting loop over all control points
            if (jz > ctrlpts) {
                break  #break if last control point was reached
            } else {
                ctrlx <- read.table(inctrl, header = FALSE)[(jz - 1) * (ntimemax + 2) + 
                    2, 3]  #control point x coordinate
                ctrly <- read.table(inctrl, header = FALSE)[(jz - 1) * (ntimemax + 2) + 
                    2, 4]  #control point y coordinate
                points(x = as.vector(t(ctrlx)), y = as.vector(t(ctrly)), col = 'yellow', 
                    pch = 3)  #plotting control point
                ctrltext <- as.character(read.table(inctrl, header = FALSE)[(jz - 1) * 
                    (ntimemax + 2) + 2, 1])

                text(x = as.numeric(as.vector(t(ctrlx))) + 0.035 * (xmax - xmin), y = as.numeric(as.vector(t(ctrly))) + 
                    0.035 * (xmax - xmin), labels = ctrltext, cex = 1, col = 'yellow', 
                    font = 1)  #control point label

                jz <- jz + 1
            }
        }
    }

    # Plotting hydrograph profiles
    if (hydrograph > 0) {

        par(new = TRUE)
        inhyd <- paste(wkdir, prefix, 'results/', prefix, 'files/', prefix, 'hydprofiles.txt', 
            sep = '')

        jz <- 1  #initializing counter for loop over all hydrographs
        repeat {
            # starting loop over all hydrographs
            if (jz > (ninhyd + nouthyd)) {
                break  #break if last hydrograph was reached
            } else {
                if (jz <= ninhyd) {
                    # color of hydrograph profile depending on type:
                    hydcol <- 'green'
                } else {
                    hydcol <- 'purple'
                }

                hydx <- read.table(inhyd, header = FALSE)[jz + 1, 2:4]  #hydrograph x coordinates
                hydy <- read.table(inhyd, header = FALSE)[jz + 1, 5:7]  #hydrograph y coordinates
                lines(x = as.vector(t(hydx)), y = as.vector(t(hydy)), col = hydcol, lwd = 1.5, 
                    lty = 3)  #profile
                points(x = as.vector(t(hydx))[2], y = as.vector(t(hydy))[2], col = hydcol, 
                    pch = 19, cex = 0.5)  #centre point
                hydtext <- as.character(read.table(inhyd, header = FALSE)[jz + 1, 1])

                hyddx1 <- as.numeric(as.vector(t(hydx))[1])  #coordinates for hydrograph label:
                hyddx3 <- as.numeric(as.vector(t(hydx))[3])
                if (hyddx1 > hyddx3) {
                    hyddx <- hyddx1
                    hyddy <- as.numeric(as.vector(t(hydy))[1])
                } else {
                    hyddx <- hyddx3
                    hyddy <- as.numeric(as.vector(t(hydy))[3])
                }

                text(x = hyddx, y = hyddy, labels = hydtext, cex = 1, col = hydcol, font = 1, 
                    pos = 4, offset = 0.3)  #hydrograph label

                jz <- jz + 1
            }
        }
    }

    # Plotting bounding box and axes
    box()
    par(cex = 0.8)
    axis(side = 1, tck = -0.01, labels = NA)  #x axis
    axis(side = 2, tck = -0.01, labels = NA)  #y axis
    axis(side = 1, lwd = 0, line = -0.6)  #x axis
    axis(side = 2, lwd = 0, line = -0.6)  #y axis
    par(cex = 1)

    # Plotting header text
    htext = vector('expression', 1)

    #Defining string for time unit
    if ( slomo == 1 ) {
        secunit0 <- 's'
    } else if ( slomo == 60 ) {
        secunit0 <- 'min'
    } else if ( slomo == 3600 ) {
        secunit0 <- 'h'
    } else if ( slomo == 86400 ) {
        secunit0 <- 'days'
    } else if ( slomo == 604800 ) {
        secunit0 <- 'weeks'
    } else if ( slomo == 2592000 ) {
        secunit0 <- 'months'
    } else if ( slomo == 31536000 ) {
        secunit0 <- 'years'
    } else if ( slomo == 315360000 ) {
        secunit0 <- 'dec'
    }else {
        secunit0 <- paste('x', as.character(slomo), 's')
    }

    if (fill != 'iscore') {
        # defining x position and label:
        if (ntimesteps <= ntimemax) {
            secpassed <- tint * ntimesteps  #time passed since start of simulation
            if (secpassed > tstop) 
                secpassed <- tstop
            htext[1] <- substitute(expression(italic(t) == secformat ~ secunit), list(secformat = format(round(secpassed, 
                1), nsmall = 1), secunit = format(secunit0)))[2]
            posx <- xmin + xdiff/5
        } else if (mstring == 'hflow') {
            htext[1] <- 'Maximum flow height'
            posx <- xmin + xdiff/2
        } else if (mstring == 'tflow') {
            htext[1] <- 'Maximum flow kinetic energy'
            posx <- xmin + xdiff/2
        } else if (mstring == 'pflow') {
            htext[1] <- 'Maximum flow pressure'
            posx <- xmin + xdiff/2
        } else if (mstring == 'basechange') {
            htext[1] <- 'Final change of basal topography'
            posx <- xmin + xdiff/2
        } else if (mstring == 'htsun') {
            htext[1] <- 'Maximum tsunami height'
            posx <- xmin + xdiff/2
        } else if (mstring == 'treach') {
            htext[1] <- 'Time of reach'
            posx <- xmin + xdiff/2
        }
    } else if (mstring == 'iii_hflow') {
        htext[1] <- 'Impact indicator index'
        posx <- xmin + xdiff/2
    } else if (mstring == 'iii_tflow') {
        htext[1] <- 'Impact indicator index'
        posx <- xmin + xdiff/2
    } else if (mstring == 'iii_pflow') {
        htext[1] <- 'Impact indicator index'
        posx <- xmin + xdiff/2
    } else if (mstring == 'dii') {
        htext[1] <- 'Deposition indicator index'
        posx <- xmin + xdiff/2
    }

    if (ydiff > xdiff * 1.2) {
        # y positions and scaling of header text and header legend
        posy1 <- ymax + ydiff/14
        posy2 <- ymax + ydiff/10
        posy3 <- ymax + ydiff/14
        posy4 <- ymax + ydiff/26
        lcex <- 1
    } else if (ydiff > xdiff/1.2) {
        posy1 <- ymax + ydiff/12
        posy2 <- ymax + ydiff/8
        posy3 <- ymax + ydiff/13
        posy4 <- ymax + ydiff/34
        lcex <- 1
    } else {
        posy1 <- ymax + ydiff/11
        posy2 <- ymax + ydiff/7
        posy3 <- ymax + ydiff/11.5
        posy4 <- ymax + ydiff/27
        lcex <- 0.8
    }

    text(x = posx, y = posy1, labels = htext[1], cex = 1.4 * lcex, col = 'black')  #printing text

    #Defining string for velocity unit
    if ( slomo == 1 ) {
        velunit0 <- 'm/s'
    } else if ( slomo == 60 ) {
        velunit0 <- 'm/min'
    } else if ( slomo == 3600 ) {
        velunit0 <- 'm/h'
    } else if ( slomo == 86400 ) {
        velunit0 <- 'm/day'
    } else if ( slomo == 604800 ) {
        velunit0 <- 'm/week'
    } else if ( slomo == 2592000 ) {
        velunit0 <- 'm/month'
    } else if ( slomo == 31536000 ) {
        velunit0 <- 'm/year'
    } else if ( slomo == 315360000 ) {
        velunit0 <- 'm/dec'
    } else {
        velunit0 <- paste('m/', as.character(slomo), 's')
    }

    # Plotting header legend (flow velocity)
    if (fill != 'iscore' && ntimesteps <= ntimemax) {

        par(new = TRUE)
        if (model <= 3) {
            arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * 
                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
                col = 'red', lwd = 0.5)  #legend for velocity
        } else if (model > 3 && model < 7) {
            arrows(x0 = xmin + 0.62 * xdiff, y0 = posy3, x1 = xmin + 0.62 * xdiff + 0.05 * 
                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
                col = 'brown', lwd = 0.5)  #legend for PHASE 1 velocity
            arrows(x0 = xmin + 0.62 * xdiff + 0.05 * min(xdiff, ydiff), y0 = posy4, x1 = xmin + 
                0.62 * xdiff, y1 = posy4, length = 0.025, code = 2, angle = 45, col = 'blue', 
                lwd = 0.5)  #legend for PHASE 2 velocity
        } else if (model == 7) {
            arrows(x0 = xmin + 0.58 * xdiff, y0 = posy3, x1 = xmin + 0.58 * xdiff + 0.05 * 
                min(xdiff, ydiff), y1 = posy3, length = 0.025, code = 2, angle = 45, 
                col = rgb(1, 0, 0), lwd = 0.5)  #legend for PHASE 1 velocity
            arrows(x0 = xmin + 0.58 * xdiff, y0 = posy4, x1 = xmin + 0.58 * xdiff + 0.05 * 
                min(xdiff, ydiff), y1 = posy4, length = 0.025, code = 2, angle = 45, 
                col = rgb(0, 1, 0), lwd = 0.5)  #legend for PHASE 2 velocity
            arrows(x0 = xmin + 0.68 * xdiff, y0 = (posy3 + posy4)/2, x1 = xmin + 0.68 * 
                xdiff + 0.05 * min(xdiff, ydiff), y1 = (posy3 + posy4)/2, length = 0.025, 
                code = 2, angle = 45, col = rgb(0, 0, 1), lwd = 0.5)  #legend for PHASE 3 velocity
        }

        textarrow = vector('expression', 1)
        textarrow[1] = substitute(expression(italic(v)[max] == vformat ~ velunit), list(vformat = format(round(maxvflow, 
            1), nsmall = 1), velunit = format(velunit0)))[2]  #defining label for flow velocity legend
        text(x = xmax - xdiff/2.15, y = posy2, labels = textarrow[1], col = 'black', 
            pos = 4, cex = lcex)

        if (model == 7) {
            text(x = xmin + 0.61 * xdiff, y = posy3, labels = 'P1', col = rgb(1, 0, 0), 
                pos = 4, cex = lcex)  #printing legend text for PHASE 1 velocity
            text(x = xmin + 0.61 * xdiff, y = posy4, labels = 'P2', col = rgb(0, 1, 0), 
                pos = 4, cex = lcex)  #printing legend text for PHASE 2 velocity
            text(x = xmin + 0.71 * xdiff, y = (posy3 + posy4)/2, labels = 'P3', col = rgb(0, 
                0, 1), pos = 4, cex = lcex)  #printing legend text for PHASE 3 velocity
        }
    }

    # Plotting footer legend (release heights, hydrographs and observed impact area
    # and deposit) x positions:
    if (releasemass > 0) {
        if (model < 7) {
            pos1 <- 0.035
        } else {
            pos1 <- 0.01
            pos1w <- 0.125
        }
        pos2 <- 0.335
    } else {
        pos2 <- 0.035
    }
    pos3 <- 0.635

    if (ydiff > xdiff * 1.2) {
        # y positions and scaling of labels
        posy1 <- 0.05
        posy2 <- 0.025
        posy3 <- 0
        lcex <- 1
    } else if (ydiff > xdiff/1.2) {
        posy1 <- 0.066
        posy2 <- 0.026
        posy3 <- -0.014
        lcex <- 1
    } else {
        posy1 <- 0.071
        posy2 <- 0.029
        posy3 <- -0.013
        lcey <- 0.8
    }

    if (releasemass > 0) {
        if (model <= 3) {
            legend('bottomleft', lwd = NA, legend = 'Release area', text.col = 'black', 
                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.1, posy1), 
                cex = lcex)
            legend('bottomleft', lty = 3, lwd = 1.5, col = 'red', legend = '', text.col = 'black', 
                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, posy2), cex = lcex)  #legend for release area
        } else if (model == 7) {
            legend('bottomleft', lwd = NA, legend = 'Release areas', text.col = 'black', 
                bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1 - 0.075, posy1), 
                cex = lcex)
            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(1, 0, 0), legend = 'P1', 
                text.col = 'red', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                    posy2), cex = lcex)  #legend for PHASE 1 release area
            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(0, 1, 0), legend = 'P2', 
                text.col = 'green', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1w, 
                    (posy2 + posy3)/2), cex = lcex)  #legend for PHASE 2 release area
            legend('bottomleft', lty = 3, lwd = 1.5, col = rgb(0, 0, 1), legend = 'P3', 
                text.col = 'blue', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos1, 
                    posy3), cex = lcex)  #legend for PHASE 3 release area
        }
    }

    if (model == 7 && mstring == 'basechange' && fill != 'iscore') {
        par(xpd = TRUE)

        legend(x = xmax + xdiff/50, y = ymax - ydiff/25, xjust = 0, title = mlabel, legend = NA, 
            bty = 'n', horiz = FALSE, cex = lcex)  #legend title

        legend(x = xmax, y = ymax - ydiff/10, xjust = 0, legend = c(ctext), lty = c(1, 
            1, 1, 1, 1, 1), lwd = 0.5, col = c('white', 'black', 'black', 'black', 'black', 
            'black'), text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #main legend

        legend(x = xmax, y = ymax - ydiff/2, xjust = 0, legend = c(ctext_neg), lty = c(1, 
            1, 1, 1, 1, 1), lwd = 0.5, col = c('gainsboro', 'gainsboro', 'gainsboro', 
            'gainsboro', 'gainsboro', 'white'), text.col = 'black', bty = 'n', horiz = FALSE, 
            cex = lcex)  #main legend for negative values

        legend(x = xmax, y = ymax - ydiff/1.1, xjust = 0, legend = c(paste('P1', phase1, 
            sep = ' '), paste('P2', phase2, sep = ' '), paste('P3', phase3, sep = ' ')), 
            fill = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), border = 'black', 
            text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #legend for fraction

        # legend(x=xmax+xdiff/50, y=ymax-ydiff/9, xjust=0, title='Max', legend=NA,
        # bty='n', horiz=FALSE, cex=lcex) #maximum legend(x=xmax+xdiff/50,
        # y=ymax-ydiff/1.3275, xjust=0, title='Min', legend=NA, bty='n', horiz=FALSE,
        # cex=lcex) #minimum

    } else if (model == 7 && mstring != 'htsun' && mstring != 'treach' && fill != 'iscore') {

        par(xpd = TRUE)

        legend(x = xmax + xdiff/50, y = ymax - ydiff/25, xjust = 0, title = mlabel, legend = NA, 
            bty = 'n', horiz = FALSE, cex = lcex)  #legend title

        legend(x = xmax, y = ymax - ydiff/10, xjust = 0, legend = c(ctext), lty = c(1, 
            1, 1, 1, 1, 1), lwd = 0.5, col = c('white', 'black', 'black', 'black', 'black', 
            'black'), text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #main legend

        legend(x = xmax, y = ymax - ydiff/1.75, xjust = 0, legend = c(paste('P1', phase1, 
            sep = ' '), paste('P2', phase2, sep = ' '), paste('P3', phase3, sep = ' ')), 
            fill = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), border = 'black', 
            text.col = 'black', bty = 'n', horiz = FALSE, cex = lcex)  #legend for fraction

        legend(x = xmax + xdiff/50, y = ymax - ydiff/9, xjust = 0, title = 'Max', legend = NA, 
            bty = 'n', horiz = FALSE, cex = lcex)  #maximum
    }

    if (hydrograph > 0) {
        legend('bottomleft', lwd = NA, legend = 'Hydrographs', text.col = 'black', bty = 'n', 
            horiz = TRUE, x.intersp = 0.25, inset = c(pos2 - 0.1, posy1), cex = lcex)
        legend('bottomleft', lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = 'green', 
            legend = 'Input', text.col = 'green', bty = 'n', horiz = TRUE, x.intersp = 0.25, 
            inset = c(pos2, posy2), cex = lcex)  #legend for input hydrograph
        legend('bottomleft', lty = 3, lwd = 1.5, pch = 19, pt.cex = 0.5, col = 'purple', 
            legend = 'Output', text.col = 'purple', bty = 'n', horiz = TRUE, x.intersp = 0.25, 
            inset = c(pos2, posy3), cex = lcex)  #legend for output hydrograph
    }
    if (depdef == 1 || impdef == 1) {
        legend('bottomleft', lwd = NA, legend = 'Observation', text.col = 'black', bty = 'n', 
            horiz = TRUE, x.intersp = 0.25, inset = c(pos3 - 0.1, posy1), cex = lcex)
        legend('bottomleft', lty = 1, lwd = 1.5, col = 'red', legend = 'Impact area', 
            text.col = 'red', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos3, 
                posy2), cex = lcex)  #legend for observed impact area
        legend('bottomleft', lty = 1, lwd = 1.5, col = 'orange', legend = 'Deposit', 
            text.col = 'orange', bty = 'n', horiz = TRUE, x.intersp = 0.25, inset = c(pos3, 
                posy3), cex = lcex)  #legend for observed deposit
    }

    # Closing plot file
    dev.off()

    cat(' completed.
')
  }
}
