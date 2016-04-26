####################################################
### Cross-link FIA plots with climatic variables ###
####################################################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

## main function:
extract_tmpppt = function(fia_db, p)
{
  require('raster')
  require('rgdal')
  tmps = raster(paste0(p,'WORLDCLIM/bio_1.bil'))
  ppts = raster(paste0(p,'WORLDCLIM/bio_12.bil'))
  coords = cbind(fia_db$LON,
                 fia_db$LAT)
  coords = SpatialPoints(coords, proj4string = CRS("+proj=longlat +datum=NAD83"))
  coords = spTransform(coords, crs(tmps))
  tmps = extract(tmps,coords)
  ppts = extract(ppts,coords)
  return(list(tmps=tmps, ppts=ppts))
}

# shows the tolerance of FIA stands in the climatic space
plot_tolerance_climatic = function(tmps, ppts, val, title_label, cutoff=NA)
{
  xr = sample(length(val))
  plot(tmps[xr]/10, ppts[xr]/10,
       pch=1,
       col=colm100[cut(val[xr],breaks=seq(0,1,le=101),include.lowest = T)],
       xlim=c(-5,30),bty='n',las=1,
       xlab='Temperature (⁰C)',ylab='Precipitation (mm/month)')
  mtext(title_label,font = 2)
  fudgeit(c(0, 1),colm100,smallplot=c(0.825,0.85,0.2,0.8))
}
