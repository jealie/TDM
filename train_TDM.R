#############################
### Drought TDM functions ###
#############################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

# computes the TDM based on the input variables
# if threshold is set to `NULL` (default), this routine computes the median model
# otherwise, it computes the thresholded model
get_model = function(tmps, ppts, val, threshold=NULL)
{
  require('raster')
  diverging.colors = colorRampPalette(c(
    '#a50026',
    '#d73027',
    '#f46d43',
    '#fdae61',
    '#fee08b',
    '#d9ef8b',
    '#a6d96a',
    '#66bd63',
    '#1a9850',
    '#006837'))
  colm100=rev(diverging.colors(100))
  
  x.dr = which(!is.na(val) & !is.na(tmps) )
  sub.x.dr = sample(x.dr)
  resolution = 50
  g = get.grid(tmps, ppts, sub.x.dr, resolution) # computes the 50x50 grid
  x = sort(unique(g$x))
  y = sort(unique(g$y))
  z = matrix(as.numeric(g$z),ncol=resolution+3,nrow=resolution+3)
  cutter = rasterFromXYZ(g)
  cutter@data@values[(cutter@data@values)==1] = 1:(sum((cutter@data@values)==1))
  cutted = extract(cutter, cbind(tmps[x.dr]/10, ppts[x.dr]/10))
  inter_model = rasterFromXYZ(g)
  if (is.null(threshold)) {
    # returns the median model
    inter_model@data@values[(inter_model@data@values)==1] = sapply(1:(sum((inter_model@data@values)==1)),
                                                                   function(i){
                                                                     subs = x.dr[which(cutted==i)]
                                                                     data = val[subs]
                                                                     # mean results in similar models
                                                                     return(median(data,na.rm=T))
                                                                   })
  } else if (threshold < 0) {
    # computes the standard deviation model
    inter_model@data@values[(inter_model@data@values)==1] = sapply(1:(sum((inter_model@data@values)==1)),
                                                                   function(i){
                                                                     subs = x.dr[which(cutted==i)]
                                                                     data = val[subs]
                                                                     return(sd(data,na.rm=T))
                                                                   })
  } else {
    # computes the likelihood model (based on threshold)
    inter_model@data@values[(inter_model@data@values)==1] = sapply(1:(sum((inter_model@data@values)==1)),
                                                                   function(i){
                                                                     subs = x.dr[which(cutted==i)]
                                                                     data = val[subs]
                                                                     return(quantile(data,na.rm=T, probs=threshold)[[1]])
                                                                   })
  }
  mask = rasterFromXYZ(g)
  inter_model@data@values[which( mask@data@values!=1)] = NA
  return(inter_model)
}

# helper function to project a variable into a different scale
do.scale = function(x,wrt=NA)
{
  if (is.na(wrt[1])) {
    M = max(x,na.rm=T)
    m = min(x,na.rm=T)
  } else {
    M = max(wrt,na.rm=T)
    m = min(wrt,na.rm=T)
  }
  return ((x-m)/(M-m))
}

# inverse operation
undo.scale = function(x,wrt)
{
  M = max(wrt,na.rm=T)
  m = min(wrt,na.rm=T)
  return (x*(M-m)+m)
}

# grids tmp/ppt data at a given resolution
# also excludes cells with too few data (default: nmin=10)
get.grid = function(tmp, ppt, sub.x.dr, resolution, nmin=10)
{
  require('FNN')
  scaled.tmp = do.scale(tmp[sub.x.dr],tmp)
  scaled.ppt = do.scale(ppt[sub.x.dr],ppt)
  x = seq(from=min(scaled.tmp, na.rm=T),to=max(scaled.tmp, na.rm=T),le=resolution+1)
  y = seq(from=min(scaled.ppt, na.rm=T),to=max(scaled.ppt, na.rm=T),le=resolution+1)
  cat(paste0('Taking the resolution: x=',abs(diff(undo.scale(x[1:2],tmp)))/10,
             ' and y=',abs(diff(undo.scale(y[1:2],ppt)))/10,'\n'))
  m_tmp = x[1] - diff(x[1:2])
  m_ppt = y[1] - diff(y[1:2])
  M_tmp = x[length(x)] + diff(x[(length(x)-1):(length(x))])
  M_ppt = y[length(y)] + diff(y[(length(y)-1):(length(y))])
  g = expand.grid(x=c(m_tmp,x,M_tmp),
                  y=c(m_ppt,y,M_ppt))
  g = cbind(g, z=get.knnx(cbind(scaled.tmp, 
                                scaled.ppt),
                          cbind(g$x,
                                g$y),
                          k=nmin)$nn.dist[,nmin] < 1/resolution)
  return(data.frame(x=undo.scale(g$x,tmp)/10,
                    y=undo.scale(g$y,ppt)/10,
                    z=g$z))
}

