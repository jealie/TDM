##########################
### Maps of tolerances ###
##########################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

# outputs the map with tolerances in the Radviz coordinate system
plot_radviz_map = function(fia_db)
{
  tolerances_noNA = fia_db[,colnames(fia_db) %in% c("ID", "DTI", "STI", "WTI", "LON", "LAT")]
  tolerances_noNA = tolerances_noNA[complete.cases(tolerances_noNA),]
  tolerances_noNA = aggregate(tolerances_noNA,by=list(tolerances_noNA$ID),mean)[,-1]
  rad = get_coords_radviz(tolerances_noNA[,c(3,2,4,1)], normalize = F)
  cat('Starting to plot the Radviz legend...')
  circledraw(3,col=grey(0.95))
  text(1+0.25,0,'Drought\nTolerance\nIndex',cex=1.5,xpd=T)
  text(-0.6+0.1,0.7+0.35,'Shade\nTolerance\nIndex',cex=1.5,xpd=T)
  text(-0.6+0.1,-0.7-0.35,'Waterlogging\nTolerance\nIndex',cex=1.5,xpd=T)
  coords_span = seq(0.01,1,le=101)
  coords_legend = expand.grid(coords_span,coords_span,coords_span)
  coords_legend = cbind(coords_legend, 1:nrow(coords_legend))
  coords_rad = get_coords_radviz(coords_legend, normalize = F)
  points(coords_rad, col=comp_hsv(coords_rad))
  cat('Done. Starting to plot the cumulative distribution...')
  circledraw(3,col=grey(0.95))
  res = 50
  xbreaks = seq(-0.5,1,le=res+1)
  ybreaks = seq(-0.86604,0.86604,le=res+1)
  xgroup = cut(rad[,1], xbreaks, include.lowest=TRUE)
  ygroup = cut(rad[,2], ybreaks, include.lowest=TRUE)
  m = table(xgroup,ygroup)
  mmax = quantile(m,probs = 0.95)[[1]]
  colormap_density = c(NA,partial.rainbow(mmax))
  end_col = colormap_density[length(colormap_density)]
  colormap_density = c(colormap_density,rep(end_col, max(m)-length(colormap_density)))
  image((xbreaks[1:res]+xbreaks[2:(res+1)])/2,
        (ybreaks[1:res]+ybreaks[2:(res+1)])/2,
        m,
        col = colormap_density,
        add=T)
  circledraw(3,col=grey(0.95),add=T)
  text(1+0.25,0,'Drought\nTolerance\nIndex',cex=1.5,xpd=T)
  text(-0.6+0.1,0.7+0.35,'Shade\nTolerance\nIndex',cex=1.5,xpd=T)
  text(-0.6+0.1,-0.7-0.35,'Waterlogging\nTolerance\nIndex',cex=1.5,xpd=T)
  legend('bottomright', xpd=T,bty='n',
         fill = c(grey(0.95),
                  partial.rainbow(mmax)[mmax/5],
                  partial.rainbow(mmax)[mmax*2/5],
                  partial.rainbow(mmax)[mmax*3/5],
                  partial.rainbow(mmax)[mmax*4/5],
                  end_col),
         legend = c('0%', '20%', '40%', '60%','80%','100%'),
         title='\n\nCumulative distribution')
  cat('Done. Starting to plot the map...')
  require('maps')
  require('mapdata')
  par(mar=c(0,0,0,0), mai=c(0,0,0,0), oma=c(0,0,0,0))
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55))
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55),fill=T,col=grey(0.95),add=T)
  cols = comp_hsv(rad)
  x = sample(1:nrow(tolerances_noNA),nrow(tolerances_noNA))
  points(tolerances_noNA$LON[x],
         tolerances_noNA$LAT[x],
         col = cols[x],
         pch=20,cex=0.5)
  cat('Done.')
}


# greatly inspired by the (now defunct) dprep package
get_coords_radviz = function(dataset, name='', normalize=T)
{
  isfact=FALSE
  n=dim(dataset)[1]
  p=dim(dataset)[2]
  classes=dataset[,p]
  if (class(classes)=="factor")
  {isfact=TRUE
  classnames= levels(classes)
  }
  classnumbers=1:length(unique(classes))
  classes=as.numeric(classes,drop=TRUE)
  varnames=colnames(dataset)
  #created projection of each observation
  if (normalize) {
    dataset=as.matrix(mmnorm(dataset))
  } else { 
    dataset = as.matrix(dataset)
  }
  dataset=dataset[,-p]
  sumrows=rowSums(dataset)
  columns=seq(0,(p-2))
  angles=(2*pi*columns)
  angles=angles/(p-1)
  cosines=cos(angles)
  sines=sin(angles)
  proj.x=(dataset%*%cosines)
  proj.x=proj.x/sumrows
  proj.y=(dataset%*%sines)
  proj.y=proj.y/sumrows
  r = cbind(proj.x, proj.y)
  rownames(r)=NULL
  return(r)
}


# greatly inspired by the (now defunct) dprep package
mmnorm = function (data,minval=0,maxval=1)
{
  #This is a function to apply min-max normalization to a matrix or dataframe.
  #Min-max normalization subtracts the minimum of an attribute from each value
  #of the attribute and then divides the difference by the range of the attribute.
  #These new values are multiplied by the given range of the attribute
  #and finally added to the given minimum value of the attribute.
  #These operations transform the data into [minval,mxval].
  #Usually minval=0 and maxval=1.
  #Uses the scale function found in the R base package.
  #Input: data= The matrix or dataframe to be scaled
  #store all attributes of the original data
  d=dim(data)
  c=class(data)
  cnames=colnames(data)
  #remove classes from dataset
  classes=data[,d[2]]
  data=data[,-d[2]]
  minvect=apply(data,2,min)
  maxvect=apply(data,2,max)
  rangevect=maxvect-minvect
  zdata=scale(data,center=minvect,scale=rangevect)
  newminvect=rep(minval,d[2]-1)
  newmaxvect=rep(maxval,d[2]-1)
  newrangevect=newmaxvect-newminvect
  zdata2=scale(zdata,center=FALSE,scale=(1/newrangevect))
  zdata3=zdata2+newminvect
  zdata3=cbind(zdata3,classes)
  if (c=="data.frame") zdata3=as.data.frame(zdata3)
  colnames(zdata3)=cnames
  return(zdata3)
}

# draw the triangle for the legend key
circledraw = function(numpts=200,radius=1,col='white',add=F)
{
  xy=rep(360/numpts, numpts)
  xy <- c(0, cumsum(xy))
  xy<-xy*pi/180
  t2xy = list(x = radius*cos(xy), y = radius*sin(xy))
  if (!add) plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1] > pin[2]) xlim <- (pin[1]/pin[2]) * xlim else ylim <- (pin[2]/pin[1]) * ylim
  plot.window(xlim, ylim, "", asp = 1)
  if (!add) {
    polygon(t2xy$x, t2xy$y, border = "black", col = col)
  } else {
    polygon(c(t2xy$x,rev(t2xy$x*2)), c(t2xy$y,rev(t2xy$y*2)), border = NA, col = 'white')
    polygon(t2xy$x, t2xy$y, border = "black", col = NULL)
  }
  invisible(NULL)
}

# computes the color gradient in the HSV space
comp_hsv = function(rad)
{
  return(hsv(h=(
    acos(abs(rad[,1])/sqrt(rowSums(rad^2)))
    * ifelse(((rad[,2]<0) & (rad[,1]>0)) | ((rad[,2]>0) & (rad[,1]<0)),-1,1)
    + ifelse(rad[,1]<0,pi,0)
    + ifelse((rad[,2]<0) & (rad[,1]>0),2*pi,0)
  )/2/pi,
  s=sqrt(rowSums(rad^2))*0.95+0.05,
  v=1-0.1*(1-sqrt(rowSums(rad^2)))^2 ))
}


# separate maps with only shade, drought and waterlogging tolerance
plot_separate_maps = function(fia_db)
{
  require('maps')
  require('mapdata')
  par(mar=c(0,0,0,0), mai=c(0,0,0.5,0), oma=c(0,0,0,0))
  x = sample(1:nrow(fia_db),nrow(fia_db))
  
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55))
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55),fill=T,col=grey(0.95),add=T)
  points(fia_db$LON[x],
         fia_db$LAT[x],
         col = colm100[cut(fia_db$DTI[x],breaks=seq(0,1,le=101),include.lowest = T)],
         pch=20,cex=0.25)
  fudgeit(c(0,1),colm100,legend.shrink=0.3)
  title('Drought tolerance index', cex.main=3)
  
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55))
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55),fill=T,col=grey(0.95),add=T)
  points(fia_db$LON[x],
         fia_db$LAT[x],
         col = colm100[cut(fia_db$STI[x],breaks=seq(0,1,le=101),include.lowest = T)],
         pch=20,cex=0.25)
  fudgeit(c(0,1),colm100,legend.shrink=0.3)
  title('Shade tolerance index', cex.main=3)
  
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55))
  map("worldHires","usa", xlim=c(-126,-65), ylim=c(20,55),fill=T,col=grey(0.95),add=T)
  points(fia_db$LON[x],
         fia_db$LAT[x],
         col = colm100[cut(fia_db$WTI[x],breaks=seq(0,1,le=101),include.lowest = T)],
         pch=20,cex=0.25)
  fudgeit(c(0,1),colm100,legend.shrink=0.3)
  title('Waterlogging tolerance index', cex.main=3)
}
