################################
### Functions to output maps ###
################################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

# load a raster layer of the USA shape
load('mask_usa.RData')

# green to red color palette with 100 steps:
colm100=rev(colorRampPalette(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837'))(100))

# The 'rainbow' color palette without the ambiguous colors at the spectrum end (blue to purple)
partial.rainbow = function(n)
{
  return(rev(rainbow(n,start=0,end=0.7)))
}

# helper function to put a colorbar on a raster plot
fudgeit = function(d, colramp, ...){
  require('fields')
  x = seq(0,1,le=length(d))
  y = seq(0,1,le=length(d))
  image.plot(list(x=x,y=y,z=d), col = colramp, legend.only = T, add=T,...)
}
