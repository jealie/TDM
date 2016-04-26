###########################################################
### Cross-link FIA plots with FUTURE climatic variables ###
###########################################################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

## main function, getting the forecasted output of the TDM `tdm_model` for:
# - projected climatic conditions (`rcp=45` or `rcp=85`)
# - current climatic conditions (`rcp=0`)
get_median_prediction = function(tdm_model, p, rcp, clim_scenarios=NULL)
{
  require('raster')
  if (rcp == 0) {
    # current tmp/ppt
    return(get_prediction_for_model(tdm_model, climatic_model_name=paste0(p,'WORLDCLIM/bio_'), climatic_model_ext='.bil'))
  } else {
    # list of climatic scenarios with both RCP4.5 and RCP8.5 models
    if (is.null(clim_scenarios)) {
      clim_scenarios = c('AC', 'BC', 'CC', 'CN', 'GF', 'GS', 'HD', 'HG', 'HE', 'IN', 'IP', 'MI', 'MR', 'MC', 'MP', 'MG', 'NO')
    }
    predictions = list()
    for (clim_model in clim_scenarios) {
      predictions[clim_model] = get_prediction_for_model(tdm_model, 
                                                         climatic_model_name=paste0(p,'WORLDCLIM/',tolower(clim_model),rcp,'bi70'), 
                                                         climatic_model_ext='.tif')
    }
    if (length(predictions) > 1) {
      # median scenario
      predictions_stack = stack(predictions)
      predictions_median = calc(predictions_stack, median_raster)
      return(predictions_median)
    } else {
      # individual scenario
      return(predictions[[1]])
    }
  }
}

# prediction for a specific raster
get_prediction_for_model = function(tdm_model, climatic_model_name, climatic_model_ext)
{
  require('raster')
  require('rgdal')
  
  ## load temperatures and precipitations from worldclim
  cat('\nProcessing scenario:',climatic_model_name,'\n')
  tmp_worldwide = raster(paste0(climatic_model_name, '1', climatic_model_ext))
  tmp = crop(tmp_worldwide, mask_usa)
  values(tmp)[is.na(mask_usa@data@values)] = NA
  ppt_worldwide = raster(paste0(climatic_model_name, '12', climatic_model_ext))
  ppt = crop(ppt_worldwide, mask_usa)
  values(ppt)[is.na(mask_usa@data@values)] = NA
  
  # get the model predictions
  current_map = mask_usa
  values(current_map) = extract(tdm_model, cbind(tmp@data@values, ppt@data@values)/10)
  return(current_map)
}

# helper function to derive the median of models predictions
median_raster = function(x) {
  if (sum(is.na(x)) > length(x)/2) {
    # more than one half of the models is NA
    return (NA)
  } else {
    return(median(x,na.rm=T))
  }
}

