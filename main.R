####################################################################
## Source Code Supplement to Liénard, Harrisson and Strigul, 2016 ##
####################################################################
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Any results obtained with this code should refer to the following publication:
# Liénard, J. and Harrisson, J. and Strigul, N., "U.S. Forest Response to Projected Climate-Related Stress: a Tolerance Perspective", Global Change Biology (2016), doi: 10.1111/gcb.13291
#####

#####
# This is the master file of the program, which goes through the majors steps of the Tolerance Distribution Model (TDM) approach
# The sub-routines are located in separate files in the same directory
#
# Note 1: you need to download the FIA and Worldclim datasets separately to run this script.
# BIL/HDR Rasters from Worldclim must be accessible in data_path/WORDLIM/bio_XXX (e.g.: data_path/WORDLIM/bio_1.bil and data_path/WORDLIM/bio_1.hdr)
# FIA database for each state must be downloaded as CSV and placed in data_path/WORDLIM/XXX_COND.CSV ; data_path/WORDLIM/XXX_PLOT.CSV ; data_path/WORDLIM/XXX_TREE.CSV
# (e.g.: data_path/WORDLIM/AL_COND.CSV, data_path/WORDLIM/AL_PLOT.CSV and data_path/WORDLIM/AL_TREE.CSV)
#
# Note 2: this script requires several third party packages, which you can install with:
# install.packages('doParallel','raster','rgdal','fields','maps','mapdata','FNN')
#####


### step 0: set the environment variables
code_path = '/media/jean/ext4/forest_code/TDM/' # change to where the code is
data_path = '/media/jean/ext4/forest_data/' # change to where the data are (FIA and WORLDCLIM directories)
setwd(code_path)
source('plot_fct.R') # loads helper functions to plot maps


### step 1: load the FIA database and compute stand-level tolerance indices
source('fia_extract.R')
fia_db = fia_extract(data_path)


### step 2: output the maps of current tolerances (radviz visualization of Panel 1 and Supplementary Results)
source('radviz.R')
plot_separate_maps(fia_db) # 3 maps, one for each tolerance
plot_radviz_map(fia_db) # the Radviz visualization of panel 1


### step 3: load climatic variables and cross-link with stand-level tolerances 
source('clim_extract.R')
tmpppt = extract_tmpppt(fia_db, data_path)
plot_tolerance_climatic(tmpppt$tmps, tmpppt$ppts, fia_db$DTI, 'Drought Tolerance Index') # output the bottom visualization of Panel 1
plot_tolerance_climatic(tmpppt$tmps, tmpppt$ppts, fia_db$STI, 'Shade Tolerance Index')
plot_tolerance_climatic(tmpppt$tmps, tmpppt$ppts, fia_db$WTI, 'Waterlogging Tolerance Index')


### step 4: compute and plot the drought TDM and its standard deviation
source('train_TDM.R')
# Computes the TDM and its standard deviation:
current_model = get_model(tmpppt$tmps, tmpppt$ppts, fia_db$DTI)
plot(current_model, main='Drought TDM', col=colm100, asp=NA)
sd_model = get_model(tmpppt$tmps, tmpppt$ppts, fia_db$DTI, -1)
plot(sd_model, main='Model standard deviation', col=colm100, asp=NA)

# computes also the models for the 90%, 95% and 99% likelihood analyses
thresh90_model = get_model(tmpppt$tmps, tmpppt$ppts, fia_db$DTI, 0.9)
thresh95_model = get_model(tmpppt$tmps, tmpppt$ppts, fia_db$DTI, 0.95)
thresh99_model = get_model(tmpppt$tmps, tmpppt$ppts, fia_db$DTI, 0.99)


### step 5: use TDM with current and projected climatic variables
source('predict_TDM.R')
# Get the current potential distribution (like Figure 2d in main text and Figure 14-top in Supp. Results)
potential_distrib = get_median_prediction(current_model, data_path, rcp=0)
plot(potential_distrib, col=colm100) 

# Forecast expected tolerance under a single projected climatic scenario (like in the Supp. Results figures 16 to 18):
scenario_AC_85 = get_median_prediction(current_model, data_path, rcp=85, clim_scenarios='AC') # in this example we use the scenario "AC"
plot(scenario_AC_85, col=colm100, main='Climatic Model AC with RCP 85')

# Computes the median predicted tolerance distribution across all climatic scenarios (like Supp. Results figure 14-bottom)
median_pred_45 = get_median_prediction(current_model, data_path, rcp=45)
plot(median_pred_45, col=colm100, main='Median predictions with RCP 45')
median_pred_85 = get_median_prediction(current_model, data_path, rcp=85)
plot(median_pred_85, col=colm100, main='Median predictions with RCP 85')

# Computes the difference between potential distribution and median predicted tolerance across all climatic scenarios (like 3c-d in main text)
diff_45 = potential_distrib
values(diff_45) = values(median_pred_45) - values(diff_45) 
plot(diff_45, main='Magnitude of projected drought tolerance (with RCP 45)')
diff_85 = potential_distrib
values(diff_85) = values(median_pred_85) - values(diff_85) 
plot(diff_85, main='Magnitude of projected drought tolerance (with RCP 85')

# Compute the vulnerability of current forests with a threshold of 90%, 95% and 99%
vuln_90percent = get_median_prediction(thresh90_model, data_path, rcp=0)
vuln_95percent = get_median_prediction(thresh95_model, data_path, rcp=0)
vuln_99percent = get_median_prediction(thresh99_model, data_path, rcp=0)

# Output the summary vulnerabilities of current forests with thresholds at 90%, 95% and 99% (like Figure 3e-f)
vuln_RCP45_all = mask_usa
values(vuln_RCP45_all) = (values(median_pred_45) > values(vuln_90percent)) + 
  2*(values(median_pred_45) > values(vuln_95percent)) + 
  4*(values(median_pred_45) > values(vuln_99percent))
plot(vuln_RCP45_all, main='Confidence in drought tolerance shift (with RCP 45)', col=c('#ebffcc','#FFD796',rep('#FC8368',2),rep('#FF0000',4)), legend=F)
legend('bottomright', c('Within expected range', '> 90%', '> 95%', '> 99%'), fill=c('#ebffcc','#FFD796','#FC8368','#FF0000'), bty='n')
vuln_RCP85_all = mask_usa
values(vuln_RCP85_all) = (values(median_pred_85) > values(vuln_90percent)) + 
  2*(values(median_pred_85) > values(vuln_95percent)) + 
  4*(values(median_pred_85) > values(vuln_99percent))
plot(vuln_RCP85_all, main='Confidence in drought tolerance shift (with RCP 85)', col=c('#ebffcc','#FFD796',rep('#FC8368',2),rep('#FF0000',4)), legend=F)
legend('bottomright', c('Within expected range', '> 90%', '> 95%', '> 99%'), fill=c('#ebffcc','#FFD796','#FC8368','#FF0000'), bty='n')

