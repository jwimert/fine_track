##!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
subroutines ported from l3a_si_finetracker_mod.f90
"""
#
import numpy as np
import scipy.linalg



##!
##!==============================================================================
##!
##
## FINE_TRACK 
##
## (Main Routine)
##
## Pass in photon heights and compute sea ice segment height
##
## Notes:
##
## - coarse_mn is needed to subtract from photon heights to get near zero,
##   if coarse_mn is not available, try using an average
##
##
##!
##!==============================================================================
##!

#
def fine_track(ph_h_in, coarse_mn, n_photons, n_shots, bin_size, lb_bin, ub_bin, wf_table, wf_bins, wf_sd, wf_mn, cal19_corr, cal19_width, cal19_strength, cal19_dt, dead_time):
#
# Computes fine_track surface height
#
# Input:
#
# ph_h_in - collection of photon heights (n = n_photons)
# coarse_mn - coarse tracker mean height
# n_photons - number of photons collected
# n_shots - shots covered to collect n_photons
# bin_size - histogram bin size
# lb_bin - histogram lower bound
# ub_bin - histogram upper bound
# wf_table - expected waveform table
# wf_bins - expected waveform bins
# wf_sd - expected waveform table stdev values
# wf_mn - expected waveform table mean values
# cal19_corr - first photon bias correction look-up table (6, 160, 498)
# cal19_width - cal19 width ancillary table (6, 498)
# cal19_strength - cal19 strength ancillary table (6, 160)
# cal19_dt -  cal19 dead-time ancillary table (6)
# dead_time - system dead-time (cal19_corr input)
#
#
# Output:
#
# h_surf - surface height
# w_gauss - stdev
#

#
# Set histogram bin values (center and edges)
#
  bins_center = np.arange(lb_bin, ub_bin+bin_size, bin_size)
  bins_edges = np.arange(lb_bin - bin_size*0.5, ub_bin + bin_size, bin_size)

#
# Compute photon_Rate
#
  photon_rate = n_photons / n_shots

#
# Call hist_full:
#
  ph_h, hist_full_out, bins_full = \
    hist_full(ph_h_in, coarse_mn, bins_edges)

#
# Call hist_trim1: 
# (+/- lb/ub_bin, but centered around histogram mode)
# Compute mean_trim1
#
  ph_h_trim1, hist_trim1_out, bins_trim1 = \
    hist_trim1(ph_h, hist_full_out, bins_center, bins_edges, bin_size, lb_bin, ub_bin)

  mean_trim1 = np.mean(ph_h_trim1)

#
# Call hist_trim2: 
# (+/- 2*stdev)
#
  ph_h_trim2, hist_trim2_out, bins_trim2 = \
    hist_trim2(ph_h_trim1, bins_edges)

#
# Call hist_trim3:
# (remove bins LT 2 before/after first/last occurance of bin GE 2)
#
  ph_h_trim3, hist_trim3_out = \
    hist_trim3(ph_h_trim2, hist_trim2_out, bins_center, bins_edges, bin_size)

##
## Trim histogram and expected waveform table
##

#
# Subtrack mean (from after trim1) from photon heights and re-create histogram
#
  ph_h_trim_fit = ph_h_trim3 - mean_trim1

#
# Construct histogram using new heights
#
  hist_trim_fit, bins_trim_fit = np.histogram(ph_h_trim_fit, bins_edges) 

#
# Count number of trimmed histograms
#
  n_photons_trim = hist_trim_fit.sum()

#
# Call Gauss Fit:
# fit waveform to expected waveform table
#

  error_surface, biquad_h, biquad_sd = \
    gauss_fit(hist_trim_fit, wf_table, wf_bins, wf_sd, wf_mn, bin_size, bins_center)

#
# Call First Photon Bias
#
  fpb_corr_m, fpb_width, fpb_strength, fpb_dt = \
    fpb_corr(cal19_corr, cal19_width, cal19_strength, cal19_dt, \
    dead_time, photon_rate, hist_full_out, bins_center)

#
# Build surface height, set w_gauss
#
  h_surf = coarse_mn + mean_trim1 + biquad_h - fpb_corr_m

  w_gauss = biquad_sd

#
# Return
#
  return h_surf, w_gauss, fpb_corr_m









##!
##!==============================================================================
##!
##
##
## HIST_FULL
##
## Construct initial histogram using all photons
##
##!
##!==============================================================================
##!

#
def hist_full(ph_h, coarse_h, bin_edges):
#
# remove coarse mean from heights
#
  ph_h_zero = ph_h - coarse_h
#
# compute histogram
#
  hist, bins = np.histogram(ph_h_zero, bin_edges) 
#
  return ph_h_zero, hist, bins 


##!
##!==============================================================================
##!
##
##
## HIST_TRIM1
##
## +/- lb/ub_bin, but centered around histogram mode
##
##
##!
##!==============================================================================
##!

#
def hist_trim1(ph_h_zero, hist_full, bins_center, bins_edges, bin_size, lb_bin, ub_bin):
#
# Find Mode
#
  mode_index = hist_full.argmax()
#
# Truncate photons outside window around histogram mode
#
  min_win = bins_center[mode_index] + lb_bin
  max_win = bins_center[mode_index] + ub_bin
#
  min_win = min_win - bin_size*0.5
  max_win = max_win + bin_size*0.5
#
# Keep photons within window
#
  n_ph_trim1 = 0
  list_temp = []
  for i in np.arange(ph_h_zero.size):
    if ((ph_h_zero[i]>min_win)&(ph_h_zero[i]<max_win)):
      n_ph_trim1 = n_ph_trim1 + 1
      list_temp.append(ph_h_zero[i])
#
  ph_h_trim1 = np.array(list_temp)
#
# compute histogram
#
  hist_trim1, bins_trim1 = np.histogram(ph_h_trim1, bins_edges) 
#
# return
#
  return ph_h_trim1, hist_trim1, bins_trim1



##!
##!==============================================================================
##!
##
##
## HIST_TRIM2
##
## +/- 2*stdev
##
##
##!
##!==============================================================================
##!

#
def hist_trim2(ph_h, bins_edges):
#
# Compute mean and standard deviation 
#
  trim1_mean = ph_h.mean()
  trim1_stdev = ph_h.std()
#
# Compute trim window top and bottom
#
  trim2_height_bot = trim1_mean - (2 * trim1_stdev)
  trim2_height_top = trim1_mean + (2 * trim1_stdev)
#
# Keep photons within window
#
  n_ph_trim2 = 0
  list_temp = []
  for i in np.arange(ph_h.size):
    if ((ph_h[i]>trim2_height_bot)&(ph_h[i]<trim2_height_top)):
      n_ph_trim2 = n_ph_trim2 + 1
      list_temp.append(ph_h[i])
#
  ph_h_trim2 = np.array(list_temp)
#
# compute histogram
#
  hist_trim2, bins_trim2 = np.histogram(ph_h_trim2, bins_edges) 

#
# return
#
  return ph_h_trim2, hist_trim2, bins_trim2




##!
##!==============================================================================
##!
##
##
## HIST_TRIM3
##
## remove bins LT 2 before/after first/last occurance of bin GE 2
##
##
##!
##!==============================================================================
##!

#
def hist_trim3(ph_h, hist_trim2, bins_center, bins_edges, bin_size):
#
# All bins with less than n_photon_trim photons before the first bin
# and after the last bin with at least n_photon_trim photons are removed.
#
  hist_trim = hist_trim2
  for i in np.arange(hist_trim.size):
    if (hist_trim[i] >= 2):
      pre_trim_h0 = bins_center[i] - bin_size*0.5
      break
    else:
      hist_trim[i] = 0

  for i in reversed(np.arange(hist_trim.size)):
    if (hist_trim[i] >= 2):
      pre_trim_h1 = bins_center[i] + bin_size*0.5
      break
    else:
      hist_trim[i] = 0
#
# Determine and save photons within populated bins after the pre-trim
#

  num_wins_ph2 = 0
  list_temp = []
  for i in np.arange(ph_h.size):
    if ((ph_h[i]>pre_trim_h0)&(ph_h[i]<pre_trim_h1)):
      num_wins_ph2 = num_wins_ph2 + 1
      list_temp.append(ph_h[i])
  
  ph_h_trim3 = np.array(list_temp)
#
# compute histogram
#
  hist_trim3, bins_trim2 = np.histogram(ph_h_trim3, bins_edges) 

#
#   return
#
  return ph_h_trim3, hist_trim3





##!
##!==============================================================================
##!
##
##
## GAUSS_FIT
##
## Use table of expected waveforms to compute value of height and stdev
##
##
##!
##!==============================================================================
##!

#
def gauss_fit(hist_trim_fit, wf_table, wf_bins, wf_sd, wf_mn, bin_size, bins_center):
#
# Find first and last non-zero bin of histogram
#
# Store trimmed hist and bins
#
  hist_gauss_pretrim = hist_trim_fit

  for i in np.arange(hist_gauss_pretrim.size):
    min_bin = i
    if (hist_gauss_pretrim[i] > 0):
      break

  for i in reversed(np.arange(hist_gauss_pretrim.size)):
    max_bin = i
    if (hist_gauss_pretrim[i] > 0):
      break
      
  hist_temp = []
  bins_temp = []
  for i in np.arange(min_bin, max_bin+1, 1, dtype=int):
    hist_temp.append(hist_gauss_pretrim[i])
    bins_temp.append(bins_center[i])
    
  hist_gauss_trim = np.array(hist_temp)
  bins_gauss_trim = np.array(bins_temp)

#
# Compute normalized trimmed histogram
#
  norm_gauss_hist = hist_gauss_trim / hist_gauss_trim.sum()

#
# Trim input waveform table
#
  wf_min_bin = round( abs(np.min(wf_bins) - bins_center[min_bin]) / bin_size)
  wf_max_bin = round( abs(np.min(wf_bins) / bin_size) + abs(bins_center[max_bin] / bin_size) )

  wf_table_trim = wf_table[:, :, wf_min_bin : wf_max_bin+1]
  wf_bins_trim = wf_bins[wf_min_bin : wf_max_bin+1]

#
# Compute error surface
#
  error_surface = np.zeros((wf_mn.size, wf_sd.size))

  for i in np.arange(wf_mn.size):
    for j in np.arange(wf_sd.size):
        #
        # Add check for bad peak
        #      
        if ( (sum(wf_table_trim[i, j, :]) > 0.0) & (np.nanargmax(wf_table_trim[i, j, :]) != 0) &  (np.nanargmax(wf_table_trim[i, j, :]) != wf_bins_trim.size) ):
            error_surface[i,j] = sum( ( (wf_table_trim[i, j, :]/sum(wf_table_trim[i, j, :])) - norm_gauss_hist[:])**2 )
        else:
            error_surface[i,j] = np.nan
  
#
# Replace Nans with Max
#
  error_surface[np.isnan(error_surface)] = np.nanmax(error_surface)

#
# Find local minimum
#
  h_min = np.unravel_index(np.nanargmin(error_surface), error_surface.shape)[0]
  sd_min = np.unravel_index(np.nanargmin(error_surface), error_surface.shape)[1]

#
# Check for edge cases
#
  if (sd_min == 0):
    sd_min = 1
  if (h_min == 0):
    h_min = 1
    
  if (sd_min == wf_sd.size - 1):
    sd_min = wf_sd.size - 2
  if (h_min == wf_mn.size - 1):
    h_min = wf_mn.size - 2

#
# setup biquad fit
#
    
  biquad_h = wf_mn[h_min-1 : h_min+2]
  biquad_sd = wf_sd[sd_min-1 : sd_min+2]
  biquad_error = np.zeros((3, 3))

  for ii in np.arange(0,3):
    for jj in np.arange(0,3):
        biquad_error[ii,jj] = error_surface[h_min-1+ii,sd_min-1+jj]

  biquad_error = np.transpose(biquad_error)

#
# Construact 21x21 mesh to fit with 2nd order quadratic curve
#
  mesh_xi, mesh_yi = np.meshgrid(biquad_h, biquad_sd)
  x_biquad = np.linspace(biquad_h[0],biquad_h[2],21)
  y_biquad = np.linspace(biquad_sd[0],biquad_sd[2],21)

  data = np.c_[mesh_xi.ravel(), mesh_yi.ravel(), biquad_error.ravel()]

  X,Y = np.meshgrid(x_biquad, y_biquad)
  XX = X.flatten()
  YY = Y.flatten()

#
# best-fit quadratic curve (2nd-order)
#
  A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
  C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
    
#
# evaluate it on a grid
#
  Z_quad = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)

#
# Locate minimum
#
  z_quad_min_0 = np.unravel_index(Z_quad.argmin(), Z_quad.shape)[0]
  z_quad_min_1 = np.unravel_index(Z_quad.argmin(), Z_quad.shape)[1]

#
# Save local height and gaussian
#
  biquad_h = x_biquad[z_quad_min_1]
  biquad_sd = y_biquad[z_quad_min_0]

#
#   return
#
  return error_surface, biquad_h, biquad_sd



##!
##!==============================================================================
##!
##
##
##
## FPB_CORR
##
## Compute first photon bias correction
##
##
##!
##!==============================================================================
##!

#
def fpb_corr(cal19_corr, cal19_width, cal19_strength, cal19_dt, dead_time, photon_rate, hist, bins):
#
# Input:
#   cal19_corr (6, 160, 498)
#   cal19_width (6, 498)
#   cal19_strength (6, 160)
#   cal19_dt (6)
#   dead_time
#   photon_rate
#   hist
#   bins
#


#
# Compute width from histogram
#

#
# Compute 10th and 90th energy percentile
#

  hist_sum = sum(hist)
  ht_10p = hist_sum * 0.1
  ht_90p = hist_sum * 0.9

#
# Construct accumulated energy array
#
  hist_accum_temp = []
  hist_accum_temp.append(hist[0])
  for ii in np.arange(1,hist.size):
    hist_accum_temp.append(hist[ii]+hist_accum_temp[ii-1])
    
  hist_accum = np.array(hist_accum_temp)

#
# Find 10th and 90th percentile bins
#
  for ii in np.arange(1,hist.size):
    bin_10p = ii
    if (hist_accum[ii]>=ht_10p):
        break

  for ii in np.arange(1,hist.size):
    bin_90p = ii
    if (hist_accum[ii]>=ht_90p):
        break

#
# Interpolate width bins to percentile height
#
  h1 = bins[bin_10p-1] + \
    (bins[bin_10p] - bins[bin_10p-1]) * \
    (ht_10p - hist_accum[bin_10p-1]) / (hist_accum[bin_10p] - hist_accum[bin_10p-1])
  h2 = bins[bin_90p-1] + \
    (bins[bin_90p] - bins[bin_90p-1]) * \
    (ht_90p - hist_accum[bin_90p-1]) / (hist_accum[bin_90p] - hist_accum[bin_90p-1])

#
# Convert width from m to ns
#
  width_ns = (h2-h1) / (3.0E8/2) * 1E9

#
# Use CAL19 ancillary arrays to compute indexes for CAL19 look-up table
#

#
# Compute dead time index
#
  diff_dead_time = (np.abs(cal19_dt - dead_time).min())
  index_dead_time = (np.abs(cal19_dt - dead_time).argmin())

#
# Compute strength index
#
  diff_strength = (np.abs(cal19_strength[index_dead_time] - photon_rate).min())
  index_strength = (np.abs(cal19_strength[index_dead_time] - photon_rate).argmin())

#
# Compute width index
#
  diff_width = (np.abs(cal19_width[index_dead_time] - width_ns).min())
  index_width = (np.abs(cal19_width[index_dead_time] - width_ns).argmin())

#
# Use cal_19 lookup table to find correction
# Convert from ps to m
#
  fpb_corr_ps = cal19_corr[index_dead_time, index_strength, index_width]
  fpb_corr_m = fpb_corr_ps * 1E-12 * (3E8/2)

#
# Set output parameters
#
  fpb_width = width_ns
  fpb_strength = photon_rate
  fpb_dt = dead_time

#
# return
#
  return fpb_corr_m, fpb_width, fpb_strength, fpb_dt




##!
##!==============================================================================
##!
##
##
## SPEC_SHOT_FILTER
##
## Determine if there are any specular shots, and set photon heights within
## specular shot to invalid (-9999.0)
##
##
##!
##!==============================================================================
##!

#
def spec_shot_filter(i0, i1, ATL03_ph_height, ATL03_mframe, ATL03_pulse, ATL03_ph_conf) :
#
# Input:
#
# i0 - starting index of ATL03 to search
# i1 - ending index of ATL03 to search
# ATL03_ph_height - photon heights
# ATL03_mframe - mainframe
# ATL03_pulse - pulse
# ATL03_ph_conf - photon signal confidence interval
#
#
# Output:
#
# n_spec_shots - number of specular shots encountered
# spec_shot_mf - value of mainframe for specular shot
# spec_shot_pulse - value of pulse for specular shot
#

#
# Initialize
#
  spec_shot_out = []
#   spec_shot_mf = []
#   spec_shot_pulse = []	
  n_spec_shots = 0

#   print('MFRAMES:',ATL03_mframe[i0],ATL03_mframe[i1])
  for ii in np.arange(ATL03_mframe[i0],ATL03_mframe[i1]+1):
#     print('MFRAME',ii)
    for jj in np.arange(0,200):
      n_photons = np.count_nonzero( (ATL03_mframe == ii) & (ATL03_pulse== jj) & (ATL03_ph_conf[:,2] >= 3 ) )
      if (n_photons > 16):
#         print('SPECULAR SHOT',ii,jj, n_photons)
        n_spec_shots = n_spec_shots + 1
#         spec_shot_mf.append(ii)
#         spec_shot_pulse.append(jj)
        ATL03_ph_height[(ATL03_mframe == ii) & (ATL03_pulse== jj) ] = -9999.0

#   print(' ')
#   print('specular shots found:')
#   print(n_spec_shots)
#   print(' ')
        spec_shot_out.append(
        	{
        		'mframe': ii,
        		'pulse': jj,
        		'n_photons': n_photons	
        	}        
        )


#
#   return
#
  return n_spec_shots, spec_shot_out
