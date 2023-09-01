##!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_atl07_segs.py
Jesse Wimert (06/2023)

This routine reads photon height data from an ATL03 file and 
computes sea ice heights using a python fine_track routine.
An ATL07 file is read to assist in collecting photons and applying
geophysical corrections. After computing heights, the new heights
are compared to the corresponding ATL07 heights.

INPUTS:
        -t0 delta time start time to process sea ice segments
        -t1 delta time end time to process sea ice segments
         (default, first 100 segments of ATL07 are processed)
        -spec_filter filter photons which are part of a specular shot (y or n)
         (default, y)
         (warning: slow)
        -output output option (default, slim, debug)
         (default, default)

OUTPUTS:
        -summary.txt summary file of sea ice height results
        
        -summary.png summary plot comparing sea ice heights
        
        -fine_track_summary_xxxx.png (optional, '-output debug') individual
         fine track summary plot for each sea ice segment

Steps:

Read ancillary data
Read data from ATL07
Read photon heights and data from ATL03
Collect photons within finetrack window
Use routines from finetrack_subroutines.py to compute finetrack segment height

Notes:
- no skew correction computed or applied 
- output files depend on output option

"""
#
import os
import argparse
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from fine_trap import fine_track
from fine_trap import spec_shot_filter
# from finetrack_subroutines import fine_track



##
## Function to parse command line arguments
##
 
def parse():
    """ Parse command line input arguments. """

    parser = argparse.ArgumentParser(description='Compute ATL07 sea ice heights')
    parser.add_argument('-t0', action='store', default='0.0', dest='t0',
                        help='start time of processing',
                        type=str)
    parser.add_argument('-t1', action='store', default='0.0', dest='t1',
                        help='plot upper bound',
                        type=str)
    parser.add_argument('-spec_filter', action='store', default='n', 
                        dest='spec_filter',
                        help="Perform specular shot filter before processing \
                        ['y' or 'n', default 'y']", type=str)
    parser.add_argument('-output', action='store', default='default', 
                        dest='output_flag',
                        help="Select group of output file(s) \
                        ['default', 'slim', or 'debug', default 'default']", type=str)
    inps = parser.parse_args()
    return inps

##
## Print welcome message
##

os.system('clear')
print(' ')
print('==============================================')
print(' ')
print('run_atl07_segs.py')
print(' ')
print(' ')
print('Read ATL07 sea ice segments and recreate')
print('sea ice heights from ATL03 photon heights')
print('using fine_trap.py routine')
print(' ')
print(' ')
print('Jesse Wimert')
print('Last update: August 2023')
print(' ')
print('==============================================')
print(' ')
print(' ')


##
## Read command line input parameters
##

inps = parse()
t0 = float(inps.t0)
t1 = float(inps.t1)
spec_filter = inps.spec_filter
output_flag = inps.output_flag


##
## Input Files:
##

#
# Set input ATL03 and ATL07 files
#
ATL03_file = '/Users/jwimert/Desktop/func_test_antarc_s2/ATL03_20190317122653_12070212_006_02.h5'
ATL07_file = '/Users/jwimert/Desktop/func_test_antarc_s2/ATL07-02_20190317111321_12070201_006_02.h5'

#
# Set input waveform table
#
wf_table_file = '/Users/jwimert/Desktop/func_test_antarc_s2/run_atl07/latest/sea_ice_waveforms_202308224.h5'

#
# Load Cal19 tables
#
cal19_in = np.load('cal19.npy')
cal19_dead_time = np.load('cal19_dead_time.npy')
cal19_strength = np.load('cal19_strength.npy')
cal19_width = np.load('cal19_width.npy')


##
## Set Constants:
##

#
# Dead time (ns)
# Note: dead time for center strong and weak beams
# For other beams/spots, need to read from ATL03
#
dead_time_strong = 2.98125
dead_time_weak = 3.0525

#
# Set histogram bins
# (lb/ub_bin : CENTER of lowest/highest bin)
#
bin_size = 0.025
lb_bin = -2.0
ub_bin = 3.5

# bins_center = np.arange(lb_bin, ub_bin+bin_size, bin_size)
# bins_edges = np.arange(lb_bin - bin_size*0.5, ub_bin + bin_size, bin_size)

#
# Set padding time for collecting photons (s)
# 
# Note:
# 
# Photons collected by creating a window around the delta_time for the ATL07 segment
# window = atl07_time +/- (atl07_seg_length + time_pad) / 2
#
time_pad = 0.5 


##
## Read Expected Waveform Table
##

#
# Open file
#
hf = h5py.File(wf_table_file, 'r')

#
# Read waveform table
#
wf_table = hf['wf_table'][:][:][:]

#
# Read ancillary axes arrays
#
wf_bins = hf['binz'][:]
wf_sd = hf['sdz'][:]
wf_mn = hf['mnz'][:]

hf.close()


##
## Read ATL07 file
##

#
# Open file
#
hf = h5py.File(ATL07_file, 'r')

#
# Set processing to center strong beam
#
orientation = hf.get('/orbit_info/sc_orient')[0]

if orientation == 0:
    beam = 'gt2l'
elif orientation == 1:
    beam = 'gt2r'

#
# Set beam
#
gtx = hf.get(beam)

#
# Read data
#
ATL07_delta_time = gtx.get('sea_ice_segments/delta_time')[:]

heightx = hf.get(beam+'/sea_ice_segments/heights')
ATL07_h_uncorr = heightx.get('height_segment_height_uncorr')[:]
ATL07_seg_length = heightx.get('height_segment_length_seg')[:]
ATL07_coarse_mn = gtx.get('sea_ice_segments/stats/height_coarse_mn')[:]

ATL07_h_corr = heightx.get('height_segment_height')[:]
ATL07_w_gauss = heightx.get('height_segment_w_gaussian')[:]
ATL07_htcorr_skew = heightx.get('height_segment_htcorr_skew')[:]

ATL07_n_pulse = heightx.get('height_segment_n_pulse_seg')[:]
ATL07_n_pulse_used = heightx.get('height_segment_n_pulse_seg_used')[:]
ATL07_quality_flag = heightx.get('height_segment_fit_quality_flag')[:]

ATL07_n_photons_actual = gtx.get('sea_ice_segments/stats/n_photons_actual')[:]
ATL07_n_photons_used = gtx.get('sea_ice_segments/stats/n_photons_used')[:]
ATL07_photon_rate = gtx.get('sea_ice_segments/stats/photon_rate')[:]

ATL07_fpb_corr = gtx.get('sea_ice_segments/stats/fpb_corr')[:]
# ATL07_fpb_width = gtx.get('sea_ice_segments/stats/fpb_corr_width')[:]
# ATL07_fpb_strength = gtx.get('sea_ice_segments/stats/fpb_strength')[:]
# ATL07_fpb_avg_dt = gtx.get('sea_ice_segments/stats/fpb_avg_dt')[:]

ATL07_lpe = gtx.get('sea_ice_segments/geophysical/height_segment_lpe')[:]
ATL07_mss = gtx.get('sea_ice_segments/geophysical/height_segment_mss')[:]
ATL07_dynib = gtx.get('sea_ice_segments/geophysical/height_segment_dynib')[:]
ATL07_oc = gtx.get('sea_ice_segments/geophysical/height_segment_ocean')[:]

ATL07_histogram = gtx.get('sea_ice_segments/stats/hist_photon_heights')[:]


hf.close()

#
# Filter invalid data
#
ATL07_h_uncorr[ATL07_h_corr > 3.0e+37] = 0.0
ATL07_h_corr[ATL07_h_corr > 3.0e+37] = 0.0
ATL07_h_corr[ATL07_h_corr < -3.0e+37] = 0.0
ATL07_w_gauss[ATL07_w_gauss > 3.0e+37] = 0.0
ATL07_htcorr_skew[ATL07_htcorr_skew > 3.0e+37] = 0.0
ATL07_htcorr_skew[ATL07_htcorr_skew < -3.0e+37] = 0.0


##
## Read ATL03 file
##

#
# Open file
#
hf = h5py.File(ATL03_file, 'r')

#
# Read data
#
ATL03_ph_delta_time = hf[beam + '/heights/delta_time'][:]
ATL03_ph_height = hf[beam + '/heights/h_ph'][:]

#
# Pulse number needed to find specular returns
#
ATL03_mframe = hf[beam + '/heights/pce_mframe_cnt'][:]
ATL03_pulse = hf[beam + '/heights/ph_id_pulse'][:]
ATL03_ph_conf = hf[beam + '/heights/signal_conf_ph'][:]

hf.close()


##
## Compute ATL07 indicies for requested time bound
##



if (t0 == 0.0) & (t1 != 0.0):
  t0 = ATL07_delta_time[0]

if (t0 != 0.0) & (t1 == 0.0):
  t1 = ATL07_delta_time[ATL07_delta_time.size - 1]

if (t0 == 0.0) & (t1 == 0.0):
  t0 = ATL07_delta_time[0]
  t1 = ATL07_delta_time[0] + 1.0

if (t0 != 0.0):
  for i in np.arange(0, ATL07_delta_time.size):
    if (ATL07_delta_time[i] > t0):
      atl07_i0 = i
      break

if (t1 != 0.0):
  for i in np.arange(0, ATL07_delta_time.size):
    if (ATL07_delta_time[i] > t1):
      atl07_i1 = i
      break

##
## Print summary of input files and sleceted options
##

print(' ')
print('Input ATL03:')
print(' ')
print(ATL03_file)
print('ATL03 time range:')
print(ATL03_ph_delta_time[0],ATL03_ph_delta_time[ATL03_ph_delta_time.size-1])
print(' ')
print('Input ATL07:')
print(ATL07_file)
print('ATL07 time range:')
print(ATL07_delta_time[0],ATL07_delta_time[ATL07_delta_time.size-1])
print(' ')
print(' ')
print('Requested Time Range:')
print(t0)
print(t1)
print(' ')
print('number of segments to compute:')
print(atl07_i1 - atl07_i0 - 1)
print(' ')
print(' ')
print('Specular return filter?:')
print(spec_filter)
print(' ')
print(' ')
print('Output Files:')
print(output_flag)
print(' ')
print(' ')



ATL07_h_uncorr[ATL07_h_corr > 3.0e+37] = 0.0
ATL07_h_corr[ATL07_h_corr > 3.0e+37] = 0.0
ATL07_h_corr[ATL07_h_corr < -3.0e+37] = 0.0
ATL07_w_gauss[ATL07_w_gauss > 3.0e+37] = 0.0
ATL07_htcorr_skew[ATL07_htcorr_skew > 3.0e+37] = 0.0
ATL07_htcorr_skew[ATL07_htcorr_skew < -3.0e+37] = 0.0

##
## Specular filter
##

if (spec_filter == 'y'):

  print(' ')
  print('Perform specular shot filter')
  print(' ')
  
  i0 = np.searchsorted(ATL03_ph_delta_time,t0)
  i1 = np.searchsorted(ATL03_ph_delta_time,t1)

  n_phots_spec_full = (ATL03_ph_height[i0:i1] > -9000.0).sum()
 
#
# Call spec_shot_filter
#
  
  n_spec_shots, spec_shot_out = \
    spec_shot_filter(i0, i1, ATL03_ph_height, ATL03_mframe, ATL03_pulse, ATL03_ph_conf) 

  n_phots_spec_filt = (ATL03_ph_height[i0:i1] > -9000.0).sum()

  print('Number of specular shots found:')
  print(n_spec_shots)
  print('Number of photons filtered')
  print(n_phots_spec_full - n_phots_spec_filt)
  print(' ')

  output_spec_df = pd.DataFrame(spec_shot_out)    
  output_spec_df.to_csv('specular_shots.csv', index=False)

##
## Prepare for looping through segments
##

#
# Initialize arrays for storing surface heights and gaussian widths
#
output_temp = []

#
# Loop through segments and compute heights
#
#
for i_test in np.arange(atl07_i0, atl07_i1+1):

#
# Select ATL07 segment, search for start/stop time of segment
# Assume orbital speed of 6.9 km/s
#
    t1 = ATL07_delta_time[i_test] - (ATL07_seg_length[i_test] + time_pad) / 2.0 / 6900.0
    t2 = ATL07_delta_time[i_test] + (ATL07_seg_length[i_test] + time_pad) / 2.0 / 6900.0

#
# Collect index of selected photons
# --Photons within time window and lb/ub_win height 
# Computing uncorrected heights, adjusting lb/ub_win using ATL07_coarse_mn
#
    photons_loc = np.where( (ATL03_ph_delta_time > t1) & (ATL03_ph_delta_time < t2) & \
      ( ATL03_ph_height > (ATL07_coarse_mn[i_test] + lb_bin) ) & \
      ( ATL03_ph_height < (ATL07_coarse_mn[i_test] + ub_bin)) )

#
# Count shots and photons
#
    running_pulse = ((ATL03_mframe[photons_loc] - ATL03_mframe[photons_loc][0])) * 200 + \
      ATL03_pulse[photons_loc]
    n_shots = np.max(running_pulse) - np.min(running_pulse) + 1
    n_photons = ATL03_ph_delta_time[photons_loc].size
    
#
# Call fine_track to compute height and gaussian width
#
    try:
        h_surf_temp, w_gauss_temp, fpb_corr = \
          fine_track(ATL03_ph_height[photons_loc], ATL07_coarse_mn[i_test], \
          n_photons, n_shots, bin_size, lb_bin, ub_bin, \
          wf_table, wf_bins, wf_sd, wf_mn, \
          cal19_in, cal19_width, cal19_strength, cal19_dead_time, dead_time_strong)

#
# Compute corrected surface height
#
        h_surf_corr = h_surf_temp - \
          ATL07_lpe[i_test] - ATL07_mss[i_test] - ATL07_dynib[i_test] - ATL07_oc[i_test]

#
#         diff = (ATL07_h_corr[i_test] - ATL07_htcorr_skew[i_test]) - h_surf_corr
#         print(i_test, diff)

#
# Save output data
#
        
        output_temp.append(
        	{
        		'seg_count': i_test,
        		'atl07_dt': ATL07_delta_time[i_test],
        		'h_surf': h_surf_corr,
        		'w_gauss': w_gauss_temp,
        		'fpb_corr': fpb_corr,
        		'atl07_h': ATL07_h_corr[i_test] - ATL07_htcorr_skew[i_test],
        		'atl07_w': ATL07_w_gauss[i_test],
        		'atl07_fpb':ATL07_fpb_corr[i_test]	
        	}        
        )
       
    except:
        print(i_test, 'ERROR RUNNING fine_track')

#
# Save output as dataframe and write file
#

output_df = pd.DataFrame(output_temp)    
output_df.to_csv('run_atl07_segs_summary.csv', index=False)


fig, (ax0, ax1, ax2, ax3) = plt.subplots(4, figsize=(12, 20))
# , height_ratios=[4, 3, 4, 6])
fig.suptitle('Fine Track Summary')
ax0.scatter(output_df['atl07_dt'], output_df['atl07_h'], label='ATL07 height')
ax0.scatter(output_df['atl07_dt'], output_df['h_surf'], label='computed height')
ax0.legend()
ax0.title.set_text('Sea ice segment heights')
ax1.scatter(output_df['atl07_dt'], output_df['atl07_h'] - output_df['h_surf'])
ax1.title.set_text('Sea ice segment height diff')
ax2.scatter(output_df['atl07_dt'], output_df['atl07_w'], label='ATL07_w_gauss')
ax2.scatter(output_df['atl07_dt'], output_df['w_gauss'], label='computed_w_gauss')
ax2.legend()
ax2.title.set_text('Sea ice segment w_gauss')
ax3.scatter(output_df['atl07_dt'], output_df['atl07_fpb'], label='ATL07 fpb')
ax3.scatter(output_df['atl07_dt'], output_df['fpb_corr'], label='computed fpb')
ax3.legend()
ax3.title.set_text('Sea ice segment first photon bias correction')

outputfile = 'run_atl07_segs_summary.png'
fig.savefig(outputfile, dpi=fig.dpi)
plt.close('all')
