# fine_track
Python interpretation of ICESat-2 sea ice algorithm fine tracker (part of atlas_l3a_si)

There are two main routines included:
1) fine_trap.py contains the main fine_track routine as well as additional subroutines including the specular shot filter and plotting routines
2) run_atl07_segs.py contains a routine to compute a range of sea ice segment heights using fine_trap.py

run_atl07_segs.py requires and ATL07 and corresponding ATL03 to process. Additionally, the expected waveform table and calibration data from ATL03 is needed. This is supplied by running finetrack_prep and linking the resulting atl03_calibrations_and_wf_tables.h5.

run_atl07_segs.py has optional input parameters:
- t0: start time for processing (delta time)
- t1: end time for processing (delta time)
- b: beam selection (default all)
- spec_filter: perform specular shot filter
- debug_plots: provide plots for each finetrack segment
- csv_file: plot summary csv file of finetrack height segments

within run_atl07_segs.py, the following file paths and names are hardcoded in:
- ATL03
- ATL07
- ATL03 calibrations file

If no times are entered, the first second of ATL07 data will be processed. If one of t0/t1 is missing, the start/end time is used.

Directory sample_output contains sample output.

A folder with test data is available on the iceproc6 machine:
/home/jwimert/finetrack_python

- ATL03_20190317122653_12070212_006_02.h5
- ATL07-02_20190317111321_12070201_006_02.h5
- atl03_calibrations_and_wf_tables.h5


Sample output:


![run_atl07_segs_summary](https://github.com/jwimert/fine_track/assets/61554992/c48c43f5-26c1-4e4a-b25a-c3e205050202)

