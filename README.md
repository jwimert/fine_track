# fine_track
Python interpretation of ICESat-2 sea ice algorithm fine tracker (part of atlas_l3a_si)

There are two main routines included:
1) fine_trap.py contains the main fine_track routine as well as additional subroutines including the specular shot filter
2) run_atl07_segs.py contains a routine to compute a range of sea ice segment heights using fine_trap.py

run_atl07_segs.py requires and ATL07 and corresponding ATL03 to process. Additionally, the expected waveform table and cal19 data is also needed. The expected waveform table can be created using the create_wf_table with the desired ATL03. Cal19 data can be read and stored from an ATL03.

run_atl07_segs.py has optional input parameters:
- t0: start time for processing (delta time)
- t1: end time for processing (delta time)
- spec_filter: perform specular shot filter

within run_atl07_segs.py, the following file paths and names are hardcoded in:
- ATL03
- ATL07
- expected waveform table
- cal19 data

If no times are entered, the first second of ATL07 data will be processed. If one of t0/t1 is missing, the start/end time is used.

A folder with test data is available on the iceproc6 machine:
/home/jwimert/finetrack_python

ATL03_20190317122653_12070212_006_02.h5
ATL07-02_20190317111321_12070201_006_02.h5
cal19_dead_time.npy
cal19.npy
cal19_strength.npy
cal19_width.npy
sea_ice_waveforms_202308224.h5
