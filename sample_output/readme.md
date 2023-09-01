Sample output for run_atl07_segs.py

- ATL03_20190317122653_12070212_006_02.h5
- ATL07-02_20190317111321_12070201_006_02.h5
- sea_ice_waveforms_202308224.h5
- t0 38060910.0 
- t1 38060915.0


python run_atl07_segs.py -t0 38060910.0 -t1 38060915.0 -spec_filter 'y'

 
==============================================
 
run_atl07_segs.py
 
 
Read ATL07 sea ice segments and recreate
sea ice heights from ATL03 photon heights
using fine_trap.py routine
 
 
Jesse Wimert
Last update: August 2023
 
==============================================
 
 
 
Input ATL03:
 
/Users/jwimert/Desktop/func_test_antarc_s2/ATL03_20190317122653_12070212_006_02.h5
ATL03 time range:
38060812.91238648 38061273.680097215
 
Input ATL07:
/Users/jwimert/Desktop/func_test_antarc_s2/ATL07-02_20190317111321_12070201_006_02.h5
ATL07 time range:
38060852.85098807 38061027.13915148
 
 
Requested Time Range:
38060910.0
38060915.0
 
number of segments to compute:
4871
 
 
Specular return filter?:
y
 
 
Output Files:
default
 
 
 
Perform specular shot filter
 
Number of specular shots found:
256
Number of photons filtered
4599
 
31256 ERROR RUNNING fine_track


