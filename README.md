# fine_track
Python interpretation of ICESat-2 sea ice algorithm fine tracker (part of atlas_l3a_si)

There are two main routines included:
1) fine_trap.py contains the main fine_track routine as well as additional subroutines including the specular shot filter
2) run_atl07_segs.py contains a routine to compute a range of sea ice segment heights using fine_trap.py
