# PISM Drainage basin delineation tool

Experimental extensions to the package originally developped for PISM.
See fork information on github for the original code and appropriate documentation.

These changes include:
- possibility of passing a vector field (`upstream_area`) instead of originally just surface elevation (`upslope_area`).
- possibility of using or not the GSL RKF45 stepping algorithm via a new `use_gsl` parameter, and of specifying a few
other parameters such as `elevation_step` and `path_length`.
