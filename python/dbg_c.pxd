# -*- mode: python -*-

cdef extern from "../src/dbg.hh":
    # bint upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bint output)
    # bint upstream_area(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, bint output)
    bint upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bint output, bint use_gsl, double elevation_step, int path_length)
    bint upstream_area(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, bint output, bint use_gsl, double elevation_step, int path_length)
    bint accumulated_flow(double *x, int Mx, double *y, int My, double *z, double *mask, int n_samples)
    bint initialize_mask(int Mx, int My, double *thickness, int *mask)

