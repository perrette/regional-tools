#ifndef _DBG_H_
#define _DBG_H_

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

int initialize_mask(int Mx, int My, double *thickness, int* mask);

int upslope_area(double *x, int Mx, double *y, int My, double *z, int *mask, bool output, bool use_gsl, double elevation_step, int path_length);

int upstream_area(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, bool output, bool use_gsl, double elevation_step, int path_length);

int accumulated_flow(double *x, int Mx, double *y, int My, double *z, double *my_mask, int n_samples);

int particle_tracking(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, 
        double x0, double y0,
        double *xl , double *yl, int *il, int *jl, int &nl, int n_max, bool use_gsl);

#endif /* _DBG_H_ */
