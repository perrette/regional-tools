// particle tracking
// aimed at comparing GSL+RKF45 scheme to a more naive approach 
// #include "dbg_extensions.hh"
#include <iostream>
#include "streamline.hh"

using namespace std;

GriddedData::GriddedData(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z, double *my_u, double *my_v, int *my_m) 
    : elevation(my_Mx, my_My), vel_x(my_Mx, my_My), vel_y(my_Mx, my_My) , mask(my_Mx, my_My)
{
    x  = my_x;
    Mx = my_Mx;

    y  = my_y;
    My = my_My;

    z   = my_z;
    u   = my_u;
    v   = my_v;
    m   = my_m;

    // We assume that the grid is uniform.
    dx = x[1] - x[0];
    dy = y[1] - y[0];
    spacing = dx > dy ? dx : dy;

    one_over_dx = 1.0 / dx;
    one_over_dy = 1.0 / dy;

    if (z != NULL) {
        elevation.wrap(z);
    }
    if (u != NULL) {
        vel_x.wrap(u);
    }
    if (v != NULL) {
        vel_y.wrap(v);
    }
    if (m != NULL) {
        mask.wrap(m);
    }
}
//
// GriddedData::~GriddedData() : {
// }

inline int GriddedData::find_cell(const double position[], int &i, int &j) {

    i = floor((position[0] - x[0]) * one_over_dx);
    j = floor((position[1] - y[0]) * one_over_dy);

  // bail if we ended up outside the grid
  if (i < 0 || i + 1 > Mx - 1 || j < 0 || j + 1 > My - 1) {
    i = j = -1;
    return 1;
  }

  return 0;
}

template<typename T>
void get_corner_values(Array2D<T> &arr, int &i, int &j, T &A, T &B, T &C, T &D) {
    A = arr(i,     j);
    B = arr(i,     j + 1);
    C = arr(i + 1, j + 1);
    D = arr(i + 1, j);
}

void GriddedData::evaluate(const double position[], double *elevation, double *gradient, double *vel) {

    int ierr, i, j;
    double A, B, C, D;

    ierr = this->find_cell(position, i, j);

    // Pretend that outside the grid the surface is perfectly flat, the elevation
    // of the sea level (0) and ice-free, with zero velocity
    if (ierr != 0) {

        if (elevation != NULL)
            *elevation = 0;

        if (gradient != NULL)
            gradient[0] = gradient[1] = 0;

        if (vel != NULL)
            vel[0] = vel[1] = 0;

        return;
    }

    double delta_x = position[0] - x[i],
           delta_y = position[1] - y[j];

    // 
    double alpha = one_over_dx * delta_x,
           beta  = one_over_dy * delta_y;

    // interpolation weights
    double aa = (1 - alpha) * (1 - beta),
           bb = (1 - alpha) *      beta,
           cc = alpha       *      beta,
           dd = alpha       * (1 - beta);

    if (vel != NULL) {

        get_corner_values(this->vel_x, i, j, A, B, C, D);
        vel[0] = aa * A + bb * B + cc * C + dd * D ;

        get_corner_values(this->vel_y, i, j, A, B, C, D);
        vel[1] = aa * A + bb * B + cc * C + dd * D ;
    }

    if (elevation != NULL or gradient != NULL) {
        get_corner_values(this->elevation, i, j, A, B, C, D);
    }

    if (elevation != NULL) {
        *elevation = aa * A + bb * B + cc * C + dd * D ;
    }

    // the gradient
    if (gradient != NULL) {
        double gamma = one_over_dx * one_over_dy * (A + C - B - D);

        gradient[0] = (D - A) * one_over_dx + delta_x * gamma;
        gradient[1] = (B - A) * one_over_dy + delta_y * gamma;
    }

} // end of GriddedData::evaluate()

inline void GriddedData::get_vector(const double position[], double f[]) {
    if (has_vectors()) {
        // use the vector field when provided
        evaluate(position, NULL, NULL, f);
    } else {
        evaluate(position, NULL, f, NULL);
        f[0] = -f[0];
        f[1] = -f[1];
    }
}


        // void get_vector(const double position[], double f[]);
inline int NaiveStepping::apply_step(double *position, const double &step_length) {

    double vector[2], vmag;

    dem->get_vector(position, vector);
    vmag = sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

    if (vmag == 0) return GSL_SUCCESS - 1;

    double scal = step_length/vmag;

    position[0] = position[0] + scal*vector[0];
    position[1] = position[1] + scal*vector[1];

    return GSL_SUCCESS;
}


inline int GSLStepping::apply_step(double *position, const double &step_length) {

    double vector[2], vmag, err[2];
    dem->get_vector(position, vector);
    vmag = sqrt(vector[0]*vector[0] + vector[1]*vector[1]);

    if (vmag == 0) return GSL_SUCCESS-1;

    double status = gsl_odeiv_step_apply(step,
            0,         // starting time (irrelevant)
            step_length / vmag, // step size (units of time)
            // step_length / gradient_magnitude, // step size (units of time)
            position, err, NULL, NULL, &system);

    // if (status != GSL_SUCCESS) {
    //     printf ("error, return value=%d\n", status);
    // }
    return status;
}

template <class T>
int particle_tracking_intern(double position[], T* stepping, 
        const double &step_length, const int &n_max, double *xl, double *yl, int *il, int *jl, int &length) {

    int i, j, status;

    GriddedData* dem = stepping->dem;

    for (int step_counter = 0; step_counter < n_max; ++step_counter) {
        status = dem->find_cell(position, i, j);

        // stop if outside domain
        if (status != 0)
            break;

        // stop if bad mask value
        if (dem->is_ice_free(i, j)) 
            break;

        xl[step_counter] = position[0];
        yl[step_counter] = position[1];
        il[step_counter] = i;
        jl[step_counter] = j;
        length++; // = step_counter + 1;

        status = stepping->apply_step(position, step_length);

        if (status != GSL_SUCCESS) {
            // printf ("error, return value=%d\n", status);
            break;
        };
    } 
    return 0;
}

int particle_tracking(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, 
        double x0, double y0,
        double *xl , double *yl, int *il, int *jl, int &nl, int n_max=-1, bool use_gsl=true) { // parameters {

    GriddedData dem(x, Mx, y, My, z, u, v, mask);

    // default value for maximum streamline length
    if (n_max <= 0) n_max = dem.Mx + dem.My;

    if ((z != NULL) and (u !=NULL or v !=NULL) ) {
        printf("streamline: provide either z or u,v. Stop.\n");
        return 1; // provide either z or u
    } else if ((u !=NULL and v ==NULL) or (u==NULL and v!=NULL) ) {
        printf("streamline: provide both u AND v (or z). Stop.\n");
        return 1; // provide either z or u
    }

    printf("streamline: use velocity field? %d \n", dem.has_vectors());

    // int status;
    int steps_per_cell = 2;

    double step_length = dem.spacing / steps_per_cell, position[2];

    position[0] = x0;
    position[1] = y0;

    // void *stepping;

    if (use_gsl) {
        // stepping = &gsl_stepping;
        GSLStepping gsl_stepping(&dem);
        particle_tracking_intern(position, &gsl_stepping, step_length, n_max, xl, yl, il, jl, nl);
    } else {
        // stepping = &naive_stepping;
        NaiveStepping naive_stepping(&dem);
        particle_tracking_intern(position, &naive_stepping, step_length, n_max, xl, yl, il, jl, nl);
    };

    return 0;
}
