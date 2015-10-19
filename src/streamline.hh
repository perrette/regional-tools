#ifndef _STREAMLINE_H_
#define _STREAMLINE_H_

// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <cmath>
#include "Array2D.hh"
#include "dbg.hh"

class GriddedData {
    public:
        int Mx, My;         // dimensions
        double *x, *y;      // grid
        double *z, *u, *v;  // data on the grid
        int *m;             // integer mask

        GriddedData(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *m);
        ~GriddedData() {};

        int find_cell(const double position[], int &i, int &j);
        void evaluate(const double position[], double *elevation, double *gradient, double *vel);

        bool has_vectors() { return (u != NULL and v != NULL); };  // is the vector field defined?
        bool has_mask() { return (m != NULL); };  // is the vector field defined?
        bool is_ice_free(int &i, int &j) { return has_mask() && mask(i,j) == ICE_FREE; };
        void get_vector(const double position[], double f[]);

        double one_over_dx, one_over_dy;
        double spacing, dx, dy;

    protected:
      // void get_corner_values(int i, int j,
      //                        double &A, double &B, double &C, double &D);
        
        Array2D<double> elevation;
        Array2D<double> vel_x;
        Array2D<double> vel_y;
        Array2D<int> mask;
};


// one step...
class NaiveStepping {
    public:
        GriddedData *dem;

        NaiveStepping(GriddedData *my_dem) : dem(my_dem) {};

        int apply_step(double *position, const double &step_length);
    // private:
        // static int right_hand_side(double t, const double y[], double f[], void* params) ; // use the same form as GSL
};

static int right_hand_side(double t, const double y[], double f[], void* params) {
    ( (GriddedData*)params )->get_vector(y, f);
    return GSL_SUCCESS;
};

class GSLStepping : public NaiveStepping {
    public:
        GSLStepping(GriddedData* dem) : NaiveStepping(dem) {
            system = {right_hand_side, NULL, 2, dem};
            step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);
            // dbg_context ctx = {system, step, 2, path_length, 0, elevation_step};
        };
        ~GSLStepping() {
            gsl_odeiv_step_free(step);
        }
        int apply_step(double *position, const double &step_length);
    private:
        gsl_odeiv_system system;
        gsl_odeiv_step *step;
};

#endif
