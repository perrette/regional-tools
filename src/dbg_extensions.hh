// Extensions to be able to use input vectors instead of simply surface elevation
#include "DEM.hh"

// Extend DEM class to use surface velocity instead of gradient
class SurfaceData: public DEM {
public:
  SurfaceData(double *x, int Mx, double *y, int My, double *z, double *u, double *v);
  ~SurfaceData() {}

  void evaluate_vel(const double *position, double *vel);

  double *u, *v;

protected:
  Array2D<double> vel_x, vel_y;
};


SurfaceData::SurfaceData(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z, double *my_u, double *my_v)
  : DEM(my_x, my_Mx, my_y, my_My, my_z), vel_x(my_Mx, my_My), vel_y(my_Mx, my_My) {

  u   = my_u;
  v   = my_v;
  vel_x.wrap(u);
  vel_y.wrap(v);
}

// like evaluate, but return velocity at a given position instead of elevation gradient
void SurfaceData::evaluate_vel(const double *position, double *vel) {

  int ierr, i, j;
  double A, B, C, D;

  ierr = this->find_cell(position, i, j);

  // Pretend that outside the grid the surface is perfectly flat, the elevation
  // of the sea level (0) and ice-free.
  if (ierr != 0) {

    if (vel != NULL)
      vel[0] = vel[1] = 0;

    return;
  }

  double
    delta_x = position[0] - this->x[i],
    delta_y = position[1] - this->y[j];

  // 
  double
      alpha = this->one_over_dx * delta_x,
      beta  = this->one_over_dy * delta_y;

  A = this->vel_x(i,     j);
  B = this->vel_x(i,     j + 1);
  C = this->vel_x(i + 1, j + 1);
  D = this->vel_x(i + 1, j);
  
  vel[0] =     ( (1 - alpha) * (1 - beta) * A +
                 (1 - alpha) *      beta  * B +
                 alpha       *      beta  * C +
                 alpha       * (1 - beta) * D );
  
  A = this->vel_y(i,     j);
  B = this->vel_y(i,     j + 1);
  C = this->vel_y(i + 1, j + 1);
  D = this->vel_y(i + 1, j);
  
  vel[1] =     ( (1 - alpha) * (1 - beta) * A +
                 (1 - alpha) *      beta  * B +
                 alpha       *      beta  * C +
                 alpha       * (1 - beta) * D );

} // end of SurfaceData::evaluate_vel()


static int right_hand_side_vel(double t, const double y[], double f[], void* params) {

  ((SurfaceData*)params)->evaluate_vel(y, f);

  return GSL_SUCCESS;

}

