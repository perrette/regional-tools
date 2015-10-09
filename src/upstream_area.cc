#include "dbg_internal.hh"
#include "DEM.hh"
#include <map>

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


static int streamline_vel(dbg_context ctx, int i_start, int j_start,
                      Array2D<int> &old_mask, Array2D<int> &new_mask) {
  SurfaceData *dem = (SurfaceData*)ctx.system.params;

  int mask_counter = 0,
    n_max = (dem->Mx + dem->My) * ctx.steps_per_cell,
    i = i_start, j = j_start,
    i_old, j_old,
    status;

  double step_length = dem->spacing / ctx.steps_per_cell,
    position[2],
    err[2],
    elevation,
    vel[2],
    vel_magnitude;

  std::map<int,int> values;

  int mask_value = old_mask(i_start, j_start);

  // stop if the current cell already has a value assigned
  if (mask_value > 0 || mask_value == ICE_FREE)
    return 0;

  position[0] = dem->x[i_start] + dem->dx * 0.5;
  position[1] = dem->y[j_start] + dem->dy * 0.5;

  // get elevation
  dem->evaluate(position, &elevation, NULL);
  
  // if there is no ice or we're below the minimum elevation, we're done
  if (elevation < ctx.min_elevation)
    return 0;

  // if we're above the maximum elevation, wait.
  if (elevation > ctx.max_elevation)
    return 1;

  for (int step_counter = 0; step_counter < n_max; ++step_counter) {

    i_old = i; j_old = j;
    status = dem->find_cell(position, i, j);

    if (status != 0)
      break;

    mask_value = old_mask(i, j);

    if (mask_value == ICE_FREE)   // ice-free
      break;

    if ((i != i_old || j != j_old) && (mask_value > 0)) {
      values[mask_value]++;
      mask_counter++;

      if (mask_counter == ctx.path_length)
        break;
    }

    // dem->evaluate(position, NULL, gradient);
    dem->evaluate_vel(position, vel);

    // gradient_magnitude = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1]);
    vel_magnitude = sqrt(vel[0]*vel[0] + vel[1]*vel[1]);

    // take a step
    
    // if (vel_magnitude > 0) {
    //     double scal = step_length/vel_magnitude; 
    //     position[0] += vel[0] * scal;
    //     position[1] += vel[1] * scal;
    // }

    status = gsl_odeiv_step_apply(ctx.step,
                                  0,         // starting time (irrelevant)
                                  step_length / vel_magnitude, // step size (units of time)
                                  // step_length / gradient_magnitude, // step size (units of time)
                                  position, err, NULL, NULL, &ctx.system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }

  } // time-stepping loop (step_counter)

  // Find the mask value that appears more often than others.
  int most_frequent_mask_value = NO_VALUE;
  int number_of_occurences = 0;

  std::map<int,int>::iterator k;
  for (k = values.begin(); k != values.end(); ++k) {
    if (k->second > number_of_occurences) {
      number_of_occurences = k->second;
      most_frequent_mask_value = k->first;
    }
  }

  new_mask(i_start, j_start) = most_frequent_mask_value;

  if (most_frequent_mask_value == NO_VALUE)
    return 1;

  return 0;
}


int upstream_area(double *x, int Mx, double *y, int My, double *z, double *u, double *v, int *mask, bool output) {
  int remaining = 0;
  const double elevation_step = 10;

  SurfaceData dem(x, Mx, y, My, z, u, v);

  Array2D<int> my_mask(Mx, My), new_mask(Mx, My);
  my_mask.wrap(mask);
  if (new_mask.allocate() != 0)
    return -1;

#pragma omp parallel default(shared)
  {
    gsl_odeiv_system system = {right_hand_side_vel, NULL, 2, &dem};
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);
    dbg_context ctx = {system, step, 2, 5, 0, elevation_step};

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
        new_mask(i, j) = my_mask(i, j);

    do {
#pragma omp for schedule(dynamic)
      for (int j = 0; j < My; j++)
        for (int i = 0; i < Mx; i++)
          my_mask(i, j) = new_mask(i, j);

#pragma omp single
      remaining = 0;

#pragma omp for schedule(dynamic) reduction(+:remaining)
      for (int j = 0; j < My; j++) { // Note: traverse in the optimal order
        for (int i = 0; i < Mx; i++) {
          remaining += streamline_vel(ctx, i, j, my_mask, new_mask);
        }
      }

      ctx.min_elevation = ctx.max_elevation;
      ctx.max_elevation += elevation_step;

    } while (remaining > 0);

    gsl_odeiv_step_free(step);

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
	my_mask(i, j) = new_mask(i, j);

  } // end of the parallel block

  return 0;
}