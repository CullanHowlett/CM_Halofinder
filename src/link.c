// This contains all the extra routines (other than general library functions) required to produce the linking length
// lookup table associated with a linking length that varies as a function of the particle's distance from the the user-
// defined origin (note that in fact we only give each cell it's own linking length as the variation within the cell is
// negligible and we have, on average, one particle per cell).
//
// The steps required are as follows:
//   - The routine Create_link_tab is called at the beginning of the main, serial fof routine. This has the final result
//     of producing a lookup table of 'Nlink' linking lengths between the minimum and maximum distances from the origin
//     on each processor. 
//   - To do this we first  calculate the maximum redshift needed for the lookup table by minimization.
//   - We then calculate the linking lengths for Nlink redshifts between 0 and zmax.
//   - This then means that for each cell we now just find the appropriate element in the lookup table based on the processors 
//     distance from the origin and the location of the cell in the processor.

#include "vars.h"

void Create_Link(void) {
 
  int i, status, iter = 0, iter_max = 2000;
  double zmax = 1.0, lower = 0.0, upper = 3.0;
  const double link_sim = 0.2;

  // We first need to find the redshift corresponding to the maximum distance 
  // calculated for each processor, we do this by minimizing over the distance redshift relation.
  gsl_function F;
  F.params = 0;
  F.function = minimizer;

  gsl_min_fminimizer * s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
  gsl_min_fminimizer_set(s,&F,zmax,lower,upper);

  do {
    iter++;
    status = gsl_min_fminimizer_iterate(s);

    lower = gsl_min_fminimizer_x_lower(s);
    upper = gsl_min_fminimizer_x_upper(s);
    zmax  = gsl_min_fminimizer_x_minimum(s);
    
    status = gsl_min_test_interval(lower, upper, 1.0e-7, 1.0e-6);
  } while (status == GSL_CONTINUE && iter < iter_max);
  gsl_min_fminimizer_free(s);

  // We then loop over the range of redshifts and for each redshift 
  // we calculate the distance and linking length
  double * disttab = (double *)malloc(Nlink*sizeof(double));
  double * linktab = (double *)malloc(Nlink*sizeof(double));

  double zstep = zmax/((double)Nlink-1.0);  
  for (i=0; i<Nlink; i++) {
    double z = (i*zstep);
    double Omegaz     = (Omega*(1.0+z)*(1.0+z)*(1.0+z))/(Omega*(1.0+z)*(1.0+z)*(1.0+z)+(1.0-Omega));
    double delta_sim  = (18*PI*PI + 82*(Omegaz-1.0) - 39*(Omegaz-1.0)*(Omegaz-1.0)) / Omegaz;
    double delta_2lpt = pow((1.0 - (1.686/3.0) - ((1.686*1.686*pow(Omegaz,(-1.0/143.0)))/21.0)), -3.0);
    disttab[i] = rz(z)*rz(z);
    linktab[i] = pow(link_sim*(delta_sim/delta_2lpt), (1.0/3.0));
  }

  link_acc    = gsl_interp_accel_alloc();
  link_spline = gsl_spline_alloc(gsl_interp_cspline, Nlink);
  gsl_spline_init(link_spline, disttab, linktab, Nlink);

  free(disttab);
  free(linktab);

}

// Here we calculate the difference between dmax and d(zmax) to find zmax. 
// For a given value of z we evaluate the integral of 1/E(z') between 0 and z
// ==========================================================================
double minimizer(double z, void * params) {
  double d = rz(z);
  return (sqrt(dmax)-d)*(sqrt(dmax)-d);
}

// This returns the distance for a given value of z
// ================================================
double rz(double z) {
  double result, error;
  gsl_function F;
  F.function = Ez;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
  gsl_integration_qags(&F,0.0,z,0.0,1e-7,1000,w,&result,&error); 
  gsl_integration_workspace_free (w);
      
  // Put in units of Mpc/h
  result *= 2997.92458;

  return ((dmax-result)*(dmax-result));
}

double Ez(double z, void * params) {
  return sqrt(Omega*(1.0+z)*(1.0+z)*(1.0+z)+(1.0-Omega));
}   
