#include "vars.h"

// MPI Variables
// =============
int ierr;                                             // The MPI call status
int NTask;                                            // The total number of tasks
int ThisTask;                                         // The number of each task
int Local_nx, Local_ny, Local_nz;                     // The grid position of each task
double rmin[3], rmax[3], rmin_buff[3], rmax_buff[3];  // The extents of each task

// Particle data
struct part_data * P;

#ifdef VARLINK
// The linking length lookup table
gsl_spline * link_spline;
gsl_interp_accel * link_acc;
#endif

// Simulation variables
// ====================
unsigned int maxparticles;      // The maximum number of particles on each processor
unsigned int nparticles_tot;    // The actual number of particles on each processor
#ifdef VARLINK
double dmin, dmax;              //  The minimum and maximum distance of each processor from the origin
#endif

// Gadget-Style header
#ifdef GADGET_STYLE
struct io_header_1 header;
#endif

// Run parameters
// ==============
char InputDir[100];        // The input directory
char OutputDir[100];       // The output directory
char InputFileBase[100];   // The base input filename
char OutputFileBase[100];  // The base output filename
int Nx, Ny, Nz;            // The number of tasks in each direction
int particles;             // The number of particles in the simulation (per side)
int nread;                 // The number of processors reading in the data at once
int ninputfiles;           // The number of input files
int starting_file;         // The first file to start from
int nwrite;                // The number of processors writing out the data at once
double boxsize;            // The edge length of the whole simulation (Mpc/h)
double boundarysize;       // The size of the boundary region beyond each processor (Mpc/h)
double linklength;         // The linking length (Mpc/h)
double nphalomin;          // The minimum number of particles per halo
#ifdef VARLINK
int Nlink;                            // The size of the linking length lookup table
double Omega;                         // The matter density at z=0
double Origin_x, Origin_y, Origin_z;  // The position of the origin
#endif 

