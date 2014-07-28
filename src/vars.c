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

// Seperated particle data
struct part_data_half * P_pos, * P_vel;

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
char InputFileBase[500];   // The base input filename
char OutputFileBase[500];  // The base output filename
#ifdef LIGHTCONE
#ifdef UNFORMATTED
char InputInfoFile[500];   // The PICOLA info file filename
#endif
#endif
int Seed;                  // The seed for the subsampling random number generator
int nread;                 // The number of processors reading in the data at once
int nwrite;                // The number of processors writing out the data at once
int Px, Py, Pz;            // The number of particles in the simulation (per side)
int Nx, Ny, Nz;            // The number of tasks in each direction
int InputStyle;            // How to work out which files to read in
int ninputfiles;           // The number of input files
int starting_file;         // The first file to start from
double Buffer;             // The amount of buffer memory to allocate to account for the inhomogeneity of the dark matter field
double SampInHalos;        // The fractional subsampling rate for particles inside halos
double SampOutHalos;       // The fractional subsmapling rate for particles outside of halos
double Lxmin, Lxmax;       // The x-boundaries of the simulation
double Lymin, Lymax;       // The y-boundaries of the simulation
double Lzmin, Lzmax;       // The z-boundaries of the simulation
double boundarysize;       // The size of the boundary region beyond each processor (Mpc/h)
double linklength;         // The linking length (Mpc/h)
double nphalomin;          // The minimum number of particles per halo
#ifdef VARLINK
double Omega;                         // The matter density at z=0
double Origin_x, Origin_y, Origin_z;  // The position of the origin
#endif 

