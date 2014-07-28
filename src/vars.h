#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// GSL libraries
#include <gsl/gsl_rng.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

// MPI and FFTW libraries
#include <mpi.h>

// Some definitions
#define  PI          3.14159265358979323846 // PI (obviously)

// MPI Variables
// =============
extern int ierr;                                             // The MPI call status
extern int NTask;                                            // The total number of tasks
extern int ThisTask;                                         // The number of each task
extern int Local_nx, Local_ny, Local_nz;                     // The grid position of each task
extern double rmin[3], rmax[3], rmin_buff[3], rmax_buff[3];  // The extents of each task

// Particle data
extern struct part_data {
#ifdef PARTICLE_ID
  unsigned long long ID;
#endif
#ifdef MEMORY_MODE
  float Pos[3];
  float Vel[3];
#else
  double Pos[3];
  double Vel[3];
#endif
} * P;

// Seperated Particle Data
#ifdef MEMORY_MODE
extern struct part_data_half {
  float Var[3];
} * P_pos, * P_vel;
#else
extern struct part_data_half {
  double Var[3];
} * P_pos, * P_vel;
#endif

#ifdef VARLINK
// The linking length lookup table
extern gsl_spline * link_spline;
extern gsl_interp_accel * link_acc;
#endif

// Simulation variables
// ====================
extern unsigned int maxparticles;      // The maximum number of particles on each processor
extern unsigned int nparticles_tot;    // The actual number of particles on each processor
#ifdef VARLINK
extern double dmin, dmax;              //  The minimum and maximum distance of each processor from the origin
#endif

// Gadget-Style header
#ifdef GADGET_STYLE
extern struct io_header_1 {
  unsigned int npart[6];      // npart[1] gives the number of particles in the file, other particle types are ignored
  double mass[6];             // mass[1] gives the particle mass
  double time;                // Cosmological scale factor of snapshot
  double redshift;            // redshift of snapshot
  int flag_sfr;               // Flags whether star formation is used
  int flag_feedback;          // Flags whether feedback from star formation is included
  unsigned int npartTotal[6]; // npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                              // the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder.
  int flag_cooling;           // Flags whether radiative cooling is included
  int num_files;              // Determines the number of files that are used for a snapshot
  double BoxSize;             // Simulation box size (in code units)
  double Omega0;              // matter density
  double OmegaLambda;         // vacuum energy density
  double HubbleParam;         // little 'h'
  int flag_stellarage;        // flags whether the age of newly formed stars is recorded and saved
  int flag_metals;            // flags whether metal enrichment is included
  int hashtabsize;            // gives the size of the hashtable belonging to this snapshot file
  char fill[84];              // Fills to 256 Bytes
}
header;
#endif

// Run parameters
// ==============
extern char InputFileBase[500];   // The base input filename
extern char OutputFileBase[500];  // The base output filename
#ifdef LIGHTCONE
#ifdef UNFORMATTED
extern char InputInfoFile[500];   // The PICOLA info file filename
#endif
#endif
extern int Seed;                  // The seed for the subsampling random number generator
extern int nread;                 // The number of processors reading in the data at once
extern int nwrite;                // The number of processors writing out the data at once
extern int Px, Py, Pz;            // The number of particles in the simulation (per side)
extern int Nx, Ny, Nz;            // The number of tasks in each direction
extern int InputStyle;            // How to work out which files to read in
extern int ninputfiles;           // The number of input files
extern int starting_file;         // The first file to start from
extern double Buffer;             // The amount of buffer memory to allocate to account for the inhomogeneity of the dark matter field
extern double SampInHalos;        // The fractional subsampling rate for particles inside halos
extern double SampOutHalos;       // The fractional subsmapling rate for particles outside of halos
extern double Lxmin, Lxmax;       // The x-boundaries of the simulation
extern double Lymin, Lymax;       // The y-boundaries of the simulation
extern double Lzmin, Lzmax;       // The z-boundaries of the simulation
extern double nphalomin;          // The minimum number of particles per halo
extern double linklength;         // The linking length (Mpc/h)
extern double boundarysize;       // The size of the boundary region beyond each processor (Mpc/h)
#ifdef VARLINK
extern double Omega;                         // The matter density at z=0
extern double Origin_x, Origin_y, Origin_z;  // The position of the origin
#endif 

// Prototypes
// ==========
void FOF(void);
void Read_Data(void);
void Set_Params(void);
void Output_Halos(void);
void Checkhalo(unsigned int i);
void Read_Parameterfile(char * fname);
void FatalError(char * filename, int linenum);
void Befriend(unsigned int ip, unsigned int ip2);
void Subsample(int in, unsigned int nsubsamp, unsigned int * subsamp);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
#ifdef VARLINK
void Create_Link(void);
double rz(double z);
double Ez(double z, void * params);
double minimizer(double z, void * params);
#endif
