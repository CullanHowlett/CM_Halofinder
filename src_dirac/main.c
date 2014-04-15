#include "vars.h"

// This is the main driver program for the Friends-of-Friends and MPI routines
int main(int argc, char **argv) {

  // Set up MPI
  // ==========
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(argc < 2) {
    if(ThisTask == 0) {
      fprintf(stdout, "\nParameters are missing.\n");
      fprintf(stdout, "Call with <ParameterFile>\n\n");
    }
    ierr = MPI_Finalize();
    exit(0);
  }
 
  // Read in run parameters and allocate memory for global variables such as the particle data from
  // the inital conditions 
  Read_Parameterfile(argv[1]);

  if(ThisTask == 0) {
    printf("\nRun Parameters\n");
    printf("==============\n");
#ifdef VARLINK
    printf("Cosmology:\n");
    printf("  Omega_Matter(z=0) = %lf\n", Omega);
#endif
    printf("Simulation:\n");
    printf("  Particles(x, y, z) = %d, %d, %d\n", Px, Py, Pz);
    printf("  Boxsize(x, y, z) = %lf, %lf, %lf\n", Lxmax-Lxmin, Lymax-Lymin, Lzmax-Lzmin);
    printf("  Boundary size = %lf\n", boundarysize);
    printf("  Linking Length = %lf\n", linklength);
    printf("  Particles per Halo = %lf\n", nphalomin);
#ifdef VARLINK
    printf("  Linking Length table size = %lf\n", Nlink);
    printf("  Origin(x, y, z) = %lf, %lf, %lf\n", Origin_x, Origin_y, Origin_z);
#endif
    printf("\n");
  }

  if (ThisTask == 0) {
    printf("Setting memory and variables\n");
    printf("============================\n\n");
  }
  Set_Params();

  // A parallel read routine that reads the inital conditions and ends with each processor containing
  // a subsample of the initial particles, where each subsample consists of particles with similar locations
  // This routine also finds if a particle is near a boundary and additionally passes a copy of the particle
  // to the processors that share the particular boundary
  if (ThisTask == 0) {
    printf("Reading in Dark Matter field\n");
    printf("============================\n");
  }
#ifdef READ_INFO
  Read_Data_Info();
#else
  Read_Data();
#endif

  ierr = MPI_Barrier(MPI_COMM_WORLD);

  // This is the actual Friends-of-Friends routine, which can be run effectively serially if each
  // processor has all the necessary particles
  if (ThisTask == 0) {
    printf("\nPerforming Friends-of-friends algorithm\n");
    printf("=======================================\n");
  }
  FOF();

  if (ThisTask == 0) printf("Done :)\n\n");
  ierr = MPI_Finalize();

  return 0;
}

// This catches I/O errors occuring for fread(). In this case we better stop.
// ==========================================================================
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("\nERROR: I/O error (fread) on task=%d has occured\n\n", ThisTask);
      fflush(stdout);
      FatalError("read_data.c", 441);
    }
  return nread;
}

// Error message
// =============
void FatalError(char * filename, int linenum) {
  printf("Fatal Error at line %d in file %s\n\n", linenum, filename);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, 1);
  exit(1);
}
