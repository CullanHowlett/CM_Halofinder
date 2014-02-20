# The executable name
# ===================
EXEC    = CM_Halofinder_test

# Choose the machine you are running on. Currently versions for SCIAMA and DARWIN are implemented
# ===============================================================================================
MACHINE = SCIAMA
#MACHINE = DARWIN

# Options for optimization
# ========================
OPTIMIZE  = -O3 -Wall

# Various C preprocessor directives that change the way CM_Halofinder is made
# ===========================================================================
#PERIODIC = -DPERIODIC	          # Periodic Box
#OPTIONS += $(PERIODIC)

#VARLINK = -DVARLINK               # Redshift dependent linking length
#OPTIONS += $(VARLINK)

MEMORY_MODE = -DMEMORY_MODE       # Uses floats for particle data rather than doubles
OPTIONS += $(MEMORY_MODE)

LIGHTCONE = -DLIGHTCONE          # Read in data from a lightcone simulation rather than a snapshot
OPTIONS += $(LIGHTCONE)          # If used with UNFORMATTED option then the file must be organised in PICOLA's LIGHTCONE output style

#PARTICLE_ID = -DPARTICLE_ID       # Reads in and remembers the unsigned long long particle ID's then outputs them if necessary
#OPTIONS += $(PARTICLE_ID)

GADGET_STYLE = -DGADGET_STYLE          # Read files with Gadget's '1' style format, with the corresponding
OPTIONS += $(GADGET_STYLE)             # header and correct velocity units

#OUTPUT_PARTICLES = -DOUTPUT_PARTICLES  # Outputs all the particles in each halo, with output file format: nhalos, then for each each halo
#OPTIONS += $(OUTPUT_PARTICLES)         # npart then for each particle the position, velocity and ID if requested. Otherwise we only output the group
                                       # properties of each halo

READ_INFO = -DREAD_INFO          # Reads in the dark matter field using an associated info file, which contains a list of files in which each slice
OPTIONS += $(READ_INFO)          # of the field can be found. This will be quicker and less memory intensive than the standard read routine

# Run some checks on option compatability
# =======================================
ifdef PARTICLE_ID
ifndef OUTPUT_PARTICLES
   $(warning WARNING: We are not outputting the particles (OUTPUT_PARTICLES not on) so keeping track of particle IDs (with option PARTICLE_ID) is unnecessary.)
endif
endif

# Setup libraries and compile the code
# ====================================
ifeq ($(MACHINE),SCIAMA)
  CC = mpiCC
  GSL_INCL  = -I/opt/gridware/libs/gcc/gsl/1.14/include/
  GSL_LIBS  = -L/opt/gridware/libs/gcc/gsl/1.14/lib/  -lgsl -lgslcblas
  MPI_INCL  = -I/opt/gridware/mpi/gcc/openmpi/1_4_3/include
  MPI_LIBS  = -L/opt/gridware/mpi/gcc/openmpi/1_4_3/lib/ -lmpi
endif

ifeq ($(MACHINE),DARWIN)
  CC = mpiicc	
  GSL_INCL  = -I/usr/local/Cluster-Apps/gsl/1.9/include/
  GSL_LIBS  = -L/usr/local/Cluster-Apps/gsl/1.9/lib/  -lgsl -lgslcblas
  MPI_INCL  = -L/usr/local/Cluster-Apps/intel/impi/3.1/include
  MPI_LIBS  = -L/usr/local/Cluster-Apps/intel/impi/3.1/lib -lmpi
endif

LIBS   =   -lm $(MPI_LIBs) $(GSL_LIBS)

CFLAGS =   $(OPTIMIZE) $(GSL_INCL) $(MPI_INCL) $(OPTIONS)

OBJS = src/main.o src/startup.o src/FOF.o src/vars.o
ifdef VARLINK
OBJS += src/link.o
endif
ifdef READ_INFO
OBJS += src/read_data_info.o
else
OBJS += src/read_data.o 
endif

INCL   = src/vars.h  Makefile

CM_Halofinder: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o  $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f src/*.o *~ src/*~ $(EXEC)
