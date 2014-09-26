# The executable name
# ===================
EXEC    = CM_Halofinder

# Choose the machine you are running on. Currently versions for SCIAMA and DARWIN are implemented
# ===============================================================================================
MACHINE = SCIAMA
#MACHINE = DARWIN

# Options for optimization
# ========================
OPTIMIZE  = -O3 -Wall

# Various C preprocessor directives that change the way CM_Halofinder is made
# ===========================================================================
#PERIODIC = -DPERIODIC	           # Periodic Box
#OPTIONS += $(PERIODIC)

#VARLINK = -DVARLINK               # Redshift dependent linking length
#OPTIONS += $(VARLINK)

MEMORY_MODE = -DMEMORY_MODE       # Uses floats for particle data rather than doubles
OPTIONS += $(MEMORY_MODE)

LIGHTCONE = -DLIGHTCONE           # Read in data from a PICOLA LIGHTCONE simulation rather than a snapshot
OPTIONS += $(LIGHTCONE)           # If data is in binary then the UNFORMATTED option must be used rather then GADGET_STYLE

#PARTICLE_ID = -DPARTICLE_ID       # Reads in and remembers the unsigned long long particle ID's then outputs them if necessary
#OPTIONS += $(PARTICLE_ID)

#GADGET_STYLE = -DGADGET_STYLE     # Read files with Gadget's '1' style format, with the corresponding header
#OPTIONS += $(GADGET_STYLE)        # This option is incompatible with PICOLA LIGHTCONE simulations, in this case we must use the UNFORMATTED option instead

UNFORMATTED = -DUNFORMATTED       # Read files in PICOLA's UNFORMATTED style. A PICOLA info-file containing the number of 
OPTIONS += $(UNFORMATTED)         # particles in each file must be specified. Within each chunk all the data for a given particle is written contiguously.

#OUTPUT_PARTICLES = -DOUTPUT_PARTICLES  # Outputs all the particles in each halo, along with the group properties of the halo.
#OPTIONS += $(OUTPUT_PARTICLES)         # Output format is ID (if requested), Position, Velocity

INERTIA = -DINERTIA		    # Calculate the unique elements of the moment of inertia tensor for each halo.
OPTIONS += $(INERTIA)		# Output is Ixx, Iyy, Izz, Ixy, Ixz, Iyz

DISPERSION = -DDISPERSION	# Calculate the unique elements of the velocity dispersion tensor for each halo.
OPTIONS += $(DISPERSION)	# Output is sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz

OUTPUT_UNFORMATTED = -DOUTPUT_UNFORMATTED # All output (halos and particle subsamples) is written in unformatted binary.
OPTIONS += $(OUTPUT_UNFORMATTED)          # Otherwise output is in ASCII

# Run some checks on option compatability
# =======================================
ifdef PARTICLE_ID
ifdef LIGHTCONE
   $(warning WARNING: PICOLA LIGHTCONE output does not output particle IDs)
endif
endif

ifdef GADGET_STYLE
ifdef LIGHTCONE
   $(error ERROR: LIGHTCONE AND GADGET_STYLE are not compatible, for binary output with PICOLA LIGHTCONE simulations please choose the UNFORMATTED option.)
endif
endif

ifdef UNFORMATTED
ifndef LIGHTCONE 
   $(error ERROR: UNFORMATTED option is currently incompatible with snapshot simulations, for binary output with snapshot simulations please choose the GADGET_STYLE option.)
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

OBJS = src_v2/main.o src_v2/startup.o src_v2/FOF.o src_v2/vars.o
ifdef VARLINK
OBJS += src_v2/link.o
endif
ifdef READ_INFO
OBJS += src_v2/read_data_info.o
else
OBJS += src_v2/read_data.o 
endif

INCL   = src_v2/vars.h  Makefile

CM_Halofinder: $(OBJS) 
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o  $(EXEC)

$(OBJS): $(INCL) 

clean:
	rm -f src_v2/*.o *~ src_v2/*~ $(EXEC)
