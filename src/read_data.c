// This file contains all subroutines associated with MPI-based file input and data organisation

#include "vars.h"

void Read_Data(void) {

  // This routine gets each processor to open a different initial conditions file and then pass the data
  // to the correct processor based on the particle position using MPI_Alltoall. Hence the greatest
  // speed is obtained when we use particle data that is spread over multiple files. Ideally we use the same
  // number of input files as the number of processors so that no processor wastes time sitting idle, however
  // this routine can support any combination of processor and input files numbers.

  // This routine is set up to read a formatted file containing 6 columns: x, y, z, v_x, v_y, v_z
  
  // NOTE: Due to a limitation in MPI we are unable to send more that 2^31 bytes of information in a single 
  // MPI_Alltoall command without a segmentation fault due to integer overflow. Whilst this seems like a lot the problem 
  // becomes clear when we realise that we send all the positions/velocities in terms of their individual bytes. 
  // Hence this means that the maximum number of particles a single processor can send/receive at once is ~170,000,000 
  // for single floating point precision and ~85,000,000 for double. This limit can be reached surprising easily if the number 
  // of files containing the full input simulation is quite small, in tests the code could NOT handle a (768 Mpc/h)^3 box with 1024^3 particles
  // spread over 32 files in double precision. In this case we could either: convert to single precision, increase the number of files
  // the input is spread over, or try and modify the particle data structures and pass over less information at a time.
  
  // I have split the passing of the positions and velocities into separate directions to as above. I can't change nparticles or cnparticles to unsigned int
  // as MPI only allows ints for these (MPI_Alltoallv complains if these are unsigned. However this shouldn't matter as we are no longer multiplying by the 
  // number of bytes being sent at this point so the maximum number of particles we can transfer is the max of an int. One further improvement that i'd like to make
  // (but later) is to create a derived MPI datatype for the particles structure so that we can pass all the particle data (position & velocities) at once 
  // but without having to multiply by the number of bytes and hence allowing us to pass 2^31 particles at a time.
  
  FILE * fp=NULL;
  char buf[500];
  char ** filelist=NULL;
  int i, imax, k, q;
  int readproc, file_exists, filenumber=ninputfiles, nreadid;
  int *readid;
  int processor_number, neighbour_x, neighbour_y, neighbour_z, processor_comp[3], subscript[3];
  unsigned int j, m, n, npartfile=0, Nout=0, Nout_glob=0; 
  int nparticles_sum, nparticles_recv_sum;
  int * data_index, * nparticles, * nparticles_recv, * cnparticles, * cnparticles_recv;
  double xmin=Lxmin, ymin=Lymin, zmin=Lzmin;
  double xmax=Lxmax, ymax=Lymax, zmax=Lzmax;
  struct part_data * P_file=NULL;
#ifdef PARTICLE_ID
  unsigned long long * P_sorted_id, * P_temp_id;
#endif
  part_datatype * P_sorted_posx, * P_sorted_posy, * P_sorted_posz;
  part_datatype * P_sorted_velx, * P_sorted_vely, * P_sorted_velz;
  part_datatype * P_temp_posx, * P_temp_posy, * P_temp_posz;
  part_datatype * P_temp_velx, * P_temp_vely, * P_temp_velz;
  //struct part_data_half * P_sorted_pos, * P_sorted_vel, * P_temp_pos, * P_temp_vel;
#ifdef LIGHTCONE
#ifdef UNFORMATTED 
  char * name, * bufa;
  int cont, ninfo=0, dummy, * filenumbers=NULL;
  unsigned int npart, * npartfiles=NULL;
  float val, xminf, yminf, zminf, xmaxf, ymaxf, zmaxf;
#ifdef PARTICLE_ID
  unsigned long long ival;
#endif
#else
  char bufcount[500];
#endif
#else
#ifdef GADGET_STYLE
  int dummy;
  float val;
#ifdef PARTICLE_ID
  unsigned long long ival;
#endif
#else
  char bufcount[500];
#endif
#endif

  nparticles_tot=0;

  // Only a small subset of the processors can read in the data. In order to balance the load
  // as much as possible, without knowing the machine configuration, we will space out the
  // 'reading' processors linearly.
  // First generate an array of which processors can read in data
  nreadid=0;
  readid = (int *)malloc(nread*sizeof(int));
  for (i=0; i<NTask; i++) {
    if((i % (int)int((double)NTask/(double)nread)) == 0) {
      readid[nreadid] = i;
      nreadid++;
      if (nreadid == nread) break;       // This is here to account for when nread is not a factor of NTask        
    }
  }

  readproc=0;
  for (k=0; k<nread; k++) {
    if (readid[k] == ThisTask) readproc=1;
  }
  
  // If Inputstyle!=0 we need to get the necessary processors to read in the list of input files and store them.
  if (InputStyle) {
    ninputfiles=0;
    if (readproc) {
      if(!(fp=fopen(InputFileBase, "rb"))) {
        printf("\nERROR: Task %d unable to file %s.\n\n", ThisTask, InputFileBase); 
        FatalError("read_data.c", 92);
      }
      while(fgets(buf,500,fp)) {
        // Check for empty/commented lines and ignore them
        char buf1[500], buf2[500];
        if(sscanf(buf, "%s%s", buf1, buf2) < 1) continue;
        if((buf1[0] == '%')) continue;
        ninputfiles++;
      }
      rewind(fp);
      filelist = (char **)malloc(ninputfiles*sizeof(char*));
      for (q=0; q<ninputfiles; q++) filelist[q] = (char*)malloc(500*sizeof(char));
      ninputfiles=0;
      while(fgets(buf,500,fp)) {
        // Check for empty/commented lines and ignore them
        char buf1[500], buf2[500];
        if(sscanf(buf, "%s%s", buf1, buf2) < 1) continue;
        if((buf1[0] == '%')) continue;
        strcpy(filelist[ninputfiles],buf1);
        ninputfiles++;
      }
      fclose(fp);
    }
    ierr = MPI_Bcast(&ninputfiles, 1, MPI_INT, readid[0], MPI_COMM_WORLD);
  }

#ifdef LIGHTCONE
#ifdef UNFORMATTED
  // For unformatted lightcones we need the PICOLA infofile to get the total number of particles in each file.
  // There are plenty of ways of reading in the lightcones that don't require the infofile but this is the quickest,
  // and as this mode is designed for PICOLA-LIGHTCONE simulations there will be an infofile available anyway
  // Read in the info file and determine which files we need to read in 
  if (readproc) {
    if(!(fp = fopen(InputInfoFile, "r"))) {
      printf("\nERROR: Can't open info file '%s'.\n\n", buf);
      FatalError("read_data.c", 124);
    }

    while(fgets(buf,500,fp)) ninfo++;
    filenumbers = (int *)malloc(ninfo*sizeof(int));
    npartfiles = (unsigned int *)malloc(ninfo*sizeof(unsigned int));
    rewind(fp);
    ninfo=0;
    while(fgets(buf,500,fp)) {
      if(strncmp(buf,"#",1)==0) continue;
      if(sscanf(buf,"%d %f %f %f %f %f %f %u\n",&filenumber,&xminf,&yminf,&zminf,&xmaxf,&ymaxf,&zmaxf,&npart)!=8) { 
        printf("Bad input format in info file '%s\n", InputInfoFile); 
        FatalError("read_data.c", 135);
      }  
      filenumbers[ninfo] = filenumber;
      npartfiles[ninfo] = npart;
      ninfo++;
    }
    fclose(fp);
  }
#endif
#endif

  if (ThisTask == 0) {
    printf("%d tasks reading in %d files...\n", nread, ninputfiles);
    fflush(stdout);
  }
        
  imax = (int)ceil((double)ninputfiles/(double)nread);
  for(i=0; i<imax; i++) {

    if (ThisTask == 0) {
      printf("Task %d performing loop over files: %d of %d\n", ThisTask, i+1, imax);
      printf("------------------------------------\n");
      fflush(stdout);
    }

    ierr = MPI_Barrier(MPI_COMM_WORLD);

    // Now find out which processors were chosen and give each of them a unique file to read in
    // If InputStyle=0 we just give them a file 'InputFileBase+filenumber'. Otherwise'filenumber' is
    // the entry in 'filelist' that contains the filename for each processor to read in.
    file_exists=0;
    for (k=0; k<nread; k++) {
      if (readid[k] == ThisTask) {
        filenumber = i*nread+k;
        if (filenumber < ninputfiles) {
          file_exists=1;
          if (InputStyle) {
            sprintf(buf, "%s", filelist[filenumber]);
          } else {
            sprintf(buf, "%s.%d", InputFileBase, filenumber+starting_file);
          }
        }
      }
    }

    nparticles = (int *)malloc(NTask*sizeof(int));
    nparticles_recv = (int *)malloc(NTask*sizeof(int));

    // Set the default number of particles to 1. This is so we that we don't have any tasks trying to send zero data
    // to any other tasks. This is removed later.
    for(k=0; k<NTask; k++) nparticles[k] = nparticles_recv[k] = 1;

// Lightcone simulations
// =====================
#ifdef LIGHTCONE

// Binary
#ifdef UNFORMATTED

    // Set the default number of particles to 1. This is so we that we don't have any tasks trying to send zero data
    // to any other tasks. This is removed later.
    for(k=0; k<NTask; k++) nparticles[k] = nparticles_recv[k] = 1;

    if(file_exists) {
    
      // Find out how many particles are in the file using the info file data
      if (InputStyle) {       
        // Isolate the filenumber
        bufa = strdup(buf);
        name = strsep(&bufa, ".");
        filenumber = atoi(strsep(&bufa, "."));
        for (k=0; k<ninfo; k++) {
          if (filenumbers[k] == filenumber) {
            npartfile = npartfiles[k];
            break;
          }
        }
      } else {
        for (k=0; k<ninfo; k++) {
          if (filenumbers[k] == (filenumber+starting_file)) {
            npartfile = npartfiles[k];
            break;
          }
        }
      }

      if(!(fp=fopen(buf, "rb"))) {
        printf("\nERROR: Task %d unable to file %s.\n\n", ThisTask, buf); 
        FatalError("read_data.c", 121);
      }

      P_file = (struct part_data *)malloc(npartfile*sizeof(struct part_data));

      npartfile=0;
      while(1) {
        cont = 0;
        cont = fread(&dummy, sizeof(dummy), 1, fp);
        if (cont) {
          my_fread(&npart, sizeof(unsigned int), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          for (j=npartfile; j<npart+npartfile; j++) {
#ifdef PARTICLE_ID
            my_fread(&(id), sizeof(unsigned int), 1, fp);
            P_file[j].ID = ival;
#endif
            for (k=0; k<3; k++) {
              my_fread(&(val), sizeof(float), 1, fp);
              P_file[j].Pos[k] = val;
            }
            for (k=0; k<3; k++) {
              my_fread(&(val), sizeof(float), 1, fp);
              P_file[j].Vel[k] = val;
            }
          }
          my_fread(&dummy, sizeof(dummy), 1, fp);
          npartfile += npart;
        } else {
          fclose(fp);
          break;
        }
      }

    /*while(1) {

      cont = 0;

      // Set the default number of particles to 1. This is so we that we don't have any tasks trying to send zero data
      // to any other tasks. This is removed later.
      for(k=0; k<NTask; k++) nparticles[k] = nparticles_recv[k] = 1;

      if(file_exists) {
        cont = fread(&dummy, sizeof(dummy), 1, fp);
        if (cont) {
          my_fread(&npartfile, sizeof(unsigned int), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          P_file = (struct part_data *)malloc(npartfile*sizeof(struct part_data));
          for (j=0; j<npartfile; j++) {
#ifdef PARTICLE_ID
            my_fread(&(id), sizeof(unsigned int), 1, fp);
            P_file[j].ID = ival;
#endif
            for (k=0; k<3; k++) {
              my_fread(&(val), sizeof(float), 1, fp);
              P_file[j].Pos[k] = val;
            }
            for (k=0; k<3; k++) {
              my_fread(&(val), sizeof(float), 1, fp);
              P_file[j].Vel[k] = val;
            }
          }
          my_fread(&dummy, sizeof(dummy), 1, fp);
        } else {
          fclose(fp);
          file_exists = 0;
        }
      }

      MPI_Allreduce(MPI_IN_PLACE, &cont, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if(!(cont)) break;

      if (file_exists) {*/

// ASCII
#else 

    if(file_exists) {
      if(!(fp=fopen(buf, "r"))) {
        printf("\nERROR: Task %d unable to find file %s.\n\n", ThisTask, buf); 
        FatalError("read_data.c", 165);
      }      
      
      npartfile = 0;
      while(fgets(bufcount,500,fp)) npartfile++;
      P_file = (struct part_data *)malloc(npartfile*sizeof(struct part_data));
      rewind(fp);
      for(j=0; j<npartfile; j++) {
#ifdef MEMORY_MODE
#ifdef PARTICLE_ID
        if((fscanf(fp, "%llu %f %f %f %f %f %f\n", &(P_file[j].ID), &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                   &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 7)) {
#else
        if((fscanf(fp, "%f %f %f %f %f %f\n", &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                              &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 6)) {
#endif
#else
#ifdef PARTICLE_ID
        if((fscanf(fp, "%llu %lf %lf %lf %lf %lf %lf\n", &(P_file[j].ID), &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                         &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 7)) {
#else
        if((fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                    &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 6)) {
#endif
#endif
          printf("\nERROR: Task %d has incorrect line format in file %s.\n\n", ThisTask, buf); 
          ierr = MPI_Abort(MPI_COMM_WORLD, 2);
        }
      }
      fclose(fp);
#endif

// Snapshot
// ========
#else

// Binary
#ifdef GADGET_STYLE

    if(file_exists) {
      if(!(fp=fopen(buf, "rb"))) {
        printf("\nERROR: Task %d unable to find file %s.\n\n", ThisTask, buf); 
        FatalError("read_data.c", 124);
      }
      
      my_fread(&dummy, sizeof(dummy), 1, fp);
      my_fread(&header, sizeof(header), 1, fp);
      my_fread(&dummy, sizeof(dummy), 1, fp);
      npartfile = header.npart[1];
      P_file = (struct part_data *)malloc(npartfile*sizeof(struct part_data));

      // Read in the particle positions
      my_fread(&dummy, sizeof(dummy), 1, fp);
      for(j=0; j<npartfile; j++) {
        for (k=0; k<3; k++) {
          my_fread(&(val), sizeof(float), 1, fp);
          P_file[j].Pos[k] = val;
        }
      }
      my_fread(&dummy, sizeof(dummy), 1, fp);

      // Read in the particle velocities
      my_fread(&dummy, sizeof(dummy), 1, fp);
      for(j=0; j<npartfile; j++) {
        for (k=0; k<3; k++) {
          my_fread(&(val), sizeof(float), 1, fp);
          P_file[j].Vel[k] = val;
        }
      }
      my_fread(&dummy, sizeof(dummy), 1, fp);

#ifdef PARTICLE_ID
      // Read in the particle ids
      my_fread(&dummy, sizeof(dummy), 1, fp);
      for(j=0; j<npartfile; j++) {
        my_fread(&(ival), sizeof(unsigned int), 1, fp);
        P_file[j].ID = ival;
      }
      my_fread(&dummy, sizeof(dummy), 1, fp);
#endif

      fclose(fp);

// ASCII
#else

    if(file_exists) {
      if(!(fp=fopen(buf, "r"))) {
        printf("\nERROR: Task %d unable to find file %s.\n\n", ThisTask, buf); 
        FatalError("read_data.c", 165);
      }      
      
      npartfile = 0;
      while(fgets(bufcount,500,fp)) npartfile++;
      P_file = (struct part_data *)malloc(npartfile*sizeof(struct part_data));
      rewind(fp);
      for(j=0; j<npartfile; j++) {
#ifdef MEMORY_MODE
#ifdef PARTICLE_ID
        if((fscanf(fp, "%llu %f %f %f %f %f %f\n", &(P_file[j].ID), &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                   &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 7)) {
#else
        if((fscanf(fp, "%f %f %f %f %f %f\n", &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                              &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 6)) {
#endif
#else
#ifdef PARTICLE_ID
        if((fscanf(fp, "%llu %lf %lf %lf %lf %lf %lf\n", &(P_file[j].ID), &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                         &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 7)) {
#else
        if((fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &(P_file[j].Pos[0]), &(P_file[j].Pos[1]), &(P_file[j].Pos[2]), 
                                                    &(P_file[j].Vel[0]), &(P_file[j].Vel[1]), &(P_file[j].Vel[2])) != 6)) {
#endif
#endif
          printf("\nERROR: Task %d has incorrect line format in file %s.\n\n", ThisTask, buf); 
          ierr = MPI_Abort(MPI_COMM_WORLD, 2);
        }
      }
      fclose(fp);
#endif
#endif

      printf("Task %d finished reading %d particles from file %s...\n", ThisTask, npartfile, buf);
      fflush(stdout);

      // This is how we identify particles as outside the simulation. For non-periodic boxes we still keep
      // particles around the simulation edge, if they are supplied, otherwise the halos will be incorrect
      xmin = Lxmin; ymin = Lymin; zmin = Lzmin;
      xmax = Lxmax; ymax = Lymax; zmax = Lzmax;
#ifndef PERIODIC
      xmin -= boundarysize; ymin -= boundarysize; zmin -= boundarysize;
      xmax += boundarysize; ymax += boundarysize; zmax += boundarysize;
#endif

      // For each particle we calculate which processor it should go to and count the number of particles
      // that are being sent to each processor
      for(j=0; j<npartfile; j++) {
        if((P_file[j].Pos[0] < xmin) || (P_file[j].Pos[0] > xmax) || (P_file[j].Pos[1] < ymin) || (P_file[j].Pos[1] > ymax) || (P_file[j].Pos[2] < zmin) || (P_file[j].Pos[2] > zmax)) {
          Nout++;
          continue;
        } 
        processor_comp[0]=(int)floor(Nx*((P_file[j].Pos[0]-Lxmin)/(Lxmax-Lxmin)));
        processor_comp[1]=(int)floor(Ny*((P_file[j].Pos[1]-Lymin)/(Lymax-Lymin)));
        processor_comp[2]=(int)floor(Nz*((P_file[j].Pos[2]-Lzmin)/(Lzmax-Lzmin)));
#ifndef PERIODIC
        if (processor_comp[0] < 0) processor_comp[0]++;
        if (processor_comp[1] < 0) processor_comp[1]++;
        if (processor_comp[2] < 0) processor_comp[2]++;
#endif     
        if (processor_comp[0] >= Nx) processor_comp[0]--;
        if (processor_comp[1] >= Ny) processor_comp[1]--;
        if (processor_comp[2] >= Nz) processor_comp[2]--;
        
        processor_number=Ny*Nx*processor_comp[2]+Nx*processor_comp[1]+processor_comp[0];
        nparticles[processor_number]++;

        // This section loops over each particle and for each particle we calculate the processor the particle would be
        // on if it had a position displaced by the boundarysize. We then produce an array of three numbers (x, y, z) that
        // indicate  which of the boundaries the particle is associated with. The number itself can be either 0, -1 or 1:
        // -1 - The particle is a boundary particle with the -x, -y or -z boundary
        //  0 - The particle is not a boundary particle with that particular boundary  
        //  1 - The particle is a boundary particle with the +x, +y or +z boundary
        // IFDEF loops ensure that if the particle were to be assigned to a neighbour on a periodic boundary then
        // if we are using isolated boundary conditions we don't increment the count to that neighbour
        // and as such don't (in the end) pass the particle.
        for (k=0; k<3; k++) subscript[k]=0;
        for (k=-1; k<=1; k+=2) {
          neighbour_x=(int)floor(Nx*((P_file[j].Pos[0]-Lxmin+k*boundarysize)/(Lxmax-Lxmin)));
          if (neighbour_x != processor_comp[0]) {
#ifndef PERIODIC
            if ((neighbour_x < 0) || (neighbour_x >= Nx)) continue;
#endif
            subscript[0]=k;
          }
        }
        for (k=-1; k<=1; k+=2) {
          neighbour_y=(int)floor(Ny*((P_file[j].Pos[1]-Lymin+k*boundarysize)/(Lymax-Lymin)));
          if (neighbour_y != processor_comp[1]) {
#ifndef PERIODIC
            if ((neighbour_y < 0) || (neighbour_y >= Ny)) continue;
#endif
            subscript[1]=k;
          }
        }
        for (k=-1; k<=1; k+=2) {
          neighbour_z=(int)floor(Nz*((P_file[j].Pos[2]-Lzmin+k*boundarysize)/(Lzmax-Lzmin)));
          if (neighbour_z != processor_comp[2]) {
#ifndef PERIODIC
            if ((neighbour_z < 0) || (neighbour_z >= Nz)) continue;
#endif
            subscript[2]=k;
          }
        }

        // This set of if loops then increments the counts based on the array of three numbers produced
        // above. The caveat is that if the 3 number array contains more than one non zero number we have
        // to increment the counts for each non zero number individually and then for the full combination
        // For example: If a particle is in the lower left corner it will have to be passed to the neighbour
        // on the left and the neighbour on the bottom as well as the neighbour on the bottom left corner.
        neighbour_x=processor_comp[0]+subscript[0];
        neighbour_y=processor_comp[1]+subscript[1];
        neighbour_z=processor_comp[2]+subscript[2];
#ifdef PERIODIC
        if (neighbour_x < 0) neighbour_x=Nx-1;
        if (neighbour_y < 0) neighbour_y=Ny-1;
        if (neighbour_z < 0) neighbour_z=Nz-1;
        if (neighbour_x >= Nx) neighbour_x=0;
        if (neighbour_y >= Ny) neighbour_y=0;
        if (neighbour_z >= Nz) neighbour_z=0;
#endif
        if (subscript[0] != 0) {
          processor_number=Ny*Nx*processor_comp[2]+Nx*processor_comp[1]+neighbour_x;
          nparticles[processor_number]++;
          if (subscript[1] != 0) {
            processor_number=Ny*Nx*processor_comp[2]+Nx*neighbour_y+neighbour_x;
            nparticles[processor_number]++;
            if (subscript[2] != 0) {
              processor_number=Ny*Nx*neighbour_z+Nx*neighbour_y+neighbour_x;
              nparticles[processor_number]++;
            }
          }
          if (subscript[2] != 0) {
            processor_number=Ny*Nx*neighbour_z+Nx*processor_comp[1]+neighbour_x;
            nparticles[processor_number]++;
          }
        }
        if (subscript[1] != 0) {
          processor_number=Ny*Nx*processor_comp[2]+Nx*neighbour_y+processor_comp[0];
          nparticles[processor_number]++;
          if (subscript[2] != 0) {
            processor_number=Ny*Nx*neighbour_z+Nx*neighbour_y+processor_comp[0];
            nparticles[processor_number]++;
          }
        }
        if (subscript[2] != 0) {
          processor_number=Ny*Nx*neighbour_z+Nx*processor_comp[1]+processor_comp[0];
          nparticles[processor_number]++;
        } 
      }
    }
      
    ierr = MPI_Alltoall(&(nparticles[0]), 1, MPI_UNSIGNED, &(nparticles_recv[0]), 1, MPI_UNSIGNED, MPI_COMM_WORLD);
 
    nparticles_sum = nparticles[0];
    nparticles_recv_sum = nparticles_recv[0];
    cnparticles = (int *)calloc(NTask,sizeof(int));
    cnparticles_recv = (int *)calloc(NTask,sizeof(int)); 

    // Produce a cumulative count of the number of particles being sent to each processor
    for (k=1; k<NTask; k++) {
      nparticles_sum      += nparticles[k];
      nparticles_recv_sum += nparticles_recv[k];
      cnparticles[k]       = cnparticles[k-1]+nparticles[k-1];
      cnparticles_recv[k]  = cnparticles_recv[k-1]+nparticles_recv[k-1];
    }

    //P_sorted = (struct part_data *)malloc(nparticles_sum*sizeof(struct part_data));
#ifdef PARTICLE_ID
    P_sorted_id = (unsigned long long *)malloc(nparticles_sum*sizeof(unsigned long long));
#endif
    //P_sorted_pos = (struct part_data_half *)malloc(nparticles_sum*sizeof(struct part_data_half));
    //P_sorted_vel = (struct part_data_half *)malloc(nparticles_sum*sizeof(struct part_data_half));
    P_sorted_posx = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));
    P_sorted_posy = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));
    P_sorted_posz = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));
    P_sorted_velx = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));
    P_sorted_vely = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));
    P_sorted_velz = (part_datatype *)malloc(nparticles_sum*sizeof(part_datatype));

    //if (file_exists) {
      //printf("Task %d sending %u particles...\n", ThisTask, nparticles_sum-NTask);
      //fflush(stdout);
    //}
    //MPI_Barrier(MPI_COMM_WORLD);
    //printf("Task %d receiving %u particles...\n", ThisTask, nparticles_recv_sum-NTask);
    //fflush(stdout);
          
    if (nparticles_tot+nparticles_recv_sum-NTask >= maxparticles) {
      printf("\nERROR: Task %d receiving more particles than allocated for.\n", ThisTask);
      printf("       Please increase the 'Buffer' size.\n\n");
      FatalError("read_data.c", 331);
      fflush(stdout);
    }
    
    if (file_exists) {

      // Again loop over each particle but this time sort the data so that each segment of the array is
      // nparticles(processor) long, begins at cnparticles(processor) and contains only the data being
      // sent to the particular processor
      data_index = (int *)calloc(NTask,sizeof(int));
      for (j=0; j<npartfile; j++) {
        if((P_file[j].Pos[0] < xmin) || (P_file[j].Pos[0] > xmax) || (P_file[j].Pos[1] < ymin) || (P_file[j].Pos[1] > ymax) || (P_file[j].Pos[2] < zmin) || (P_file[j].Pos[2] > zmax)) continue;
        processor_comp[0]=(int)floor(Nx*((P_file[j].Pos[0]-Lxmin)/(Lxmax-Lxmin)));
        processor_comp[1]=(int)floor(Ny*((P_file[j].Pos[1]-Lymin)/(Lymax-Lymin)));
        processor_comp[2]=(int)floor(Nz*((P_file[j].Pos[2]-Lzmin)/(Lzmax-Lzmin)));
#ifndef PERIODIC
        if (processor_comp[0] < 0) processor_comp[0]++;
        if (processor_comp[1] < 0) processor_comp[1]++;
        if (processor_comp[2] < 0) processor_comp[2]++;
#endif 
        if (processor_comp[0] >= Nx) processor_comp[0]--;
        if (processor_comp[1] >= Ny) processor_comp[1]--;
        if (processor_comp[2] >= Nz) processor_comp[2]--;

        processor_number=Ny*Nx*processor_comp[2]+Nx*processor_comp[1]+processor_comp[0];
#ifdef PARTICLE_ID
        P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
        //for (k=0; k<3; k++) {
        //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
        //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
        //}
        P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
        P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
        P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
        P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
        P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
        P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
        data_index[processor_number]++;

        // This is a repeat of the above process, looping over each particle  and finding which boundaries it
        // is near. This time however we store the particle data itself in a single array taking care to note that,
        // as with the particle counts, a single particles may have to be given to more than one neighbour and as such have
        // to be put in the array multiple times in different locations.
        for (k=0; k<3; k++) subscript[k]=0;
        for (k=-1; k<=1; k+=2) {
          neighbour_x=(int)floor(Nx*((P_file[j].Pos[0]-Lxmin+k*boundarysize)/(Lxmax-Lxmin)));
          if (neighbour_x != processor_comp[0]) {
#ifndef PERIODIC
            if ((neighbour_x < 0) || (neighbour_x >= Nx)) continue;
#endif
            subscript[0]=k;
          }
        }
        for (k=-1; k<=1; k+=2) {
          neighbour_y=(int)floor(Ny*((P_file[j].Pos[1]-Lymin+k*boundarysize)/(Lymax-Lymin)));
          if (neighbour_y != processor_comp[1]) {
#ifndef PERIODIC
            if ((neighbour_y < 0) || (neighbour_y >= Ny)) continue;
#endif
            subscript[1]=k;
          }
        }
        for (k=-1; k<=1; k+=2) {
          neighbour_z=(int)floor(Nz*((P_file[j].Pos[2]-Lzmin+k*boundarysize)/(Lzmax-Lzmin)));
          if (neighbour_z != processor_comp[2]) {
#ifndef PERIODIC
            if ((neighbour_z < 0) || (neighbour_z >= Nz)) continue;
#endif
            subscript[2]=k;
          }
        }

        // This rather confusing section is very similar to that above and differs only in the fact that the
        // subscripts for where the data is saved have to be very convoluted to store the particle data
        // in the correct place. 
        neighbour_x=processor_comp[0]+subscript[0];
        neighbour_y=processor_comp[1]+subscript[1];
        neighbour_z=processor_comp[2]+subscript[2];
#ifdef PERIODIC
        if (neighbour_x < 0) neighbour_x=Nx-1;
        if (neighbour_y < 0) neighbour_y=Ny-1;
        if (neighbour_z < 0) neighbour_z=Nz-1;
        if (neighbour_x >= Nx) neighbour_x=0;
        if (neighbour_y >= Ny) neighbour_y=0;
        if (neighbour_z >= Nz) neighbour_z=0;
#endif     
        if (subscript[0] != 0) {
          processor_number=Ny*Nx*processor_comp[2]+Nx*processor_comp[1]+neighbour_x;
#ifdef PARTICLE_ID
          P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
          //for (k=0; k<3; k++) {
          //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
          //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
          //}
          P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
          P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
          P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
          P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
          P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
          P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
          data_index[processor_number]++;
          if (subscript[1] != 0) {
            processor_number=Ny*Nx*processor_comp[2]+Nx*neighbour_y+neighbour_x;
#ifdef PARTICLE_ID
            P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif 
            //for (k=0; k<3; k++) {
            //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
            //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
            //}
            P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
            P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
            P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
            P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
            P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
            P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
            data_index[processor_number]++;
            if (subscript[2] != 0) {
              processor_number=Ny*Nx*neighbour_z+Nx*neighbour_y+neighbour_x;
#ifdef PARTICLE_ID
              P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
              //for (k=0; k<3; k++) {
              //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
              //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
              //}
              P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
              P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
              P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
              P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
              P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
              P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
              data_index[processor_number]++;
            }
          }
          if (subscript[2] != 0) {
            processor_number=Ny*Nx*neighbour_z+Nx*processor_comp[1]+neighbour_x;
#ifdef PARTICLE_ID
            P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
            //for (k=0; k<3; k++) {
            //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
            //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
            //}
            P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
            P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
            P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
            P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
            P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
            P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
            data_index[processor_number]++;
          }
        }
        if (subscript[1] != 0) {
          processor_number=Ny*Nx*processor_comp[2]+Nx*neighbour_y+processor_comp[0];
#ifdef PARTICLE_ID
          P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
          //for (k=0; k<3; k++) {
          //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
          //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
          //}
          P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
          P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
          P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
          P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
          P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
          P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
          data_index[processor_number]++;
          if (subscript[2] != 0) {
            processor_number=Ny*Nx*neighbour_z+Nx*neighbour_y+processor_comp[0];
#ifdef PARTICLE_ID
            P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
            //for (k=0; k<3; k++) {
            //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
            //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
            //}
            P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
            P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
            P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
            P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
            P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
            P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
            data_index[processor_number]++;
          }
        }
        if (subscript[2] != 0) {
          processor_number=Ny*Nx*neighbour_z+Nx*processor_comp[1]+processor_comp[0];
#ifdef PARTICLE_ID
          P_sorted_id[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].ID;
#endif
          //for (k=0; k<3; k++) {
          //  P_sorted_pos[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Pos[k];
          //  P_sorted_vel[cnparticles[processor_number]+data_index[processor_number]].Var[k] = P_file[j].Vel[k];
          //}
          P_sorted_posx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[0];
          P_sorted_posy[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[1];
          P_sorted_posz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Pos[2];
          P_sorted_velx[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[0];
          P_sorted_vely[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[1];
          P_sorted_velz[cnparticles[processor_number]+data_index[processor_number]] = P_file[j].Vel[2];
          data_index[processor_number]++;
        }
      }
      free(data_index);
      free(P_file);
    }

#ifdef PARTICLE_ID
    P_temp_id = (unsigned long long *)malloc(nparticles_recv_sum*sizeof(unsigned long long));
#endif
    //P_temp_pos = (struct part_data_half *)malloc(nparticles_recv_sum*sizeof(struct part_data_half));
    //P_temp_vel = (struct part_data_half *)malloc(nparticles_recv_sum*sizeof(struct part_data_half));
    P_temp_posx = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));
    P_temp_posy = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));
    P_temp_posz = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));
    P_temp_velx = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));
    P_temp_vely = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));
    P_temp_velz = (part_datatype *)malloc(nparticles_recv_sum*sizeof(part_datatype));

    if (ThisTask == 0) {
      printf("Transferring data...\n");
      fflush(stdout);
    }

#ifdef PARTICLE_ID

    for (k=0; k<NTask; k++) {
      nparticles[k]       *= sizeof(unsigned long long);
      nparticles_recv[k]  *= sizeof(unsigned long long);
      cnparticles[k]      *= sizeof(unsigned long long);
      cnparticles_recv[k] *= sizeof(unsigned long long);
    }

    ierr = MPI_Alltoallv(&(P_sorted_id[0]), &(nparticles[0]), &(cnparticles[0]), MPI_BYTE, &(P_temp_id[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_BYTE, MPI_COMM_WORLD);

    for (k=0; k<NTask; k++) {
      nparticles[k]       /= sizeof(unsigned long long);
      nparticles_recv[k]  /= sizeof(unsigned long long);
      cnparticles[k]      /= sizeof(unsigned long long);
      cnparticles_recv[k] /= sizeof(unsigned long long);
    }

#endif

    //for (k=0; k<NTask; k++) {
    //  nparticles[k]       *= sizeof(struct part_data_half);
    //  nparticles_recv[k]  *= sizeof(struct part_data_half);
    //  cnparticles[k]      *= sizeof(struct part_data_half);
    //  cnparticles_recv[k] *= sizeof(struct part_data_half);
    //}
    
    //ierr = MPI_Alltoallv(&(P_sorted_pos[0]), &(nparticles[0]), &(cnparticles[0]), MPI_PARTDATA, &(P_temp_pos[0]),
    //                     &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_PARTDATA, MPI_COMM_WORLD);

    //ierr = MPI_Alltoallv(&(P_sorted_vel[0]), &(nparticles[0]), &(cnparticles[0]), MPI_PARTDATA, &(P_temp_vel[0]),
    //                      &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_PARTDATA, MPI_COMM_WORLD);
  
#ifdef MEMORY_MODE
    
    ierr = MPI_Alltoallv(&(P_sorted_posx[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_posx[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
                         
    ierr = MPI_Alltoallv(&(P_sorted_posy[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_posy[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
                         
    ierr = MPI_Alltoallv(&(P_sorted_posz[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_posz[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
                         
    ierr = MPI_Alltoallv(&(P_sorted_velx[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_velx[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_vely[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_vely[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_velz[0]), &(nparticles[0]), &(cnparticles[0]), MPI_FLOAT, &(P_temp_velz[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_FLOAT, MPI_COMM_WORLD);
                         
#else

    ierr = MPI_Alltoallv(&(P_sorted_posx[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_posx[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_posy[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_posy[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_posz[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_posz[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_velx[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_velx[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
          
    ierr = MPI_Alltoallv(&(P_sorted_vely[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_vely[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
        
    ierr = MPI_Alltoallv(&(P_sorted_velz[0]), &(nparticles[0]), &(cnparticles[0]), MPI_DOUBLE, &(P_temp_velz[0]),
                         &(nparticles_recv[0]), &(cnparticles_recv[0]), MPI_DOUBLE, MPI_COMM_WORLD);
                               
#endif

//#ifdef LIGHTCONE
//#ifndef UNFORMATTED
//    free(nparticles);
//#endif
//#else
    free(nparticles);
//#endif
    free(cnparticles);
#ifdef PARTICLE_ID
    free(P_sorted_id);
#endif
    //free(P_sorted_pos);
    //free(P_sorted_vel);
    free(P_sorted_posx);
    free(P_sorted_posy);
    free(P_sorted_posz);
    free(P_sorted_velx);
    free(P_sorted_vely);
    free(P_sorted_velz);

    //for (k=0; k<NTask; k++) {
    //  nparticles_recv[k]  /= sizeof(struct part_data_half);
    //  cnparticles_recv[k] /= sizeof(struct part_data_half);
    //}

#ifdef PERIODIC
    if (ThisTask == 0) {
      printf("\nWrapping particles over periodic boundaries...\n");
      fflush(stdout);
    }
#endif

    for (k=0; k<NTask; k++) {
      for (q=0; q<nparticles_recv[k]-1; q++) {
        m=cnparticles_recv[k]+q;
        n=nparticles_tot+(unsigned int)q;
#ifdef PERIODIC
        // If we are using periodic boundary conditions then we have to wrap the particles 
        // passed over the periodic boundaries based on the boxsize
        /*if (P_temp_pos[m].Var[0] > rmax_buff[0]) P_temp_pos[m].Var[0]=rmin[0]-((Lxmax-Lxmin)-P_temp_pos[m].Var[0]);
        if (P_temp_pos[m].Var[0] < rmin_buff[0]) P_temp_pos[m].Var[0]=rmax[0]+P_temp_pos[m].Var[0];
        if (P_temp_pos[m].Var[1] > rmax_buff[1]) P_temp_pos[m].Var[1]=rmin[1]-((Lymax-Lymin)-P_temp_pos[m].Var[1]);
        if (P_temp_pos[m].Var[1] < rmin_buff[1]) P_temp_pos[m].Var[1]=rmax[1]+P_temp_pos[m].Var[1];
        if (P_temp_pos[m].Var[2] > rmax_buff[2]) P_temp_pos[m].Var[2]=rmin[2]-((Lzmax-Lzmin)-P_temp_pos[m].Var[2]);
        if (P_temp_pos[m].Var[2] < rmin_buff[2]) P_temp_pos[m].Var[2]=rmax[2]+P_temp_pos[m].Var[2];*/
        if (P_temp_posx[m] > rmax_buff[0]) P_temp_posx[m]=rmin[0]-((Lxmax-Lxmin)-P_temp_posx[m]);
        if (P_temp_posx[m] < rmin_buff[0]) P_temp_posx[m]=rmax[0]+P_temp_posx[m];
        if (P_temp_posy[m] > rmax_buff[1]) P_temp_posy[m]=rmin[1]-((Lymax-Lymin)-P_temp_posy[m]);
        if (P_temp_posy[m] < rmin_buff[1]) P_temp_posy[m]=rmax[1]+P_temp_posy[m];
        if (P_temp_posz[m] > rmax_buff[2]) P_temp_posz[m]=rmin[2]-((Lzmax-Lzmin)-P_temp_posz[m]);
        if (P_temp_posz[m] < rmin_buff[2]) P_temp_posz[m]=rmax[2]+P_temp_posz[m];
        
#endif
#ifdef PARTICLE_ID
        P[n].ID = P_temp_id[m];
#endif
        //for (j=0; j<3; j++) {
        //  P[n].Pos[j] = P_temp_pos[m].Var[j];
        //  P[n].Vel[j] = P_temp_vel[m].Var[j];
        //}
        P[n].Pos[0] = P_temp_posx[m];
        P[n].Pos[1] = P_temp_posy[m];
        P[n].Pos[2] = P_temp_posz[m];
        P[n].Vel[0] = P_temp_velx[m];
        P[n].Vel[1] = P_temp_vely[m];
        P[n].Vel[2] = P_temp_velz[m];
   
        // This is a simple check that each processor does to make sure that every particle on the processor 
        // is within the previously calculated boundaries for that processor
        if((P[n].Pos[0] > rmax_buff[0]) || (P[n].Pos[0] < rmin_buff[0])) { printf("\nERROR: Task %d has incorrect particle in the x-direction: %lf, %lf, %lf\n\n", ThisTask, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]); FatalError("read_data.c", 408);}
        if((P[n].Pos[1] > rmax_buff[1]) || (P[n].Pos[1] < rmin_buff[1])) { printf("\nERROR: Task %d has incorrect particle in the y-direction: %lf, %lf, %lf\n\n", ThisTask, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]); FatalError("read_data.c", 409);}
        if((P[n].Pos[2] > rmax_buff[2]) || (P[n].Pos[2] < rmin_buff[2])) { printf("\nERROR: Task %d has incorrect particle in the z-direction: %lf, %lf, %lf\n\n", ThisTask, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2]); FatalError("read_data.c", 410);}

      }
      nparticles_tot += (unsigned int)(nparticles_recv[k]-1);
    }

//#ifdef LIGHTCONE
//#ifndef UNFORMATTED
//    free(nparticles_recv);
//#endif
//#else
    free(nparticles_recv);
//#endif
    free(cnparticles_recv);
#ifdef PARTICLE_ID
    free(P_temp_id);
#endif
    //free(P_temp_pos);
    //free(P_temp_vel);
    free(P_temp_posx);
    free(P_temp_posy);
    free(P_temp_posz);
    free(P_temp_velx);
    free(P_temp_vely);
    free(P_temp_velz);

    ierr = MPI_Barrier(MPI_COMM_WORLD);
    if (ThisTask == 0) { 
      printf("------------------------------------\n\n");
      fflush(stdout);
    }

//#ifdef LIGHTCONE
//#ifdef UNFORMATTED
//  }
//  free(nparticles);
//  free(nparticles_recv);
//#endif
//#endif
    
  }

  // If Inputstyle!=0 we need to deallocate the 'filelist'.
  if (InputStyle) {
    for (k=0; k<nread; k++) {
      if (readid[k] == ThisTask) {
        for (q=0; q<ninputfiles; q++) free(filelist[q]);
        free(filelist);
      }
    }
  }
  free(readid);

#ifdef LIGHTCONE
#ifdef UNFORMATTED
  // Free the info file data if we used it
  free(npartfiles);
  free(filenumbers);
#endif
#endif
      
  ierr = MPI_Barrier(MPI_COMM_WORLD);
  printf("Task %d has %d particles in total...\n", ThisTask, nparticles_tot);
  fflush(stdout);
  ierr = MPI_Reduce(&Nout, &Nout_glob, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ThisTask == 0) {
    printf("There were %u particles outside the simulation box...\n", Nout_glob);
    fflush(stdout);
  }

  return;
}
