// This file contains all subroutines associated with MPI-based file input and data organisation

#include "vars.h"

void Read_Data_Info(void) {

  // This routine uses an info file containing the particle extents on each file to determine which files
  // each processor needs to read in. The info file must contain a list with each row corresponding to a different
  // slice and columns filenum, xmin, ymin, zmin, xmax, ymax, zmax, number of particles. The info file must
  // be named in the same way as the input dark matter files, but appended with .info rather than a processor
  // number, .i.e., InputDir/InputFileBase.info
  
  FILE * fp;
  char infofile[500], buf[500];
  int filenumber, * filenumbers, nfiles, readflag_glob;
  int i, nprocgroup, groupTask, masterTask, * readflag;
  unsigned int j, npart, * npartfile, * Nout, Nout_glob, Nout_tot=0;
  double xmin, xmax, ymin, ymax, zmin, zmax;
  struct part_data P_file; 
#ifdef PERIODIC
  struct part_data P_wrap;
#endif
#ifdef LIGHTCONE
#ifdef UNFORMATTED
  int k, dummy;
  unsigned int nchunk;
#endif
#else
#ifdef GADGET_STYLE
  int k, dummy;
  struct part_data * P_tmp;
#endif
#endif

  nfiles=0;
  nparticles_tot=0;
  if (ThisTask == 0) printf("%d tasks reading in %d files...\n", nread, ninputfiles);

  readflag = (int *)malloc(ninputfiles*sizeof(int));
  Nout = (unsigned int*)malloc(ninputfiles*sizeof(unsigned int));

  for (i=0; i<ninputfiles; i++) {
    readflag[i] = 0;
    Nout[i] = 0;
  }

  // Master loop to ensure that we only have a certain number of processors reading in at any one point
  nprocgroup = NTask / nread;
  if (NTask % nread) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {

      // Read in the info file and determine which files we need to read in 
      sprintf(infofile, "%s/%s.info", InputDir, InputFileBase);
      if(!(fp = fopen(infofile, "r"))) {
        printf("\nERROR: Can't open info file '%s'.\n\n", buf);
        FatalError("read_data_info.c", 33);
      }

      filenumbers = (int *)malloc(ninputfiles*sizeof(int));
      npartfile   = (unsigned int *)calloc(ninputfiles,sizeof(unsigned int));

      while(fgets(buf,500,fp)) {

        if(strncmp(buf,"#",1)==0) continue;
      
        if(sscanf(buf,"%d %lf %lf %lf %lf %lf %lf %u\n",&filenumber,&xmin,&ymin,&zmin,&xmax,&ymax,&zmax,&npart)!=8) { printf("Bad input format in info file '%s\n", infofile); FatalError("read_data_info.c", 41);}  
        npartfile[filenumber] += npart;
        if((rmin_buff[0] <= xmax) && (rmax_buff[0] >= xmin) && (rmin_buff[1] <= ymax) && (rmax_buff[1] >= ymin) && (rmin_buff[2] <= zmax) && (rmax_buff[2] >= zmin)) {
          filenumbers[filenumber] = 1;
#ifdef PERIODIC
        } else {
          if (Local_nx == 0) {
            xmin -= (Lxmax-Lxmin);
            xmax -= (Lxmax-Lxmin);
          } else if (Local_nx == Nx-1) {
            xmin += (Lxmax-Lxmin);
            xmax += (Lxmax-Lxmin);
          }
          if (Local_ny == 0) {
            ymin -= (Lymax-Lymin);
            ymax -= (Lymax-Lymin);
          } else if (Local_ny == Ny-1) {
            ymin += (Lymax-Lymin);
            ymax += (Lymax-Lymin);
          }
          if (Local_nz == 0) {
            zmin -= (Lzmax-Lzmin);
            zmax -= (Lzmax-Lzmin);
          } else if (Local_nz == Nz-1) {
            zmin += (Lzmax-Lzmin);
            zmax += (Lzmax-Lzmin);
          }

          if((rmin_buff[0] <= xmax) && (rmax_buff[0] >= xmin) && (rmin_buff[1] <= ymax) && (rmax_buff[1] >= ymin) && (rmin_buff[2] <= zmax) && (rmax_buff[2] >= zmin)) filenumbers[filenumber] = 1;         
#endif
        }
      }
      fclose(fp);

      // Now loop over all the necessary input files and save the required particles
      for (i=0; i<ninputfiles; i++) {
 
        if (!(filenumbers[i])) continue;

        readflag[i]++;
        sprintf(buf, "%s/%s.%d", InputDir, InputFileBase, i);

// LIGHTCONE simulations
// =====================
#ifdef LIGHTCONE

// Binary
#ifdef UNFORMATTED
       
        if(!(fp = fopen(buf, "rb"))) {
          printf("\nERROR: Can't open input file '%s'.\n\n", buf);
          FatalError("read_data_info.c", 80);
        }        

        npart = 0;
        while(npart != npartfile[i]) {
          my_fread(&dummy, sizeof(dummy), 1, fp);
          my_fread(&nchunk, sizeof(unsigned int), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          my_fread(&dummy, sizeof(dummy), 1, fp);
          for (j=0; j<nchunk; j++) {
#ifdef PARTICLE_ID
            my_fread(&(P_file.ID), sizeof(unsigned long long), 1, fp);
#endif
            for (k=0; k<3; k++) my_fread(&(P_file.Pos[k]), sizeof(float), 1, fp);
            for (k=0; k<3; k++) my_fread(&(P_file.Vel[k]), sizeof(float), 1, fp);
 
// ASCII
#else

        if(!(fp = fopen(buf, "r"))) {
          printf("\nERROR: Can't open input file '%s'.\n\n", buf);
          FatalError("read_data_info.c", 97);
        }
        for(j=0; j<npartfile[i]; j++) {
#ifdef MEMORY_MODE
#ifdef PARTICLE_ID
          if((fscanf(fp, "%llu %f %f %f %f %f %f\n", &(P_file.ID), &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                     &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 7)) {
#else
          if((fscanf(fp, "%f %f %f %f %f %f\n", &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 6)) {
#endif
#else
#ifdef PARTICLE_ID
          if((fscanf(fp, "%llu %lf %lf %lf %lf %lf %lf\n", &(P_file.ID), &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                           &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 7)) {
#else
          if((fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                      &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 6)) {
#endif
#endif
            printf("Task %d has error reading file %s\n", ThisTask, buf); 
            FatalError("read_data_info.c", 326);
          }

#endif

// Snapshot simulations
// ====================
#else

// Binary
#ifdef GADGET_STYLE

        P_tmp = (struct part_data *)malloc(npartfile[i]*sizeof(struct part_data));

        if(!(fp = fopen(buf, "rb"))) {
          printf("\nERROR: Can't open input file '%s'.\n\n", buf);
          FatalError("read_data_info.c", 90);
        }
        my_fread(&dummy, sizeof(dummy), 1, fp);
        my_fread(&header, sizeof(header), 1, fp);
        my_fread(&dummy, sizeof(dummy), 1, fp);

        // Read in the particle positions
        my_fread(&dummy, sizeof(dummy), 1, fp);
        for(j=0; j<npartfile[i]; j++) {
          for (k=0; k<3; k++) my_fread(&(P_tmp[j].Pos[k]), sizeof(float), 1, fp);
        }
        my_fread(&dummy, sizeof(dummy), 1, fp);

        // Read in the particle velocities
        my_fread(&dummy, sizeof(dummy), 1, fp);
        for(j=0; j<npartfile[i]; j++) {
          for (k=0; k<3; k++) my_fread(&(P_tmp[j].Vel[k]), sizeof(float), 1, fp);
        }
        my_fread(&dummy, sizeof(dummy), 1, fp);

#ifdef PARTICLE_ID
        // Read in the particle ids
        my_fread(&dummy, sizeof(dummy), 1, fp);
        for(j=0; j<npartfile[i]; j++) my_fread(&(P_tmp[j].ID), sizeof(unsigned long long), 1, fp);
        my_fread(&dummy, sizeof(dummy), 1, fp);
#endif

        for (j=0; j<npartfile[i]; j++) {
          P_file = P_tmp[j];

// ASCII
#else

        if(!(fp = fopen(buf, "r"))) {
          printf("\nERROR: Can't open input file '%s'.\n\n", buf);
          FatalError("read_data_info.c", 97);
        }
        for(j=0; j<npartfile[i]; j++) {
#ifdef MEMORY_MODE
#ifdef PARTICLE_ID
          if((fscanf(fp, "%llu %f %f %f %f %f %f\n", &(P_file.ID), &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                     &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 7)) {
#else
          if((fscanf(fp, "%f %f %f %f %f %f\n", &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 6)) {
#endif
#else
#ifdef PARTICLE_ID
          if((fscanf(fp, "%llu %lf %lf %lf %lf %lf %lf\n", &(P_file.ID), &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                           &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 7)) {
#else
          if((fscanf(fp, "%lf %lf %lf %lf %lf %lf\n", &(P_file.Pos[0]), &(P_file.Pos[1]), &(P_file.Pos[2]), 
                                                      &(P_file.Vel[0]), &(P_file.Vel[1]), &(P_file.Vel[2])) != 6)) {
#endif
#endif
            printf("Task %d has error reading file %s\n", ThisTask, buf); 
            FatalError("read_data_info.c", 326);
          }

#endif
#endif

          if((P_file.Pos[0] < Lxmin) || (P_file.Pos[0] > Lxmax) || (P_file.Pos[1] < Lymin) || (P_file.Pos[1] > Lymax) || (P_file.Pos[2] < Lzmin) || (P_file.Pos[2] > Lzmax)) {
            Nout[i]++;
            continue;
          }

          // Only keep those particles that belong on this processor
          if((P_file.Pos[0] >= rmin_buff[0]) && (P_file.Pos[0] <= rmax_buff[0]) && 
             (P_file.Pos[1] >= rmin_buff[1]) && (P_file.Pos[1] <= rmax_buff[1]) && 
             (P_file.Pos[2] >= rmin_buff[2]) && (P_file.Pos[2] <= rmax_buff[2])) {
            P[nparticles_tot] = P_file;
            nparticles_tot++;
            if (nparticles_tot >= maxparticles) {
              printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
              printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
              FatalError("read_data_info.c", 138);
            }
          }
#ifdef PERIODIC
          if (Local_nx == 0) {
            P_wrap = P_file;
            P_wrap.Pos[0] = P_file.Pos[0]-(Lxmax-Lxmin);
            if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
               (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
               (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
              P[nparticles_tot] = P_wrap;
              nparticles_tot++;
              if (nparticles_tot >= maxparticles) {
                printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                FatalError("read_data_info.c", 138);
              }
            }
            if (Local_ny == 0) {
              P_wrap.Pos[1] = P_file.Pos[1]-(Lymax-Lymin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
              if (Local_nz == 0) {
                P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              } else if (Local_nz == Nz-1) {
                P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              }
            } else if (Local_ny == Ny-1) {
              P_wrap.Pos[1] = P_file.Pos[1]+(Lymax-Lymin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
              if (Local_nz == 0) {
                P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot > maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              } else if (Local_nz == Nz-1) {
                P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              }
            }
            P_wrap.Pos[1] = P_file.Pos[1];
            if (Local_nz == 0) {
              P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            } else if (Local_nz == Nz-1) {
              P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            }
          } else if (Local_nx == Nx-1) {
            P_wrap = P_file;
            P_wrap.Pos[0] = P_file.Pos[0]+(Lxmax-Lxmin);
            if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
               (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
               (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
              P[nparticles_tot] = P_wrap;
              nparticles_tot++;
              if (nparticles_tot >= maxparticles) {
                printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                FatalError("read_data_info.c", 138);
              }
            }
            if (Local_ny == 0) {
              P_wrap.Pos[1] = P_file.Pos[1]-(Lymax-Lymin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
              if (Local_nz == 0) {
                P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              } else if (Local_nz == Nz-1) {
                P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              }
            } else if (Local_ny == Ny-1) {
              P_wrap.Pos[1] = P_file.Pos[1]+(Lymax-Lymin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
              if (Local_nz == 0) {
                P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              } else if (Local_nz == Nz-1) {
                P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              }
            }
            P_wrap.Pos[1] = P_file.Pos[1];
            if (Local_nz == 0) {
              P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            } else if (Local_nz == Nz-1) {
              P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            }
          }
          P_wrap = P_file;
          if (Local_ny == 0) {
            P_wrap.Pos[1] = P_file.Pos[1]-(Lymax-Lymin);
            if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
               (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
               (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
              P[nparticles_tot] = P_wrap;
              nparticles_tot++;
              if (nparticles_tot >= maxparticles) {
                printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                FatalError("read_data_info.c", 138);
              }
            }
            if (Local_nz == 0) {
              P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            } else if (Local_nz == Nz-1) {
              P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            }
          } else if (Local_ny == Ny-1) {
              P_wrap.Pos[1] = P_file.Pos[1]+(Lymax-Lymin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
              if (Local_nz == 0) {
                P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              } else if (Local_nz == Nz-1) {
                P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
                if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                   (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                   (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                  P[nparticles_tot] = P_wrap;
                  nparticles_tot++;
                  if (nparticles_tot >= maxparticles) {
                    printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                    printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                    FatalError("read_data_info.c", 138);
                  }
                }
              }
            }
            P_wrap = P_file;
            if (Local_nz == 0) {
              P_wrap.Pos[2] = P_file.Pos[2]-(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            } else if (Local_nz == Nz-1) {
              P_wrap.Pos[2] = P_file.Pos[2]+(Lzmax-Lzmin);
              if((P_wrap.Pos[0] >= rmin_buff[0]) && (P_wrap.Pos[0] <= rmax_buff[0]) && 
                 (P_wrap.Pos[1] >= rmin_buff[1]) && (P_wrap.Pos[1] <= rmax_buff[1]) && 
                 (P_wrap.Pos[2] >= rmin_buff[2]) && (P_wrap.Pos[2] <= rmax_buff[2])) {
                P[nparticles_tot] = P_wrap;
                nparticles_tot++;
                if (nparticles_tot >= maxparticles) {
                  printf("\nERROR: Number of particles of Task %d greater than maximum number of particles, \n", ThisTask);      
                  printf("       please check Px, Py and Pz are correct then increase the Buffer\n");      
                  FatalError("read_data_info.c", 138);
                }
              }
            }
#endif

#ifdef LIGHTCONE
#ifdef UNFORMATTED
          }
          npart += nchunk;
          my_fread(&dummy, sizeof(dummy), 1, fp);
#endif
        }
#else
        }
#ifdef GADGET_STYLE
        free(P_tmp);
#endif
#endif

        fclose(fp);
      }
      free(npartfile);
      free(filenumbers);

      printf("Task %d has %d particles in total...\n", ThisTask, nparticles_tot);

    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (i=0; i<ninputfiles; i++) {
    ierr = MPI_Reduce(&(readflag[i]), &readflag_glob, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    ierr = MPI_Reduce(&(Nout[i]), &Nout_glob, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    if (ThisTask == 0) Nout_tot += Nout_glob/readflag_glob;
  }
  if (ThisTask == 0) printf("There were %u particles outside the simulation box...\n", Nout_tot);

  free(Nout);
  free(readflag);

  return;
}
