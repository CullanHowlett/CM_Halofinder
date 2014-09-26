// This is the main FoF routine. The Friends-of-Friends algorithm works by chaining together particles that are 
// seperated by less than the theoretically motivated linking length. Each group of particles is then specified as
// a halo so long as the number of constituent particles in the particle chain is greater than some number. In our case
// this is specified as an run parameter, but typical values include 10 or 20 particles. To speed this algorithm up we
// divide the region on the processor into cells, where each cell contains, on average, 1 particle. 

// In the algorithm each particle has an associated 'pointer' (not an actual pointer), which it starts by pointing to
// itself. We then loop over all particles and, inside that, we loop over all the other particles in the cell and then 
// over all forward cells. If two particles are within the linking length, we point each particle's pointer to where the
// other particle was pointing, hence creating a chain/loop of particles. This also serves to connect two chains together 
// if any of the particles in one chain are within a linking length of any particles in the other chain. Once we have 
// chained all the particles together we than identify the halos and give each halo a head particle. We then output the
// halos by calculating the position of the halo's centre of mass and if this is within the processor boundaries we
// follow the particle chain, outputting the data for each constituent particle.

#include "vars.h"

int nhalos;
int *nparthalo, *ihalo;
unsigned int * next;
double linksq;

void FOF(void) {

  gsl_rng * rgen;
  int nx, ny, nz;
  int nphalo, nhalos;
  int ix, iy, iz, ix1, iy1, ix2, iy2, iz2;
  unsigned int i, ip, ip2, ind;
  unsigned int ninhalo=0, nouthalo=0;
  unsigned int * head, *chain, * first;
  unsigned int * inhalo=NULL, * outhalo=NULL;
  double lcell[3];
  double sampbuffer=1.5;
#ifdef VARLINK
  double distance;
#endif 

  if (ThisTask == 0) printf("Creating subcells...\n");
  
  rgen = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(rgen, Seed);

  // set grid of cells
#ifdef PERIODIC 
  lcell[0]=(Lxmax-Lxmin)/Px;
  lcell[1]=(Lymax-Lymin)/Py;
  lcell[2]=(Lzmax-Lzmin)/Pz;
#else
  lcell[0]=(Lxmax-Lxmin+2.0*boundarysize)/Px;
  lcell[1]=(Lymax-Lymin+2.0*boundarysize)/Py;
  lcell[2]=(Lzmax-Lzmin+2.0*boundarysize)/Pz;
#endif

#ifdef VARLINK
  // Creates the linking length lookup table
  if (ThisTask == 0) printf("Creating linking length lookup table...\n");
  Create_Link();
#else
  linksq = linklength*linklength;
#endif

  nx=(int)ceil((rmax_buff[0]-rmin_buff[0])/lcell[0])+1;
  ny=(int)ceil((rmax_buff[1]-rmin_buff[1])/lcell[1])+1;
  nz=(int)ceil((rmax_buff[2]-rmin_buff[2])/lcell[2])+1;
  
  next  = (unsigned int *)malloc(nparticles_tot*sizeof(unsigned int));  // halo loop
  chain = (unsigned int *)calloc(nparticles_tot,sizeof(unsigned int));  // chain in cell
  first = (unsigned int *)calloc(nx*ny*nz,sizeof(unsigned int));        // first particle in cell
 
  // These now start at 1 (allows us to use unsigned ints)
  next--;
  chain--;

  // set chain
  for (i=1; i<=nparticles_tot; i++) {
    next[i] = i;

    ix=(int)floor((P[i-1].Pos[0]-rmin_buff[0])/lcell[0]);
    iy=(int)floor((P[i-1].Pos[1]-rmin_buff[1])/lcell[1]);
    iz=(int)floor((P[i-1].Pos[2]-rmin_buff[2])/lcell[2]);

    ind = (unsigned int)iz*(unsigned int)ny*(unsigned int)nx+(unsigned int)iy*(unsigned int)nx+(unsigned int)ix;

    ip=first[ind];
    if (ip == 0) {
      first[ind]=i;
    } else {
      chain[i]=ip;
      first[ind]=i;
    }
  } 

  // HALO FINDER
  // look at particles cell by cell
  for (iz=0; iz<nz; iz++) {
    for (iy=0; iy<ny; iy++) {
      for (ix=0; ix<nx; ix++) {
        ip=first[(unsigned int)iz*(unsigned int)ny*(unsigned int)nx+(unsigned int)iy*(unsigned int)nx+(unsigned int)ix];
        if (ip == 0) continue;

#ifdef VARLINK
        // Calculates the distance of the cell from the origin and uses the lookup table to return the 
        // square of the linking length
        distance=(rmin_buff[0]-Origin_x+(ix-0.5)*lcell[0])*(rmin_buff[0]-Origin_x+(ix-0.5)*lcell[0]) +
                 (rmin_buff[1]-Origin_y+(iy+0.5)*lcell[1])*(rmin_buff[1]-Origin_y+(iy+0.5)*lcell[1]) +
                 (rmin_buff[2]-Origin_z+(iz+0.5)*lcell[2])*(rmin_buff[2]-Origin_z+(iz+0.5)*lcell[2]);
        gsl_spline_eval(link_spline, distance, link_acc);
#endif  

        //follow this chain
        do {
          //scan this cell for friends
          ip2=ip;
          do {
            ip2=chain[ip2];
            if (ip2 == 0) break; 
            Befriend(ip, ip2);
          } while(1);

          //scan forward cells for friends     
          for (iz2=iz; iz2<=iz+1; iz2++) {
            if (iz2 >= nz) continue;
            if (iz2 == iz) {
              iy1=iy;
            } else {
              iy1=iy-1;
            }             

            for (iy2=iy1; iy2<=iy+1; iy2++) {
              if ((iy2 >= ny) || (iy2 < 0)) continue;
              if ((iz2 == iz) && (iy2 == iy)) { 
                ix1 = ix;
              } else {
                ix1 = ix-1;
              }

              for (ix2=ix1; ix2<=ix+1; ix2++) {
                if ((ix2 >= nx) || (ix2 < 0)) continue;
                if ((ix2 == ix) && (iy2 == iy) && (iz2 == iz)) continue;
                     
                // do the neighbour cell
                ip2=first[(unsigned int)iz2*(unsigned int)ny*(unsigned int)nx+(unsigned int)iy2*(unsigned int)nx+(unsigned int)ix2];
                if (ip2 == 0) continue;
                Befriend(ip, ip2);                    
                do { 
                  ip2=chain[ip2];
                  if (ip2 == 0) break;    
                  Befriend(ip, ip2);
                } while(1);  
              }
            }
          }

          ip=chain[ip];
          if (ip == 0) break; 
        } while(1);

      }
    }
  }
  
  // RETRIEVE HALOS
  chain++;
  free(chain);
  free(first);
#ifdef VARLINK
  gsl_spline_free(link_spline);
  gsl_interp_accel_free(link_acc);
#endif
  
  
  if (SampInHalos>0)  {
    unsigned int dummy = (unsigned int)rint(sampbuffer*nparticles_tot*SampInHalos);
    if (dummy > nparticles_tot) dummy = nparticles_tot;
    inhalo  = (unsigned int *)malloc(dummy*sizeof(int));
  }
  if (SampOutHalos>0) {
    unsigned int dummy = (unsigned int)rint(sampbuffer*nparticles_tot*SampOutHalos);
    if (dummy > nparticles_tot) dummy = nparticles_tot;
    outhalo = (unsigned int *)malloc(dummy*sizeof(int));
  }
  head = (unsigned int *)calloc(nparticles_tot,sizeof(unsigned int));
  head--;

  if (ThisTask == 0) printf("Retrieving halos...\n");

  // first count heads
  nhalos=0;
  for (i=1; i<=nparticles_tot; i++) {
    if (head[i] == 0) { 
      ip=i;
      nphalo=0;

      do {
        ip=next[ip];
        nphalo++;
        if (ip == i) break;
        head[ip]=1;
      } while(1);

      if (nphalo < nphalomin) {
        head[ip]=1;
        
        // Subsample particles not in halos
        if(SampOutHalos>0) {
          do {
            ip=next[ip];
            if (gsl_rng_uniform(rgen) < SampOutHalos) {
              // We only want to output particles that are within the processor boundaries 
#ifdef PERIODIC
              if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0] >= rmax[0])) break;
              if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1] >= rmax[1])) break;
              if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2] >= rmax[2])) break;
#else
              if (Local_nx == Nx-1) {
                if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0]  > rmax[0])) break;
              } else {
                if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0] >= rmax[0])) break;
              }
              if (Local_ny == Ny-1) {
                if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1]  > rmax[1])) break;
              } else {
                if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1] >= rmax[1])) break;
              }
              if (Local_nz == Nz-1) {
                if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2]  > rmax[2])) break;
              } else {
                if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2] >= rmax[2])) break;
              }
#endif
              outhalo[nouthalo]=ip;
              nouthalo++;
              if((nouthalo>=sampbuffer*nparticles_tot*SampOutHalos) || (nouthalo>=nparticles_tot)) {
                printf("\nERROR: Task %d has subsampled more particles than it has memory for.\n", ThisTask);
                printf("       This is unexpected, but can be avoided by increasing the parameter sampbuffer in FOF.c, line 27.\n\n");
                FatalError((char *)"FOF.c", 237); 
              }
            }
            if (ip == i) break;
          } while(1);
        }
        
      } else {
        nhalos++;
        
        // Subsample particles in halos
        if(SampInHalos>0) {
          do {
            ip=next[ip];
            if (gsl_rng_uniform(rgen) < SampInHalos) {
              // We only want to output particles that are within the processor boundaries 
#ifdef PERIODIC
              if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0] >= rmax[0])) break;
              if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1] >= rmax[1])) break;
              if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2] >= rmax[2])) break;
#else
              if (Local_nx == Nx-1) {
                if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0]  > rmax[0])) break;
              } else {
                if ((P[ip-1].Pos[0] < rmin[0]) || (P[ip-1].Pos[0] >= rmax[0])) break;
              }
              if (Local_ny == Ny-1) {
                if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1]  > rmax[1])) break;
              } else {
                if ((P[ip-1].Pos[1] < rmin[1]) || (P[ip-1].Pos[1] >= rmax[1])) break;
              }
              if (Local_nz == Nz-1) {
                if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2]  > rmax[2])) break;
              } else {
                if ((P[ip-1].Pos[2] < rmin[2]) || (P[ip-1].Pos[2] >= rmax[2])) break;
              }
#endif
              inhalo[ninhalo]=ip;
              ninhalo++;
              if((ninhalo>=sampbuffer*nparticles_tot*SampInHalos) || (ninhalo>=nparticles_tot)) {
                printf("\nERROR: Task %d has subsampled more particles than it has memory for.\n", ThisTask);
                printf("       This is unexpected, but can be avoided by increasing the parameter sampbuffer in FOF.c, line 27.\n\n");
                FatalError((char *)"FOF.c", 279); 
              }
            }
            if (ip == i) break;
          } while(1);
        }
        
      }          
    } 
  }
     
  nparthalo = (int *)malloc(nhalos*sizeof(int));
  ihalo     = (int *)malloc(nhalos*sizeof(int));

  if (ThisTask == 0) printf("Checking the halos...\n\n");

  nhalos = 0;
  for (i=1; i<=nparticles_tot; i++) {
    if (head[i] == 0) Checkhalo(i);
  }

  head++;
  free(head);

  // Write out the subsampled data
  if (SampInHalos > 0) {
    if (ThisTask == 0) {
      printf("Outputting the particles in halos\n");
      printf("=================================\n\n");
    }
    Subsample(1, ninhalo, inhalo);
  }
  if (SampOutHalos > 0) {
    if (ThisTask == 0) {
      printf("Outputting the particles outside halos\n");
      printf("======================================\n\n");
    }
    Subsample(0, nouthalo, outhalo);
  }

  if (SampInHalos>0)  free(inhalo);
  if (SampOutHalos>0) free(outhalo);

  // Write out the halo data
  if (ThisTask == 0) {
    printf("Outputting the halos\n");
    printf("====================\n\n");
  }
  Output_Halos();

  // Write out the info file
  Output_Info(ninhalo, nouthalo);  

  next++;
  free(P);
  free(next);
  free(ihalo);
  free(nparthalo);

  return;
}

void Befriend(unsigned int ip, unsigned int ip2) {

  //befriend, each one points to where the other was pointing (true symmetry)

  int unjoined;
  unsigned int ip3;
  double dx, dy, dz, distsq;

  unjoined=1;  

  dx = P[ip-1].Pos[0]-P[ip2-1].Pos[0];
  dy = P[ip-1].Pos[1]-P[ip2-1].Pos[1];
  dz = P[ip-1].Pos[2]-P[ip2-1].Pos[2];

  distsq = dx*dx+dy*dy+dz*dz;

  if (distsq > linksq) return;

  ip3=ip;
  do {
    ip3=next[ip3];
    if (ip3 == ip) break;
    if (ip3 == ip2) { 
      unjoined=0;
      break;
    }
  } while(1);    

  if (unjoined) {
    ip3 = next[ip];
    next[ip] = next[ip2];
    next[ip2] = ip3;
  } 

  return;
}

void Checkhalo(unsigned int i) {

  int nphalo=0;
  unsigned int j;
  double xh=0, yh=0, zh=0;
      
  j=i;
  do {
    nphalo++;
    xh += P[j-1].Pos[0];
    yh += P[j-1].Pos[1];
    zh += P[j-1].Pos[2];
    j=next[j];
    if (j == i) break;
  } while(1);

  xh /= (double)nphalo;
  yh /= (double)nphalo;
  zh /= (double)nphalo;

  // We only want to output halos that have a centre of mass within the processor boundaries 
  if ((xh < rmin[0]) || (xh >= rmax[0])) return;
  if ((yh < rmin[1]) || (yh >= rmax[1])) return;
  if ((zh < rmin[2]) || (zh >= rmax[2])) return;

  // This checks that for any halos with a centre of mass on the processor none of the 
  // constituent particles extend beyond the boundary, otherwise the halos aren't
  // complete and the boundary regions must be made larger.
  j=i;
  do {
    if ((P[j-1].Pos[0]-rmin_buff[0] < linklength) || (rmax_buff[0]-P[j-1].Pos[0] < linklength)) {
      printf("\nWARNING: Particle very near boundary in x-direction. Halos may be incorrect\n");
      printf("         We highly recommend repeating the run with increased boundary size.\n\n");
    }
    if ((P[j-1].Pos[1]-rmin_buff[1] < linklength) || (rmax_buff[1]-P[j-1].Pos[1] < linklength)) {
      printf("\nWARNING: Particle very near boundary in y-direction. Halos may be incorrect\n");
      printf("         We highly recommend repeating the run with increased boundary size.\n\n");
    }
    if ((P[j-1].Pos[2]-rmin_buff[2] < linklength) || (rmax_buff[2]-P[j-1].Pos[2] < linklength)) {
      printf("\nWARNING: Particle very near boundary in z-direction. Halos may be incorrect\n");
      printf("         We highly recommend repeating the run with increased boundary size.\n\n");
    }
    j=next[j];
    if (j == i) break;
  } while(1);

  ihalo[nhalos]     = i;
  nparthalo[nhalos] = nphalo;

  nhalos++;

  return;
}

void Subsample(int in, unsigned int nsubsamp, unsigned int * subsamp) {
  
  FILE *fp; 
  char buf[300];
  int nprocgroup, groupTask, masterTask;
  unsigned int i, j;
#ifdef OUTPUT_UNFORMATTED
  int k, dummy;
  unsigned int pc, blockmaxlen;
  float * block;
#ifdef PARTICLE_ID
  unsigned long long * blockid;
#endif
#endif

  nprocgroup = NTask / nwrite;
  if (NTask % nwrite) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;

#ifdef OUTPUT_UNFORMATTED
  blockmaxlen = (unsigned int)(nparticles_tot/3);
#endif
  
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if (nsubsamp > 0) {
        if(in) {
          sprintf(buf, "%s_inhalo.%d", OutputFileBase, ThisTask);
        } else {
          sprintf(buf, "%s_outhalo.%d", OutputFileBase, ThisTask);
        }
        if(!(fp = fopen(buf, "w"))) {
          printf("\nError: Unable to open output file %s.\n\n", buf);
          FatalError((char *)"FOF.c", 465);
        }
        
        // If unformatted, output all positions then all velocities then all ID's, otherwise output one particle at a time.
#ifdef OUTPUT_UNFORMATTED
        dummy = sizeof(unsigned int);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        my_fwrite(&nsubsamp, sizeof(unsigned int), 1, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        block = (float *)malloc(nparticles_tot*sizeof(float));

        // Write coordinates
        dummy = sizeof(float)*3*nsubsamp;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(i=0, pc=0; i<nsubsamp; i++) {
          j = subsamp[i];
          for(k = 0; k < 3; k++) block[3*pc+k] = (float)(P[j-1].Pos[k]);
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3*pc, fp);
	        pc = 0;
	      }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3*pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // Write velocities
        dummy = sizeof(float)*3*nsubsamp;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(i=0, pc=0; i<nsubsamp; i++) {
          j = subsamp[i];
          for(k = 0; k < 3; k++) block[3*pc+k] = (float)(P[j-1].Vel[k]);
          pc++;
          if(pc == blockmaxlen) {
            my_fwrite(block, sizeof(float), 3*pc, fp);
	        pc = 0;
	      }
        }
        if(pc > 0) my_fwrite(block, sizeof(float), 3*pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // Write ID's (if necessary)
#ifdef PARTICLE_ID
        blockid = (unsigned long long *)block;
        blockmaxlen = (unsigned int)((nparticles_tot*sizeof(float)) / sizeof(unsigned long));

        dummy = sizeof(unsigned long long)*nsubsamp;
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        for(i=0, pc=0; i<nsubsamp; i++) {
          j = subsamp[i];
          blockid[pc] = P[j-1].ID;
          pc++;
          if(pc == blockmaxlen) {
	        my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
	        pc = 0;
	      }
        }
        if(pc > 0) my_fwrite(blockid, sizeof(unsigned long long), pc, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#endif

        free(block);   
#else
        for(i=0; i<nsubsamp; i++) {
          j=subsamp[i];
#ifdef PARTICLE_ID
          fprintf(fp,"%15llu %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P[j-1].ID, (float)P[j-1].Pos[0], (float)P[j-1].Pos[1], (float)P[j-1].Pos[2], (float)P[j-1].Vel[0], (float)P[j-1].Vel[1], (float)P[j-1].Vel[2]);
#else
          fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", (float)P[j-1].Pos[0], (float)P[j-1].Pos[1], (float)P[j-1].Pos[2], (float)P[j-1].Vel[0], (float)P[j-1].Vel[1], (float)P[j-1].Vel[2]);
#endif
        }
#endif
        fclose(fp);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }   
  return;
}

void Output_Halos(void) {

  FILE *fp; 
  char buf[300];
  int i, j, nprocgroup, groupTask, masterTask, halosize;
  unsigned int ind;
  double invnparthalo, xh, yh, zh, vxh, vyh, vzh;
  float * haloblock;
#ifdef INERTIA
  double xh2, yh2, zh2, xhyh, xhzh, yhzh;
#endif
#ifdef DISPERSION
  double vxh2, vyh2, vzh2, vxhvyh, vxhvzh, vyhvzh;
#endif
#ifdef OUTPUT_UNFORMATTED
  int dummy;
#endif
#ifdef OUTPUT_PARTICLES
  float P_out[6];
#endif
#ifndef OUTPUT_PARTICLES
  unsigned int pc;
#ifdef OUTPUT_UNFORMATTED
  unsigned int blockmaxlen;
#endif
#endif

  // How many variables are we writing out for each halo?
  halosize = 7;
#ifdef INERTIA
  halosize += 6;
#endif
#ifdef DISPERSION
  halosize += 6;
#endif

  nprocgroup = NTask / nwrite;
  if (NTask % nwrite) nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if (nhalos > 0) {
        sprintf(buf, "%s.%d", OutputFileBase, ThisTask);
        if(!(fp = fopen(buf, "w"))) {
          printf("\nError: Unable to open output file %s.\n\n", buf);
          FatalError((char *)"FOF.c", 592);
        }

        // We split the 'ifdefs' into outputting particles or just group properties. This maskes for some redundant
        // code and could be shorter, but this is confusing enough as it is without trying to compress it more.

        // Outputting each constituent particle of the halo (this is almost the same for formatted and unformatted)
        // =======================================================================================================
#ifdef OUTPUT_PARTICLES

        // Output the number of halos
#ifdef OUTPUT_UNFORMATTED
        dummy = sizeof(nhalos); 
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        my_fwrite(&nhalos, sizeof(int), 1, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
        fprintf(fp, "%15d\n", nhalos);
#endif
        // Output each halo separately with all the particles following each individual halo
        // For each halo calculate the group properties
        haloblock = (float *)malloc(halosize*sizeof(float));
        for(i=0; i<nhalos; i++) {
	      invnparthalo = 1.0/(double)nparthalo[i];
          xh=0; yh=0; zh=0; vxh=0; vyh=0; vzh=0;
#ifdef INERTIA
          xh2=0; yh2=0; zh2=0; xhyh=0; xhzh=0; yhzh=0;
#endif
#ifdef DISPERSION
          vxh2=0; vyh2=0; vzh2=0; vxhvyh=0; vxhvzh=0; vyhvzh=0;
#endif
          j = ihalo[i];
          do {
            xh  += (double)P[j-1].Pos[0];
            yh  += (double)P[j-1].Pos[1];
            zh  += (double)P[j-1].Pos[2];
            vxh += (double)P[j-1].Vel[0];
            vyh += (double)P[j-1].Vel[1];
            vzh += (double)P[j-1].Vel[2];
#ifdef INERTIA
            xh2  += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[0];
            yh2  += (double)P[j-1].Pos[1]*(double)P[j-1].Pos[1];
            zh2  += (double)P[j-1].Pos[2]*(double)P[j-1].Pos[2];
	        xhyh += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[1];
	        xhzh += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[2];
	        yhzh += (double)P[j-1].Pos[1]*(double)P[j-1].Pos[2];
#endif
#ifdef DISPERSION
            vxh2   += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[0];
            vyh2   += (double)P[j-1].Vel[1]*(double)P[j-1].Vel[1];
            vzh2   += (double)P[j-1].Vel[2]*(double)P[j-1].Vel[2];
            vxhvyh += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[1];
            vxhvzh += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[2];
            vyhvzh += (double)P[j-1].Vel[1]*(double)P[j-1].Vel[2];
#endif
            j=next[j];
            if (j == ihalo[i]) break;
          } while(1);

          // Put whatever we need to output into the output array
#ifdef OUTPUT_UNFORMATTED
          ind = 0;
          haloblock[ind] = (float)nparthalo[i];
#else
          ind = -1;
#endif
          haloblock[ind+1] = xh*invnparthalo;
	      haloblock[ind+2] = yh*invnparthalo;
	      haloblock[ind+3] = zh*invnparthalo;
  	      haloblock[ind+4] = vxh*invnparthalo;
	      haloblock[ind+5] = vyh*invnparthalo;
	      haloblock[ind+6] = vzh*invnparthalo;
#ifdef INERTIA
	      haloblock[ind+7] = yh2+zh2-invnparthalo*(yh*yh+zh*zh);
 	      haloblock[ind+8] = xh2+zh2-invnparthalo*(xh*xh+zh*zh);
	      haloblock[ind+9] = xh2+yh2-invnparthalo*(xh*xh+yh*yh);
	      haloblock[ind+10] = invnparthalo*xh*yh-xhyh;	
	      haloblock[ind+11] = invnparthalo*xh*zh-xhzh;
	      haloblock[ind+12] = invnparthalo*yh*zh-yhzh;
#ifdef DISPERSION
	      haloblock[ind+13] = (vxh2-vxh*vxh*invnparthalo)*invnparthalo;
	      haloblock[ind+14] = (vyh2-vyh*vyh*invnparthalo)*invnparthalo;
 	      haloblock[ind+15] = (vzh2-vzh*vzh*invnparthalo)*invnparthalo;
  	      haloblock[ind+16] = (vxhvyh-vxh*vyh*invnparthalo)*invnparthalo;
	      haloblock[ind+17] = (vxhvzh-vxh*vzh*invnparthalo)*invnparthalo;
	      haloblock[ind+18] = (vyhvzh-vyh*vzh*invnparthalo)*invnparthalo;
#endif
#else
#ifdef DISPERSION
	      haloblock[ind+7] = (vxh2-vxh*vxh*invnparthalo)*invnparthalo;
	      haloblock[ind+8] = (vyh2-vyh*vyh*invnparthalo)*invnparthalo;
 	      haloblock[ind+9] = (vzh2-vzh*vzh*invnparthalo)*invnparthalo;
  	      haloblock[ind+10] = (vxhvyh-vxh*vyh*invnparthalo)*invnparthalo;
	      haloblock[ind+11] = (vxhvzh-vxh*vzh*invnparthalo)*invnparthalo;
	      haloblock[ind+12] = (vyhvzh-vyh*vzh*invnparthalo)*invnparthalo;
#endif
#endif
          // Output the halo properties
#ifdef OUTPUT_UNFORMATTED
          dummy = sizeof(float)*halosize;
          my_fwrite(&dummy, sizeof(dummy), 1, fp);
          my_fwrite(&(haloblock[0]), sizeof(float), halosize, fp);
          my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
#ifdef INERTIA
#ifdef DISPERSION
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[6], haloblock[7], haloblock[8], haloblock[9], 
                                haloblock[10], haloblock[11], haloblock[12], haloblock[13], haloblock[14], haloblock[15], haloblock[16], haloblock[17]);
#else  
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[6], haloblock[7], haloblock[8], haloblock[9], 
                                haloblock[10], haloblock[11]);
#endif
#else
#ifdef DISPERSION
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[12], haloblock[13], haloblock[14], haloblock[15], 
                                haloblock[16], haloblock[17]);
#else
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5]);
#endif
#endif
#endif
          // Output the constituent particles
          j = ihalo[i];
          P_out[0] = (float)P[j-1].Pos[0];
          P_out[1] = (float)P[j-1].Pos[1];
          P_out[2] = (float)P[j-1].Pos[2];
          P_out[3] = (float)P[j-1].Vel[0];
          P_out[4] = (float)P[j-1].Vel[1];
          P_out[5] = (float)P[j-1].Vel[2];
#ifdef OUTPUT_UNFORMATTED
#ifdef PARTICLE_ID
          dummy = sizeof(P_out)+sizeof(unsigned long long);
#else
          dummy = sizeof(P_out);
#endif
          my_fwrite(&dummy, sizeof(dummy), 1, fp);
#ifdef PARTICLE_ID
          my_fwrite(&(P[j-1].ID), sizeof(unsigned long long), 1, fp);
#endif
          my_fwrite(&(P_out[0]), sizeof(P_out), 1, fp);
          my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
#ifdef PARTICLE_ID
          fprintf(fp,"%15llu %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P[j-1].ID, P_out[0], P_out[1], P_out[2], P_out[3], P_out[4], P_out[5]);
#else
          fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P_out[0], P_out[1], P_out[2], P_out[3], P_out[4], P_out[5]);
#endif
#endif
          do {
            j=next[j];
            if (j == ihalo[i]) break; 
            P_out[0] = (float)P[j-1].Pos[0];
            P_out[1] = (float)P[j-1].Pos[1];
            P_out[2] = (float)P[j-1].Pos[2];
            P_out[3] = (float)P[j-1].Vel[0];
            P_out[4] = (float)P[j-1].Vel[1];
            P_out[5] = (float)P[j-1].Vel[2];
#ifdef OUTPUT_UNFORMATTED
#ifdef PARTICLE_ID
            dummy = sizeof(P_out)+sizeof(unsigned long long);
#else
            dummy = sizeof(P_out);
#endif
            my_fwrite(&dummy, sizeof(dummy), 1, fp);
#ifdef PARTICLE_ID
            my_fwrite(&(P[j-1].ID), sizeof(unsigned long long), 1, fp);
#endif
            my_fwrite(&(P_out[0]), sizeof(P_out), 1, fp);
            my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
#ifdef PARTICLE_ID
            fprintf(fp,"%15llu %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P[j-1].ID, P_out[0], P_out[1], P_out[2], P_out[3], P_out[4], P_out[5]);
#else
            fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P_out[0], P_out[1], P_out[2], P_out[3], P_out[4], P_out[5]);
#endif
#endif
          } while(1);
        }
#else

        // Otherwise we only care about the halo group properties (formatted outputs line by line whereas unformatted outputs in chunks)
        // =============================================================================================================================

#ifdef OUTPUT_UNFORMATTED
        // Output the number of halos
        dummy = sizeof(nhalos); 
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
        my_fwrite(&nhalos, sizeof(int), 1, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);

        // Allocate memory to store the chunk of group properties (we have 4.0*nparticles_tot bytes free from deallocating head)
        haloblock = (float *)malloc(nparticles_tot*sizeof(float));
        blockmaxlen = (unsigned int)(nparticles_tot/halosize);
        dummy = nhalos*halosize*sizeof(float);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
        // Allocate enough for a single halo at a time
        haloblock = (float *)malloc((halosize-1)*sizeof(float));
#endif
        // Calculate the group properties for each halo and output when block is full
        // Loop over all the halos
        for(i=0, pc=0; i<nhalos; i++) {
	      invnparthalo = 1.0/(double)nparthalo[i];
          xh=0; yh=0; zh=0; vxh=0; vyh=0; vzh=0;
#ifdef INERTIA
          xh2=0; yh2=0; zh2=0; xhyh=0; xhzh=0; yhzh=0;
#endif
#ifdef DISPERSION
          vxh2=0; vyh2=0; vzh2=0; vxhvyh=0; vxhvzh=0; vyhvzh=0;
#endif
          j = ihalo[i];
          do {
            xh  += (double)P[j-1].Pos[0];
            yh  += (double)P[j-1].Pos[1];
            zh  += (double)P[j-1].Pos[2];
            vxh += (double)P[j-1].Vel[0];
            vyh += (double)P[j-1].Vel[1];
            vzh += (double)P[j-1].Vel[2];
#ifdef INERTIA
            xh2  += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[0];
            yh2  += (double)P[j-1].Pos[1]*(double)P[j-1].Pos[1];
            zh2  += (double)P[j-1].Pos[2]*(double)P[j-1].Pos[2];
	        xhyh += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[1];
	        xhzh += (double)P[j-1].Pos[0]*(double)P[j-1].Pos[2];
	        yhzh += (double)P[j-1].Pos[1]*(double)P[j-1].Pos[2];
#endif
#ifdef DISPERSION
            vxh2   += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[0];
            vyh2   += (double)P[j-1].Vel[1]*(double)P[j-1].Vel[1];
            vzh2   += (double)P[j-1].Vel[2]*(double)P[j-1].Vel[2];
            vxhvyh += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[1];
            vxhvzh += (double)P[j-1].Vel[0]*(double)P[j-1].Vel[2];
            vyhvzh += (double)P[j-1].Vel[1]*(double)P[j-1].Vel[2];
#endif
            j=next[j];
            if (j == ihalo[i]) break;
          } while(1);

          // Put whatever we need to output into the output array
          // Convert the number of particles per halo to a float for ease of writing out.
#ifdef OUTPUT_UNFORMATTED
          ind = pc*halosize;
          haloblock[ind] = (float)nparthalo[i];
#else
          ind = -1;
#endif
          haloblock[ind+1] = xh*invnparthalo;
	      haloblock[ind+2] = yh*invnparthalo;
	      haloblock[ind+3] = zh*invnparthalo;
  	      haloblock[ind+4] = vxh*invnparthalo;
	      haloblock[ind+5] = vyh*invnparthalo;
	      haloblock[ind+6] = vzh*invnparthalo;
#ifdef INERTIA
	      haloblock[ind+7] = yh2+zh2-invnparthalo*(yh*yh+zh*zh);
 	      haloblock[ind+8] = xh2+zh2-invnparthalo*(xh*xh+zh*zh);
	      haloblock[ind+9] = xh2+yh2-invnparthalo*(xh*xh+yh*yh);
	      haloblock[ind+10] = invnparthalo*xh*yh-xhyh;	
	      haloblock[ind+11] = invnparthalo*xh*zh-xhzh;
	      haloblock[ind+12] = invnparthalo*yh*zh-yhzh;
#ifdef DISPERSION
	      haloblock[ind+13] = (vxh2-vxh*vxh*invnparthalo)*invnparthalo;
	      haloblock[ind+14] = (vyh2-vyh*vyh*invnparthalo)*invnparthalo;
 	      haloblock[ind+15] = (vzh2-vzh*vzh*invnparthalo)*invnparthalo;
  	      haloblock[ind+16] = (vxhvyh-vxh*vyh*invnparthalo)*invnparthalo;
	      haloblock[ind+17] = (vxhvzh-vxh*vzh*invnparthalo)*invnparthalo;
	      haloblock[ind+18] = (vyhvzh-vyh*vzh*invnparthalo)*invnparthalo;
#endif
#else
#ifdef DISPERSION
	      haloblock[ind+7] = (vxh2-vxh*vxh*invnparthalo)*invnparthalo;
	      haloblock[ind+8] = (vyh2-vyh*vyh*invnparthalo)*invnparthalo;
 	      haloblock[ind+9] = (vzh2-vzh*vzh*invnparthalo)*invnparthalo;
  	      haloblock[ind+10] = (vxhvyh-vxh*vyh*invnparthalo)*invnparthalo;
	      haloblock[ind+11] = (vxhvzh-vxh*vzh*invnparthalo)*invnparthalo;
	      haloblock[ind+12] = (vyhvzh-vyh*vzh*invnparthalo)*invnparthalo;
#endif
#endif
#ifdef OUTPUT_UNFORMATTED
          pc++;
          if (pc == blockmaxlen) {
            my_fwrite(&(haloblock[0]), sizeof(float), pc*halosize, fp);
            pc = 0;
          }        
        }
        if (pc > 0) my_fwrite(&(haloblock[0]), sizeof(float), pc*halosize, fp);
        my_fwrite(&dummy, sizeof(dummy), 1, fp);
#else
#ifdef INERTIA
#ifdef DISPERSION
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[6], haloblock[7], haloblock[8], haloblock[9], 
                                haloblock[10], haloblock[11], haloblock[12], haloblock[13], haloblock[14], haloblock[15], haloblock[16], haloblock[17]);
#else  
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[6], haloblock[7], haloblock[8], haloblock[9], 
                                haloblock[10], haloblock[11]);
#endif
#else
#ifdef DISPERSION
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5], haloblock[12], haloblock[13], haloblock[14], haloblock[15], 
                                haloblock[16], haloblock[17]);
#else
          fprintf(fp,"%15d %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", 
                  nparthalo[i], haloblock[0], haloblock[1], haloblock[2], haloblock[3], haloblock[4], haloblock[5]);
#endif
#endif
        }
#endif
#endif
        fclose(fp);
        free(haloblock);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }   
  return;
}

// Generate the info file which contains a list of all the output files, the 8 corners of the region in those files and the number of halos and subsampled particles in the region
// ===============================================================================================================================================================================
void Output_Info(unsigned int ninhalo, unsigned int nouthalo) {

  FILE * fp; 
  char buf[300];
  int i;

  int * nhalos_table = (int *)malloc(sizeof(int) * NTask);
  unsigned int * ninhalo_table = (unsigned int *)malloc(sizeof(unsigned int) * NTask);
  unsigned int * nouthalo_table = (unsigned int *)malloc(sizeof(unsigned int) * NTask);
  MPI_Allgather(&nhalos, 1, MPI_INT, nhalos_table, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&ninhalo, 1, MPI_UNSIGNED, ninhalo_table, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
  MPI_Allgather(&nouthalo, 1, MPI_UNSIGNED, nouthalo_table, 1, MPI_UNSIGNED, MPI_COMM_WORLD);

  if (ThisTask == 0) {
    sprintf(buf, "%s.info", OutputFileBase);
    if(!(fp = fopen(buf, "w"))) {
      printf("\nERROR: Can't write in file '%s'.\n\n", buf);
      FatalError((char *)"FOF.c", 934);
    }
    fprintf(fp, "#    FILENUM      XMIN         YMIN        ZMIN         XMAX         YMAX         ZMAX         NHALOS     NSUBSAMPIN     NSUBSAMPOUT    \n");
    for (i=0; i<NTask; i++) {
      int Task_nz = i/(Nx*Ny);
      int Task_ny = (i-Nx*Ny*Task_nz)/Nx;
      int Task_nx = i-Nx*Task_ny-Nx*Ny*Task_nz; 
      double x0 = Task_nx*((Lxmax-Lxmin)/(double)Nx)+Lxmin;
      double y0 = Task_ny*((Lymax-Lymin)/(double)Ny)+Lymin;
      double z0 = Task_nz*((Lzmax-Lzmin)/(double)Nz)+Lzmin;
      double x1 = x0+((Lxmax-Lxmin)/(double)Nx);
      double y1 = y0+((Lymax-Lymin)/(double)Ny);
      double z1 = z0+((Lzmax-Lzmin)/(double)Nz);
      fprintf(fp, "%12d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12d %12u %12u\n", i, x0, y0, z0, x1, y1, z1, nhalos_table[i], ninhalo_table[i], nouthalo_table[i]);
    }
    fclose(fp);
  }

  free(nhalos_table);
  free(ninhalo_table);
  free(nouthalo_table);

  return;
}

