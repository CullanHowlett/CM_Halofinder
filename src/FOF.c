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
                FatalError("FOF.c", 231); 
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
                FatalError("FOF.c", 268); 
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

  // Write out the halo data
  if (ThisTask == 0) {
    printf("Outputting the halos\n");
    printf("====================\n\n");
  }
  Output_Halos();
  
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

  next++;
  free(P);
  free(next);
  free(ihalo);
  free(nparthalo);
  if (SampInHalos>0)  free(inhalo);
  if (SampOutHalos>0) free(outhalo);

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
  
  nprocgroup = NTask / nwrite;
  
  if (NTask % nwrite) nprocgroup++;
  
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  
  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      if(in) {
        sprintf(buf, "%s_inhalo.%d", OutputFileBase, ThisTask);
      } else {
        sprintf(buf, "%s_outhalo.%d", OutputFileBase, ThisTask);
      }
      if(!(fp = fopen(buf, "w"))) {
        printf("\nError: Unable to open output file %s.\n\n", buf);
        FatalError("FOF.c", 342);
      }
      for(i=0; i<nsubsamp; i++) {
        j=subsamp[i];
#ifdef PARTICLE_ID
        fprintf(fp,"%15llu %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", P[j-1].ID, (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
#else
        fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));   
#endif
      }
      fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }   
  return;
}

void Output_Halos(void) {

  FILE *fp; 
  char buf[300];
  int i, j, nprocgroup, groupTask, masterTask;

  nprocgroup = NTask / nwrite;

  if (NTask % nwrite) nprocgroup++;

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++) {
    if (ThisTask == (masterTask + groupTask)) {
      sprintf(buf, "%s.%d", OutputFileBase, ThisTask);
      if(!(fp = fopen(buf, "w"))) {
        printf("\nError: Unable to open output file %s.\n\n", buf);
        FatalError("FOF.c", 342);
      }
#ifdef OUTPUT_PARTICLES
      fprintf(fp, "%15d\n", nhalos);
#endif
      for(i=0; i<nhalos; i++) {
	double invnparthalo = 1.0/(double)nparthalo[i];
        double xh=0, yh=0, zh=0;
        double vxh=0, vyh=0, vzh=0;
#ifdef INERTIA
        double xh2=0, yh2=0, zh2=0;
	double xhyh=0, xhzh=0, yhzh=0;
#endif
#ifdef DISPERSION
        double vxh2=0, vyh2=0, vzh2=0;
	double vxhvyh=0, vxhvzh=0, vyhvzh=0;
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
	xh *= invnparthalo;
	yh *= invnparthalo;
	zh *= invnparthalo;
	vxh *= invnparthalo;
	vyh *= invnparthalo;
	vzh *= invnparthalo;
        fprintf(fp,"%12d %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf %12.6lf ", nparthalo[i], xh, yh, zh, vxh, vyh, vzh);
#ifdef INERTIA
	double Ixx = yh2+zh2-nparthalo[i]*(yh*yh+zh*zh);
 	double Iyy = xh2+zh2-nparthalo[i]*(xh*xh+zh*zh);
	double Izz = xh2+yh2-nparthalo[i]*(xh*xh+yh*yh);
	double Ixy = nparthalo[i]*xh*yh-xhyh;	
	double Ixz = nparthalo[i]*xh*zh-xhzh;
	double Iyz = nparthalo[i]*yh*zh-yhzh;
        fprintf(fp,"%15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf ", Ixx, Iyy, Izz, Ixy, Ixz, Iyz);
#endif
#ifdef DISPERSION
	double Sigmaxx = vxh2*invnparthalo-vxh*vxh;
	double Sigmayy = vyh2*invnparthalo-vyh*vyh;
	double Sigmazz = vzh2*invnparthalo-vzh*vzh;
	double Sigmaxy = vxhvyh*invnparthalo-vxh*vyh;
	double Sigmaxz = vxhvzh*invnparthalo-vxh*vzh;
	double Sigmayz = vyhvzh*invnparthalo-vyh*vzh;
        fprintf(fp,"%15.6lf %15.6lf %15.6lf %15.6lf %15.6lf %15.6lf ", Sigmaxx, Sigmayy, Sigmazz, Sigmaxy, Sigmaxz, Sigmayz);
#endif
	fprintf(fp, "\n");
#ifdef OUTPUT_PARTICLES
        j = ihalo[i];
#ifdef PARTICLE_ID
        fprintf(fp,"%15llu ", P[j-1].ID);
#endif
        fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
        do {
          j=next[j];
          if (j == ihalo[i]) break;
#ifdef PARTICLE_ID
        fprintf(fp,"%15llu ", P[j-1].ID);
#endif
        fprintf(fp,"%15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n", (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
        } while(1);
#endif
      }
      fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }   
  return;
}
