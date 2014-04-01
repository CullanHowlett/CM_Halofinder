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
int * nparthalo, *ihalo;
unsigned int * next;
double linksq;

void FOF(void) {

  int nx, ny, nz;
  int nphalo, nhalos;
  int ix, iy, iz, ix1, iy1, ix2, iy2, iz2;
  unsigned int i, ip, ip2, ind;
  unsigned int * head, *chain, * first;
  double lcell[3];
#ifdef VARLINK
  double distance;
#endif 

  if (ThisTask == 0) printf("Creating subcells...\n");

  // set grid of cells 
  lcell[0]=(Lxmax-Lxmin)/Px;
  lcell[1]=(Lymax-Lymin)/Py;
  lcell[2]=(Lzmax-Lzmin)/Pz;

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
        nphalo=nphalo+1;
        if (ip == i) break;
        head[ip]=1;
      } while(1);

      if (nphalo < nphalomin) {
        head[ip]=1;   
      } else {
        nhalos=nhalos+1;
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
#ifdef PERIODIC
  if ((xh < rmin[0]) || (xh >= rmax[0])) return;
  if ((yh < rmin[1]) || (yh >= rmax[1])) return;
  if ((zh < rmin[2]) || (zh >= rmax[2])) return;

  // This checks that for any halos with a centre of mass on the processor none of the 
  // constituent particles extend beyond the boundary, otherwise the halos aren't
  // complete and the boundary regions must be made larger.
  j=i;
  do {
    if ((P[j-1].Pos[0]-rmin_buff[0] < linklength) || (rmax_buff[0]-P[j-1].Pos[0] < linklength)) {
      printf("Particle very near boundary in x-direction. Repeat run with increased boundary size\n");
    }
    if ((P[j-1].Pos[1]-rmin_buff[1] < linklength) || (rmax_buff[1]-P[j-1].Pos[1] < linklength)) {
      printf("Particle very near boundary in y-direction. Repeat run with increased boundary size\n");
    }
    if ((P[j-1].Pos[2]-rmin_buff[2] < linklength) || (rmax_buff[2]-P[j-1].Pos[2] < linklength)) {
      printf("Particle very near boundary in z-direction. Repeat run with increased boundary size\n");
    }
    j=next[j];
    if (j == i) break;
  } while(1);

#else
  if (Local_nx == Nx-1) {
    if ((xh < rmin[0]) || (xh  > rmax[0])) return;
  } else {
    if ((xh < rmin[0]) || (xh >= rmax[0])) return;
  }
  if (Local_ny == Ny-1) {
    if ((yh < rmin[1]) || (yh  > rmax[1])) return;
  } else {
    if ((yh < rmin[1]) || (yh >= rmax[1])) return;
  }
  if (Local_nz == Nz-1) {
    if ((zh < rmin[2]) || (zh  > rmax[2])) return;
  } else {
    if ((zh < rmin[2]) || (zh >= rmax[2])) return;
  }
#endif

  ihalo[nhalos]     = i;
  nparthalo[nhalos] = nphalo;

  nhalos++;

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
      sprintf(buf, "%s/%s.%d", OutputDir, OutputFileBase, ThisTask);
      if(!(fp = fopen(buf, "w"))) {
        printf("Error. Can't write in file '%s'\n", buf);
        MPI_Abort(MPI_COMM_WORLD, 10);
      }
#ifdef OUTPUT_PARTICLES
      fprintf(fp, "%12d\n", nhalos);
      for(i=0; i<nhalos; i++) {
        j = ihalo[i];
        double xh=0, yh=0, zh=0;
        double vxh=0, vyh=0, vzh=0;
        do {
          xh += P[j-1].Pos[0];
          yh += P[j-1].Pos[1];
          zh += P[j-1].Pos[2];
          vxh += P[j-1].Vel[0];
          vyh += P[j-1].Vel[1];
          vzh += P[j-1].Vel[2];
          j=next[j];
          if (j == ihalo[i]) break;
        } while(1);
        fprintf(fp,"%12d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", nparthalo[i], xh/nparthalo[i], yh/nparthalo[i], zh/nparthalo[i], vxh/nparthalo[i], vyh/nparthalo[i], vzh/nparthalo[i]);   
        j = ihalo[i];
#ifdef PARTICLE_ID
        fprintf(fp,"%12llu %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", P[j-1].ID, (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
#else
        fprintf(fp,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
#endif
        do {
          j=next[j];
          if (j == ihalo[i]) break;
#ifdef PARTICLE_ID
          fprintf(fp,"%12llu %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", P[j-1].ID, (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
#else
          fprintf(fp,"%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", (float)(P[j-1].Pos[0]),(float)(P[j-1].Pos[1]),(float)(P[j-1].Pos[2]),(float)(P[j-1].Vel[0]),(float)(P[j-1].Vel[1]),(float)(P[j-1].Vel[2]));
#endif
        } while(1);
      }
#else
      for(i=0; i<nhalos; i++) {
        j = ihalo[i];
        double xh=0, yh=0, zh=0;
        double vxh=0, vyh=0, vzh=0;
        do {
          xh += P[j-1].Pos[0];
          yh += P[j-1].Pos[1];
          zh += P[j-1].Pos[2];
          vxh += P[j-1].Vel[0];
          vyh += P[j-1].Vel[1];
          vzh += P[j-1].Vel[2];
          j=next[j];
          if (j == ihalo[i]) break;
        } while(1);
        fprintf(fp,"%12d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n", nparthalo[i], xh/nparthalo[i], yh/nparthalo[i], zh/nparthalo[i], vxh/nparthalo[i], vyh/nparthalo[i], vzh/nparthalo[i]);   
      }
#endif
      fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }   
  return;
}
