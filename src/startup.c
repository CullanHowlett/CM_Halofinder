#include "vars.h"

void Set_Params(void) {

  double extra;

#ifdef VARLINK
  int i, origin_processor, origin_processor_comp[3];
  double dmin_buff[3], dmin_nobuff[3], dmax_buff[3], dmax_nobuff[3], dmin_comp[3], dmax_comp[3];
#endif

  // Check that the geometry matches the number of processors
  if (NTask != Nx*Ny*Nz) {
    if (ThisTask == 0) printf("\nERROR: The number of cells does not match the number of processors... Whoops!!\n\n");
    FatalError("startup.c", 18);
  }

  // Check that the number of processors reading and writing is not greater than the number
  // of processors otherwise we might not read in or write out all the data
  if (nread > NTask) {
    if (ThisTask == 0) printf("\nWARNING: Tried to read in with more processors than are being used...setting nread = NTask...\n\n");
    nread=NTask;
  }
  if (nwrite > NTask) {
    if (ThisTask == 0) printf("\nWARNING: Tried to write halo data with more processors than are being used...setting nwrite = NTask...\n\n");
    nwrite=NTask;
  }

  // Here we define the x, y and z coordinates of each processor as well as the maximum and minimum
  // extents of the processor, with and without including the boundary regions. We don't care if
  // the boundaries are periodic or not as the particle positions themselves are wrapped to fit in these
  // limits
  Local_nz = ThisTask/(Nx*Ny);
  Local_ny = (ThisTask-Nx*Ny*Local_nz)/Nx;
  Local_nx = ThisTask-Nx*Local_ny-Nx*Ny*Local_nz;
  rmin[0]=Local_nx*((Lxmax-Lxmin)/(double)Nx)+Lxmin;
  rmin[1]=Local_ny*((Lymax-Lymin)/(double)Ny)+Lymin;
  rmin[2]=Local_nz*((Lzmax-Lzmin)/(double)Nz)+Lzmin;
  rmax[0]=rmin[0]+((Lxmax-Lxmin)/(double)Nx);
  rmin_buff[0]=rmin[0]-boundarysize;
  rmax_buff[0]=rmax[0]+boundarysize;
  rmax[1]=rmin[1]+((Lymax-Lymin)/(double)Ny);
  rmin_buff[1]=rmin[1]-boundarysize;
  rmax_buff[1]=rmax[1]+boundarysize;
  rmax[2]=rmin[2]+((Lzmax-Lzmin)/(double)Nz);
  rmin_buff[2]=rmin[2]-boundarysize;
  rmax_buff[2]=rmax[2]+boundarysize;
//#ifndef PERIODIC
//  if (Local_nx == 0) rmin_buff[0]=rmin[0];
//  if (Local_nx == Nx-1) rmax_buff[0]=rmax[0];
//  if (Local_ny == 0) rmin_buff[1]=rmin[1];
//  if (Local_ny == Ny-1) rmax_buff[1]=rmax[1];
//  if (Local_nz == 0) rmin_buff[2]=rmin[2];
//  if (Local_nz == Nz-1) rmax_buff[2]=rmax[2];
//#endif

#ifdef VARLINK
  // This finds the maximum and minimum distances from the origin for each processor. To make it general
  // , i.e., so that it holds for both periodic and non-periodic boundaries, and for any origin, is suprising
  // non-trivial and requires us to calculate each of the three directions individually before finding the
  // maximum and minimum of each and combining them. We also have to identify on which processor, if any,
  // the origin is located and set the minimum distance to the origin for that processor to 0 explicitly. 
  dmin = 0.0;
  dmax = 0.0;
  origin_processor_comp[0]=int(floor(Nx*(Origin_x/(Lxmax-Lxmin))));
  origin_processor_comp[1]=int(floor(Ny*(Origin_y/(Lymax-Lymin))));
  origin_processor_comp[2]=int(floor(Nz*(Origin_z/(Lzmax-Lzmin))));
  if (ThisTask == 0) {
    if ((origin_processor_comp[0] >= Nx) || (origin_processor_comp[0] < 0) ||
        (origin_processor_comp[1] >= Ny) || (origin_processor_comp[1] < 0) ||
        (origin_processor_comp[2] >= Nz) || (origin_processor_comp[2] < 0)) 
        printf("\nWARNING: origin out of simulation bounds (is this intentional?)\n\n");
  }
  dmin_buff[0] = (rmin_buff[0]-Origin_x)*(rmin_buff[0]-Origin_x);
  dmin_buff[1] = (rmin_buff[1]-Origin_y)*(rmin_buff[1]-Origin_y);
  dmin_buff[2] = (rmin_buff[2]-Origin_z)*(rmin_buff[2]-Origin_z);
  dmin_nobuff[0] = (rmin[0]-Origin_x)*(rmin[0]-Origin_x);
  dmin_nobuff[1] = (rmin[1]-Origin_y)*(rmin[1]-Origin_y);
  dmin_nobuff[2] = (rmin[2]-Origin_z)*(rmin[2]-Origin_z);
  dmax_buff[0] = (rmax_buff[0]-Origin_x)*(rmax_buff[0]-Origin_x);
  dmax_buff[1] = (rmax_buff[1]-Origin_y)*(rmax_buff[1]-Origin_y);
  dmax_buff[2] = (rmax_buff[2]-Origin_z)*(rmax_buff[2]-Origin_z);
  dmax_nobuff[0] = (rmax[0]-Origin_x)*(rmax[0]-Origin_x);
  dmax_nobuff[1] = (rmax[1]-Origin_y)*(rmax[1]-Origin_y);
  dmax_nobuff[2] = (rmax[2]-Origin_z)*(rmax[2]-Origin_z);
  for (i=0; i<3; i++) {
    dmin_comp[i] = dmin_buff[i];
    if (dmax_buff[i]   < dmin_comp[i]) dmin_comp[i] = dmax_buff[i];
    if (dmax_nobuff[i] < dmin_comp[i]) dmin_comp[i] = dmax_nobuff[i];
    if (dmin_nobuff[i] < dmin_comp[i]) dmin_comp[i] = dmin_nobuff[i];
    dmax_comp[i] = dmin_buff[i];
    if (dmax_buff[i]   > dmax_comp[i]) dmax_comp[i] = dmax_buff[i];
    if (dmax_nobuff[i] > dmax_comp[i]) dmax_comp[i] = dmax_nobuff[i];
    if (dmin_nobuff[i] > dmax_comp[i]) dmax_comp[i] = dmin_nobuff[i];
    dmin += dmin_comp[i];
    dmax += dmax_comp[i];
  }
  origin_processor=Ny*Nx*origin_processor_comp[2]+Nx*origin_processor_comp[1]+origin_processor_comp[0];
  if (ThisTask == origin_processor) dmin = 0.0;
#endif

  // Creates space for the particle data (including extra space for the boundary particles based on the
  // percentage of particles found in the boundaries)
  extra=2*(Nx*(boundarysize/(Lxmax-Lxmin))+Ny*(boundarysize/(Lymax-Lymin))+Nz*(boundarysize/(Lzmax-Lzmin)));
  extra=extra+4*Nx*Ny*(boundarysize/(Lxmax-Lxmin))*(boundarysize/(Lymax-Lymin));
  extra=extra+4*Nx*Nz*(boundarysize/(Lxmax-Lxmin))*(boundarysize/(Lzmax-Lzmin));
  extra=extra+4*Ny*Nz*(boundarysize/(Lymax-Lymin))*(boundarysize/(Lzmax-Lzmin));
  extra=extra+8*Nx*Ny*Nz*(boundarysize/(Lxmax-Lxmin))*(boundarysize/(Lymax-Lymin))*(boundarysize/(Lzmax-Lzmin));

  maxparticles=(unsigned int)ceil(Buffer*(1.0+extra)*((unsigned int)Px/Nx)*((unsigned int)Py/Ny)*((unsigned int)Pz/Nz));

  P = (struct part_data *)malloc(maxparticles*sizeof(struct part_data));

  return;
}

void Read_Parameterfile(char * fname) {

#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  char buf[500],buf1[500],buf2[500],buf3[500];
  int i,j,nt;
  int id[MAXTAGS];
  int errorFlag = 0;

  // read parameter file on all processes for simplicity

  nt = 0;

  strcpy(tag[nt], "Nx");
  addr[nt] = &Nx;
  id[nt++] = INT;

  strcpy(tag[nt], "Ny");
  addr[nt] = &Ny;
  id[nt++] = INT;

  strcpy(tag[nt], "Nz");
  addr[nt] = &Nz;
  id[nt++] = INT;

  strcpy(tag[nt], "Px");
  addr[nt] = &Px;
  id[nt++] = INT;

  strcpy(tag[nt], "Py");
  addr[nt] = &Py;
  id[nt++] = INT;

  strcpy(tag[nt], "Pz");
  addr[nt] = &Pz;
  id[nt++] = INT;
  
  strcpy(tag[nt], "Lxmin");
  addr[nt] = &Lxmin;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Lymin");
  addr[nt] = &Lymin;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Lzmin");
  addr[nt] = &Lzmin;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Lxmax");
  addr[nt] = &Lxmax;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Lymax");
  addr[nt] = &Lymax;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Lzmax");
  addr[nt] = &Lzmax;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Seed");
  addr[nt] = &Seed;
  id[nt++] = INT;
  
  strcpy(tag[nt], "Buffer");
  addr[nt] = &Buffer;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Boundarysize");
  addr[nt] = &boundarysize;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Linklength");
  addr[nt] = &linklength;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "Nphalomin");
  addr[nt] = &nphalomin;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "SampInHalos");
  addr[nt] = &SampInHalos;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "SampOutHalos");
  addr[nt] = &SampOutHalos;
  id[nt++] = FLOAT;
  
  strcpy(tag[nt], "InputStyle");
  addr[nt] = &InputStyle;
  id[nt++] = INT; 
  
  strcpy(tag[nt], "Ninputfiles");
  addr[nt] = &ninputfiles;
  id[nt++] = INT;
  
  strcpy(tag[nt], "Starting_file");
  addr[nt] = &starting_file;
  id[nt++] = INT;

  strcpy(tag[nt], "Nread");
  addr[nt] = &nread;
  id[nt++] = INT;

  strcpy(tag[nt], "Nwrite");
  addr[nt] = &nwrite;
  id[nt++] = INT;
  
  strcpy(tag[nt], "InputFileBase");
  addr[nt] = InputFileBase;
  id[nt++] = STRING;
  
  strcpy(tag[nt], "OutputFileBase");
  addr[nt] = OutputFileBase;
  id[nt++] = STRING;

#ifdef VARLINK
  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Origin_x");
  addr[nt] = &Origin_x;
  id[nt++] = INT;

  strcpy(tag[nt], "Origin_y");
  addr[nt] = &Origin_y;
  id[nt++] = INT;

  strcpy(tag[nt], "Origin_z");
  addr[nt] = &Origin_z;
  id[nt++] = INT;
#endif

  if((fd = fopen(fname, "r"))) {
    while(!feof(fd)) {
      buf[0] = 0;
      fgets(buf, 500, fd);

      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2) continue;
      if(buf1[0] == '%') continue;

      for(i = 0, j = -1; i < nt; i++) {
        if(strcmp(buf1, tag[i]) == 0)  {
          j = i;
	        tag[i][0] = 0;
	        break;
        }
      }
      
      if(j >= 0) {
	      switch (id[j]) {
	        case FLOAT:
	          *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
	          strcpy((char *)addr[j], buf2);
	          break;
	        case INT:
	          *((int *) addr[j]) = atoi(buf2);
	          break;
	      }
      } else {
        if(ThisTask == 0) fprintf(stdout,"\nError: in run parameters file %s:  Tag '%s' not allowed or multiple defined.\n",fname,buf1);
        errorFlag = 1;
      }
    }
    fclose(fd);
    for(i = 0; i < nt; i++) {
      if(*tag[i]) {
        if(ThisTask == 0) fprintf(stdout, "\nError: I miss a value for tag '%s' in parameter file '%s'.\n\n", tag[i], fname);
        errorFlag = 1;
      }
    }
  } else {
    if(ThisTask == 0) fprintf(stdout,"\ERROR: Parameter file %s not found.\n\n",fname);
    errorFlag = 1;
  }

  if(errorFlag) {
    MPI_Finalize();
    exit(0);
  }

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS

  return;
}
