import math
import numpy as np
import scipy as sp
import os

# Compilation options
MEMORY_MODE = 1
PARTICLE_ID = 1

# The parameters
Nx = 8
Ny = 8
Nz = 8

Px = 2048
Py = 2048
Pz = 2048

Lxmin = 0.0
Lymin = 0.0
Lzmin = 0.0
Lxmax = 1200.0
Lymax = 1200.0
Lzmax = 1200.0

Buffer = 1.25
Boundarysize = 35.0

ninputfiles = 512

# Perform the calculation. We assume that all the particles are spread uniformly over the input files
extra=2*(Nx*(Boundarysize/(Lxmax-Lxmin))+Ny*(Boundarysize/(Lymax-Lymin))+Nz*(Boundarysize/(Lzmax-Lzmin)));
extra=extra+4*Nx*Ny*(Boundarysize/(Lxmax-Lxmin))*(Boundarysize/(Lymax-Lymin));
extra=extra+4*Nx*Nz*(Boundarysize/(Lxmax-Lxmin))*(Boundarysize/(Lzmax-Lzmin));
extra=extra+4*Ny*Nz*(Boundarysize/(Lymax-Lymin))*(Boundarysize/(Lzmax-Lzmin));
extra=extra+8*Nx*Ny*Nz*(Boundarysize/(Lxmax-Lxmin))*(Boundarysize/(Lymax-Lymin))*(Boundarysize/(Lzmax-Lzmin));

maxparticles=math.ceil(Buffer*(1.0+extra)*(float(Px)/float(Nx))*(float(Py)/float(Ny))*(float(Pz)/float(Nz)));

# The amount of memory per particle
if (MEMORY_MODE):
    if (PARTICLE_ID):
        part = 32.0
    else:
        part = 24.0 
else:
    if (PARTICLE_ID):
        part = 56.0
    else:
        part = 48.0

fac_1 = part*maxparticles                # Particle structure
fac_2 = (part*Px*Py*Pz)/ninputfiles      # Opening and reading input data
fac_3 = fac_1/ninputfiles                # Temporary memory for sorted/transferred data
fac_4 = 16.0*maxparticles                # FOF algorithm

MEM_FOF  = fac_1+fac_4
MEM_READ = fac_1+fac_2+fac_3
MEM_TRANSFER = fac_1+2.0*fac_3

print MEM_FOF, MEM_READ, MEM_TRANSFER
      
MEM_TOT = (max([MEM_FOF,MEM_READ,MEM_TRANSFER]))/(1024.0*1024.0*1024.0)
print MEM_TOT, 'GB'


