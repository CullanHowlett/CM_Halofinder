import math
import numpy as np
import scipy as sp
import os

# Set the limits for the files
Lxmin = -3072.0
Lymin = 0.0
Lzmin = 0.0
Lxmax = -1501.0
Lymax = 3072.0
Lzmax = 3072.0
Infofile = "/users/cullanh/mock_catalogues/codes/CM_Halofinder/files/PICOLA_MICE_L3072_N4096_M4096_R101_lightcone.info"
Infile = "/scratch/mm855/L3072_N4096_R101_MICE/PICOLA_MICE_L3072_N4096_M4096_R101_lightcone"
Outfile = "/users/cullanh/mock_catalogues/codes/CM_Halofinder/files/CM_HALOFINDER_MICE_L3072_N4096_M4096_R101_LC1.input"

# Read in the info file
info = open(Infofile,'r')
output = open(Outfile, 'w')
info.readline()  # Remove the header
for line in info:
  ln = line.split()
  filenum = int(ln[0])
  xmin = float(ln[1])
  ymin = float(ln[2])
  zmin = float(ln[3])
  xmax = float(ln[4])
  ymax = float(ln[5])
  zmax = float(ln[6])
  npart = int(ln[7])
  if npart > 0:
      if ((Lxmin <= xmax) and (Lxmax >= xmin) and (Lymin <= ymax) and (Lymax >= ymin) and (Lzmin <= zmax) and (Lzmax >= zmin)):
          strout = str(Infile+".%d\n" % filenum)
          output.write(strout)
info.close()
output.close()


