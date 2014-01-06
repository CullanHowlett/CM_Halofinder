#!/bin/bash
#PBS -l nodes=32:ppn=8
#PBS -l walltime=10:00:00
#PBS -N CM_Halofinder
#PBS -j oe
#PBS -m abe
#PBS -q cluster.q
#PBS -V

cd $PBS_O_HOME
. /etc/profile.d/modules.sh
. $HOME/useful_stuff/modules_2lpt

module load torque
module load mpi/gcc/openmpi/1.4.3

#mpirun /users/cullanh/mock_catalogues/CM_Halofinder/CM_Halofinder /mnt/astro4/cullanh/PICOLA/L750_N250_R2008_4mesh/CM_HALOS_2LPT_L750_N250_R2008_4mesh.param
mpirun /users/cullanh/mock_catalogues/CM_Halofinder/CM_Halofinder /mnt/lustre/cullanh/PICOLA/ashleymocks/L1280_N1536_R9001/CM_HALOFINDER_L1280_N1536_R9001.param
#mpirun /users/cullanh/mock_catalogues/CM_Halofinder/CM_Halofinder /users/cullanh/mock_catalogues/CM_Halofinder/files/run_parameters.dat
#mpirun /users/cullanh/mock_catalogues/CM_Halofinder/CM_Halofinder /mnt/lustre/cullanh/GADGET/L1280_N1536_R9001/CM_HALOFINDER_GADGET_L1280_N1536_R9001_z0p15.param




 
