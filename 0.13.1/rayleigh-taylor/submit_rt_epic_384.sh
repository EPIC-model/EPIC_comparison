#!/bin/bash --login

#SBATCH --job-name=EPIC-RT
  # %x gives job-name (SLURM_JOB_NAME)
  # %j gives jobid (individual SLURM_JOB_ID)
  # %A gives jobid (master     SLURM_ARRAY_JOB_ID)
  # %a gives array task id number
  #  https://slurm.schedmd.com/sbatch.html
#SBATCH --nodes=8
#SBATCH --ntasks=1024
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --account=e710
#SBATCH --partition=standard
#SBATCH --qos=long

# refer to ARCHER2 documentation for more information
#     https://docs.archer2.ac.uk/user-guide/

# set number of threads per process
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export FI_OFI_RXM_SAR_LIMIT=64K
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
export MPI_DIR=$MPICH_DIR
export NETCDF_C_DIR=$NETCDF_DIR
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
export FC=ftn

# assumes MPI-enabled build

# adjust directory to top-level epic install directory (prefix)
EPIC_DIR="/work/e710/e710/mf248/gcc/11.2.0"

srun ../../xthi/src/xthi_mpi_mp

# run from time 0 to 6
#srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_384.config

# run from time 6 to 7 (first restart)
#srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_384.config --restart epic_rt_384_0000000002_parcels.nc

# run from time 7 to 8 (second restart)
#srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_384.config --restart epic_rt_384_0000000003_parcels.nc

# run from time 8 to 9 (third restart)
#srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_384.config --restart epic_rt_384_0000000004_parcels.nc

# run from time 9 to 10 (fourth restart)
#srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_384.config --restart epic_rt_384_0000000005_parcels.nc
