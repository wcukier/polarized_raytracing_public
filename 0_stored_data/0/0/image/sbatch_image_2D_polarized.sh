#!/bin/bash

#####SBATCH -J IM2
######SBATCH --constraint=HSW24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --time=03:59:00
#SBATCH --output diag_image_$1.output
#SBATCH --mail-user=wcukier.notify@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail


timedatectl | cat

module purge
module load anaconda3/2021.11
module load gcc/8
module load intel/2021.1.2
module load intel-mpi/intel/2021.3.1
module load hdf5/intel-2021.1/intel-mpi/1.10.6


export OMP_NUM_THREADS=1
dir=$1

# rm ../../../data/image/im_E_*
# rm ../../../data/image/im_Q_*
# rm ../../../data/image/im_U_*
# rm ../../../data/image/*.gif
# rm log.log
# rm log2.log
# rm im

num=$(ls ../../../${dir}/geodesics_sync/geodesics_input_* | wc -l)
files=(../../../${dir}/geodesics_sync/geodesics_input_*)
first=${files[0]}
# sub=${first:45:6}


sub=$(echo $first | grep -o '_[0-9][0-9]*' | grep -Eio '[0-9]*')
echo
echo "${sub}"
echo "${num}"
echo

gfortran -g -cpp -fcheck=all -O -fbounds-check -mcmodel=large loop_image_2D_polarized.f90 -o im${dir} -fopenmp -Dnxy=200 -DFRAC=1 -DDIRECTORY=\"$dir\"

srun ./im${dir} ${num} ${sub}

rm im${dir}

gfortran -g -cpp -fcheck=all -O -fbounds-check -mcmodel=large loop_image_2D_direct.f90 -o imd${dir} -fopenmp -Dnxy=200 -DFRAC=1 -DDIRECTORY=\"$dir\"

srun ./imd${dir} ${num} ${sub}

rm imd${dir}

# rm fort.*
