#!/bin/bash

#SBATCH -J IM3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --time=6-23:59:00
#SBATCH --output diag_image_3d.output
#SBATCH --mail-user=wcukier.notify@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-type=fail


date
echo "************************"
cat sbatch_image_3D_polarized.sh
echo "************************"
echo ""
echo ""
echo ""

export OMP_NUM_THREADS=1

rm ../../../data/image/im_*_phase*

num=$(ls ../../../data/geodesics_sync/geodesics_input_* | wc -l)
files=(../../../data/geodesics_sync/geodesics_input_*)
first=${files[0]}
sub=${first:45:6}

gfortran -g -cpp -fcheck=all -O3 -fbounds-check -mcmodel=large loop_image_3D_polarized.f90 -o im -fopenmp

srun ./im ${num} ${sub}


rm fort.*

echo "************************"
date