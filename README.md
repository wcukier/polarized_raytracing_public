This is the codebase I used for my Fall 2022 JP on black hole polarized raytracing, advised by Benjamin Crinquand and Anatoly Spitkovsky.  There are a lot of files here, most of them I did not contribute to.  
- The directory `Polar` should contain the `GRZeltron` code along with the `GEOKERR` code that I used for the simple magnetic geometries.  If you use these please cite  Parfrey et al. (2019) and  Dexter and Agol (2009) respectivly
- The directory `0_stored_data/0/0/image` contains everything relevant I worked on.  I know it is a mess, I might get around to making this more usable someday. In that directory the following files are relevant:
  - `loop_image_2D_polarized.f90` and `loop_image_3D_polarized.f90`  which can be used to get the Stokes parameters for the 2d and 3d simulations respectivly
  - `loop_image_2D_direct.f90` which for the 2d simulation can extract the direct and indirect images separately
  - `sbatch_image_2d_polarized.sh`, `sbatch_image_2d_direct.sh`, and `sbatch_image_3d_polarized.sh` which run the above codes on a `slurm` HPC
  - `Visualize_v2.jl` which is used for creating the images
  - `azimuthal_fourier.jl` which does the azimuthal fourier transform into the β modes
  - `gen_figures.jl` which specifically creates figures at 20˚ and 70˚ for all the runs I did here

If you want to cite this code, please reach out to me at wcukier@princeton.edu --we don't have a paper together yet.
