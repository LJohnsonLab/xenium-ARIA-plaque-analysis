#!/bin/bash
#SBATCH --job-name=r_job
#SBATCH --output=r_job_%j.out
#SBATCH --error=r_job_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --account=coa_lajo247_uksr
#SBATCH --partition=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jm.arbones@uky.edu

# Start from a clean module environment to avoid conflicts with previously loaded software
module purge

# Load the working R module that provides R and Rscript
module load r-4.4.0-gcc-9.3.0-7nl5q44

# Load the system curl library/module.
# This is needed so the R package "curl" can compile and link correctly.
module load curl-7.76.1-gcc-9.3.0-37a6bcy

# Load the GMP library/module.
# This provides gmp.h and related libraries needed by packages such as leidenAlg.
module load gmp-6.2.1-gcc-9.3.0-iqbexj7

# Load ImageMagick.
# This is needed by the R package "magick".
module load imagemagick-7.0.8-7-gcc-9.3.0-zygjkan

# Add CMake to PATH.
# Some packages, such as fs and therefore Seurat dependencies, need cmake during compilation.
export PATH=/mnt/gpfs3_amd/share/apps/spack/0.23.0/opt/spack/linux-centos8-zen2/gcc-13.3.0/cmake-3.23.5-xtcjpemtsajntmn77vtuii3stexynmft/bin:$PATH

# Define the root directory of the GMP installation.
# We use this below to expose the correct include and library directories to the compiler/linker.
export GMP_ROOT=/mnt/gpfs3_amd/share/apps/spack/0.16.2-3602-af806a8c1e/opt/spack/linux-centos8-zen2/gcc-9.3.0/gmp-6.2.1-iqbexj764rd7nanxakt33capjvtwhds2

# Add GMP header files to the compiler include search path.
# This is what lets compilation find gmp.h.
export CPATH=$GMP_ROOT/include:$CPATH

# Add GMP libraries to the linker search path at compile time.
export LIBRARY_PATH=$GMP_ROOT/lib:$LIBRARY_PATH

# Add GMP libraries to the runtime shared-library search path.
# This helps R load compiled packages that depend on GMP.
export LD_LIBRARY_PATH=$GMP_ROOT/lib:$LD_LIBRARY_PATH


Rscript banksy_parameter_sweep.R

