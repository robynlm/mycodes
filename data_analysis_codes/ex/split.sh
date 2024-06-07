#!/bin/sh
#SBATCH --job-name=split
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=15
#SBATCH --time=72:00:00
#SBATCH --mem=60G
#SBATCH --partition=sciama2.q

now=$(date +"%T")
echo "Start time : $now"

module purge
module load anaconda3/2019.03
module load ffmpeg/4.1.4

module load system/intel64
module load intel_comp/2019.2
module load automake/1.16 pkg-config/0.29 jpeg-turbo/2.0.3
module load gsl/2.5 papi/5.7.0 hwloc/2.1.0 hdf5/1.8.17 openssl/1.1.1
module load curl/7.54.0 libz/1.2.11 perl/5.26
module load fftw/3.3.8


for SIMNAME in pflrw_d3e2_L1821_t1_N64_LCDM_withlapse_1plog pflrw_d3e2_L2500_t2_N64_LCDM 
do
 python split_files.py $SIMNAME
done

now=$(date +"%T")
echo "End time : $now"

