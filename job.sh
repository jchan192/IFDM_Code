#PBS -l nodes=10
#PBS -q large-md
#PBS -N FDM1e9
#PBS -l walltime=24:00:00
#PBS -M jchan@astr.tohoku.ac.jp

source /opt/modules/default/init/bash
module switch PrgEnv-cray PrgEnv-gnu
module load gsl/2.4_gnu-7.3
module load fftw3
echo ${FFTW_DIR}
echo ${GSL_DIR}

cd "/work/chanhi/IFDM_Code"

rm BECMPI.o
CC src/CosmoBEC_MPI_stan.cpp -o BECMPI.o --std=c++11 -I${FFTW_DIR}/include  -L${FFTW_DIR}/lib -lfftw3 -lfftw3_mpi -lm -g -lfftw3_threads -fopenmp -I${GSL_DIR}/include -L${GSL_DIR}/lib -lgsl -lgslcblas -ldl -O3

ofile="output"
icfile="IC"
 
# mpirun -n cores ./BECMPI.o    mfdm            g        xmax   res    t_write   t_stop    ti          opath                 icpath                           rerun
aprun -n 256    ./BECMPI.o   7e-23          0         0.2    512    30001     30001     0       ${ofile}/binary/       ${ofile}/binary/                 1   > log.txt
# mpirun -n 1024    ./BECMPI.o   1e-22          0         0.1    1024    26001     26001    26000       ${ofile}/binary/       ${icfile}/sorted_icN8z50L10.dat  0   > log.txt
