ifdm: CosmoBEC_MPI_stan.cpp
	mpic++ CosmoBEC_MPI_stan.cpp -o BECMPI.o --std=c++11 -lfftw3 -lfftw3_mpi -lm -g -lfftw3_threads  -lgsl -lgslcblas -ldl -O3

readsnap: readsnap.cpp
	g++ readsnap.cpp -o ../readsnap --std=c++11

