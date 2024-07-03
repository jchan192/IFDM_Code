			------------------------
			|     User Guide       |
			------------------------

(1) Running IFDM
(2) Read output
--------------------------------------------------------------------------------------------

How to run the IFDM?
./job.sh 

There are 10 parameters:
mfdm    = particle mass of mfdm in eV 
g       = interaction term in eV cm^3
xmax    = Boxlength / 2
res     = number of grids in 1 axis. Total grid = res**3 
t_write = write output at every t_write step
t_stop  = stop at step t_stop
ti      = starting step. If rerun = 0, ti = 0. If rerun = 1, ti = preferred snapshot number. 
opath   = path to output directory. Example:: output/binary/
icpath  = path to input file.       Example:: IC/sorted_icN8z50L10.dat
rerun   = 0 if it is a new simulation. = 1 If you want to rerun a simulation.


Example:
#                          mfdm           g        xmax   res    t_write   t_stop    ti          opath                 icpath                           rerun
 mpirun -n 16 ./BECMPI.o   2.5e-24        0        5.     128    100        1000     0        ${ofile}/binary3/     ${icfile}/sorted_icN8z50L10.dat      0    

Example of reruning a simulation from ti = 1000
#                          mfdm           g        xmax   res    t_write   t_stop    ti          opath                  icpath                           rerun
 mpirun -n 16 ./BECMPI.o   2.5e-24        0        5.     128    100        2000     1000        ${ofile}/binary3/     ${ofile}/binary3/                 1


Notes: Reruning a simulation needs to change the icpath to the directory that contains the starting snapshot.


---------------------------------------------------------------------------------------------

How to read the outputs?
In the src directory, make readsnap.
./read res zslice snapnum ifile ofile

There are 5 parameters:
res = number of grids in 1 axis
zslice = the position of the slice along the z-axis 
snapnum = the number of the output. Example: output999999.dat has snapnum=999999
ifile = input directory  
ofile = output directory

Example:
./readsnap 128 64 999998 output/binary3/ output/slice3/
It reads a snapshot at output/binary3/output999998.dat and output the slice to output/binary3/output999998.dat.z64

Output file contains 4 columns: 
x [Mpc] y[Mpc] density[solarMass/Mpc**3] phase of wavefunction. 
