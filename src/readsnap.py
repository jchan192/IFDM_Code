import numpy as np

N=512 # resolution
kg2eV = 1.783e-36; # J
kg2solar = 1.989e30; # solar 

with open('../output/binary6/output1500.dat', 'rb') as f:
    mfdm   = np.fromfile(f,dtype=np.double, count=1)
    lambdda = np.fromfile(f,dtype=np.double, count=1) # ignore
    a = np.fromfile(f,dtype=np.double, count=1)
    x = np.fromfile(f,dtype=np.double, count=N) #coordinates along one axis
    psi = np.fromfile(f,dtype=np.cdouble, count=N*N*N) # wavefunction

rho = np.abs(psi)**2*mfdm*kg2eV / kg2solar # solarmass * Mpc^3
print('m = ',mfdm)
print('a = ',a)
