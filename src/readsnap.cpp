#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>

using complex = std::complex<double>; 
using vector_real = std::vector<double>;
using vector_complex = std::vector<std::complex<double>>;

int main(int argc, char **argv){
  
  std::string inp_opath;
  std::string inp_icpath;
  int N, zslice, t;

  if (argc != 6) {
    std::cout << "ERROR: Arguments "<< argc << "!= 5 \n";
    std::exit(0);
  }
  else{
    N = std::atoi(argv[1]);
    zslice = std::atoi(argv[2]);
    t = std::atoi(argv[3]);
    inp_opath = argv[4];
    inp_icpath = argv[5];    
  }

  vector_real x(N);
  vector_complex psi(N*N*N);
  double ai;
  double mfdm, lambda;
  const double kg2eV = 1.783e-36; //J
  const double kg2solar = 1.989e30; // solar 
  std::stringstream ifilename;
  std::stringstream ofilename;

  ifilename << inp_opath << "/output"<< t <<".dat";
  ofilename << inp_icpath << "/output" << t <<".dat.z" << zslice;

  std::ifstream rfile (ifilename.str(), std::ios::in | std::ios::binary);
  std::ofstream ofile (ofilename.str());
  
  // read
  if (rfile){
    rfile.read(reinterpret_cast<char*> (&mfdm), sizeof(double)); 
    rfile.read(reinterpret_cast<char*> (&lambda), sizeof(double));
    rfile.read(reinterpret_cast<char*> (&ai), sizeof(double)); // y and z are the same anyway
    rfile.read(reinterpret_cast<char*> (&x[0]), N*sizeof(double)); // y and z are the same anyway
    rfile.read(reinterpret_cast<char*> (&psi[0]), 2*N*N*N * sizeof(double));
  }else{
    std::cout << "ERROR: cant read binary file !!\n";
  }
  
  // for(int i =0;i<N;++i){
  //   std::cout << x[i] << "\n"; 
  // }
  
  
  std::cout << "writing snapshot mfdm = " << mfdm<< " and g = " << lambda <<" at a = " << ai <<"\n";
  std::cout << "input file = " << ifilename.str() << "\n output file = " << ofilename.str() << "\n";
  mfdm *= kg2eV / kg2solar;
  
  // // write
  if (ofile){
    for (int i = 0; i < N; ++i){
      for (int j = 0; j < N; ++j){
	std::stringstream data_fdm;
	data_fdm << x[i] << " " << x[j] << " "<< pow(abs(psi[(i*N+j)*N+zslice]),2)*mfdm  << " " << std::arg(psi[(i*N+j)*N+zslice])<<"\n"; 
	ofile.write(data_fdm.str().c_str(), data_fdm.str().length());
      }
    }
  }
  
  ifilename.str(std::string());
  ofilename.str(std::string());
  ofile.close();
  return 0;
}
