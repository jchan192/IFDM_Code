#include <mpi.h>
extern const double cm2Mpc;
extern double mfdm, g, mfdm_solar;

using vector_realf = std::vector<float>;

void write_binary_params(MPI_Params &mpi, Params &par, Operators &opr, int t){

  std::stringstream filename_fdm;
  filename_fdm << par.opath.str() << "output" << t << ".dat";
  double gg = g * (pow(cm2Mpc,3));

  std::ofstream ffdm = std::ofstream(filename_fdm.str(), std::ios::out | std::ios::binary);
  
  if (ffdm){
    ffdm.write((char *) &mfdm, sizeof(double));
    ffdm.write((char *) &gg, sizeof(double));
    ffdm.write((char *) &par.ai, sizeof(double));
    ffdm.write((char *) &par.x[0], opr.size * sizeof(double));
  }
  ffdm.close();
  
}
void write_binary_vg (MPI_Params &mpi, Params &par, Operators &opr, vector_realf &vxg){
    // Writing data into a file in the format of:
    // index, density, real potential.

  int N = par.res;

  // file path
  std::stringstream filename_fdm;
  filename_fdm << par.opath.str() << "voutput.dat";

//  std::cout << "writing at" << par.opath.str() << "\n";
  // file name
  MPI_File file;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD,
		filename_fdm.str().c_str(),
		MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY,
		MPI_INFO_NULL,
		&file);
  
  MPI_File_set_view(file,
		    mpi.my_rank * (int) N*N*mpi.local_n0 * sizeof(float),
		    MPI_FLOAT, 
		    MPI_FLOAT,
		    "native",
		    MPI_INFO_NULL);
  
  MPI_File_write_all(file, &vxg[0], N*N*mpi.local_n0, MPI_FLOAT, &status);
  MPI_File_close(&file);
}

void write_binary_wfc (MPI_Params &mpi, Params &par, Operators &opr, int t){
    // Writing data into a file in the format of:
    // index, density, real potential.

  int N = par.res;

  // file path
  std::stringstream filename_fdm;
  filename_fdm << par.opath.str() << "output" << t << ".dat";

//  std::cout << "writing at" << par.opath.str() << "\n";
  // file name
  MPI_File file;
  MPI_Status status;
  MPI_File_open(MPI_COMM_WORLD,
		filename_fdm.str().c_str(),
		MPI_MODE_CREATE|MPI_MODE_APPEND|MPI_MODE_WRONLY,
		MPI_INFO_NULL,
		&file);
  
  MPI_File_write_ordered(file, &opr.wfc[0], 2*N*N*mpi.local_n0, MPI_DOUBLE, &status);
  MPI_File_close(&file);
}

void writelog( std::stringstream &text){

  std::ofstream log_file("output/log_file.txt", std::ios_base::out | std::ios_base::app );
  log_file << text.str() << std::endl;
  text.str(std::string());
} 

void writeall(MPI_Params &mpi, Params &par, Operators &opr, int t){

  std::stringstream logtxt;
  
  if (mpi.my_rank==0){
    logtxt << "Running:" << (int)((double)t/((double)par.timesteps)*100) << "%"; 
//    writelog(logtxt);
    write_binary_params(mpi,par,opr,t);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  write_binary_wfc(mpi,par,opr,t);
}
