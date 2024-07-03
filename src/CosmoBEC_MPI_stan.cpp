#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include "matlib.h"
#include "struct3D.h"
#include "write3D.h"
#include "Cosmoinput3D.h"
#include <cmath>
#include "gadget_IC.h"
#include <gsl/gsl_integration.h>
#include <fftw3-mpi.h>
#include <mpi.h>


void Timestep_criteria(MPI_Params &mpi, Params &par, Operators &opr, int tn,
		       gsl_function &dt_KE_integrator, gsl_function &dt_PE_integrator, gsl_function &dt_QI_integrator){

  double dt_KE, dt_PE, dt_QI, maxrho, maxphi, max_dt, error, allmaxrho, allmaxphi;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

//  if (mpi.my_rank==0) std::cout << "Working on Timestep criteria... \n";
  maxphi = std::abs(*max_element(opr.phi.begin(), opr.phi.end(), [](const double &a, const double &b) { return std::abs(a) < std::abs(b); }));
  maxrho = pow(abs(*max_element(opr.wfc.begin(), opr.wfc.end(), [](const complex &a, const complex &b) { return abs(a) < abs(b); })),2);

  MPI_Allreduce(&maxphi, &allmaxphi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&maxrho, &allmaxrho, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    dt_KE = 4./ (3.*M_PI) * (mfdm/pow(c,2))/hbar * pow(par.dx,2) * pow(par.ai,2);
  //  dt_KE = 4./ (9.*M_PI) * (mfdm/pow(c,2))/hbar * pow(par.dx,2) * pow(par.ai,2);
  //  dt_KE = 4. * (mfdm/pow(c,2))/hbar * pow(par.dx,2) * pow(par.ai,2) * 10;
  //  dt_KE = 4. * (mfdm/pow(c,2))/hbar * pow(par.dx,2) * pow(par.ai,2);

  //  if (mfdm <= 1e-24)   dt_KE = 4. * (mfdm/pow(c,2))/hbar * pow(par.dx,2) * pow(par.ai,2) * 10;
  dt_PE = 2.*M_PI * hbar/(mfdm/pow(c,2)) / allmaxphi * par.ai;
  //  dt_PE = 2.*M_PI * hbar/(mfdm/pow(c,2)) / allmaxphi * par.ai /5.;

  dt_QI = 2.*M_PI * hbar/ fabs(g) / allmaxrho * pow(par.ai,3);
  if (g == 0) dt_QI = 0.0;

  if (mpi.my_rank == 0) std::cout << " max dt_KE = " << dt_KE 
  				  << " max dt_PE = " << dt_PE
  				  << " max dt_QI = " << dt_QI  << "\n";

  // find min dt
  max_dt = dt_KE;
  int n = 3;
  if (max_dt > dt_PE) {
    max_dt = dt_PE;
    n = 2;
  }
  if ((max_dt > dt_QI) and (g!=0)) {
    max_dt = dt_QI;
    n = 4;
  }

  // compute a_next
  double H =  H_0 * sqrt(Omega_m0/pow(par.ai,3) + (1.0 - Omega_m0));
  double da = H * par.ai * max_dt;
  double a_next = par.ai + da;
  if (mpi.my_rank == 0 ) std::cout << "ai = " << par.ai << " a_next = " << a_next <<" da = " << da << " n = " << n <<" ";
  
  gsl_integration_qags(&dt_KE_integrator, par.ai, par.ai + da, 0, 1e-10, 1000, w, &par.dt_KE, &error); 
  gsl_integration_qags(&dt_PE_integrator, par.ai, par.ai + da, 0, 1e-10, 1000, w, &par.dt_PE, &error); 
  gsl_integration_qags(&dt_QI_integrator, par.ai, par.ai + da, 0, 1e-10, 1000, w, &par.dt_QI, &error); 
  //  gsl_integration_qags(&dt_PE_integrator, par.ai, par.ai + da * 0.5, 0, 1e-10, 1000, w, &par.dt_PE, &error); 
  // gsl_integration_qags(&dt_QI_integrator, par.ai, par.ai + da * 0.5, 0, 1e-10, 1000, w, &par.dt_QI, &error); 

  par.ai = a_next;

  // if (mpi.my_rank == 0) std::cout << " dt_KE = " << par.dt_KE 
  // 				  << " dt_PE = " << par.dt_PE
  // 				  << " dt_QI = " << par.dt_QI << "\n\n";
  
}

void poisson(MPI_Params &mpi, Params &par, Operators &opr, fftw_plan phik_fplan, fftw_plan phik_bplan)
{  
  int N = opr.size;
  int iglobal;
  
  for (int i=0; i < mpi.local_n0; ++i){
    for (int j=0; j < N; ++j){
      for (int k=0; k < N; ++k){
        opr.phi[k+(2*(N/2+1))*(j+i*N)] = mfdm_solar * (pow(abs(opr.wfc[k+N*(j+i*N)]), 2)) - par.rhoc ; // for poisson solver
      }}}  

  // // Poisson solver in Fourier method

 
  fftw_execute(phik_fplan);

  for (int i = 0; i < mpi.local_n0; ++i){
    for (int j = 0; j < N; ++j){
      for (int k = 0; k < N/2+1; ++k){
        iglobal = i + mpi.my_rank * mpi.local_n0;

        if (pow(par.kx[iglobal],2)+pow(par.ky[j],2)+pow(par.kz[k],2) == 0){
          opr.phikk[k+(N/2+1)*(j+i*N)] = complex(0.0,0.0);  
        }
        else{
          opr.phikk[k+(N/2+1)*(j+i*N)] = - 4.0  * M_PI * bigG * opr.phikk[k+(N/2+1)*(j+i*N)]/(pow(par.kx[iglobal],2) + pow(par.ky[j],2) + pow(par.kz[k],2));
        }
      }
    }
  }

  fftw_execute(phik_bplan);

  for (int i=0; i < mpi.local_n0; ++i){
    for (int j=0; j < N; ++j){
      for (int k=0; k < 2*(N/2+1); ++k){
        if (k < N) opr.phi[k+(2*(N/2+1))*(j+i*N)] *= par.fftnorm;
        else       opr.phi[k+(2*(N/2+1))*(j+i*N)] = 0.0; // pad zero
      }}}
}

void Split_op_MPI(MPI_Params &mpi, Params &par, Operators &opr, fftw_plan wfc_fplan, fftw_plan wfc_bplan, fftw_plan phik_fplan, fftw_plan phik_bplan, int tn)
{ 
  
  int t0 = time(NULL);
  int N = opr.size;
  std::stringstream logtxt;
  
  if (mpi.my_rank == 0) std::cout << " tn = "<< tn <<"\n";

  for (int i = 0; i < mpi.local_n0; ++i){
  for (int j = 0; j < N; ++j){    
  for (int k = 0; k < N; ++k){
    opr.wfc[k+N*(j+N*i)] *= exp(- 0.5 * ( par.dt_PE * opr.phi[k+(2*(N/2+1))*(j+i*N)] * mfdm/pow(c,2)/hbar + par.dt_QI * g /hbar * pow(abs(opr.wfc[k+N*(j+N*i)]),2)) * complex(0.0,1.0));
  }}}

  fftw_execute(wfc_fplan);

  for (int i = 0; i < mpi.local_n0; ++i){
  for (int j = 0; j < N; ++j){    
  for (int k = 0; k < N; ++k){
        opr.wfc[k+N*(j+N*i)] *=  exp(- opr.KE_phase[k+N*(j+N*i)] * par.dt_KE * complex(0.0, 1.0)) * par.fftnorm; 
  }}}

  fftw_execute(wfc_bplan);

  poisson(mpi,par,opr, phik_fplan, phik_bplan);  

  for (int i = 0; i < mpi.local_n0; ++i){
  for (int j = 0; j < N; ++j){    
  for (int k = 0; k < N; ++k){
    opr.wfc[k+N*(j+N*i)] *= exp(- 0.5 * ( par.dt_PE * opr.phi[k+(2*(N/2+1))*(j+i*N)] * mfdm/pow(c,2)/hbar + par.dt_QI * g /hbar * pow(abs(opr.wfc[k+N*(j+N*i)]),2)) * complex(0.0,1.0));
  }}}

  // if (mpi.my_rank==0){
  //   logtxt << "Splitop takes "<<time(NULL)-t0<<"\t"<<std::endl;
  //   writelog(logtxt);
  // }

}
void Pre_FFTW(MPI_Params &mpi, Operators &opr, fftw_plan &phik_fplan, fftw_plan &phik_bplan, fftw_plan &wfc_fplan, fftw_plan &wfc_bplan){

  int N = opr.size;
  vector_complex backup(N*N*mpi.local_n0);

  // for some reason, calling fftw_plan will turn opr.wfc to 0 sometimes. Example N=480
  for (int ijk=0; ijk < N*N*mpi.local_n0; ijk++){
    backup[ijk] = opr.wfc[ijk];
  }

  fftw_complex *outp = reinterpret_cast<fftw_complex*>(&opr.phikk[0]);  // standard declaration of output fttw
  fftw_complex *inw  = reinterpret_cast<fftw_complex*>(&opr.wfc[0]);   // standard declaration of input fttw
  fftw_complex *outw = reinterpret_cast<fftw_complex*>(&opr.wfc[0]);  // standard declaration of output fttw

  phik_fplan = fftw_mpi_plan_dft_r2c_3d(N, N, N, &opr.phi[0], outp, MPI_COMM_WORLD, FFTW_MEASURE);
  phik_bplan = fftw_mpi_plan_dft_c2r_3d(N, N, N, outp, &opr.phi[0], MPI_COMM_WORLD, FFTW_MEASURE);
  wfc_fplan = fftw_mpi_plan_dft_3d(N, N, N, inw, outw, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
  wfc_bplan = fftw_mpi_plan_dft_3d(N, N, N, inw, outw, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);     

  for (int ijk=0; ijk < N*N*mpi.local_n0; ijk++){
    opr.wfc[ijk] = backup[ijk];
  }
}
void Simulation(MPI_Params &mpi, Params &par, Operators &opr) 
{
  int N = opr.size;
  std::stringstream logtxt;

  fftw_plan phik_fplan, phik_bplan, wfc_fplan, wfc_bplan;
  Pre_FFTW(mpi, opr, phik_fplan, phik_bplan, wfc_fplan, wfc_bplan);
  std::cout << "in Simulation : ";
  std::cout << opr.wfc[0] << " " << opr.wfc[1] << " " << opr.wfc[2]<<" from"<< mpi.my_rank<<"\n";

  // initialize scalefactor integrator
  
  double dt_params[] = { H_0, Omega_m0 };
  gsl_function dt_KE_integrator, dt_PE_integrator, dt_QI_integrator;
  dt_QI_integrator.function = &Ha4_preintegrat;
  dt_QI_integrator.params = dt_params;
  dt_KE_integrator.function = &Ha3_preintegrat;
  dt_KE_integrator.params = dt_params;
  dt_PE_integrator.function = &Ha2_preintegrat;
  dt_PE_integrator.params = dt_params;
  //  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  //  std::cout << "in Simulation : ";
  // std::cout << opr.wfc[0] << " " << opr.wfc[1] << " " << opr.wfc[2]<<" from"<< mpi.my_rank<<"\n";

  int write_num = 0;
  int sten_counter = 0;
  int write_sten_counter = 2;

  for (int tn = par.ti; tn < par.timesteps; ++tn){

    //    std::cout << "in Simulation : ";
    //std::cout << opr.wfc[0] << " " << opr.wfc[1] << " " << opr.wfc[2]<<" from"<< mpi.my_rank<<"\n";
    int tstart = time(NULL);

    poisson(mpi, par, opr, phik_fplan, phik_bplan);
    Timestep_criteria(mpi, par, opr, tn, dt_KE_integrator, dt_PE_integrator, dt_QI_integrator);

    // if (par.ai > 0.25 and write_num == 0) {
    //   //  writeall(mpi, par, opr, 9999972);
    //   if (mpi.my_rank == 0) std::cout << "write z = 3 at step" << tn+1 << "\n"; 
    //   write_num++;
    // }
    // else if (par.ai > 0.33 and write_num == 1) {
    //   //      writeall(mpi, par, opr, 9999973); 
    //   if (mpi.my_rank == 0) std::cout << "write z = 2 at step" << tn+1 << "\n"; 
    //   write_num++;
    //   //      break;
    // }
    // else if (par.ai > 0.5 and write_num == 2) {
    //   //      writeall(mpi, par, opr, 999998); 
    //   if (mpi.my_rank == 0) std::cout << "write z = 1 at step" << tn+1 << "\n"; 
    //   write_num++;
    //   //      break;
    // }
    if (par.ai > 1) {
      writeall(mpi, par, opr, 999999); 
      if (mpi.my_rank == 0) std::cout << "write z = 0 at step" << tn+1 << "\n"; 
      write_num++;
      break;
    }
    else if (par.ai > 0.714 and write_num == 4) {
      writeall(mpi, par, opr, 999998); 
      if (mpi.my_rank == 0) std::cout << "write z = 0.4 at step" << tn+1 << "\n"; 
      write_num++;
      //      break;
    }
    else if (par.ai > 0.5 and write_num == 2) {
      writeall(mpi, par, opr, 999997); 
      if (mpi.my_rank == 0) std::cout << "write z = 1 at step" << tn+1 << "\n"; 
      write_num++;
      break;
    }
    // else if (par.ai > 0.333 and write_num == 2) {
    //   writeall(mpi, par, opr, 999996); 
    //   if (mpi.my_rank == 0) std::cout << "write z = 2 at step" << tn+1 << "\n"; 
    //   write_num++;
    //   //      break;
    // }
    else if (par.ai > 0.25 and write_num == 2) {
      writeall(mpi, par, opr, write_sten_counter); 
      if (mpi.my_rank == 0) std::cout << "write z < 3  at step" << write_sten_counter << "\n"; 
      //      write_num++;
      write_sten_counter++;
      //      break;
    }
    else if (par.ai > 0.25 and write_num == 1) {
      writeall(mpi, par, opr, 999995); 
      if (mpi.my_rank == 0) std::cout << "write z = 3 at step" << tn+1 << "\n"; 
      write_num++;
      sten_counter=0;
      //      break;
    }
    else if (par.ai > 0.105 and write_num == 0) {
      writeall(mpi, par, opr, 999994); 
      if (mpi.my_rank == 0) std::cout << "write z = 8.5 at step" << tn+1 << "\n"; 
      write_num++;
      //      break;
    }
    else if (( (tn % par.twrite == 0 )) or tn == (int) par.timesteps-1 ) { 
      writeall(mpi, par, opr, tn+1) ;
    }
    if (mpi.my_rank==0) std::cout << "before write " << write_num << " " << (0.01 * write_num+1) <<"\n";
    sten_counter++;
    // if (par.ai > (0.01 * (write_num+1))){
    //   if (mpi.my_rank==0) std::cout << "writing at i=" << write_num << "\n";
    //   writeall(mpi, par, opr, write_num);
    //   write_num++;
    //   if (par.ai > 1) break;
    // }

    // if (par.ai > 1) {
    //   writeall(mpi, par, opr, 999999); 
    //   if (mpi.my_rank == 0) std::cout << "write z = 0 at step" << tn+1 << "\n"; 
    //   write_num++;
    //   break;
    // }
    // else if (( (tn % par.twrite == 0 )) or tn == (int) par.timesteps-1 ) { 
    //   writeall(mpi, par, opr, tn+1) ;
    // }

     Split_op_MPI(mpi, par, opr, wfc_fplan, wfc_bplan, phik_fplan, phik_bplan, tn);
     if (mpi.my_rank==0){
       std::cout << "simulation takes:"<< time(NULL)-tstart << " seconds"<<"\n";
     }

  }
}


int main(int argc, char **argv)
{

  int tstart;
  int my_rank, total_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_rank);
  fftw_mpi_init();
  std::stringstream logtxt;

  if (argc != 11) {
    if(my_rank==0) std::cout << "ERROR: Arguments "<< argc << "!= 11 \n";
    std::exit(0);
  }
  else{
    mfdm = std::atof(argv[1]);
    mfdm_solar = mfdm * kg2eV / kg2solar; // solarmass
    g = std::atof(argv[2]) / (pow(cm2Mpc,3)) ; // from eV cm^3 to eV Mpc^3
    inp_xmax = std::atof(argv[3]) / h_0; // Mpc
    inp_res = std::atoi(argv[4]);
    inp_twrite = std::atof(argv[5]);
    inp_timestep = std::atof(argv[6]);
    inp_ti = std::atof(argv[7]);
    inp_opath = argv[8];
    inp_icpath = argv[9];
    rerun = (std::atof(argv[10]) == 1)? true:false;
  }

  if(my_rank==0){
    std::cout << " Input parameters: \n"
	      << "mfdm=" << mfdm <<"\n"
              << "g=" << g * pow(cm2Mpc,3) <<"\n"
	      << "xmax = " << inp_xmax <<" Mpc \n"
	      << "Boxsize = " << inp_xmax*2 <<" Mpc \n"
	      << "res = " << inp_res <<"\n"
	      << "twrite = " << inp_twrite <<"\n"
	      << "timestep = " << inp_timestep <<"\n"
	      << "ti = "<< inp_ti <<"\n"
	      << "opath = "<< inp_opath <<"\n\n";
  }

  MPI_Params mpi = MPI_Params(my_rank, total_rank, inp_res);
  Params par = Params(inp_xmax, inp_res, inp_twrite, inp_timestep, inp_np, inp_im, inp_offset, inp_zoomin, inp_icpath, inp_opath, inp_ti); // xmax, res, dt, timestep, np, im
  Operators opr = Operators(mpi, par, 0.0, inp_offset);

  if (my_rank==0){    
    tstart = time(NULL);
    std::cout << " start IFDM... " << "\n"
	      << " running with " << total_rank<< " cores \n\n";
  }

  if (!rerun) Gadget_IC_Generator(mpi,par,opr);
  else read_rerun(mpi, par, opr, par.ti);
  
  //std::cout << "in Simulation : ";
  //std::cout << opr.wfc[0] << " " << opr.wfc[1] << " " << opr.wfc[2]<<" from"<< mpi.my_rank<<"\n";

  Simulation(mpi, par, opr);
  
  if (my_rank==0){
    std::cout << "simulation takes:"<< time(NULL)-tstart << " seconds"<<"\n";
  }
  
  MPI_Finalize();
  return 0;
}
