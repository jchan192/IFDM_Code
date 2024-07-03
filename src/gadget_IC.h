#include <complex>
#include <vector>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <fftw3-mpi.h>
#include <iomanip>

extern const double cm2Mpc;
extern double mfdm, g, mfdm_solar;

using vector_realf = std::vector<float>;
using vector_int = std::vector<int>;

struct io_header
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};

struct gadget_particle_data
{
  float Pos[3];
  float Vel[3];
  float Mass;
  int Type;

  float Rho, U, Temp, Ne;
};

int CIC_scheme(MPI_Params &mpi, io_header &Gadget_header, Params &par, vector_complex &wfc, 
	       double &Lg, vector_real &xg,
	       float xp,float yp,float zp, int &pn,
	       vector_realf &vxg, vector_realf &vyg, vector_realf &vzg,
	       float vxp, float vyp, float vzp){

  int Ng = par.res;
  double pmass = Gadget_header.mass[1]; // Msolar
  int gi, gj, gk, gii, gjj, gkk;
  double dx,dy,dz,tx,ty,tz;
  bool inside_zoom;
  
//if (mpi.my_rank==0) {std::cout << mpi.my_rank << " " << xp <<  " " << yp << " " << zp << " "
//				 << gi << " " << gj << " " << gk << " " << Lg << "\n" ;}

  gi = (int) floor((xp - 0.5*Lg)/Lg) ;
  gj = (int) floor((yp - 0.5*Lg)/Lg) ;
  gk = (int) floor((zp - 0.5*Lg)/Lg) ;

//  if (mpi.my_rank==1) {std::cout << mpi.my_rank << " " << xp <<  " " << yp << " " << zp << " "
//  				 << gi << " " << gj << " " << gk  << "\n" ;}

  //periodical boundary
  if(gi == -1){
    gi = Ng-1;
    xp += Gadget_header.BoxSize;
  }

  if (gj < 0) {
    gj = Ng-1;
    yp += Gadget_header.BoxSize;
  }

  if (gk < 0) {
    gk = Ng-1;
    zp += Gadget_header.BoxSize;
  }

//   if (mpi.my_rank==0) std::cout <<"\n"<< gi << " "<< gj << " "<< gk << "\n"; 

/*   // CIC */
  dx = fabs((double) xp - xg[gi])/Lg;
  dy = fabs((double) yp - xg[gj])/Lg;
  dz = fabs((double) zp - xg[gk])/Lg;
  tx = 1. - dx;
  ty = 1. - dy;
  tz = 1. - dz;
  gi -= mpi.local_n0*mpi.my_rank;
  gii = gi + 1.;
  gjj = gj + 1.;
  gkk = gk + 1.;

  if (gii == Ng) {gii = 0; }
  if (gjj == Ng) {gjj = 0; } 
  if (gkk == Ng) {gkk = 0; } 

  
  if ((gi >= 0) && (gi < mpi.local_n0)){
//  if (mpi.my_rank==0) std::cout << "for (gi>=0) " << gi << " "<< gj << " "<< gk << " "<< gii << " "<< gjj << " "<< gkk << std::endl; 

    wfc[gk + Ng * (gj + gi * Ng)] += complex(tx * ty * tz * pmass / pow(Lg,3), 0.0);
    wfc[gk + Ng * (gjj + gi * Ng)] += complex(tx * dy * tz * pmass / pow(Lg,3), 0.0);
    wfc[gkk + Ng * (gj + gi * Ng)] += complex(tx * ty * dz * pmass / pow(Lg,3), 0.0);
    wfc[gkk + Ng * (gjj + gi * Ng)] += complex(tx * dy * dz * pmass / pow(Lg,3), 0.0);

//    if (mpi.my_rank==1) std::cout << vxg[gk + Ng * (gj + gi * Ng)] << " "<<  gk + Ng * (gj + gi * Ng)<<"\n";

   vxg[gk + Ng * (gj + gi * Ng)] += tx * ty * tz * vxp / pow(Lg,3);
   vxg[gk + Ng * (gjj + gi * Ng)] += tx * dy * tz * vxp / pow(Lg,3);
   vxg[gkk + Ng * (gj + gi * Ng)] += tx * ty * dz * vxp/ pow(Lg,3);
   vxg[gkk + Ng * (gjj + gi * Ng)] += tx * dy * dz * vxp / pow(Lg,3);

    vyg[gk + Ng * (gj + gi * Ng)] += tx * ty * tz * vyp / pow(Lg,3);
    vyg[gk + Ng * (gjj + gi * Ng)] += tx * dy * tz * vyp / pow(Lg,3);
    vyg[gkk + Ng * (gj + gi * Ng)] += tx * ty * dz * vyp/ pow(Lg,3);
    vyg[gkk + Ng * (gjj + gi * Ng)] += tx * dy * dz * vyp / pow(Lg,3);

    vzg[gk + Ng * (gj + gi * Ng)] += tx * ty * tz * vzp / pow(Lg,3);
    vzg[gk + Ng * (gjj + gi * Ng)] += tx * dy * tz * vzp / pow(Lg,3);
    vzg[gkk + Ng * (gj + gi * Ng)] += tx * ty * dz * vzp/ pow(Lg,3);
    vzg[gkk + Ng * (gjj + gi * Ng)] += tx * dy * dz * vzp / pow(Lg,3);

  }

  if ((gii >= 0) && (gii < mpi.local_n0)){
//    if (mpi.my_rank==0) std::cout << "for (0 <= gii < n0) " << gi << " "<< gj << " "<< gk << " "<< gii << " "<< gjj << " "<< gkk << std::endl; 

    wfc[gk + Ng * (gjj + gii * Ng)] += complex(dx * dy * tz * pmass / pow(Lg,3), 0.0);
    wfc[gk + Ng * (gj + gii * Ng)] += complex(dx * ty * tz * pmass/ pow(Lg,3), 0.0);
    wfc[gkk + Ng * (gj + gii * Ng)] += complex(dx * ty * dz * pmass / pow(Lg,3), 0.0);
    wfc[gkk + Ng * (gjj + gii * Ng)] += complex(dx * dy * dz * pmass / pow(Lg,3), 0.0);

    vxg[gk + Ng * (gjj + gii * Ng)] += dx * dy * tz * vxp / pow(Lg,3);
    vxg[gk + Ng * (gj + gii * Ng)] += dx * ty * tz * vxp / pow(Lg,3);
    vxg[gkk + Ng * (gj + gii * Ng)] += dx * ty * dz * vxp / pow(Lg,3);
    vxg[gkk + Ng * (gjj + gii * Ng)] += dx * dy * dz * vxp / pow(Lg,3);

    vyg[gk + Ng * (gjj + gii * Ng)] += dx * dy * tz * vyp / pow(Lg,3);
    vyg[gk + Ng * (gj + gii * Ng)] += dx * ty * tz * vyp / pow(Lg,3);
    vyg[gkk + Ng * (gj + gii * Ng)] += dx * ty * dz * vyp / pow(Lg,3);
    vyg[gkk + Ng * (gjj + gii * Ng)] += dx * dy * dz * vyp / pow(Lg,3);

    vzg[gk + Ng * (gjj + gii * Ng)] += dx * dy * tz * vzp / pow(Lg,3);
    vzg[gk + Ng * (gj + gii * Ng)] += dx * ty * tz * vzp / pow(Lg,3);
    vzg[gkk + Ng * (gj + gii * Ng)] += dx * ty * dz * vzp / pow(Lg,3);
    vzg[gkk + Ng * (gjj + gii * Ng)] += dx * dy * dz * vzp / pow(Lg,3);
    
  }

//  std::cout << mpi.my_rank << " " << gi <<  " " << gj << " " << gk << " " << wfc[gk + Ng * (gj + gi * Ng)]<< "\n";  
  return 0;
}
void initial_velocity_Poisson_Solver(MPI_Params &mpi,
				     Params &par,
				     Operators &opr,
				     vector_realf &velx,
				     vector_realf &vely,
				     vector_realf &velz,
				     vector_realf &Lvelx,
				     vector_realf &Rvelx
				     ){

  int N = opr.size;
  int n = mpi.local_n0;
  int ip,im, jp,jm, kp,km, iglobal;
  float dvx, dvy, dvz;
  double unit_conversion = 1e3 / m2Mpc * s2Myr; // convert gadget velocity (km/s) to Mpc/Myr

  // for(int i=0; i<n; i++){
  //   for(int j=0; j<N; j++){
  //     for(int k=0; k<N; k++){
  // 	double ip = i +1;
  // 	if (mpi.my_rank==0) std::cout << velx[k+N*(j+ip*N)] << " " << vely[k+N*(j+i*N)] << " "<< velz[k+N*(j+i*N)] <<"\n";
  // 	//    if (mpi.my_rank==0) std::cout << Lvelx[i] << " " << Rvelx[i]  <<"\n";
  //     }}}

  if (mpi.my_rank == 0 ) std::cout << "poisson solver...\n" ;  
//  if (mpi.my_rank == 0) std::cout << velx[0]  <<" " << velx[1]<< " " << velx[2] <<"\n";
//  if (mpi.my_rank == 0) std::cout << Rvelx[0]  <<" " << Rvelx[1]<<" " << Rvelx[2] <<"\n";

  fftw_complex *inp = reinterpret_cast<fftw_complex*>(opr.phik.data());   // standard declaration of input fttw
  fftw_complex *outp = reinterpret_cast<fftw_complex*>(opr.phik.data());  // standard declaration of output fttw
  fftw_plan phik_fplan = fftw_mpi_plan_dft_3d(N, N, N, inp, outp, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);
  fftw_plan phik_bplan = fftw_mpi_plan_dft_3d(N, N, N, inp, outp, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);

  for(int i=0; i<n; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){

	im = i-1;
	jm = j-1;
	km = k-1;
	ip = i+1;
	jp = j+1;
	kp = k+1;

//	if (i==0) {im = n-1;}
	if (j==0) {jm = N-1;}
	if (k==0) {km = N-1;}
//	if (i==N-1) {ip = 0;}
	if (j==N-1) {jp = 0;}
	if (k==N-1) {kp = 0;}
	
	if (i==0 && n==1)     dvx = Rvelx[k+N*j] - Lvelx[k+N*j] ;
	else if (i==0)        dvx = velx[k+N*(j+ip*N)] - Lvelx[k+N*j] ;
	else if (i==n-1) dvx = Rvelx[k+N*j] - velx[k+N*(j+im*N)] ;
	else             dvx = velx[k+N*(j+ip*N)] - velx[k+N*(j+im*N)] ;


	dvy = vely[k+N*(jp+i*N)] - vely[k+N*(jm+i*N)] ;
	dvz = velz[kp+N*(j+i*N)] - velz[km+N*(j+i*N)] ;
	opr.phik[k+N*(j+i*N)] = complex( 0.5 * (dvx + dvy + dvz) * unit_conversion / (par.dx)  , 0.0);
      }
    }
  }
  //  if (mpi.my_rank == 0) std::cout << opr.phik[0] << " "<< opr.phik[1] << " "<< opr.phik[2] << "\n";

  fftw_execute(phik_fplan);

  for(int i=0; i<n; i++){
    for(int j=0; j<N; j++){
      for(int k=0; k<N; k++){
  	iglobal = i + mpi.my_rank * n;
  	if ( pow(par.kx[iglobal],2) + pow(par.ky[j],2) + pow(par.kz[k],2) == 0 ){
  	  opr.phik[k+N*(j+i*N)] = complex(0.0,0.0);
  	}
  	else{
	  //	  opr.phik[k+N*(j+i*N) ] = - pow(par.ai,3./2) * opr.phik[k+N*(j+i*N)]/(pow(par.kx[iglobal]/h_0,2) + pow(par.ky[j]/h_0,2) + pow(par.kz[k]/h_0,2));
	    	  opr.phik[k+N*(j+i*N) ] = - pow(par.ai,3./2) * opr.phik[k+N*(j+i*N)]/(pow(par.kx[iglobal],2) + pow(par.ky[j],2) + pow(par.kz[k],2));
  	}
      }
    }
  }
  fftw_execute(phik_bplan);
//  if (mpi.my_rank == 0) std::cout << opr.phik[0] << " "<< opr.phik[1] << " "<<opr.phik[2] << "\n";

  for (int i=0;i<N*N*n;++i){
    opr.phik[i] /= static_cast<double>(N*N*N);
//    if (mpi.my_rank == 0 ) std::cout << opr.phik[i]<<"\n";
  }
//  if (mpi.my_rank == 0) std::cout << opr.phik[0] << " "<< opr.phik[1] << " "<< opr.phik[2] << "\n\n";

}

bool equal_boxsize(double a, double b, double epsilon = 0.001)
{
  if (fabs(a - b) < epsilon){
    return true;
  }
  else{
    return false;
  }
}

void read_PosVel(io_header &Gadget_header,
		 MPI_Params &mpi,
		 Params &par,
		 gadget_particle_data &dummyP, 
		 int &Ng,
		 double &extent,
		 double &Lg,
		 vector_real &xg,
		 vector_realf &pPos, 
		 vector_realf &pVel,
		 int &chunk,
		 int &tot_chunks){
  
  pPos.clear();
  pVel.clear();
  
  int dummy;
//  par.icpath.str() = par.icpath.str() + "." + std::to_string(chunk);
  std::string buffer;

  buffer =  par.icpath.str() + "." + std::to_string(chunk);
  std::ifstream rfile (buffer, std::ios::in | std::ios::binary);  
  std::ifstream rfile2 (buffer, std::ios::in | std::ios::binary);  
  if (rfile and rfile2){

    if (mpi.my_rank == 0) std::cout << "read file from " <<  buffer <<"\n";
    
#define SKIP rfile.read(reinterpret_cast<char*> (&dummy), sizeof(dummy));
    
    SKIP; 
    rfile.read(reinterpret_cast<char*> (&Gadget_header), sizeof(Gadget_header));
    if (mpi.my_rank==0){
      std::cout << "npart = " << Gadget_header.npart[1] << "\n"
		<< "mass = "  << Gadget_header.mass[1]  << "\n"
		<< "a = " << Gadget_header.time << "\n"
		<< "redshift? = " << Gadget_header.redshift << "\n"
		<< "BoxSize = " << Gadget_header.BoxSize << "\n"
		<< "h = " << Gadget_header.HubbleParam << " " 
		<< std::endl;
    }
    SKIP;
    
    // if zoomin
    if(par.zoomin){
      Gadget_header.BoxSize *= extent ;
      std::cout << "zoomin BoxSize is now = " << Gadget_header.BoxSize << "\n";
    }else{
      extent = 1.0;
    }
    
    if (!equal_boxsize(Gadget_header.BoxSize, 2 * par.xmax * 1e3 * h_0)){
      std::cout << "ERROR: Boxsize of IC is different to xmax in Simulation !" 
		<< 2 * par.xmax * 1e3 / h_0 << " vs " << Gadget_header.BoxSize << "\n";
//      MPI_Finalize();
      std::exit(0);
    }
      
    //intialize CIC
    Gadget_header.mass[1] *= 1e10; // convert to solar mass / h 
    par.ai = Gadget_header.time; // save ai     
    par.rhoc = tot_chunks * Gadget_header.npart[1]  * Gadget_header.mass[1] / pow(Gadget_header.BoxSize ,3) * pow(1e3,3) * pow(h_0,2); // units = M_solar / Mpc**3
    
    if (mpi.my_rank == 0 ) std::cout << "par.rhoc = " << par.rhoc << " M_solar / Mpc**3\n";
    rfile2.ignore(sizeof(dummy) + 
		  sizeof(Gadget_header) +
		  sizeof(dummy) +
		  sizeof(float) +
		  sizeof(float) * 3 *Gadget_header.npart[1] +
		  sizeof(float) * 2);
    
    Lg = Gadget_header.BoxSize/Ng; // kpc
    if (mpi.my_rank == 0 ) std::cout << "Lg = " << Lg << " kpc \n";
    
    for(int i=0; i<Ng; i++){
      xg[i] = i*Lg + 0.5*Lg;
    }
    
    if (mpi.my_rank == 0 ) std::cout << "start reading IC... \n";      
    rfile.read(reinterpret_cast<char*> (&dummyP.Pos[0]), sizeof(float));
        
    for(int pn = 0; pn < Gadget_header.npart[1] * 3; pn++){
      pPos.push_back(0.0);
      pVel.push_back(0.0);
    }
    
    rfile.read(reinterpret_cast<char*> (&pPos[0]), Gadget_header.npart[1]  * 3 * sizeof(float));
    rfile2.read(reinterpret_cast<char*> (&pVel[0]), Gadget_header.npart[1]  * 3 * sizeof(float));
  }
  else{
    if (mpi.my_rank == 0 ) std::cout << " cannot read IC file !!\n" ;
    std::exit(0);
  }
  rfile.close();
  rfile2.close();
}
void reorder(std::vector<float> &vA, std::vector<int> &vI)  
{
  // reorder vA = [x0,y0,z0,x1,y1,z1,x2,y2,z2,... ] based on vI = [id0, id1, id2....]
  // such that newvA = [ vA[id0*3], vA[id0*3+1], vA[id1*3+2], vA[id1*3]... ]
  
  size_t i, j, k;
  float tx, ty, tz;

  std::cout << "Begin Reordering ... \n";
  if ((int) vI.size()*3 != vA.size()){
    std::cout << "Error: Index_vector.size/3 != Pos.size !!\n";
    MPI_Finalize();
  }

  for(i = 0; i < vI.size(); i++){
    if(i != vI[i]){
      tx = vA[3*i];
      ty = vA[3*i+1];
      tz = vA[3*i+2];
      k = i;
      while(i != (j = vI[k])){
	// every move places a value in it's final location
	vA[3*k] = vA[3*j];
	vA[3*k+1] = vA[3*j+1];
	vA[3*k+2] = vA[3*j+2];
	vI[k] = k;
	k = j;
      }
      vA[3*k] = tx;
      vA[3*k+1] = ty;
      vA[3*k+2] = tz;
      vI[k] = k;
    }
  }
}

void SortParticle(vector_realf &pPos, int &np){

  float xp, buff;
  int ngi, sortPos;
  vector_int index(np); 

  std::cout << "Begin Sorting ...  ";
  for (int pn =0; pn < np; pn++){
    index[pn] = pn;
  }

  //sort the x-axis and get indices
  std::sort(index.begin(), 
	    index.end(),
	    [&](const int& a, const int& b) {return (pPos[a*3] < pPos[b*3]);}
	   );

  // reorder particles from small to large x-axis  
  reorder(pPos, index);

}

void Pre_DistributeParticles_withGhost(MPI_Params &mpi, vector_realf &pPos, vector_real &xg, 
				       vector_int &numpart_per_rank, vector_int &displs, vector_int &Lgpr, vector_int &Ldispls){

  // find 

  int Ng = xg.size();
  int dxg = Ng/mpi.total_rank;
  int rank = 0;
  int numpart = pPos.size()/3;
  
  /* for (int node = 0; node < mpi.total_rank; node++){ */
  /*   std::cout << "Before computation: Nodes " << node << " have particles: " << numpart_per_rank[node] */
  /* 	      << " " << Lgpr[node] << " "  << displs[node] << " " << Ldispls[node] << "\n"; */
  /* } */

  // get npr & Lgpr
  Lgpr[0] = 0;
  for(int pn = 0; pn < numpart; pn++){
    if(rank+1 == mpi.total_rank){
      numpart_per_rank[rank] = numpart - std::accumulate(numpart_per_rank.begin(), 
							 numpart_per_rank.end(), 0) + 1;
      break;
    }else if(xg[(rank+1)*dxg] < pPos[pn*3]){
      rank++;
    }else{
      numpart_per_rank[rank]++;
    }

    if( (xg[(rank+1)*dxg-1] < pPos[pn*3]) && (pPos[pn*3] < xg[(rank+1)*dxg] )) Lgpr[rank+1]++;
    //    if((xg[(rank)*dxg] < pPos[pn*3]) && (pPos[pn*3] < xg[(rank+1)*dxg])) numpart_per_rank[rank]++; 
  }
  
  /* // get Lgpr for rank = 0   */
  for(int pn = numpart - 1; pn > 0; pn--){
   if (xg[Ng-1] < pPos[pn*3]) Lgpr[0]++ ;
   else break;
  }

  // get displs
  displs[0] = 0.0;
  Ldispls[0] = numpart - Lgpr[0];
  for(int rank=1; rank < mpi.total_rank; rank++){
    displs[rank]   = std::accumulate(numpart_per_rank.begin(), numpart_per_rank.begin()+rank, 0) ;
   Ldispls[rank]   = displs[rank] - Lgpr[rank];
  }
  
  /* for (int node = 0; node < mpi.total_rank; node++){ */
  /*   std::cout << "Nodes " << node << " have particles: " << numpart_per_rank[node] */
  /* 	      << " " << Lgpr[node] << " "  << displs[node] << " " << Ldispls[node] << "\n"; */
  /* } */

  /* std::cout << "all parts " << numpart <<"\n"; */
  /* std::cout << "all parts in all nodes " << std::accumulate(numpart_per_rank.begin(), numpart_per_rank.end(), 0) << "\n";   */

  /* // consider x,y,z positions */
  for(int rank = 0; rank < mpi.total_rank; rank++){
    numpart_per_rank[rank]*= 3;
    Lgpr[rank] *= 3;
    displs[rank] *= 3;
    Ldispls[rank] *= 3;
  }
}

void DistributeParticles(MPI_Params &mpi, Params &par,
			 vector_realf &pPos, vector_realf &ppPos, 
			 vector_realf &pVel, vector_realf &ppVel, 
			 vector_real &xg, 
			 vector_int &numpart_per_rank, vector_int &displs, 
			 vector_int &Lgpr, vector_int &Ldispls, double &Lg){

  if (mpi.my_rank==0) std::cout << "begin distributParts ...\n";

  MPI_Bcast(&Lg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);      
  MPI_Bcast(&xg[0], par.res, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numpart_per_rank[0], mpi.total_rank, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Lgpr[0], mpi.total_rank, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&displs[0], mpi.total_rank, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ldispls[0], mpi.total_rank, MPI_INT, 0, MPI_COMM_WORLD); 
    
  for(int pn = 0; pn < Lgpr[mpi.my_rank]+numpart_per_rank[mpi.my_rank]; pn++){
    ppPos.push_back(0.0);
    ppVel.push_back(0.0);
  }
  
  MPI_Scatterv(&pPos[0], &Lgpr[0], &Ldispls[0], MPI_FLOAT,
	       &ppPos[0], Lgpr[mpi.my_rank], MPI_FLOAT,
	       0, MPI_COMM_WORLD);
  
  MPI_Scatterv(&pPos[0], &numpart_per_rank[0], &displs[0], MPI_FLOAT,
	       &ppPos[Lgpr[mpi.my_rank]], numpart_per_rank[mpi.my_rank], MPI_FLOAT,
	       0, MPI_COMM_WORLD);
  
  MPI_Scatterv(&pVel[0], &Lgpr[0], &Ldispls[0], MPI_FLOAT,
	       &ppVel[0], Lgpr[mpi.my_rank], MPI_FLOAT,
	       0, MPI_COMM_WORLD);
  
  MPI_Scatterv(&pVel[0], &numpart_per_rank[0], &displs[0], MPI_FLOAT,
	       &ppVel[Lgpr[mpi.my_rank]], numpart_per_rank[mpi.my_rank], MPI_FLOAT,
	       0, MPI_COMM_WORLD);
  
  /* if(mpi.my_rank==1) { */
    /*   std::cout << ppPos.size() <<  " " << Lgpr[mpi.my_rank]<<"\n"; */
    /* for(int pn = 0; pn < (int) Lgpr[mpi.my_rank]/3; pn++){ */
    /*   std::cout << ppVel[pn*3] << " "<< ppVel[pn*3+1] << " "<< ppVel[pn*3+2] << "\n"; */
    /* }} */

}

void Distribute_GhostLRvxg(MPI_Params &mpi, Params &par,
		   vector_realf &vxg, vector_realf &Lvxg, vector_realf &Rvxg){
  

  int Ng = par.res;
  MPI_Request reqs[4]; 

  int Ltag = 0;
  int Rtag = 1;
  int prev_rank, next_rank;

  prev_rank = mpi.my_rank - 1 ;
  next_rank = mpi.my_rank + 1 ;

  if (mpi.my_rank == 0)                 prev_rank = mpi.total_rank - 1;
  if (mpi.my_rank == mpi.total_rank -1) next_rank = 0;

  if (mpi.my_rank == 0 ) std::cout << "Begin Distributing Ghost vxg cells.... \n"; 
  MPI_Irecv(&Lvxg[0], Ng*Ng, MPI_FLOAT, prev_rank, Ltag, MPI_COMM_WORLD, &reqs[0]);
  MPI_Irecv(&Rvxg[0], Ng*Ng, MPI_FLOAT, next_rank, Rtag, MPI_COMM_WORLD, &reqs[1]);

  MPI_Isend(&vxg[Ng*Ng*(mpi.local_n0-1)], Ng*Ng, MPI_FLOAT, next_rank, Ltag, MPI_COMM_WORLD, &reqs[2]);
  MPI_Isend(&vxg[0], Ng*Ng, MPI_FLOAT, prev_rank, Rtag, MPI_COMM_WORLD, &reqs[3]);

  MPI_Waitall(4, reqs, MPI_STATUS_IGNORE);
}
void Renormalization(MPI_Params &mpi,
		     Params &par,
		     Operators &opr){
  
  double oldrhoc, alloldrhoc;
  int Ng = opr.size;
  double unit_conversion = pow(h_0,2) * pow(1e3,3); // Note : Gadget density unit = h**2 * Msolar / kpc**3  --> Msolar / Mpc**3

  oldrhoc = 0;
  for(int i=0; i < mpi.local_n0; i++){
    for(int j=0; j < Ng; j++){
      for(int k=0; k < Ng; k++){
	oldrhoc += opr.wfc[k+Ng*(j+i*Ng)].real() * unit_conversion;
  }}}

  MPI_Allreduce(&oldrhoc, &alloldrhoc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(mpi.my_rank==0) std::cout << "oldrhoc = " << alloldrhoc << "\n";

  double norm = par.rhoc / alloldrhoc * pow(Ng,3);

  for(int i=0; i < mpi.local_n0; i++){
    for(int j=0; j < Ng; j++){
      for(int k=0; k < Ng; k++){
	opr.wfc[k+Ng*(j+i*Ng)] *= norm;
  }}}  
}

void read_gadget(io_header &Gadget_header, 
		 gadget_particle_data &dummyP, 
		 MPI_Params &mpi,
		 Params &par,
		 Operators &opr){

  int ilocal;
  int Ng = opr.size;
  int ng = mpi.local_n0;
  vector_real xg(Ng);
  vector_realf vxg, vyg, vzg;
  vector_realf Lvxg, Rvxg;
  int npart_chunk, nchunk;
  vector_realf pPos, pVel;
  vector_realf ppPos, ppVel;
  double extent = 0.2;  // same as MUSIC extent
  double Lg ; // kpc 
  vector_int numpart_per_rank(mpi.total_rank), displs(mpi.total_rank);
  vector_int Lgpr(mpi.total_rank), Ldispls(mpi.total_rank);  // Left ghost particle per rank
  vector_int Rgpr(mpi.total_rank), Rdispls(mpi.total_rank); // Right ghost particle per rank

//  int c = 0;
  for(int i=0; i<Ng*Ng*mpi.local_n0; i++){
    vxg.push_back(0.0);
    vyg.push_back(0.0);
    vzg.push_back(0.0);
  }

  for(int i=0; i<Ng*Ng; i++){
    Lvxg.push_back(0.0);
    Rvxg.push_back(0.0);
//    if (mpi.my_rank==0) std::cout << Lvxg[i] << " " << Rvxg[i]  <<"\n";
  }


  // start reading

  int tot_chunks = 4;
  for (int chunk = 0; chunk < tot_chunks; chunk++){
    if (mpi.my_rank == 0 ){
      read_PosVel(Gadget_header, mpi, par, dummyP,
		  Ng, extent, Lg,
		  xg, pPos, pVel, chunk, tot_chunks);
      
      Pre_DistributeParticles_withGhost(mpi, pPos, xg, numpart_per_rank, displs, 
					Lgpr, Ldispls);
    }
    DistributeParticles(mpi, par, pPos, ppPos, pVel, ppVel, xg, 
			numpart_per_rank, displs, Lgpr, Ldispls, Lg); 
  
  /*     //  begin CIC */
    MPI_Bcast(&Gadget_header.mass[1], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&par.rhoc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&par.ai, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Gadget_header.BoxSize, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (mpi.my_rank == 0) std::cout << "start  CIC... \n";
    for(int pn = 0; pn < (numpart_per_rank[mpi.my_rank] + Lgpr[mpi.my_rank])/3; pn++){
      CIC_scheme(mpi, Gadget_header, par, opr.wfc, Lg, xg,
    		 ppPos[pn*3], ppPos[pn*3 + 1], ppPos[pn*3+2], pn,
    		 vxg, vyg, vzg, ppVel[pn*3], ppVel[pn*3 + 1], ppVel[pn*3+2]);
    }

    // reset for next ic.dat
    ppPos.clear();
    ppVel.clear();
    std::fill(Lgpr.begin(), Lgpr.end(), 0);
    std::fill(numpart_per_rank.begin(), numpart_per_rank.end(), 0);
    if (mpi.my_rank == 0) std::cout << "\n";
  }
    
/*   //volume weighted velocity */
  for(int i=0; i<Ng*Ng*mpi.local_n0; i++){
    if (opr.wfc[i].real()!=0){
      vxg[i] *= Gadget_header.mass[1] / opr.wfc[i].real();
      vyg[i] *= Gadget_header.mass[1] / opr.wfc[i].real();
      vzg[i] *= Gadget_header.mass[1] / opr.wfc[i].real();
    }
  }

  //  if (mpi.my_rank == 0) std::cout << vxg[0]  <<"\n"
  Distribute_GhostLRvxg(mpi, par, vxg, Lvxg, Rvxg);
  initial_velocity_Poisson_Solver(mpi, par, opr, vxg, vyg, vzg, Lvxg, Rvxg);
  write_binary_vg(mpi, par, opr, vxg);

/* 	// load it to wavefunction */
  double unit_conversion = pow(h_0,2) * pow(1e3,3); // Note : Gadget density unit = h**2 * Msolar / kpc**3  --> Msolar / Mpc**3

  //  Renormalization(mpi, par, opr);

  for(int i=0; i < mpi.local_n0; i++){
    for(int j=0; j < Ng; j++){
      for(int k=0; k < Ng; k++){
	 opr.wfc[k+Ng*(j+i*Ng)] = 1.225 * sqrt(opr.wfc[k+Ng*(j+i*Ng)] * unit_conversion / mfdm_solar) * exp( mfdm/pow(c,2)/hbar * opr.phik[k+Ng*(j+i*Ng)] * complex(0.0,1.0)) ; // Msolar/Mpc^3
//	opr.wfc[k+Ng*(j+i*Ng)] = sqrt(opr.wfc[k+Ng*(j+i*Ng)] * unit_conversion / mfdm_solar)  ; // Msolar/Mpc^3
//	 if (mpi.my_rank == 0) std::cout << opr.wfc[k+Ng*(j+i*Ng)] << "\n";
//	 if (mpi.my_rank==0) std::cout << opr.wfc[k+Ng*(j+i*Ng)] << " " << vxg[k+Ng*(j+i*Ng)] << " " << vyg[k+Ng*(j+i*Ng)] << " "<< vzg[k+Ng*(j+i*Ng)] << " " << opr.phik[k+Ng*(j+i*Ng)] <<"\n";
  }}}
  
}
void scale_factor_generator(Params &par){
 
  /* double logda = (log10(1.0)-log10(par.ai)) / par.timesteps;

   for (int i = 0; i <= par.timesteps ; i++){
   par.a[i] = pow(10,log10(par.ai) + i * logda);
   }
  */
  for (int i = 0; i<=par.timesteps; i++){
    par.a[i] = par.ai + (1.0 - par.ai)/par.timesteps * i;
//    std::cout << par.a[i] << "\n";
  }  
}
void Gadget_IC_Generator(MPI_Params &mpi,
			 Params &par,
			 Operators &opr){
  
  gadget_particle_data dummyP;
  io_header Gadget_header;
  int tstart;
  std::stringstream logtxt;

  if (mpi.my_rank == 0 ) tstart = time(NULL);

  //  std::cout<< "rank = "<<mpi.my_rank << "; total rank =" << mpi.total_rank << "; local size = "<< mpi.local_n0<<"\n\n";
  read_gadget(Gadget_header,dummyP,mpi,par,opr);
  scale_factor_generator(par);

  if (mpi.my_rank==0){
    std::cout << "read gadget IC takes:"<< time(NULL)-tstart << " seconds"<<"\n\n";
  }
  //  writeall(mpi, par, opr, 0);
//  std::cout << opr.wfc[0] << " " << opr.wfc[1] << " " << opr.wfc[2]<<" from"<< mpi.my_rank<<"\n";
 
}

void read_rerun(MPI_Params &mpi, Params &par, Operators &opr, int t){
  
  std::stringstream ifilename;
  int N = opr.size;
  double dummy;


  std::string tmp =  par.icpath.str() + "output" + std::to_string(t) + ".dat";
  const char* cstr = tmp.c_str();


  MPI_File fh;
  MPI_Status status;
  int failed;

  failed = MPI_File_open(MPI_COMM_WORLD, cstr, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  
  if (failed){
    if (mpi.my_rank == 0 ) std::cout << "ERROR: cannot find icpath  !!!\n";
    if (mpi.my_rank == 0 ) std::cout << tmp << "\n";
    exit(EXIT_FAILURE);
  }
  else{
    MPI_File_read(fh, &dummy, 1, MPI_DOUBLE, &status);
    MPI_File_read(fh, &dummy, 1, MPI_DOUBLE, &status);
    MPI_File_read(fh, &par.ai, 1, MPI_DOUBLE, &status);
    MPI_File_read(fh, &par.x[0], N, MPI_DOUBLE, &status);
    dummy  /= pow(cm2Mpc,3);
    
    unsigned long long buffer = (N+3) * sizeof(double) + mpi.my_rank * 2*N*N*mpi.local_n0*sizeof(double);
    MPI_File_seek(fh, buffer, MPI_SEEK_SET);
    MPI_File_read(fh, &opr.wfc[0],  2*N*N*mpi.local_n0, MPI_DOUBLE, &status);
    MPI_File_close(&fh);
  }
  //  writeall(mpi,par,opr,t+1);
  par.rhoc = rhoc;

  /* for(int i=0; i < N*N*mpi.local_n0; i++){ */
  /* 	 opr.wfc[i] *= 2.0 ; // Msolar/Mpc^3 */
  /* } */

}
