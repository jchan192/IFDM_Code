extern const double bigG, ivR, ivL, Mh1, Mh2, hbar, c;
extern double mfdm, mfdm_solar;
struct MPI_Params
{
  int my_rank, total_rank, local_n0;
  std::vector<int> np;
  MPI_Params(int _my_rank, int _total_rank, int _res){
    my_rank=_my_rank;
    total_rank = _total_rank;
    local_n0 = _res/total_rank;
    for (int i=0 ; i < total_rank; i++){
      np.push_back(0);
    }
  }
};

struct Params
{      
  double xmax;
  unsigned int res;
  double dt_KE, dt_PE, dt_QI;
  unsigned int timesteps;
  double dx, dy, dz;
  vector_real x, y, z;
  double dkx, dky, dkz;
  vector_real kx, ky, kz;
  bool im_time;
  int np;
  vector_real px, py;
  vector_real pvx, pvy;
  double ai;
  int ti;
  vector_real a;
  double rhoc;
  int output_counter;
  bool zoomin;
  int twrite;
  double fftnorm;
  std::stringstream icpath;
  std::stringstream opath;

  // initialization
  Params(double _xmax, unsigned int _res, int _twrite, unsigned int _timesteps, int _np, bool im, double offset[2], bool _zoomin, std::string &_icpath, std::string &_opath, int _ti)
  {

    xmax = _xmax;
    res = _res;
    twrite = _twrite;
    ti = _ti;
    timesteps = _timesteps;
    a.reserve((int) timesteps + 1.0);
    dx = 2.0 * xmax / res;
    dy = 2.0 * xmax / res;
    dz = 2.0 * xmax / res;
    x.reserve(res);
    y.reserve(res);
    z.reserve(res);
    dkx = M_PI / xmax;
    dky = M_PI / xmax;
    dkz = M_PI / xmax;
    kx.reserve(res);
    ky.reserve(res);
    kz.reserve(res);
    im_time = im;
    np = _np;
    zoomin = _zoomin;
    icpath <<  _icpath;
    opath <<  _opath;
    fftnorm = 1./static_cast<double>(pow(res,3));

    for (size_t i = 0; i < res; ++i){
      x.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));
      y.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));
      z.emplace_back(xmax / res - xmax + i * (2.0 * xmax / res));

      if (i < res / 2){
	  kx.push_back(i * M_PI / xmax);
	  ky.push_back(i * M_PI / xmax);
	  kz.push_back(i * M_PI / xmax);
      }
      else{
	  kx.push_back((static_cast<double>(i) - res) * M_PI / xmax);
	  ky.push_back((static_cast<double>(i) - res) * M_PI / xmax);
	  kz.push_back((static_cast<double>(i) - res) * M_PI / xmax);
      }
    }
  }

};

double SolitonDensity (double r, double M )
{
  double rc = 1.6 * pow(M/1e9, -1./3) * (1e-22/mfdm) ; // kpc
  double rho0 = 3.1e15 * pow(2.5e-22/mfdm,2) * pow(rc,-4); // Msolar/Mpc**3
  return rho0 * pow(1+0.091*pow(r/rc, 2), -8); // Msolar/Mpc**3
}

struct Operators
{
  size_t size;
  vector_real phi, KE_phase;  // potential
  vector_complex wfc, phik, phikk; // wavefunction

  public:
  Operators(MPI_Params &mpi, Params &par, double voffset, double wfcoffset[3])
  {
    size = par.res;
    phik.reserve(size * size* mpi.local_n0);
    phikk.reserve(size * mpi.local_n0 * (size/2+1));
    phi.reserve(size * mpi.local_n0 * 2 * (size/2+1));
    wfc.reserve(size*size * mpi.local_n0);
    KE_phase.reserve(size*size * mpi.local_n0);

    for (size_t i = 0; i < mpi.local_n0; ++i){
      for (size_t j = 0; j < size; ++j){
        for (size_t k = 0; k < size; ++k){
          phi.push_back(0.0);
          if (k < (size/2+1)) phikk.push_back(complex(0.0,0.0));
          phik.push_back(complex(0.0,0.0)); // for velocity poisson solver
          wfc.push_back(complex(0.0,0.0));
          KE_phase.push_back((0.5 * hbar/(mfdm/pow(c,2)) * (pow(par.kx[i + mpi.my_rank * mpi.local_n0], 2) + pow(par.ky[j], 2) + pow(par.kz[k], 2)))); }
      }
    }
  };
};

double Ha2_preintegrat (double x, void *p) {
  double H_0 = ((double *) p)[0];
  double Omega_m0 = ((double *) p)[1];
  double H = H_0 * sqrt(Omega_m0/pow(x,3) + (1.0 - Omega_m0));
  return 1./(H * pow(x,2));
}

double Ha3_preintegrat (double x, void *p) {
  double H_0 = ((double *) p)[0];
  double Omega_m0 = ((double *) p)[1];
  double H = H_0 * sqrt(Omega_m0/pow(x,3) + (1.0 - Omega_m0));
  return 1./(H * pow(x,3));
}

double Ha4_preintegrat (double x, void *p) {
  double H_0 = ((double *) p)[0];
  double Omega_m0 = ((double *) p)[1];
  double H = H_0 * sqrt(Omega_m0/pow(x,3) + (1.0 - Omega_m0));
  return 1./(H * pow(x,4));
}
