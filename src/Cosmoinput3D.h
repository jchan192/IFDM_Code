const int THREAD_NUM = 40;
const double m2kpc = 3.086e19; // m
const double m2Mpc = 3.086e22; // m
const double cm2Mpc = 3.086e24; // cm
const double s2Myr = 60.0 * 60.0 * 24.0 * 365.0 * 1e6; // s
const double kg2eV = 1.783e-36; //kg
const double kg2solar = 1.989e30; // solar 

//const double rhoc = 1.88e-29 * 1.e-3 * 1.e+6 * pow(3.085677581e+22,3)/(2e+30); // in solar mass/(Mpc**3 h**-2) (DOUBLE checked!)
const double rhoc = 3.51059e+10;
const double h_0 = 0.677;
const double H_0 = h_0 * 100 * s2Myr / m2Mpc * 1e3 ; // Myr^-1
const double hbar = 6.58211928e-16 / s2Myr; // eV * Myr

const double Omega_m0 = 0.276;
const double bigG = 6.674e-11 * kg2solar / pow(m2Mpc,3) * pow(s2Myr,2) ;// Mpc**3/Msolar/Myr**2
const double c = 2.99792458e8 / m2Mpc * s2Myr; // Mpc/Myr

const double ivR = 0.0; 
const double ivL = 0.0;
const double Mh1 = 3e10;
const double Mh2 = 1e11;

int inp_np = 500;
bool inp_im = false;
double inp_offset[3] = {0.0025,0.0025,0.0};
//double inp_offset[3] = {0.0,0.0,0.0};
bool inp_zoomin = false;

double mfdm, mfdm_solar, g, inp_xmax, inp_res, inp_timestep, inp_ti;
int inp_twrite;
bool rerun;
std::string inp_opath, inp_icpath;

