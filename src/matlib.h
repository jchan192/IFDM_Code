#include <complex>
#include <vector>
#include <iostream>
#include <fftw3.h>
#include <random>
#include <functional>

using complex = std::complex<double>; 
using vector_real = std::vector<double>;
using vector_complex = std::vector<std::complex<double>>;

// void fft(vector_complex &x, bool inverse)
// {
//   vector_complex y(x.size(), complex(0.0, 0.0));
//   fftw_plan p;

//   fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());   // standard declaration of input fttw
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw

//   p = fftw_plan_dft_1d(x.size(), in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

//   fftw_execute(p);
//   fftw_destroy_plan(p);

//   for (size_t i = 0; i < x.size(); ++i)
//   {
//     x[i] = y[i] / sqrt(static_cast<double>(x.size()));   // fftw give unnormalized output, needs renormalization here.
//   }
// }  

// void fft2d(vector_complex &x, bool inverse, int threadn)
// {
//   vector_complex y(x.size(), complex(0.0, 0.0));
//   fftw_plan p;

//   fftw_plan_with_nthreads(threadn);  
//   fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());   // standard declaration of input fttw
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw

//   p = fftw_plan_dft_2d(sqrt(x.size()), sqrt(x.size()), in, out, (inverse ? 1 : -1), FFTW_ESTIMATE);

//   fftw_execute(p);
//   fftw_destroy_plan(p);
//   for (size_t i = 0; i < x.size(); ++i)
//   {
//     x[i] = y[i] / sqrt(static_cast<double>(x.size()));   // fftw give unnormalized output, needs renormalization here.
//   }

// }
void fft3d(vector_complex &x, bool inverse, int threadn)
{
  vector_complex y(x.size(), complex(0.0, 0.0));
  fftw_plan p;
  int n = (int)pow(x.size(),1./3)+1;

  fftw_plan_with_nthreads(threadn);  
  fftw_complex *in = reinterpret_cast<fftw_complex*>(x.data());   // standard declaration of input fttw
  fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw

  p = fftw_plan_dft_3d(n, n, n, in, out, (inverse ? FFTW_BACKWARD : FFTW_FORWARD), FFTW_ESTIMATE);

  fftw_execute(p);
  fftw_destroy_plan(p);
  for (size_t i = 0; i < x.size(); ++i)
  {
    x[i] = y[i] / sqrt(static_cast<double>(pow(n,3)));   // fftw give unnormalized output, needs renormalization here.
  }
}

// void fft_r2c(vector_real &x, vector_complex &y, bool forward)
// {
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(y.data());  // standard declaration of output fttw
//   fftw_plan p;

//   if (forward)
//   {
//     p = fftw_plan_dft_r2c_1d(x.size(), x.data(), out, FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < x.size(); ++i)
//     {
//       y[i] = y[i]/sqrt(x.size());
//     }
//   }

//   else 
//   {
//     p = fftw_plan_dft_c2r_1d(x.size(), out, x.data(), FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < x.size(); ++i)
//     {
//       x[i] = x[i]/sqrt(x.size());
//     }
//   }
// }  
// void fft2d_r2c(vector_real &ifxy, vector_complex &fxy, bool forward)
// {
//   fftw_complex *out = reinterpret_cast<fftw_complex*>(fxy.data());  // standard declaration of output fttw
//   fftw_plan p;

//   int n = sqrt(ifxy.size());
//   if (forward)
//   {
//     p = fftw_plan_dft_r2c_2d(n, n, ifxy.data(), out, FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < fxy.size(); ++i)
//     {
//       fxy[i] = fxy[i]/sqrt(ifxy.size());
//     }
//   }

//   else 
//   {
//     p = fftw_plan_dft_c2r_2d(n, n, out, ifxy.data(), FFTW_ESTIMATE);
//     fftw_execute(p);
//     fftw_destroy_plan(p);
//     for (size_t i = 0; i < ifxy.size(); ++i)
//     {
//       ifxy[i] = ifxy[i]/sqrt(ifxy.size());
//     }
//   }
// }  


vector_real der(vector_real x, vector_real y, int n)
{
  // Centered difference for periodic function y(x)

  vector_real dy(x.size());
  for (int order = 0; order < n ; ++order)
  {
    dy[0] = (y[1]-y[x.size()-1])/(x[2]-x[0]); //periodical 
    dy[x.size()-1] = (y[x.size()-2]-y[0])/(x[2]-x[0]);
  
    for (int i = 1; i < x.size(); i++)
      {
	dy[i] = (y[i-1]-y[i+1])/(x[i-1]-x[i+1]);
      }
    y = dy;
  }
  return dy;
}

void pder2d(vector_real x, vector_real f, vector_real &df, int axis)
{
  // Centered difference for periodic function y(x)

  int n = x.size();
  double h = x[2] - x[1];
  
  if (axis==0)
  {
    // df/dx
    for (int j = 0; j<n; j++)
    {
      df[j] = (f[j+n]-f[j+n*(n-1)])/(2*h);
      df[j+n*(n-1)] = (f[j]-f[j+n*(n-2)])/(2*h);    
    }
    for (int i=1; i<(n-1); i++)
    {
	for (int j=0; j<n; j++) 
	{
	  df[j+i*n] = (f[j+(i+1)*n]-f[j+(i-1)*n])/(2*h);
	}
//        if (i==10) {std::cout << df[i] << "\t";}
    }
  }

  if (axis==1)
  {
    // df/dy
    for (int i = 0; i<n; i++)
    {
      df[i*n] = (f[i*n+1]-f[i*n-1])/(2*h);
      df[(i+1)*n-1] = (f[i*n]-f[(i+1)*n-2])/(2*h);    
    }
    for (int i=0; i<n; i++)
    {
	for (int j=1; j<(n-1); j++) 
	{
	  df[j+i*n] = (f[(j+1)+i*n]-f[(j-1)+i*n])/(2*h);
	}
//        if (i==10) {std::cout << df[i] << "\n";}
    }
  }

}

double interpolate( std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate )
{
   int size = xData.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] )                                                 // special case: beyond right end
   {
      i = size - 2;
   }
   else
   {
      while ( x > xData[i+1] ) i++;
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}

float Interp2d(vector_real &x, vector_real &y, vector_real &q, double min, vector_real &p)
{
  // Bilinear interpolation: 
  // input are row major data, so q[n*i+j] = q[i][j] & q.size() = x.size()**2
  // p = target point = [xp,yp]
  int  ix1, ix2, iy1, iy2;
  float q11, q21, q12, q22, h, s;
  float x2x1, y2y1, x2x, y2y, yy1, xx1;

  s = x.size();
  h = x[1]-x[0];
  ix1 = floor((p[0]-min)/h);
  iy1 = floor((p[1]-min)/h);
  ix2 = ix1+1;
  iy2 = iy1+1;

//  std::cout<< p[0] <<"\t"<<p[1]<<"\t"<<ix1<<"\t"<<ix2<<"\t"<<iy1<<"\t"<<iy2<<"\n";

  q11 = q[ix1*s+iy1];
  q21 = q[ix2*s+iy1];
  q12 = q[ix1*s+iy2];
  q22 = q[ix2*s+iy2];

//  std::cout<< q11 <<"\t"<<q21<<"\t"<<q12<<"\t"<<q22<<"\t"<<"\n";

  x2x1 = x[ix2] - x[ix1];
  y2y1 = y[iy2] - y[iy1];
  x2x = x[ix2] - p[0];
  y2y = y[iy2] - p[1];
  yy1 = p[1] - y[iy1];
  xx1 = p[0] - x[ix1];

  return 1.0 / (x2x1 * y2y1) * (
				q11 * x2x * y2y +
				q21 * xx1 * y2y +
				q12 * x2x * yy1 +
				q22 * xx1 * yy1
				);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


double Triple_Trapazoidal(vector_real &x, vector_real &y, vector_real &z, vector_real &f){

  int N = x.size();
  int i,j,k;
  double dx = x[1]-x[0];
  vector_real tempf(N*N), tempf2(N); 
  double result;

  for(i=0;i<N;i++){
  for(j=0;j<N;j++){
    tempf[j+i*N] = (f[0+N*(j+i*N)] + f[(N-1)+N*(j+i*N)]) * 0.5 * dx;
  for(k=0;k<N;k++){
    tempf[j+i*N] += (f[k+N*(j+i*N)] + f[(k-1)+N*(j+i*N)]) * 0.5 * dx;
  }}}

  for(i=0;i<N;i++){
    tempf2[i] = (tempf[0+i*N] + tempf[(N-1)+i*N]) * 0.5 * dx;
  for(j=0;j<N;j++){
    tempf2[i] += (tempf[j+i*N] + tempf[(j-1)+i*N]) * 0.5 * dx;
  }}

  result = (tempf2[0] + tempf2[N-1]) * 0.5 * dx;
  for(i=0;i<N;i++){
    result += (tempf2[i] + tempf2[i-1]) * 0.5 * dx;
  }
  return result;
}

