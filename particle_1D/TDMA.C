#include "TDMA.h"

std::vector<double> solveTriDiag(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
{

int i = 0;
int n = a.size();
std::vector <double> x(n,0.0);
std::vector <double> f(n, 0.0);
double z;
  x[0] = c[0]/b[0];
  f[0] = d[0]/b[0];
  // forward sweep  (2:n)
  for(i=1;i<n;i++){
     z = 1/( b[i] - a[i]*x[i-1] );
     x[i] = c[i]*z;
     f[i] = ( d[i] - a[i]*f[i-1] )*z;
  }
  // backward sweep (n-1:1)
  for(i=n-2;i>=0;i--){
     f[i] = f[i]- x[i]*f[i+1];
  }

  x.clear();
  return f;

}
