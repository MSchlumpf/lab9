#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u, const double dx, const double xmin,
                const int N);
void stepUW(const int N, const double& dx, const double& dt, const double& V, double* const u);
void stepFTCS(const int N, const double& dx, const double& dt, const double& V, double* const u);
//---------------------------------------
int main(){

  const double tEnd = 5;
  const double V = 1;

  const int N  = 512;
  const double xmin = -10;
  const double xmax =  10;
  const double dx = (xmax-xmin)/(N-1);
  double dt = dx/V/1000;
  const int Na = 10; // Number of output files up to tEnd
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  cout << Na*Nk*dt << endl;
  stringstream strm;

  initialize(u0,dx, xmin,N);
  initialize(u1,dx, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);
  writeToFile(u1, "v_0", dx, xmin, N);
  
  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      // Put call to step function here
      stepUW(N, dx, dt, V, u0);
      stepFTCS(N, dx, dt, V, u1);
      // swap arrays u0 <-> u1,
      // however do not copy values, be more clever ;)
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
   strm.str("");
   strm << "v_" << i;
   writeToFile(u1, strm.str(), dx, xmin, N);
  }


  delete[] u0;
  delete[] u1;
  return 0;
}
//-----------------------------------------------

//-----------------------------------------------
void initialize(double* const u, const double dx, const double xmin,
                const int N)
{
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     if (abs(x)<=1.0)
       u[i] = 1;
     else
      u[i] =0;
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}
//-----------------------------------------------
void stepUW(const int N, const double& dx, const double& dt, const double& V, double* const u){
  double temp1 = u[0], temp2;

  u[0] -= V*dt/dx*(temp1-u[N-1]);

  for(int i=1; i<N; i++){
    temp2 = u[i];
    u[i] -= V*dt/dx*(u[i]-temp1);
    temp1 = temp2;    
  }
  
}
//-----------------------------------------------
void stepFTCS(const int N, const double& dx, const double& dt, const double& V, double* const u){
  double temp0 = u[0], temp1 = u[0], temp2;

  u[0] -= V*dt/2/dx*(u[1]-u[N-1]);

  for(int i=1; i<(N-1); i++){
      temp2 = u[i];
      u[i] -= V*dt/2/dx*(u[i+1]-temp1);
      temp1 = temp2;
  }
  u[N-1] -= V*dt/2/dx*(temp0-temp1);  
}
