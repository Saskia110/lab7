#include<iostream>
#include<cmath>
#include<fstream>

using namespace std;

void func(double* f, double* y,const double eta);

int main (){
  
  double k1[2];
  double k2[2];
  double k3[2];
  double y[2];
  const double T = 100; 
  const double dt= 0.001;
  const int N = T/dt;
  const double eta =0.5;
  double tempx;
  double tempy;
  double ytemp[2];
  double ytemp1[2];
  
  y[0]=0.0001;
  y[1]=sqrt(eta)*y[0];
  
  for(int i=0; i<N; i++){ // über k werte , stecken die dann in die unterfkt
  
  func(k1, y, eta);
      
  ytemp[0] = y[0]+dt*0.5*k1[0];
  ytemp[1] = y[1]+dt*0.5*k1[1];
  
  func(k2, ytemp, eta);

  ytemp[0] = y[0]+dt*(-k1[0]+2*k2[0]);
  ytemp[1] = y[1]+dt*(-k1[1]+2*k2[1]);
  
  func(k3, ytemp, eta );
  
  y[0] = y[0] + (dt/6)*(k1[0]+4*k2[0]+k3[0]);
  y[1] = y[1] + (dt/6)*(k1[1]+4*k2[1]+k3[1]);
  
  
  
  cout << (i+1) * dt << "\t" << y[0] << "\t" << y[1] << endl;
    
  }
  
  return 0;
}

void func(double* f, double* y, const double eta){
     
    f[0]=y[1]; // 1. Eintrag von f(y)
    f[1]=(eta-y[0]*y[0])*y[0]; // 2. Eintrag von f(y)
  
}