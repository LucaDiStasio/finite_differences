#include "finitediff.h"

//=====================  BODY  ==========================

finitediff::finitediff(){}

finitediff::finitediff(int Minp, double zinp,  vector<double> xsinp){
  M = Minp;
  N = xsinp.size();
  z = zinp;
  xs = xsinp;
}
  

finitediff::finitediff(int Minp, vector<double> fsinp, vector<double> weightsinp){
  M = Minp;
  N = weights.size();
  fs = fsinp;
  weights = weightsinp;
}
  

finitediff::finitediff(int Minp, double zinp,  vector<double> xsinp, vector<double> fsinp){
  M = Minp;
  N = xsinp.size();
  z = zinp;
  xs = xsinp;
  fs = fsinp;
}

finitediff::~finitediff(){}

void finitediff::compute_weights(){
  //-----------------------------------------------------------------------------------------------//
  // Pseudo-code from
  // "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", Bengt Fornberg.
  // Mathematics of Computation; Volume 51 (1988), Number 184, Pages 699-706.
  // see also
  // Calculation of Weights in Finite Difference Formulas, Bengt Fornberg.
  // SIAM Reviews; Volume 40 (1998), Number 3, Pages 685-691
  // for a FORTRAN version of the algorithm.
  //-----------------------------------------------------------------------------------------------//
  // Initialize variables
  vector<vector<double> >  delta;    // delta_i,j ==> delta[i][j] 0<=i<N, 0<=j<M
  delta.resize(N);
  for(unsigned int i=0; i<delta.size(); i++){
    delta[i].resize(M+1);
    for(unsigned int j=0; j<delta[i].size(); j++){
      delta[i][j] = 0.0;
    }
  }
  double c1, c2, c3, c4, c5;
  int mn;
  // Algorithm starts
  delta[0][0] = 1.0;
  c1 = 1.0;
  c4 = xs[0] - z;
  for(int i=1; i<N; i++){  
    c2 = 1.0;
    c5 = c4;
    c4 = xs[i] - z;
    mn = min(i,M);
    for(int j=0; j<i; j++){
      c3 = xs[i] - xs[j];
      c2 = c2*c3;
      if(j==i-1){
        for(int k=mn; k>0; k--){
          delta[i][k] = c1*(k*delta[i-1][k-1]-c5*delta[i-1][k])/c2;
        }
        delta[i][0] = -c1*c5*delta[i-1][0]/c2;  
      }
      for(int k=mn; k>0; k--){
        delta[j][k] = (c4*delta[j][k]-k*delta[j][k-1])/c3;
      }
      delta[j][0] = c4*delta[j][0]/c3;
    }
    c1 = c2;
  }
  //Get the final result
  weights.resize(N);
  for(unsigned int i=0; i<weights.size(); i++){
    weights[i] = delta[i][M];
  }
}

void finitediff::compute_derivative(){
  derivative = 0;
  for(unsigned int i=0; i<fs.size(); i++){
    derivative += weights[i]*fs[i];
  }
}
       
vector<double> finitediff::provide_weights(){
  return weights;
}
  
double finitediff::provide_derivative(){
  return derivative;
}
