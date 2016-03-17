#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
#include <cmath>
#include <stdlib.h>
#include <cstdlib>     //to get current directory
#include <unistd.h>    //to get home directory
#include <sys/types.h> //to get home directory
#include <pwd.h>       //to get home directory 
#include <sys/stat.h>     //create new directory
#include <sys/types.h>    //create new directory
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <boost/algorithm/string/join.hpp>
#include <exception>
#include <typeinfo>
#include <omp.h>

using namespace std;
using namespace boost;
using namespace boost::filesystem;

/*
  A class to implement explicit and implicit Finite Difference Method 
  for derivatives of every order m over arbitrarily spaced grids
*/


//===================================================
//==================  HEADER  =======================
//===================================================

class finitediff {

//===================================================  
//                  Variables
//===================================================
private:

// Control parameters
  unsigned int float_operations;                   // Number of floating point operations performed

// Output quantities
  int M;                                           // Order of the highest derivative to approximate
  int N;                                           // Number of grid points
  double z;                                        // Location at which derivative is approximated
  vector<double> xs;                               // Vector of N grid points forming the stencil for the finite difference procedure
  vector<double> fs;                               // Vector of N functional evaluation at the corresponding grid points
  vector<double> weights;                 // Vector of finite difference weights
  double derivative;                               // Approximated value of derivative at location z

  
//===================================================  
//                      Methods
//===================================================  
public:
  
  // Constructor (default)
  finitediff();
  
  //Destructor
  ~finitediff();
  
  // Constructor (init M,z,xs ==> get weights)
  finitediff(int Minp, double zinp,  vector<double> xsinp);
  
  // Constructor (init M,fs,weights ==> get derivative)
  finitediff(int Minp, vector<double> fsinp, vector<double> weightsinp);
  
  // Constructor (init M,z,xs,fs ==> get weights and derivative)
  finitediff(int Minp, double zinp,  vector<double> xsinp, vector<double> fsinp);
  
  void compute_weights(),                                               // Compute stencil weights
       compute_derivative();                                            // Compute derivative
       
  vector<double> provide_weights();                                     // Provide weights to external program
  
  double provide_derivative();                                          // Provide derivative to external program
};