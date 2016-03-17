/**
        Project: Mechanics of Extreme Thin Composite Layers for Aerospace applications
         Author: Luca Di Stasio
    Institution: Université de Lorraine & Luleå University of Technology
        Version: 10.2015
  Creation date: October 20th, 2015
    Last update: February 1st, 2016

    Description: 
          Input: 
         Output:
         
    @author Luca Di Stasio
    @version 10.2015
*/

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "finitediff.h"

using namespace std;
using namespace boost;
using namespace boost::program_options;

//==================  MAIN  =======================
int main(int argc, char** argv)
{
  /*string filename = "finite_differences_weights.csv";
  
  double arr[] = {-7.0,-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
  vector<double> xs(arr, arr+15);
  
  
  ofstream outfile(filename.c_str());
  outfile << "Order of Derivative, Order of Accuracy, ,Stencil" << endl;
  outfile << ",,";
  for(int i = 0; i<xs.size(); i++){
    outfile << xs[i];
    if(i==xs.size()-1){
      outfile << "\n";
    }else{
      outfile << ",";
    }
  }
  
  int dev_eq;
  int num_nodes;
  vector<double> print_weights;
  vector<double> x;
  vector<double> weights;
  
  double z = 0.0;
  
  for(int dev = 0; dev<5; dev++){
    if(dev%2==0 && dev!=0){
      dev_eq = dev-1;
    }else{
      dev_eq = dev;
    }
    for(int acc = 2; acc<9; acc+=2){
      if(dev!=0){
        num_nodes = dev_eq + acc;
      }else{
        num_nodes = xs.size();
      }
      
      for(int i=(xs.size()-num_nodes)/2; i<(xs.size()-num_nodes)/2+num_nodes; i++){
        x.push_back(xs[i]);
      }
      
      finitediff differences(dev,z,x);
      differences.compute_weights();
      weights = differences.provide_weights();
      
      x.clear();
      
      print_weights.resize(xs.size());
      
      for(int i=0; i<weights.size(); i++){
        print_weights[(xs.size()-num_nodes)/2+i] = weights[i];
      }
      outfile << dev << "," << acc << ",";
      for(int i = 0; i<print_weights.size(); i++){
        outfile << print_weights[i];
        if(i==print_weights.size()-1){
          outfile << "\n";
        }else{
          outfile << ",";
        }
      }
      
      weights.clear();
      print_weights.clear();
    }
  }*/
  
  int M = 1;
  
  double z = 0.0;
  
  double left[] = {-2.0,-1.0,0.0};
  vector<double> left_xs(left, left+3);
  vector<double> left_fs;
  left_fs.resize(left_xs.size());
  for(int i=0; i<left_xs.size(); i++){
    left_fs[i] = -3.0*left_xs[i] + 7.0;
  }
  finitediff left_weight_diff(M,z,left_xs);
  left_weight_diff.compute_weights();
  vector<double> left_weights = left_weight_diff.provide_weights();
  finitediff left_func_diff(M,left_fs,left_weights);
  left_func_diff.compute_derivative();
  double left_dev = left_func_diff.provide_derivative();
  
  cout << endl;
  cout << "LEFT-SIDE FIRST ORDER DERIVATIVE" << endl;
  cout << endl;
  cout << "Nodes :";
  for(int i=0; i<left_xs.size(); i++){
    cout << left_xs[i];
    if(i==left_xs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Function :";
  for(int i=0; i<left_fs.size(); i++){
    cout << left_fs[i];
    if(i==left_fs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Weights :";
  for(int i=0; i<left_weights.size(); i++){
    cout << left_weights[i];
    if(i==left_weights.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Derivative :" << left_dev << endl;
      
  double right[] = {0.0,1.0,2.0};
  vector<double> right_xs(right, right+3);
  vector<double> right_fs;
  right_fs.resize(right_xs.size());
  for(int i=0; i<right_xs.size(); i++){
    right_fs[i] = -3.0*right_xs[i] + 7.0;
  }
  finitediff right_weight_diff(M,z,right_xs);
  right_weight_diff.compute_weights();
  vector<double> right_weights = right_weight_diff.provide_weights();
  finitediff right_func_diff(M,right_fs,right_weights);
  right_func_diff.compute_derivative();
  double right_dev = right_func_diff.provide_derivative();
  
  cout << endl;
  cout << "RIGHT-SIDE FIRST ORDER DERIVATIVE" << endl;
  cout << endl;
  cout << "Nodes :";
  for(int i=0; i<right_xs.size(); i++){
    cout << right_xs[i];
    if(i==right_xs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Function :";
  for(int i=0; i<right_fs.size(); i++){
    cout << right_fs[i];
    if(i==right_fs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Weights :";
  for(int i=0; i<right_weights.size(); i++){
    cout << right_weights[i];
    if(i==right_weights.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Derivative :" << right_dev << endl;
  
  double centered[] = {-1.0,0.0,1.0};
  vector<double> centered_xs(centered, centered+3);
  vector<double> centered_fs;
  centered_fs.resize(centered_xs.size());
  for(int i=0; i<centered_xs.size(); i++){
    centered_fs[i] = -3.0*centered_xs[i] + 7.0;
  }
  finitediff centered_weight_diff(M,z,centered_xs);
  centered_weight_diff.compute_weights();
  vector<double> centered_weights = centered_weight_diff.provide_weights();
  finitediff centered_func_diff(M,centered_fs,centered_weights);
  centered_func_diff.compute_derivative();
  double centered_dev = centered_func_diff.provide_derivative();
  
  cout << endl;
  cout << "CENTERED-SIDE FIRST ORDER DERIVATIVE" << endl;
  cout << endl;
  cout << "Nodes :";
  for(int i=0; i<centered_xs.size(); i++){
    cout << centered_xs[i];
    if(i==centered_xs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Function :";
  for(int i=0; i<centered_fs.size(); i++){
    cout << centered_fs[i];
    if(i==centered_fs.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Weights :";
  for(int i=0; i<centered_weights.size(); i++){
    cout << centered_weights[i];
    if(i==centered_weights.size()-1){
      cout << endl;
    }else{
      cout << ", ";
    }
  }
  cout << "Derivative :" << centered_dev << endl;
  
  
  return 0;
}