#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <string>

using namespace std;

double error (double AV, double AV2, int n = 1);

void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile, double shift = 0.);

void read (vector<double>& r, string infile);

int main (int argc, char *argv[]){
  // ----------------------------------------- Exercise 04.2 -----------------------------------------
  int M, N, L = 2;
  ifstream in("../MolecularDynamics_NVE/output_temp.dat");
  in >> M;
  N = M/L;
  in.close();

  string fileroot = "31_out";
  vector<double> r(M), ave(N), sum_prog(N), su2_prog(N), err_prog(N);

  // Temperatures
  read(r, "../MolecularDynamics_NVE/output_temp.dat");
  // Calculating averages and averages squared
  for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         sum += r[k];
      }
      ave[i] = sum/((double)L);
   }

  // Calculating sum_prog and su2_prog
  for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += ave[j]*ave[j];
      }
      sum_prog[i] = sum_prog[i]/(i+1);
      su2_prog[i] = su2_prog[i]/(double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }
  // Printing on "out1.txt"
  print(M, N, sum_prog, err_prog, fileroot+"_temp.txt");



  // Potential energy
  for (int i = 0; i<N; ++i) {
    sum_prog[i] = 0.;
    su2_prog[i] = 0.;
  }
  read(r, "../MolecularDynamics_NVE/output_epot.dat");
  // Calculating averages and averages squared
  for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         sum += r[k];
      }
      ave[i] = sum/((double)L);
   }

  // Calculating sum_prog and su2_prog
  for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += ave[j]*ave[j];
      }
      sum_prog[i] = sum_prog[i]/(i+1);
      su2_prog[i] = su2_prog[i]/(double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }
  // Printing on "out1.txt"
  print(M, N, sum_prog, err_prog, fileroot+"_epot.txt");



  // Kinetic energy
  for (int i = 0; i<N; ++i) {
    sum_prog[i] = 0.;
    su2_prog[i] = 0.;
  }
  read(r, "../MolecularDynamics_NVE/output_ekin.dat");
  // Calculating averages and averages squared
  for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         sum += r[k];
      }
      ave[i] = sum/((double)L);
   }

  // Calculating sum_prog and su2_prog
  for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += ave[j]*ave[j];
      }
      sum_prog[i] = sum_prog[i]/(i+1);
      su2_prog[i] = su2_prog[i]/(double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }
  // Printing on "out1.txt"
  print(M, N, sum_prog, err_prog, fileroot+"_ekin.txt");

  // Total energy
  for (int i = 0; i<N; ++i) {
    sum_prog[i] = 0.;
    su2_prog[i] = 0.;
  }
  read(r, "../MolecularDynamics_NVE/output_etot.dat");
  // Calculating averages and averages squared
  for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         sum += r[k];
      }
      ave[i] = sum/((double)L);
   }

  // Calculating sum_prog and su2_prog
  for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += ave[j]*ave[j];
      }
      sum_prog[i] = sum_prog[i]/(i+1);
      su2_prog[i] = su2_prog[i]/(double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }
  // Printing on "out1.txt"
  print(M, N, sum_prog, err_prog, fileroot+"_etot.txt");
  return 0;
}

double error (double AV, double AV2, int n) {
   return n == 0 ? 0. : sqrt((AV2 - AV*AV)/(double)n);
}

void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile, double shift) {
   ofstream out(outfile);
   out << M << endl << N << endl;
   for (auto i = sum.begin(); i != sum.end(); ++i)
      out << *i-shift << endl;
   for (auto i = err.begin(); i != err.end(); ++i)
      out << *i << endl;
   out.close();
}

void read (vector<double>& r, string infile) {
  ifstream in(infile);
  cout << "Calculating averages and errors for " << infile << endl;
  int M;
  in >> M;
  for (int i = 0; i<M; ++i) {
    in >> r[i];
  }
  in.close();
  return;
}
