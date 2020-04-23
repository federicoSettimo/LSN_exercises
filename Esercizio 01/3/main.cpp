#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

// Calculates error
double error (double AV, double AV2, int n = 1);

// Prints M, N, sum_prog shifted by 'shift' and err_prog on the file 'outfile'
void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile, double shift = 0.);

int main (int argc, char *argv[]){
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // ----------------------------------------- Exercise 01.1 -----------------------------------------
   // Declaring variables
   int M = 1000000, N = 100, L = M/N;
   vector<double> r(M), ave(N), sum_prog(N), su2_prog(N), err_prog(N);
   vector<double> av2(N);
   double LL = 10., d = 50., p = 2*LL/(d*M_PI);

   // Calculating averages and averages squared
   for (int i = 0; i<N; ++i) {
      int N_in = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         if (rnd.Rannyu() <= p)
          N_in += 1;
      }
      double p_i = (double)N_in/L;
      ave[i] = 2*LL/(p_i*d);
      av2[i] = ave[i]*ave[i];
   }

   // Calculating sum_prog and su2_prog
   for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += av2[j];
      }
      sum_prog[i] = sum_prog[i]/(i+1);
      su2_prog[i] = su2_prog[i]/(double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }

  // Printing on "out.txt"
  print(M, N, sum_prog, err_prog, "out.txt", M_PI);

  rnd.SaveSeed();
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
