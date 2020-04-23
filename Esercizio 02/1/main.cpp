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

   // ----------------------------------------- Exercise 02.1 -----------------------------------------
   // Declaring variables
   int M = 1000000, N = 100, L = M/N;
   vector<double> r(M), ave(N), sum_prog(N), su2_prog(N), err_prog(N);

   // Calculating averages and averages squared
   for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         sum += M_PI*0.5*cos(M_PI*0.5*rnd.Rannyu());
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
  print(M, N, sum_prog, err_prog, "out1.txt", 1);


  // Part 2
  for (int i = 0; i< N; ++i) {
    sum_prog[i] = 0;
    su2_prog[i] = 0;
  }

  // Calculating averages and averages squared
  // random numbers non uniform but with a distribution p(x) = 2*(1-x) (Taylor centered in x=1 and normalized)
  for (int i = 0; i<N; ++i) {
     double sum = 0;
     for (int j = 0; j<L; ++j) {
        //double rand = 2.*asin(rnd.Rannyu())/M_PI;
        //sum += rand;
        double rand = sqrt(-rnd.Rannyu() + 1.) + 1.;
        sum += M_PI/(4.-4.*rand)*cos(M_PI*0.5*rand);
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

 // Printing on "out2.txt"
 print(M, N, sum_prog, err_prog, "out2.txt", 1);

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
