#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <string>
#include "random.h"

using namespace std;

const double S0 = 100., T = 1., K = 100., r = 0.1, sigma = 0.25;

double C(double St, double t);
double P(double St, double t);
double d1(double St, double t);
double d2(double d_1, double t);
double N(double x);

// Prints M, N, sum_prog and err_prog on the file 'outfile'
void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile);
double error (double AV, double AV2, int n = 1);

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
   vector<double> ave_C(N), sum_prog_C(N), su2_prog_C(N), err_prog_C(N);
   vector<double> ave_P(N), sum_prog_P(N), su2_prog_P(N), err_prog_P(N);

   for (int i = 0; i<N; ++i) {
      double sumC = 0., sumP = 0.;
      for (int j = 0; j<L; ++j) {
         //double St = S0*exp((r - 0.5*sigma*sigma)*T + sigma*rnd.Rannyu());      // GBM with mu = r = 0.1 and sigma = 0.25 at t=T
         double St = S0*exp((r - 0.5*sigma*sigma)*T + sigma*rnd.Gauss(0,T));
         //double St = S0*exp(sigma*rnd.Rannyu());
         sumP += P(St, 0.);
         sumC += C(St, 0.);
      }
      ave_P[i] = sumP/((double)L);
      ave_C[i] = sumC/((double)L);
   }

   for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog_P[i] += ave_P[j];
         su2_prog_P[i] += ave_P[j]*ave_P[j];
         sum_prog_C[i] += ave_C[j];
         su2_prog_C[i] += ave_C[j]*ave_C[j];
      }
      sum_prog_P[i] = sum_prog_P[i]/(i+1);
      su2_prog_P[i] = su2_prog_P[i]/(double)(i+1);
      err_prog_P[i] = error(sum_prog_P[i], su2_prog_P[i], i);
      sum_prog_C[i] = sum_prog_C[i]/(i+1);
      su2_prog_C[i] = su2_prog_C[i]/(double)(i+1);
      err_prog_C[i] = error(sum_prog_C[i], su2_prog_C[i], i);
   }

  // Printing
  print(M, N, sum_prog_P, err_prog_P, "out1_P.txt");
  print(M, N, sum_prog_C, err_prog_C, "out1_C.txt");


  // Part 2
  sum_prog_C.assign(N,0.);
  sum_prog_P.assign(N,0.);
  su2_prog_C.assign(N,0.);
  su2_prog_P.assign(N,0.);
  // dividing T in intervals of length dt
  int t_intervals = 100;
  double dt = T/(double)t_intervals;


  for (int i = 0; i<N; ++i) {
     double sumC = 0., sumP = 0.;
     for (int j = 0; j<L; ++j) {
       double S_0 = S0, St = S0;
        for (int ti = 0; ti<t_intervals; ++ti) {
          S_0 = St;                                                                 // S_0 is the position at t_{i-1}
          St = S_0*exp((r - 0.5*sigma*sigma)*dt + sigma*rnd.Gauss(0,1)*sqrt(dt));     // calculating S at t = t_i
        }
        cout << St << endl;
        sumP += P(St, 0.);                                                          // and using the last value of S (S_T) to calculate C and P
        sumC += C(St, 0.);
     }
     ave_P[i] = sumP/((double)L);
     ave_C[i] = sumC/((double)L);
  }

  // Calculating sum_prog and su2_prog
  for  (int i = 0; i<N; ++i) {
     for (int j = 0; j<i+1; ++j) {
        sum_prog_P[i] += ave_P[j];
        su2_prog_P[i] += ave_P[j]*ave_P[j];
        sum_prog_C[i] += ave_C[j];
        su2_prog_C[i] += ave_C[j]*ave_C[j];
     }
     sum_prog_P[i] = sum_prog_P[i]/(i+1);
     su2_prog_P[i] = su2_prog_P[i]/(double)(i+1);
     err_prog_P[i] = error(sum_prog_P[i], su2_prog_P[i], i);
     sum_prog_C[i] = sum_prog_C[i]/(i+1);
     su2_prog_C[i] = su2_prog_C[i]/(double)(i+1);
     err_prog_C[i] = error(sum_prog_C[i], su2_prog_C[i], i);
  }

  print(M, N, sum_prog_P, err_prog_P, "out2_P.txt");
  print(M, N, sum_prog_C, err_prog_C, "out2_C.txt");

   rnd.SaveSeed();
   return 0;
}

double error (double AV, double AV2, int n) {
   return n == 0 ? 0. : sqrt((AV2 - AV*AV)/(double)n);
}

void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile) {
   ofstream out(outfile);
   out << M << endl << N << endl;
   for (auto i = sum.begin(); i != sum.end(); ++i)
      out << *i << endl;
   for (auto i = err.begin(); i != err.end(); ++i)
      out << *i << endl;
   out.close();
}

double C(double St, double t) {
  if (t<T) {
    double d = d1(St,t);
    return St*N(d) - K*exp(-r*(T-t))*N(d2(d,t));
  }
  else return St - K;
}

double P(double St, double t) {
  if (t<T) {
    double d = d1(St,t);
    return St*(N(d)-1.) - K*exp(-r*(T-t))*(N(d2(d,t))-1.);
  }
  else return 0.;
}

double d1(double St, double t) {
  return 1./(sigma*sqrt(T-t))*(log(S0/K) + r + sigma*sigma*0.5*(T-t));
}

double d2(double d_1, double t) {
  return d_1 - sigma*sqrt(T-t);
}

double N(double x) {
  return 0.5*(1. + erf(x/sqrt(2)));
}
