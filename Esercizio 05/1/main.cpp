#include "point.h"
#include "random.h"
#include "functions.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

double error (double AV, double AV2, int n = 1);

void print (int M, int N, const vector<double>& sum, const vector<double>& err, const string outfile, double shift = 0.);

int main (){
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



  // M steps for the simulation
  int M = 100000000;
  // r_max calculated in rmax.cpp
  double rmax = 2.58;
  double rmax2 = 4.84;
  // storing all the r
  vector<double> r(M);
  // starting in (0,0,0): reasonable for l = 0
  // ------------------------------- to do: try different starting points ---------------------------
  point p(1,1,1);
  point pnew;
  cout << "n = 1, l = 0, m = 0\n";
  for (int i = 0; i<M; ++i) {
    if (i%1000000 == 0) cout << "Calculating trajectory... step " << i << " of " << M << endl;
    move(p, pnew, rnd, rmax);
    double alpha = A(pnew, p);
    if (rnd.Rannyu() <= alpha) {
      p.setcoord(pnew);
    }
    r[i] = p.getr();
  }

  // calculating block averages with N steps per block
  int N = 100000, L = M/N;
  vector<double> ave(N), sum_prog(N), su2_prog(N), err_prog(N);
  cout << "\nCalculating block averages with " << N << " steps per block...\n";
  p.setcoord(1,1,4);

  for (int i = 0; i<N; ++i) {
    double sum = 0;
    for (int j = 0; j<L; ++j) {
      int k = j + i*L;
      sum += r[k];
    }
    ave[i] = sum/((double)L);
  }

  for  (int i = 0; i<N; ++i) {
    for (int j = 0; j<i+1; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += ave[j]*ave[j];
    }
    sum_prog[i] = sum_prog[i]/(i+1);
    su2_prog[i] = su2_prog[i]/(double)(i+1);
    err_prog[i] = error(sum_prog[i], su2_prog[i], i);
  }

  print(M, N, sum_prog, err_prog, "out_r.txt", 1.5);



  //------------------------------- to do: same for 2,1,0 ---------------------------
  cout << "\n\nn = 2, l = 1, m = 0\n";
  for (int i = 0; i<M; ++i) {
    if (i%1000000 == 0) cout << "Calculating trajectory... step " << i << " of " << M << endl;
    move(p, pnew, rnd, rmax2);
    double alpha = A2(pnew, p);
    if (rnd.Rannyu() <= alpha) {
      p.setcoord(pnew);
    }
    r[i] = p.getr();
  }

  // calculating block averages with N steps per block
  for (int i = 0; i<N; ++i) {
     ave[i] = 0.;
     sum_prog[i] = 0.;
     su2_prog[i] = 0.;
     err_prog[i] = 0.;
  }
  cout << "\nCalculating block averages with " << N << " steps per block...\n";

  for (int i = 0; i<N; ++i) {
    double sum = 0;
    for (int j = 0; j<L; ++j) {
      int k = j + i*L;
      sum += r[k];
    }
    ave[i] = sum/((double)L);
  }

  for  (int i = 0; i<N; ++i) {
    for (int j = 0; j<i+1; ++j) {
      sum_prog[i] += ave[j];
      su2_prog[i] += ave[j]*ave[j];
    }
    sum_prog[i] = sum_prog[i]/(i+1);
    su2_prog[i] = su2_prog[i]/(double)(i+1);
    err_prog[i] = error(sum_prog[i], su2_prog[i], i);
  }

  print(M, N, sum_prog, err_prog, "out_r2.txt", 5.);
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
