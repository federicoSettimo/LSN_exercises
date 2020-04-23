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
void print (int N, const vector<double>& r2, const vector<double>& err, const string outfile = "out.txt", double shift = 0.);

template<typename T>
double norm(const vector<T>& x);

int sgn(double x) {return x>=0 ? 1 : -1;}

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

  // ----------------------------------------- Exercise 02.2 -----------------------------------------
  int N = 100;            // number of steps per walk
  int M = 10000;          // number of walk
  vector<int> pos(3);     // x, y and z coordinate of the RW
  vector<double> r2(N);   // average value of r2 for each step
  vector<double> er(N);   // average of (r^2)^2 for each step

  for (int i = 0; i<M; ++i) {      // M random walks
    pos[0] = 0;
    pos[1] = 0;
    pos[2] = 0;
    for (int j = 0; j<N; ++j) {    // each of N steps
      pos[(int)rnd.Rannyu(0,3)] += sgn(rnd.Rannyu(-1.,1.)); // 1 passo in direzione e verso random di lunghezza 1
      double r_i =norm(pos);
      r2[j] += r_i/M;         // avg of all r^2 at step j
      er[j] += r_i*r_i/(M);       // avg of all (r^2)^2 at step j
    }
  }

  // <r^2> --> sqrt( <r^2> ) and error on r from the error on r^2
  for (int i = 0; i<N; ++i) {
    r2[i] = sqrt(r2[i]);
    if (r2[i] != 0) {
      er[i] = error(r2[i]*r2[i], er[i])/(2.*r2[i]);  // now r2 actually contains r and er contains the error
    }
    else er[i] = 0.;
  }

  // Printing avgs
  print(N, r2, er, "out1.txt");

  // Part 2
  // init variables to 0 and pos as a double
  r2.assign(N, 0.);
  er.assign(N, 0.);
  pos.clear();
  vector<double> pos2(N);

  // Like before but random numbers as theta and phi. pos still in cartesian coordinates
  for (int i = 0; i<M; ++i) {
    pos2[0] = 0;
    pos2[1] = 0;
    pos2[2] = 0;
    for (int j = 0; j<N; ++j) {
      double theta = 2 * M_PI * rnd.Rannyu();
      double phi = acos(1 - 2 * rnd.Rannyu());
      pos2[0] += sin(phi) * cos(theta);
      pos2[1] += sin(phi) * sin(theta);
      pos2[2] += cos(phi);
      double r_i =norm(pos2);
      r2[j] += r_i/M;
      er[j] += r_i*r_i/(M);
    }
  }

  // Like before
  for (int i = 0; i<N; ++i) {
    r2[i] = sqrt(r2[i]);
    if (r2[i] != 0) {
      er[i] = error(r2[i]*r2[i], er[i])/(2.*r2[i]);
    }
    else er[i] = 0.;
  }

  // Printing avgs
  print(N, r2, er, "out2.txt");

  rnd.SaveSeed();
  return 0;
}

double error (double AV, double AV2, int n) {
   return n == 0 ? 0. : sqrt((AV2 - AV*AV)/(double)n);
}

void print (int N, const vector<double>& sum, const vector<double>& err, const string outfile, double shift) {
   ofstream out(outfile);
   out << N << endl;
   for (auto i = sum.begin(); i != sum.end(); ++i)
      out << *i-shift << endl;
   for (auto i = err.begin(); i != err.end(); ++i)
      out << *i << endl;
   out.close();
}

template<typename T>
double norm(const vector<T>& x) {
  double n = 0.;
  for (auto i = x.begin(); i != x.end(); ++i)
    n+=(double)(*i)*(*i);
  return n;
}
