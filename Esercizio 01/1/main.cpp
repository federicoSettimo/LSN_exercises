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

// For part 3: prints avg_chi, err_chi, N, chi on the file 'outfile' ("out3.txt")
void print (double avg_chi, double err_chi, int N, const vector<double>& chi, string outfile = "out3.txt");

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

   // Initializing r with M random numbers
   for (int i = 0; i<M; ++i)
      r[i] = rnd.Rannyu();

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
  print(M, N, sum_prog, err_prog, "out1.txt", 0.5);


   // Part 2 of Exercise 01.1
   // Initializing all to 0
   for (int i = 0; i<N; ++i) {
      ave[i] = 0.;
      sum_prog[i] = 0.;
      su2_prog[i] = 0.;
      err_prog[i] = 0.;
   }

   // Calculating averages
   for (int i = 0; i<N; ++i) {
      double sum = 0;
      for (int j = 0; j<L; ++j) {
         int k = j + i*L;
         sum += (r[k] - 0.5)*(r[k] - 0.5);
      }
      ave[i] = sum/(double)L;
   }

   // Calculating sum_prog and su2_prog
   for  (int i = 0; i<N; ++i) {
      for (int j = 0; j<i+1; ++j) {
         sum_prog[i] += ave[j];
         su2_prog[i] += ave[j]*ave[j];
      }
      sum_prog[i] /= (double)(i+1);
      su2_prog[i] /= (double)(i+1);
      err_prog[i] = error(sum_prog[i], su2_prog[i], i);
   }

   // Printing on "out1.txt"
   print(M, N, sum_prog, err_prog, "out2.txt", 1./12.);

   //Deleting allocated vectors
   ave.clear();
   sum_prog.clear();
   su2_prog.clear();
   err_prog.clear();


   // Part 3 of Exercise 01.1
   int n = M/N;
   vector<int> hits(N);
   double sum_chi = 0, avg_chi = 0, su2_chi = 0, err_chi = 0;
   vector<double> chi(N);

   for (int j = 0; j<N; ++j) {
      for (int i = 0; i<n; ++i) {
         // If the considered random number is in the correct intervall (i-th), one hit is added to the i-th interval
         if (r[j*n + i] >= (double)i/N && r[j*n + i] < (double)(i+1)/N)
            hits[j] += 1;
      }
      // Calculates chi for the j-th intervall and adds it to the total chi
      chi[j] = (hits[j] - (double)n/N)*(hits[j] - (double)n/N)/((double)n/N);
      sum_chi += chi[j];
      su2_chi += chi[j]*chi[j];
   }
   // Average and error of chi
   avg_chi = sum_chi/N;
   err_chi = error(avg_chi, su2_chi/N);

   // Prints chi
   print(avg_chi, err_chi, N, chi);

   //Deleting allocated vectors
   r.clear();
   chi.clear();

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

void print (double avg_chi, double err_chi, int N, const vector<double>& chi, string outfile) {
   ofstream out(outfile);
   out << avg_chi << endl << err_chi << endl << N << endl;
   for (auto i = chi.begin(); i != chi.end(); ++i)
      out << *i << endl;
   out.close();
}
