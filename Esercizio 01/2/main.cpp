#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

double sums(const vector<double>& r);

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

   // ----------------------------------------- Exercise 01.2 -----------------------------------------
   int N = 10000;

   for (int j = 0; j<4; ++j) {                    // cicle for the 4 values of n = 1, 2, 10, 100
     vector<double> s_exp, s_cl, s;               // N = 10^4 realizations of S_n
     int n;
     if (j==0) n=1;
     else if (j==1) n=2;
     else if (j==3) n=10;
     else n=100;

     for (int i = 0; i< N; ++i) {                 // cicle on the N realizations
       vector<double> r_exp(n), r_cl(n), r(n);    // Vectors for exponential, Cauchy-Lorentz and uniformly distributed random numbers

       for (int k = 0; k<n; ++k) {                // Generateing the n random numbers
         r_exp[k] = rnd.Exponential();            // Exponential and CauchyLorentz implemented in random.cpp
         r_cl[k] = rnd.CauchyLorentz();
         r[k] = rnd.Rannyu();
       }

       // Appending the values of S_n for the i-th iteration in the vectors s
       s_exp.push_back(sums(r_exp)/n);
       s_cl.push_back(sums(r_cl)/n);
       s.push_back(sums(r)/n);
     }

     // Generating 4 files 'outj.txt', j = 1,...,4 containing N, n, s_exp, s_cl, s
     string filename = "file" + to_string(j+1) + ".txt";
     ofstream out(filename);
     out << N << endl << n << endl;
     for (auto i = s_exp.begin(); i!=s_exp.end(); ++i)
      out << *i << endl;
    for (auto i = s_cl.begin(); i!=s_cl.end(); ++i)
      out << *i << endl;
    for (auto i = s.begin(); i!=s.end(); ++i)
      out << *i << endl;
   }

   rnd.SaveSeed();
   return 0;
}


double sums(const vector<double>& r) {
  double sum = 0.;
  for (auto i = r.begin(); i!=r.end(); ++i)
    sum += *i;
  return sum;
}
