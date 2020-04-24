#include "point.h"
#include "random.h"
#include "functions.h"
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

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


  // Calculating optimal values of r_max

  // Testing which value of rmax gives <A(x|y)> = 50%
  ofstream out("rmax.txt");
  int Nsteps = 10000;
  int Nwalks = 50;
  cout << "n = 1, l = 0, m = 0:" << endl;
  for (double rmax = 1.; rmax <= 3.; rmax+=0.05) {
    point p(1,1,1), pnew;
    double avgA = 0.;
    cout << "Calculating <A(x|y)> for rmax = " << rmax << " with " << Nsteps << " steps and " << Nwalks << " walks:";
    for (int j = 0; j < Nwalks; ++j) {
      double avgAj = 0.;
      for (int i = 0; i < Nsteps; ++i) {
        move(p, pnew, rnd, rmax);
        double alpha = A(pnew, p);
        avgAj += alpha;
        if (rnd.Rannyu() <= alpha) {
          p.setcoord(pnew);
        }
      }
      avgAj /= (double)Nsteps;
      avgA += avgAj;
    }
    avgA /= Nwalks;
    cout << " <A(x|y)> = " << avgA << endl;
    out << rmax << endl << avgA << endl;
  }
  out.close();

  cout << endl << "Closeup on the region with <A(x|y)> = 0.5:" << endl;
  out.open("rmax_closeup.txt");
  for (double rmax = 2.25; rmax <= 2.75; rmax+=0.01) {
    point p, pnew;
    double avgA = 0.;
    cout << "Calculating <A(x|y)> for rmax = " << rmax << " with " << Nsteps << " steps and " << Nwalks << " walks:";
    for (int j = 0; j < Nwalks; ++j) {
      double avgAj = 0.;
      for (int i = 0; i < Nsteps; ++i) {
        move(p, pnew, rnd, rmax);
        double alpha = A(pnew, p);
        avgAj += alpha;
        if (rnd.Rannyu() <= alpha) {
          p.setcoord(pnew);
        }
      }
      avgAj /= (double)Nsteps;
      avgA += avgAj;
    }
    avgA /= Nwalks;
    cout << " <A(x|y)> = " << avgA << endl;
    out << rmax << endl << avgA << endl;
  }
  out.close();


  cout << endl << endl << "n = 2, l = 1, m = 0:" << endl;
  out.open("rmax2.txt");
  for (double rmax = 3.; rmax <= 7.; rmax+=0.1) {
    point p(1,1,5), pnew;
    double avgA = 0.;
    cout << "Calculating <A(x|y)> for rmax = " << rmax << " with " << Nsteps << " steps and " << Nwalks << " walks:";
    for (int j = 0; j < Nwalks; ++j) {
      double avgAj = 0.;
      for (int i = 0; i < Nsteps; ++i) {
        move(p, pnew, rnd, rmax);
        double alpha = A2(pnew, p);
        avgAj += alpha;
        if (rnd.Rannyu() <= alpha) {
          p.setcoord(pnew);
        }
      }
      avgAj /= (double)Nsteps;
      avgA += avgAj;
    }
    avgA /= Nwalks;
    cout << " <A(x|y)> = " << avgA << endl;
    if (avgA < 10 && avgA > 0) out << rmax << endl << avgA << endl;
    else out << rmax << endl << 0. << endl;
  }
  out.close();

  cout << endl << "Closeup on the region with <A(x|y)> = 0.5:" << endl;
  out.open("rmax_closeup2.txt");
  for (double rmax = 4.5; rmax <= 5; rmax+=0.01) {
    point p(1,1,5), pnew;
    double avgA = 0.;
    cout << "Calculating <A(x|y)> for rmax = " << rmax << " with " << Nsteps << " steps and " << Nwalks << " walks:";
    for (int j = 0; j < Nwalks; ++j) {
      double avgAj = 0.;
      for (int i = 0; i < Nsteps; ++i) {
        move(p, pnew, rnd, rmax);
        double alpha = A2(pnew, p);
        avgAj += alpha;
        if (rnd.Rannyu() <= alpha) {
          p.setcoord(pnew);
        }
      }
      avgAj /= (double)Nsteps;
      avgA += avgAj;
    }
    avgA /= Nwalks;
    cout << " <A(x|y)> = " << avgA << endl;
    if (avgA < 10  && avgA > 0) out << rmax << endl << avgA << endl;
    else out << rmax << endl << 0. << endl;
  }
  out.close();

  return 0;
}
