#include "functions.h"

void move (const point& p, point& pnew, Random& rnd, double rmax) {
  double x, y, z;
  double r, theta, phi;
  r = rnd.Rannyu(0, rmax);
  theta = rnd.Theta();
  phi = rnd.Phi();
  x = r*cos(phi)*sin(theta);
  y = r*sin(phi)*sin(theta);
  z = r*cos(theta);
  pnew.setcoord(p.getx()+x, p.gety()+y, p.getz()+z);
}

double prob (const point& p) {
  return exp(-2.*p.getr())/M_PI;
}

double prob2 (const point& p) {
  double r = p.getr();
  return exp(-1.*r)*cos(p.gettheta())*r*r*2./(8.*8.*M_PI);
}

double min (double a, double b) {return a<b? a : b;}

double A (const point& pnew, const point& p) {return min(1., prob(pnew)/prob(p));}

double A2 (const point& pnew, const point& p) {return min(1., prob2(pnew)/prob2(p));}
