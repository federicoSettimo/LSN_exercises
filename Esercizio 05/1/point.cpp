#include "point.h"

point::point () {
  x = 0.;
  y = 0.;
  z = 0.;
  theta = 0.;
  phi = 0.;
  r = 0.;
}

point::point (double X, double Y, double Z) {
  setcoord (X,Y,Z);
}

void point::setx (double X) {
  x = X;
  setr();
  settheta();
  setphi();
  return;
}

void point::sety (double Y) {
  y = Y;
  setr();
  settheta();
  setphi();
  return;
}

void point::setz (double Z) {
  z = Z;
  setr();
  settheta();
  setphi();
  return;
}

void point::setcoord (double X, double Y, double Z) {
  x = X;
  y = Y;
  z = Z;
  setr();
  settheta();
  setphi();
}

void point::setr () {
  r = sqrt(x*x + y*y + z*z);
}

void point::setphi () {
  if (r == 0.) phi = 0.;
  else phi = acos(z/r);
}

void point::settheta () {
  if ( x == 0 && y == 0) theta = 0.;
  else theta = atan2(x,y);
}

ostream& operator << (ostream& os, const point& x) {
  os << x.getx() << endl << x.gety() << endl << x.getz() << endl;
  return os;
}
