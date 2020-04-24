#ifndef __point__
#define __point__

#include <cmath>
#include <string>
#include <iostream>

using namespace std;

class point{
private:
  double x, y, z;
  double r, theta, phi;

  void settheta ();
  void setphi ();
  void setr ();
public:
  point ();
  point (double X, double Y = 0., double Z =0.);

  double getx () const {return x;};
  double gety () const {return y;};
  double getz () const {return z;};
  double getr () const {return r;};
  double gettheta () const {return theta;};
  double getphi () const {return phi;};

  void setx (double X);
  void sety (double Y);
  void setz (double Z);
  void setcoord (double X, double Y = 0., double Z = 0.);
  void setcoord (const point& p) {setcoord(p.getx(), p.gety(), p.getz());}

  friend ostream& operator << (ostream& os, const point& x);
};

#endif //__point__
