#ifndef __functions__
#define __functions__

#include "point.h"
#include "random.h"
#include <cmath>
// Moves from point p randomly and saves the result in pnew
void move (const point& p, point& pnew, Random& rnd, double rmax = 1.);

double min (double a, double b);
double prob (const point& p);
double prob2 (const point& p);
double A (const point& pnew, const point& p);
double A2 (const point& pnew, const point& p);

#endif //__functions__
