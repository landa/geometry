#ifndef _GEOM_CODE_H_
#define _GEOM_CODE_H_

#include <math.h>

// --- Forward Declarations ----------------------------------------------------

bool find_intersection(double x1, double y1, double x2, double y2,
                       double x3, double y3, double x4, double y4,
                       double& x_out, double& y_out);
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
inline double dot(const double* a, double* b);
inline double* cross(const double* a, double* b);
inline double distance_3d(const double* start, double* end);
inline double find_slope(double x1, double y1, double x2, double y2);
inline double find_intercept(double x, double y, double slope);

// -----------------------------------------------------------------------------

/**
 * Finds the intersection of two lines given two pairs of points
 *
 * Line 1 is specified by { (x1, y1), (x2, y2) }
 * Line 2 is specified by { (x3, y3), (x4, y4) }
 * The result is written into (x_out, y_out)
 *
 * Modified from:
 * http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry2
 */
bool find_intersection(double x1, double y1, double x2, double y2,
                       double x3, double y3, double x4, double y4,
                       double& x_out, double& y_out) {
  double A1 = y2 - y1;
  double A2 = y4 - y3;
  double B1 = x1 - x2;
  double B2 = x3 - x4;
  double C1 = A1 * x1 + B1 * y1;
  double C2 = A2 * x3 + B2 * y3;
  
  double det = A1 * B2 - A2 * B1;
  
  if (fabs(det) < 0.00001) return false;
  
  x_out = (B2 * C1 - B1 * C2)/det;
  y_out = (A1 * C2 - A2 * C1)/det;
  
  return true;
}

/**
 * Determines whether (testx, testy) is in the polygon specified by vertx, verty.
 * 
 * nvert is the number of vertices
 * *vertx is the array of x values
 * *verty is the array of y values
 *
 * Sourced from:
 * http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
 */
int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy) {
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
   (testx<(vertx[j]-vertx[i])*(testy-verty[i])/(verty[j]-verty[i])+vertx[i]) )
       c = !c;
  }
  return c;
}

/**
 * Computes the dot product of two 3-dimensional points
 * 
 * *start and *end are double triples
 */
inline double dot(const double *a, double *b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/**
 * Computes the cross product of two 3-dimensional points
 * 
 * *start and *end are double triples and a double triplet is returned
 */
inline double* cross(const double *a, double *b) {
  double result[3];
  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[2] - a[1]*b[0];
}

/**
 * Computes the Euclidean distance of two 3-dimensional points
 * 
 * *start and *end are double triples
 */
inline double distance_3d(const double *start, double *end) {
  return pow(
    (start[0]-end[0])*(start[0]-end[0])+
    (start[1]-end[1])*(start[1]-end[1])+
    (start[2]-end[2])*(start[2]-end[2]),
    0.5);
}

/**
 * Computes the slope of a line, provided two points (x1,y1) and (x2,y2)
 */
inline double find_slope(double x1, double y1, double x2, double y2) {
  return (y2-y1)/(x2-x1);
}

/**
 * Computes the y-intercept of a line, provided two points (x1,y1) and (x2,y2)
 */
inline double find_intercept(double x, double y, double slope) {
  return y-slope*x;
}

#endif
