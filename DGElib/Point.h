/*
 *
 * Point.h
 *
 */

#ifndef POINT_H
#define POINT_H

class Point : public SDGData
{
  double x,y;
 
 public:

  Point(double i, double j): x(i), y(j) {};
  double getX() {return x;};
  double getY() {return y;};
  
};
#endif
