/*
 * coords_3d.h
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */

#ifndef COORDS_3D_H_
#define COORDS_3D_H_
#include <string>
class coords_3d{
  double x,y,z;
 public:
  double getX();
  void setX(double x);
  double getY();
  void setY(double y);
  double getZ();
  void setZ(double z);
  coords_3d(double,double,double);
  coords_3d();
};

#endif /* COORDS_3D_H_ */
