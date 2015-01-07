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
  float x,y,z;
 public:
  float getX();
  void setX(float x);
  float getY();
  void setY(float y);
  float getZ();
  void setZ(float z);
  coords_3d(float,float,float);
  coords_3d();
};

#endif /* COORDS_3D_H_ */
