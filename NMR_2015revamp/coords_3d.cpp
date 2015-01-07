/*
 * coords_3d.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */
#include "coords_3d.h"

float coords_3d::getX() {
  return x;
}
void coords_3d::setX(float x_given) {
  x = x_given;
}
float coords_3d::getY() {
  return y;
}
void coords_3d::setY(float y_param) {
  y = y_param;
}
float coords_3d::getZ() {
  return z;
}
void coords_3d::setZ(float z_param) {
  z = z_param;
}
coords_3d::coords_3d(float x_given, float y_given, float z_given){
  x=x_given;
  y=y_given;
  z=z_given;
}
coords_3d::coords_3d(){
  x=0.0;
  y=0.0;
  z=0.0;
}

