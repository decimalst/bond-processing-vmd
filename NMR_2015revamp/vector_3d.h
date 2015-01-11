/*
 * 3dvector.h
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */

#ifndef VECTOR_3D_H_
#define VECTOR_3D_H_
#include "coords_3d.h"
#include <cmath>

class vector_3d {
  double x_coord, y_coord, z_coord;
 public:
  double get_x();
  double get_y();
  double get_z();
  void set_x(double);
  void set_y(double);
  void set_z(double);
  double get_dot_product(vector_3d *);
  double get_magnitude();
  double get_angle_from(vector_3d *);
  vector_3d(double, double, double);
  vector_3d();
  vector_3d(coords_3d, coords_3d);
};

#endif /* 3DVECTOR_H_ */
