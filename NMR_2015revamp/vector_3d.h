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
      float x_coord,y_coord,z_coord;
	public:
      float get_x();
      float get_y();
      float get_z();
      void set_x(float);
      void set_y(float);
      void set_z(float);
      float get_dot_product(vector_3d *);
      float get_magnitude();
      float get_angle_from(vector_3d *);
      vector_3d(float, float,float);
      vector_3d();
      vector_3d(coords_3d,coords_3d);
};




#endif /* 3DVECTOR_H_ */
