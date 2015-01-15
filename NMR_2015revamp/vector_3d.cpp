/*
 * vector_3d.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */
#include "vector_3d.h"

double vector_3d::get_x(){
  return x_coord;
}
double vector_3d::get_y(){
  return y_coord;
}
double vector_3d::get_z(){
  return z_coord;
}
void vector_3d::set_x(double input){
  x_coord=input;
}
void vector_3d::set_y(double input){
  y_coord=input;
}
void vector_3d::set_z(double input){
  z_coord=input;
}

double vector_3d::get_dot_product(vector_3d* input){
  //Dot product is = x1*x2 + y1*y2 + z1*z2
  return (get_x() * input->get_x()) + get_y() * input->get_y() + get_z() * input->get_z();
}

double vector_3d::get_magnitude(){
  return sqrt(pow(x_coord,2)+pow(y_coord,2) + pow(z_coord,2));
}

double vector_3d::get_angle_from(vector_3d* input){
//Get angle from a given vector, using the dot product
//Given in degrees
  return acos(get_dot_product(input)/(get_magnitude()*input->get_magnitude()));
}

vector_3d::vector_3d(double x, double y, double z){
  x_coord=x;
  y_coord=y;
  z_coord=z;
}
vector_3d::vector_3d(){
  x_coord=0.01;
  y_coord=0.01;
  z_coord=0.01;
}
vector_3d::vector_3d(coords_3d source, coords_3d dest){
  double temp_x;
  double temp_y;
  double temp_z;
  double magnitude;
  temp_x=dest.getX() - source.getX();
  temp_y=dest.getY() - source.getY();
  temp_z=dest.getZ() - source.getZ();
  magnitude= sqrt(pow(temp_x,2) + pow(temp_y,2) + pow(temp_z,2));
  temp_x=temp_x*(1.0/magnitude);
  temp_y=temp_y*(1.0/magnitude);
  temp_z=temp_z*(1.0/magnitude);
  x_coord=temp_x;
  y_coord=temp_y;
  z_coord=temp_z;
  //TODO: finish constructor to give a unit vector output
}
