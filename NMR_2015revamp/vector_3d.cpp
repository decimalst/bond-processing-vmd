/*
 * vector_3d.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */
#include "vector_3d.h"

float vector_3d::get_x(){
  return x_coord;
}
float vector_3d::get_y(){
  return y_coord;
}
float vector_3d::get_z(){
  return z_coord;
}
void vector_3d::set_x(float input){
  x_coord=input;
}
void vector_3d::set_y(float input){
  y_coord=input;
}
void vector_3d::set_z(float input){
  z_coord=input;
}

float vector_3d::get_dot_product(vector_3d* input){
  //Dot product is = x1*x2 + y1*y2 + z1*z2
  return (get_x() * input->get_x()) + get_y() * input->get_y() + get_z() * input->get_z();
}

float vector_3d::get_magnitude(){
  return sqrt(pow(x_coord,2)+pow(y_coord,2) + pow(z_coord,2));
}

float vector_3d::get_angle_from(vector_3d* input){
//Get angle from a given vector, using the dot product
//Given in degrees
  return acos(get_dot_product(input)/(get_magnitude()*input->get_magnitude()));
}

vector_3d::vector_3d(float x, float y, float z){
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
  float temp_x;
  float temp_y;
  float temp_z;
  float magnitude;
  temp_x=fabs(source.getX() - dest.getX());
  temp_y=fabs(source.getY() - dest.getY());
  temp_z=fabs(source.getZ() - dest.getZ());
  magnitude= sqrt(pow(temp_x,2) + pow(temp_y,2) + pow(temp_z,2));
  temp_x=temp_x*(1/magnitude);
  temp_y=temp_y*(1/magnitude);
  temp_z=temp_z*(1/magnitude);
  x_coord=temp_x;
  y_coord=temp_y;
  z_coord=temp_z;
  //TODO: finish constructor to give a unit vector output
}
