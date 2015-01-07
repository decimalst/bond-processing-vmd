/*
 * trajectory.h
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */

#ifndef TRAJECTORY_H_
#define TRAJECTORY_H_

#include <vector>
#include <iostream>
#include "vector_3d.h"

class trajectory{
  int number_of_frames;
  std::vector<vector_3d> trajectory_array;
 public:
  vector_3d get_Frame(int);
  std::vector<vector_3d>* get_trajectory_array();
  trajectory(int,std::vector<vector_3d>);
  trajectory();
};



#endif /* TRAJECTORY_H_ */
