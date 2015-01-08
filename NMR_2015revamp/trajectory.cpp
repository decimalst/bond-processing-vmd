/*
 * trajectory.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */
#include "trajectory.h"
vector_3d trajectory::get_Frame(int frame){
  if(frame<number_of_frames && frame>=0){
    return trajectory_array[frame];
  }
  else{
    std::cout << "Error, frame number given not legal, returning zero vector_3d" << std::endl;
    return vector_3d();
  }
}
void trajectory::set_number_of_frames(int input){
  number_of_frames=input;
}
trajectory::trajectory(int numb_frames, std::vector<vector_3d> data){
   number_of_frames=numb_frames;
    trajectory_array=data;
}
std::vector<vector_3d>* trajectory::get_trajectory_array(){
  return &trajectory_array;
}
trajectory::trajectory(){
  number_of_frames=0;
  trajectory_array=std::vector<vector_3d>();
}


