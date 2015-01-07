/*
 * simulation.h
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_
#include "molecule.h"
#include "trajectory.h"
class simulation{
  std::vector<molecule_bond>* molecules;
  float frame_time_value;

  int number_of_molecules, number_of_frames;
 public:
  std::vector<molecule_bond> * get_molecules();
  float get_time_value();
  void set_time_value(float);
  int get_number_of_frames();
  void set_number_of_frames(int);
  int get_number_of_molecules();
  void calculate_deuterium_order_parameter();
  void calculate_first_legendre();
  void calculate_second_legendre();
  void calculate_first_and_second_legendre();
  simulation();
  ~simulation();
  simulation(std::vector<molecule_bond>*);
  simulation(std::vector<molecule_bond>*,float,int,int);
};




#endif /* SIMULATION_H_ */
