/*
 * simulation.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */
#include "simulation.h"
std::vector<molecule_bond>* simulation::get_molecules() {
  return molecules;
}
float simulation::get_time_value() {
  return frame_time_value;
}
int simulation::get_number_of_frames() {
  return number_of_frames;
}
int simulation::get_number_of_molecules() {
  return number_of_molecules;
}
void simulation::set_time_value(float set_to) {
  frame_time_value = set_to;
}
void simulation::set_number_of_frames(int set_to) {
  number_of_frames = set_to;
}
simulation::simulation(std::vector<molecule_bond>* molecules_in) {
  molecules = molecules_in;
  number_of_molecules = molecules_in->size();
  frame_time_value = 0.0;
  number_of_frames = molecules_in->at(
      1).getMoleculeTrajectory()->get_trajectory_array()->size();
}
simulation::simulation() {
  number_of_molecules = 0;
  number_of_frames = 0;
  frame_time_value = 0;
  molecules = new std::vector<molecule_bond>();
}
simulation::~simulation() {
  delete molecules;
}
simulation::simulation(std::vector<molecule_bond>* molecules_in,
                       float time_value, int molecules_number, int frames) {
  molecules = molecules_in;
  frame_time_value = time_value;
  number_of_molecules = molecules_number;
  number_of_frames = frames;
}
void simulation::calculate_deuterium_order_parameter() {
//Calculates deuterium order parameter and writes it to a text file
//Deuterium order parameter is defined as the average angle the bond makes
//with the Z-axis,

//Assumes that the z axis is in the upwards z direction: this may need to be
//Changed depending on the simulation orientation.
  vector_3d zaxis = vector_3d(
      0.0, 0.0, 1.0);

  float average_over_ensemble = 0.0;
  for (int i = 0; i < 152; i++) {
    std::cout << i << " " << molecules->at(
        i).getMoleculeTrajectory()->get_trajectory_array()->size()
              << std::endl;
  }
//Declare vector of size equal to the number of bonds in the simulation
//In this vector we will store the angle for each bond with z axis
//Once the vector has the angle of the entire simulation added, we
//average using the number of frames in the simulation.
  std::vector<float> angleholder(
      number_of_molecules);
  std::cout << angleholder.size() << " " << molecules->size() << std::endl;
  for (int i = 0; i < number_of_molecules; i++) {
    std::cout << molecules->at(
        i).getMoleculeTrajectory()->get_trajectory_array()->size()
              << std::endl;
    std::cout << i << std::endl;
    for (int j = 0; j < number_of_frames; j++) {
      //std::cout << "i " << i << " j " << j << " num_of_frames " << number_of_frames << std::endl;
      angleholder[i] += molecules->at(
          i).getMoleculeTrajectory()->get_trajectory_array()->at(
          j).get_angle_from(
          &zaxis);
    }
    angleholder[i] = angleholder[i] / number_of_frames;
    average_over_ensemble += angleholder[i];
  }
  average_over_ensemble = average_over_ensemble / number_of_molecules;
  std::cout << average_over_ensemble << std::endl;

}
void simulation::calculate_first_legendre() {
//Calculates first legendre polynomial and writes it to a text file

}
void simulation::calculate_second_legendre() {
//Calculates second legendre polynomial and writes it to a text file

}
void simulation::calculate_first_and_second_legendre() {
//Calculates both first and second legendre polynomials and write them both to a textfile
//Given "m" molecular bonds, and "n" frames


  //Step 1: For some molecular bond in simulation
  for (int i = 0; i < number_of_molecules; i++) {
    //Step 2: For all values of t-"spacing between frames" from t=1 to t=(n-1)
    for (int j = 0; j < number_of_frames; j++) {
      //Step 3: For all pairs of frames in simulation with spacing t.
      for (int k = 0; k < (number_of_frames - j); k++) {
        //Compute x(0)x(t), y(0)y(t), etc. for all involved pairs for C1, C2
        //this will be x[k]*x[k+j], y[k]y[k+j], etc.

        //Store in float sums

      }
      //Average float sums of products by (1/number_of_frames-spacing(j) )
      //Store into results_array(t+=) as a result associated with the value of spacing
      //Set float sum values back to zero.
    }
  }
}
