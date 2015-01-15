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
  number_of_frames = molecules_in->at(1).getMoleculeTrajectory()
      ->get_trajectory_array()->size();
  for (int i = 0; i < molecules_in->size(); i++) {
    molecules->at(i).getMoleculeTrajectory()->set_number_of_frames(
        number_of_frames);
  }
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
  vector_3d zaxis = vector_3d(0.0, 0.0, 1.0);

  float average_over_ensemble = 0.0;
//Declare vector of size equal to the number of bonds in the simulation
//In this vector we will store the angle for each bond with z axis
//Once the vector has the angle of the entire simulation added, we
//average using the number of frames in the simulation.
  std::vector<float> angleholder(number_of_molecules);
  std::cout << angleholder.size() << " " << molecules->size() << std::endl;
  for (int i = 0; i < number_of_molecules; i++) {
    //std::cout << i << std::endl;
    for (int j = 0; j < number_of_frames; j++) {
      //std::cout << "i " << i << " j " << j << " num_of_frames " << number_of_frames << std::endl;
      angleholder[i] += fabs(
          ((3 / 2)
              * pow(
                  cos(molecules->at(i).getMoleculeTrajectory()
                      ->get_trajectory_array()->at(j).get_angle_from(&zaxis)),
                  2)) - 0.5);
    }
    angleholder[i] = angleholder[i] / (float) number_of_frames;
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
//First, declare float values to use as sums
//First legendre sums
  int number_of_frames = get_number_of_frames();
  int bonds_per_carbon = get_molecules()->at(0).getNumberOfBonds();
  std::cout << "bonds per carbon " << bonds_per_carbon << std::endl;
  int count = 0;

  double x0xt_sum = 0.0, y0yt_sum = 0.0, z0zt_sum = 0.0;
  //Second legendre sums for left
  double x0xt_squared_sum = 0.0, y0yt_squared_sum = 0.0, z0zt_squared_sum = 0.0;
  double x0y0xtyt_sum = 0.0, x0z0xtzt_sum = 0.0, y0z0ytzt_sum = 0.0;
  double p1 = 0.0, p2 = 0.0;

  //String operations to make output filename
  std::string filename;
  filename = molecules->at(0).getName() + "_"
      + molecules->at(0).getMoleculeType() + "_legendre";

  std::ofstream write_to_file(filename);
  write_to_file << "# t" << " " << "<x(0)x(t)>" << " " << "<y(0)y(t)>" << " "
                << "<z(0)z(t)>" << " " << "<x^2(0)x^2(t)>" << " "
                << "<y^2(0)y^2(t)>" << " " << "<z^2(0)z^2(t)>" << " "
                << "<x(0)y(0)x(t)y(t)>" << " " << "<x(0)z(0)x(t)z(t)>" << " "
                << "<y(0)z(0)y(t)z(t)>" << " " << "C1" << " " << "C2"
                << std::endl;

  //Then, declare arrays for the output results vs t
  std::vector<double> first_x1(number_of_molecules), first_y1(
      number_of_molecules), first_z1(number_of_molecules);
  std::vector<double> second_xsquare1(number_of_molecules), second_ysquare1(
      number_of_molecules), second_zsquare1(number_of_molecules);
  std::vector<double> second_xy1(number_of_molecules), second_xz1(
      number_of_molecules), second_yz1(number_of_molecules);
  std::vector<int> spacing_array;

  std::vector<double> first_x2(number_of_molecules), first_y2(
      number_of_molecules), first_z2(number_of_molecules);
  std::vector<double> second_xsquare2(number_of_molecules), second_ysquare2(
      number_of_molecules), second_zsquare2(number_of_molecules);
  std::vector<double> second_xy2(number_of_molecules), second_xz2(
      number_of_molecules), second_yz2(number_of_molecules);

  std::vector<double> first_x3(number_of_molecules), first_y3(
      number_of_molecules), first_z3(number_of_molecules);
  std::vector<double> second_xsquare3(number_of_molecules), second_ysquare3(
      number_of_molecules), second_zsquare3(number_of_molecules);
  std::vector<double> second_xy3(number_of_molecules), second_xz3(
      number_of_molecules), second_yz3(number_of_molecules);

  //Temp vectors to hold vectors we select
  vector_3d vector1, vector2;

  for (int spacing = 0; spacing < number_of_frames; spacing++) {
    x0xt_sum = 0.0;
    y0yt_sum = 0.0;
    z0zt_sum = 0.0;
    x0xt_squared_sum = 0.0;
    y0yt_squared_sum = 0.0;
    z0zt_squared_sum = 0.0;
    x0y0xtyt_sum = 0.0;
    x0z0xtzt_sum = 0.0;
    y0z0ytzt_sum = 0.0;
    std::cout << "spacing= " << spacing << " out of " << number_of_frames
              << std::endl;
    for (int i = 0; i < number_of_molecules; i = i + bonds_per_carbon) {
//      std::cout << "i = " << i << " i+1= "<< i+1 << std::endl;
      count = 0;
//      vector1 =
//          molecules->at(i).getMoleculeTrajectory()->get_trajectory_array()->at(
//              0);
//      vector2 =
//          molecules->at(i+1).getMoleculeTrajectory()->get_trajectory_array()->at(
//              0);
//      std::cout << "x " << vector1.get_x() << " y " << vector1.get_y() << " z " << vector1.get_z() << std::endl;
//      std::cout << "x " << vector2.get_x() << " y " << vector2.get_y() << " z " << vector2.get_z() << std::endl;
      for (int i = 0; i < number_of_molecules; i++) {
        first_x1[i] = 0.0;
        first_y1[i] = 0.0;
        first_z1[i] = 0.0;
        second_xsquare1[i] = 0.0;
        second_ysquare1[i] = 0.0;
        second_zsquare1[i] = 0.0;
        second_xy1[i] = 0.0;
        second_xz1[i] = 0.0;
        second_yz1[i] = 0.0;
        first_x2[i] = 0.0;
        first_y2[i] = 0.0;
        first_z2[i] = 0.0;
        second_xsquare2[i] = 0.0;
        second_ysquare2[i] = 0.0;
        second_zsquare2[i] = 0.0;
        second_xy2[i] = 0.0;
        second_xz2[i] = 0.0;
        second_yz2[i] = 0.0;
        first_x3[i] = 0.0;
        first_y3[i] = 0.0;
        first_z3[i] = 0.0;
        second_xsquare3[i] = 0.0;
        second_ysquare3[i] = 0.0;
        second_zsquare3[i] = 0.0;
        second_xy3[i] = 0.0;
        second_xz3[i] = 0.0;
        second_yz3[i] = 0.0;
      }

      for (int j = 0; j < number_of_frames - spacing; j++) {
        count++;

        if (bonds_per_carbon == 1) {
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x1[i] += vector1.get_x() * vector2.get_x();
          first_y1[i] += vector1.get_y() * vector2.get_y();
          first_z1[i] += vector1.get_z() * vector2.get_z();
          second_xsquare1[i] += vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare1[i] += vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare1[i] += vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy1[i] += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz1[i] += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz1[i] += vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
        }
        if (bonds_per_carbon == 2) {
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x1[i] += vector1.get_x() * vector2.get_x();
          first_y1[i] += vector1.get_y() * vector2.get_y();
          first_z1[i] += vector1.get_z() * vector2.get_z();
          second_xsquare1[i] += vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare1[i] += vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare1[i] += vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy1[i] += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz1[i] += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz1[i] += vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
          vector1 = molecules->at(i + 1).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i + 1).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x2[i] += vector1.get_x() * vector2.get_x();
          first_y2[i] += vector1.get_y() * vector2.get_y();
          first_z2[i] += vector1.get_z() * vector2.get_z();
          second_xsquare2[i] += vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare2[i] += vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare2[i] += vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy2[i] += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz2[i] += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz2[i] += vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
        }
        if (bonds_per_carbon == 3) {
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x1[i] = vector1.get_x() * vector2.get_x();
          first_y1[i] = vector1.get_y() * vector2.get_y();
          first_z1[i] = vector1.get_z() * vector2.get_z();
          second_xsquare1[i] = vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare1[i] = vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare1[i] = vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy1[i] = vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz1[i] = vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz1[i] = vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
          vector1 = molecules->at(i + 1).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i + 1).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x2[i] = vector1.get_x() * vector2.get_x();
          first_y2[i] = vector1.get_y() * vector2.get_y();
          first_z2[i] = vector1.get_z() * vector2.get_z();
          second_xsquare2[i] = vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare2[i] = vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare2[i] = vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy2[i] = vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz2[i] = vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz2[i] = vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
          vector1 = molecules->at(i + 2).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j);
          vector2 = molecules->at(i + 2).getMoleculeTrajectory()
              ->get_trajectory_array()->at(j + spacing);
          first_x3[i] = vector1.get_x() * vector2.get_x();
          first_y3[i] = vector1.get_y() * vector2.get_y();
          first_z3[i] = vector1.get_z() * vector2.get_z();
          second_xsquare3[i] = vector1.get_x() * vector2.get_x()
              * vector1.get_x() * vector2.get_x();
          second_ysquare3[i] = vector1.get_y() * vector2.get_y()
              * vector1.get_y() * vector2.get_y();
          second_zsquare3[i] = vector1.get_z() * vector2.get_z()
              * vector1.get_z() * vector2.get_z();
          second_xy3[i] = vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          second_xz3[i] = vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          second_yz3[i] = vector1.get_z() * vector1.get_y() * vector2.get_z()
              * vector2.get_y();
        }
      }
      x0xt_sum += (first_x1[i] / (double) count)
          + (first_x2[i] / (double) count) + (first_x3[i] / (double) count);
      y0yt_sum += (first_y1[i] / (double) count)
          + (first_y2[i] / (double) count) + (first_y3[i] / (double) count);
      z0zt_sum += (first_z1[i] / (double) count)
          + (first_z2[i] / (double) count) + (first_z3[i] / (double) count);
      x0xt_squared_sum += (second_xsquare1[i] / (double) count)
          + (second_xsquare2[i] / (double) count)
          + (second_xsquare3[i] / (double) count);
      y0yt_squared_sum += (second_ysquare1[i] / (double) count)
          + (second_ysquare2[i] / (double) count)
          + (second_ysquare3[i] / (double) count);
      z0zt_squared_sum += (second_zsquare1[i] / (double) count)
          + (second_zsquare2[i] / (double) count)
          + (second_zsquare3[i] / (double) count);
      x0y0xtyt_sum += (second_xy1[i] / (double) count)
          + (second_xy2[i] / (double) count) + (second_xy3[i] / (double) count);
      x0z0xtzt_sum += (second_xz1[i] / (double) count)
          + (second_xz2[i] / (double) count) + (second_xz3[i] / (double) count);
      y0z0ytzt_sum += (second_yz1[i] / (double) count)
          + (second_yz2[i] / (double) count) + (second_yz3[i] / (double) count);
    }
    x0xt_sum = x0xt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    y0yt_sum = y0yt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    z0zt_sum = z0zt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    x0xt_squared_sum = x0xt_squared_sum
        / (number_of_molecules / bonds_per_carbon) * 0.5;
    y0yt_squared_sum = y0yt_squared_sum
        / (number_of_molecules / bonds_per_carbon) * 0.5;
    z0zt_squared_sum = z0zt_squared_sum
        / (number_of_molecules / bonds_per_carbon) * 0.5;
    x0y0xtyt_sum = x0y0xtyt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    x0z0xtzt_sum = x0z0xtzt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    y0z0ytzt_sum = y0z0ytzt_sum / (number_of_molecules / bonds_per_carbon) * 0.5;
    p1 = x0xt_sum + y0yt_sum + z0zt_sum;
    p2 = 1.5
        * (x0xt_squared_sum + y0yt_squared_sum + z0zt_squared_sum
            + 2.0 * x0y0xtyt_sum + 2.0 * x0z0xtzt_sum + 2.0 * y0z0ytzt_sum)
        - 0.5;
    write_to_file << spacing << " " << x0xt_sum << " " << y0yt_sum << " "
                  << z0zt_sum << " " << x0xt_squared_sum << " "
                  << y0yt_squared_sum << " " << z0zt_squared_sum << " "
                  << x0y0xtyt_sum << " " << x0z0xtzt_sum << " " << y0z0ytzt_sum
                  << " " << p1 << " " << p2 << std::endl;
  }

}
