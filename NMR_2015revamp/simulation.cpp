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

  double x0xt_sum1 = 0.0, y0yt_sum1 = 0.0, z0zt_sum1 = 0.0;
  //Second legendre sums for left
  double x0xt_squared_sum1 = 0.0, y0yt_squared_sum1 = 0.0, z0zt_squared_sum1 =
      0.0;
  double x0y0xtyt_sum1 = 0.0, x0z0xtzt_sum1 = 0.0, y0z0ytzt_sum1 = 0.0;
  //First legendre sums for right
  double x0xt_sum2 = 0.0, y0yt_sum2 = 0.0, z0zt_sum2 = 0.0;
  //Second legendre sums for right
  double x0xt_squared_sum2 = 0.0, y0yt_squared_sum2 = 0.0, z0zt_squared_sum2 =
      0.0;
  double x0y0xtyt_sum2 = 0.0, x0z0xtzt_sum2 = 0.0, y0z0ytzt_sum2 = 0.0;
  //First legendre sums for center
  double x0xt_sum3 = 0.0, y0yt_sum3 = 0.0, z0zt_sum3 = 0.0;
  //Second legendre sums for center
  double x0xt_squared_sum3 = 0.0, y0yt_squared_sum3 = 0.0, z0zt_squared_sum3 =
      0.0;
  double x0y0xtyt_sum3 = 0.0, x0z0xtzt_sum3 = 0.0, y0z0ytzt_sum3 = 0.0;

  //Then, declare arrays for the output results vs t
  std::vector<double> first_x1, first_y1, first_z1;
  std::vector<double> second_xsquare1, second_ysquare1, second_zsquare1;
  std::vector<double> second_xy1, second_xz1, second_yz1;
  std::vector<int> spacing_array;

  std::vector<double> first_x2, first_y2, first_z2;
  std::vector<double> second_xsquare2, second_ysquare2, second_zsquare2;
  std::vector<double> second_xy2, second_xz2, second_yz2;

  std::vector<double> first_x3, first_y3, first_z3;
  std::vector<double> second_xsquare3, second_ysquare3, second_zsquare3;
  std::vector<double> second_xy3, second_xz3, second_yz3;

  std::vector<double> first_x_tot, first_y_tot, first_z_tot;
  std::vector<double> second_xsquare_tot, second_ysquare_tot,
      second_zsquare_tot;
  std::vector<double> second_xy_tot, second_xz_tot, second_yz_tot;
  //Temp vectors to hold vectors we select
  vector_3d vector1, vector2;

  //String operations to make output filename
  std::string filename;
  filename = molecules->at(0).getName() + "_"
      + molecules->at(0).getMoleculeType() + "_legendre";

  //Declare initial results arrays as empty
  for (int i = 0; i < number_of_frames; i++) {
    spacing_array.push_back(i);
    first_x1.push_back(0.0);
    first_y1.push_back(0.0);
    first_z1.push_back(0.0);
    second_xsquare1.push_back(0.0);
    second_ysquare1.push_back(0.0);
    second_zsquare1.push_back(0.0);
    second_xy1.push_back(0.0);
    second_xz1.push_back(0.0);
    second_yz1.push_back(0.0);
    first_x2.push_back(0.0);
    first_y2.push_back(0.0);
    first_z2.push_back(0.0);
    second_xsquare2.push_back(0.0);
    second_ysquare2.push_back(0.0);
    second_zsquare2.push_back(0.0);
    second_xy2.push_back(0.0);
    second_xz2.push_back(0.0);
    second_yz2.push_back(0.0);
    first_x3.push_back(0.0);
    first_y3.push_back(0.0);
    first_z3.push_back(0.0);
    second_xsquare3.push_back(0.0);
    second_ysquare3.push_back(0.0);
    second_zsquare3.push_back(0.0);
    second_xy3.push_back(0.0);
    second_xz3.push_back(0.0);
    second_yz3.push_back(0.0);
    first_x_tot.push_back(0.0);
    first_y_tot.push_back(0.0);
    first_z_tot.push_back(0.0);
    second_xsquare_tot.push_back(0.0);
    second_ysquare_tot.push_back(0.0);
    second_zsquare_tot.push_back(0.0);
    second_xy_tot.push_back(0.0);
    second_xz_tot.push_back(0.0);
    second_yz_tot.push_back(0.0);

  }
  for (int i = 0; i < number_of_molecules; i++) {
    for (int j = 0; j < number_of_frames; j++) {
      for (int k = 0; k < (number_of_frames - j); k++) {
        if (i % (bonds_per_carbon) == 0) {
          //Then we know it is a left bond
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k + j);
          x0xt_sum1 += vector1.get_x() * vector2.get_x();
          y0yt_sum1 += vector1.get_y() * vector2.get_y();
          z0zt_sum1 += vector1.get_z() * vector2.get_z();
          x0xt_squared_sum1 += pow(vector1.get_x(), 2)
              * pow(vector2.get_x(), 2);
          y0yt_squared_sum1 += pow(vector1.get_y(), 2)
              * pow(vector2.get_y(), 2);
          z0zt_squared_sum1 += pow(vector1.get_z(), 2)
              * pow(vector2.get_z(), 2);
          x0y0xtyt_sum1 += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          x0z0xtzt_sum1 += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          y0z0ytzt_sum1 += vector1.get_y() * vector1.get_z() * vector2.get_y()
              * vector2.get_z();
        }
        if (i % (bonds_per_carbon) == 1) {
          //Then we know it is a right bond
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k + j);
          x0xt_sum2 += vector1.get_x() * vector2.get_x();
          y0yt_sum2 += vector1.get_y() * vector2.get_y();
          z0zt_sum2 += vector1.get_z() * vector2.get_z();
          x0xt_squared_sum2 += pow(vector1.get_x(), 2)
              * pow(vector2.get_x(), 2);
          y0yt_squared_sum2 += pow(vector1.get_y(), 2)
              * pow(vector2.get_y(), 2);
          z0zt_squared_sum2 += pow(vector1.get_z(), 2)
              * pow(vector2.get_z(), 2);
          x0y0xtyt_sum2 += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          x0z0xtzt_sum2 += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          y0z0ytzt_sum2 += vector1.get_y() * vector1.get_z() * vector2.get_y()
              * vector2.get_z();
        }
        if (i % (bonds_per_carbon) == 2) {
          //Then we know it is a center bond
          vector1 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k);
          vector2 = molecules->at(i).getMoleculeTrajectory()
              ->get_trajectory_array()->at(k + j);
          x0xt_sum3 += vector1.get_x() * vector2.get_x();
          y0yt_sum3 += vector1.get_y() * vector2.get_y();
          z0zt_sum3 += vector1.get_z() * vector2.get_z();
          x0xt_squared_sum3 += pow(vector1.get_x(), 2)
              * pow(vector2.get_x(), 2);
          y0yt_squared_sum3 += pow(vector1.get_y(), 2)
              * pow(vector2.get_y(), 2);
          z0zt_squared_sum3 += pow(vector1.get_z(), 2)
              * pow(vector2.get_z(), 2);
          x0y0xtyt_sum3 += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          x0z0xtzt_sum3 += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          y0z0ytzt_sum3 += vector1.get_y() * vector1.get_z() * vector2.get_y()
              * vector2.get_z();
        }
      }
      if (i % (bonds_per_carbon) == 0) {
        x0xt_sum1 = x0xt_sum1 * (1.0 / ((double) (number_of_frames - j)));
        //      std::cout << "j= " << j << " " << "x0xt_sum" << x0xt_sum << std::endl;
        y0yt_sum1 = y0yt_sum1 * (1.0 / ((double) (number_of_frames - j)));
        z0zt_sum1 = z0zt_sum1 * (1.0 / ((double) (number_of_frames - j)));
        x0xt_squared_sum1 = x0xt_squared_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        y0yt_squared_sum1 = y0yt_squared_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        z0zt_squared_sum1 = z0zt_squared_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        x0y0xtyt_sum1 = x0y0xtyt_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        x0z0xtzt_sum1 = x0z0xtzt_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        y0z0ytzt_sum1 = y0z0ytzt_sum1
            * (1.0 / ((double) (number_of_frames - j)));
        first_x1.at(j) += (x0xt_sum1);
        first_y1.at(j) += (y0yt_sum1);
        first_z1.at(j) += (z0zt_sum1);
        second_xsquare1.at(j) += (x0xt_squared_sum1);
        second_ysquare1.at(j) += (y0yt_squared_sum1);
        second_zsquare1.at(j) += (z0zt_squared_sum1);
        second_xy1.at(j) += (x0y0xtyt_sum1);
        second_xz1.at(j) += (x0z0xtzt_sum1);
        second_yz1.at(j) += (y0z0ytzt_sum1);
        x0xt_sum1 = 0.0;
        y0yt_sum1 = 0.0;
        z0zt_sum1 = 0.0;
        x0xt_squared_sum1 = 0.0;
        y0yt_squared_sum1 = 0.0;
        z0zt_squared_sum1 = 0.0;
        x0y0xtyt_sum1 = 0.0;
        x0z0xtzt_sum1 = 0.0;
        y0z0ytzt_sum1 = 0.0;
      } else if (i % (bonds_per_carbon) == 1) {
        x0xt_sum2 = x0xt_sum2 * (1.0 / ((double) (number_of_frames - j)));
        //      std::cout << "j= " << j << " " << "x0xt_sum" << x0xt_sum << std::endl;
        y0yt_sum2 = y0yt_sum2 * (1.0 / ((double) (number_of_frames - j)));
        z0zt_sum2 = z0zt_sum2 * (1.0 / ((double) (number_of_frames - j)));
        x0xt_squared_sum2 = x0xt_squared_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        y0yt_squared_sum2 = y0yt_squared_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        z0zt_squared_sum2 = z0zt_squared_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        x0y0xtyt_sum2 = x0y0xtyt_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        x0z0xtzt_sum2 = x0z0xtzt_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        y0z0ytzt_sum2 = y0z0ytzt_sum2
            * (1.0 / ((double) (number_of_frames - j)));
        first_x2.at(j) += (x0xt_sum2);
        first_y2.at(j) += (y0yt_sum2);
        first_z2.at(j) += (z0zt_sum2);
        second_xsquare2.at(j) += (x0xt_squared_sum2);
        second_ysquare2.at(j) += (y0yt_squared_sum2);
        second_zsquare2.at(j) += (z0zt_squared_sum2);
        second_xy2.at(j) += (x0y0xtyt_sum2);
        second_xz2.at(j) += (x0z0xtzt_sum2);
        second_yz2.at(j) += (y0z0ytzt_sum2);
        x0xt_sum2 = 0.0;
        y0yt_sum2 = 0.0;
        z0zt_sum2 = 0.0;
        x0xt_squared_sum2 = 0.0;
        y0yt_squared_sum2 = 0.0;
        z0zt_squared_sum2 = 0.0;
        x0y0xtyt_sum2 = 0.0;
        x0z0xtzt_sum2 = 0.0;
        y0z0ytzt_sum2 = 0.0;

      } else if (i % (bonds_per_carbon) == 2) {
        x0xt_sum3 = x0xt_sum3 * (1.0 / ((double) (number_of_frames - j)));
        //      std::cout << "j= " << j << " " << "x0xt_sum" << x0xt_sum << std::endl;
        y0yt_sum3 = y0yt_sum3 * (1.0 / ((double) (number_of_frames - j)));
        z0zt_sum3 = z0zt_sum1 * (1.0 / ((double) (number_of_frames - j)));
        x0xt_squared_sum3 = x0xt_squared_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        y0yt_squared_sum3 = y0yt_squared_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        z0zt_squared_sum3 = z0zt_squared_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        x0y0xtyt_sum3 = x0y0xtyt_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        x0z0xtzt_sum3 = x0z0xtzt_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        y0z0ytzt_sum3 = y0z0ytzt_sum3
            * (1.0 / ((double) (number_of_frames - j)));
        first_x3.at(j) += (x0xt_sum3);
        first_y3.at(j) += (y0yt_sum3);
        first_z3.at(j) += (z0zt_sum3);
        second_xsquare3.at(j) += (x0xt_squared_sum3);
        second_ysquare3.at(j) += (y0yt_squared_sum3);
        second_zsquare3.at(j) += (z0zt_squared_sum3);
        second_xy3.at(j) += (x0y0xtyt_sum3);
        second_xz3.at(j) += (x0z0xtzt_sum3);
        second_yz3.at(j) += (y0z0ytzt_sum3);
        x0xt_sum3 = 0.0;
        y0yt_sum3 = 0.0;
        z0zt_sum3 = 0.0;
        x0xt_squared_sum3 = 0.0;
        y0yt_squared_sum3 = 0.0;
        z0zt_squared_sum3 = 0.0;
        x0y0xtyt_sum3 = 0.0;
        x0z0xtzt_sum3 = 0.0;
        y0z0ytzt_sum3 = 0.0;
      }
    }
  }
  if (bonds_per_carbon == 1) {
    for (int i = 0; i < number_of_frames; i++) {
      first_x1.at(i) = first_x1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y1.at(i) = first_y1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z1.at(i) = first_z1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare1.at(i) = second_xsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare1.at(i) = second_ysquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare1.at(i) = second_zsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy1.at(i) = second_xy1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz1.at(i) = second_xz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz1.at(i) = second_yz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
  }
  if (bonds_per_carbon == 2) {
    for (int i = 0; i < number_of_frames; i++) {
      first_x1.at(i) = first_x1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y1.at(i) = first_y1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z1.at(i) = first_z1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare1.at(i) = second_xsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare1.at(i) = second_ysquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare1.at(i) = second_zsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy1.at(i) = second_xy1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz1.at(i) = second_xz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz1.at(i) = second_yz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
    for (int i = 0; i < number_of_frames; i++) {
      first_x2.at(i) = first_x2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y2.at(i) = first_y2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z2.at(i) = first_z2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare2.at(i) = second_xsquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare2.at(i) = second_ysquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare2.at(i) = second_zsquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy2.at(i) = second_xy2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz2.at(i) = second_xz2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz2.at(i) = second_yz2.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
  }
  if (bonds_per_carbon == 3) {
    for (int i = 0; i < number_of_frames; i++) {
      first_x1.at(i) = first_x1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y1.at(i) = first_y1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z1.at(i) = first_z1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare1.at(i) = second_xsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare1.at(i) = second_ysquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare1.at(i) = second_zsquare1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy1.at(i) = second_xy1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz1.at(i) = second_xz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz1.at(i) = second_yz1.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
    for (int i = 0; i < number_of_frames; i++) {
      first_x2.at(i) = first_x2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y2.at(i) = first_y2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z2.at(i) = first_z2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare2.at(i) = second_xsquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare2.at(i) = second_ysquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare2.at(i) = second_zsquare2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy2.at(i) = second_xy2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz2.at(i) = second_xz2.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz2.at(i) = second_yz2.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
    for (int i = 0; i < number_of_frames; i++) {
      first_x3.at(i) = first_x3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_y3.at(i) = first_y3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      first_z3.at(i) = first_z3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xsquare3.at(i) = second_xsquare3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_ysquare3.at(i) = second_ysquare3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_zsquare3.at(i) = second_zsquare3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xy3.at(i) = second_xy3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_xz3.at(i) = second_xz3.at(i)
          * (1.0 / (((double) number_of_molecules)));
      second_yz3.at(i) = second_yz3.at(i)
          * (1.0 / (((double) number_of_molecules)));
    }
  }
  for (int i = 0; i < number_of_frames; i++) {
    if (bonds_per_carbon == 1) {
      first_x_tot.at(i) += first_x1.at(i);
      first_y_tot.at(i) += first_y1.at(i);
      first_z_tot.at(i) += first_z1.at(i);
      second_xsquare_tot.at(i) += second_xsquare1.at(i);
      second_ysquare_tot.at(i) += second_ysquare1.at(i);
      second_zsquare_tot.at(i) += second_zsquare1.at(i);
      second_xy_tot.at(i) += second_xy1.at(i);
      second_xz_tot.at(i) += second_xz1.at(i);
      second_yz_tot.at(i) += second_yz1.at(i);
    }
    if (bonds_per_carbon >= 2) {
      first_x_tot.at(i) += first_x1.at(i);
      first_y_tot.at(i) += first_y1.at(i);
      first_z_tot.at(i) += first_z1.at(i);
      second_xsquare_tot.at(i) += second_xsquare1.at(i);
      second_ysquare_tot.at(i) += second_ysquare1.at(i);
      second_zsquare_tot.at(i) += second_zsquare1.at(i);
      second_xy_tot.at(i) += second_xy1.at(i);
      second_xz_tot.at(i) += second_xz1.at(i);
      second_yz_tot.at(i) += second_yz1.at(i);
      first_x_tot.at(i) += first_x2.at(i);
      first_y_tot.at(i) += first_y2.at(i);
      first_z_tot.at(i) += first_z2.at(i);
      second_xsquare_tot.at(i) += second_xsquare2.at(i);
      second_ysquare_tot.at(i) += second_ysquare2.at(i);
      second_zsquare_tot.at(i) += second_zsquare2.at(i);
      second_xy_tot.at(i) += second_xy2.at(i);
      second_xz_tot.at(i) += second_xz2.at(i);
      second_yz_tot.at(i) += second_yz2.at(i);
    }
    if (bonds_per_carbon >= 3) {
      first_x_tot.at(i) += first_x1.at(i);
      first_y_tot.at(i) += first_y1.at(i);
      first_z_tot.at(i) += first_z1.at(i);
      second_xsquare_tot.at(i) += second_xsquare1.at(i);
      second_ysquare_tot.at(i) += second_ysquare1.at(i);
      second_zsquare_tot.at(i) += second_zsquare1.at(i);
      second_xy_tot.at(i) += second_xy1.at(i);
      second_xz_tot.at(i) += second_xz1.at(i);
      second_yz_tot.at(i) += second_yz1.at(i);
      first_x_tot.at(i) += first_x2.at(i);
      first_y_tot.at(i) += first_y2.at(i);
      first_z_tot.at(i) += first_z2.at(i);
      second_xsquare_tot.at(i) += second_xsquare2.at(i);
      second_ysquare_tot.at(i) += second_ysquare2.at(i);
      second_zsquare_tot.at(i) += second_zsquare2.at(i);
      second_xy_tot.at(i) += second_xy2.at(i);
      second_xz_tot.at(i) += second_xz2.at(i);
      second_yz_tot.at(i) += second_yz2.at(i);
      first_x_tot.at(i) += first_x3.at(i);
      first_y_tot.at(i) += first_y3.at(i);
      first_z_tot.at(i) += first_z3.at(i);
      second_xsquare_tot.at(i) += second_xsquare3.at(i);
      second_ysquare_tot.at(i) += second_ysquare3.at(i);
      second_zsquare_tot.at(i) += second_zsquare3.at(i);
      second_xy_tot.at(i) += second_xy3.at(i);
      second_xz_tot.at(i) += second_xz3.at(i);
      second_yz_tot.at(i) += second_yz3.at(i);
    }
  }
  std::ofstream write_to_file(filename);
  write_to_file << "# t" << " " << "<x(0)x(t)>" << " " << "<y(0)y(t)>" << " "
                << "<z(0)z(t)>" << " " << "<x^2(0)x^2(t)>" << " "
                << "<y^2(0)y^2(t)>" << " " << "<z^2(0)z^2(t)>" << " "
                << "<x(0)y(0)x(t)y(t)>" << " " << "<x(0)z(0)x(t)z(t)>" << " "
                << "<y(0)z(0)y(t)z(t)>" << " " << "C1" << " " << "C2"
                << std::endl;
  double first_legendre, second_legendre;
  for (int i = 0; i < number_of_frames; i++) {
    first_legendre = first_x_tot.at(i) + first_y_tot.at(i) + first_z_tot.at(i);
    second_legendre = (1.5
        * (second_xsquare_tot.at(i) + second_ysquare_tot.at(i)
            + second_zsquare_tot.at(i) + (2.0 * second_xy_tot.at(i))
            + (2.0 * second_xz_tot.at(i)) + (2.0 * second_yz_tot.at(i)))) - 0.5;
    write_to_file << spacing_array.at(i) << " " << first_x_tot.at(i) << " "
                  << first_y_tot.at(i) << " " << first_z_tot.at(i) << " "
                  << second_xsquare_tot.at(i) << " " << second_ysquare_tot.at(i)
                  << " " << second_zsquare_tot.at(i) << " "
                  << second_xy_tot.at(i) << " " << second_xz_tot.at(i) << " "
                  << second_yz_tot.at(i) << " " << first_legendre << " "
                  << second_legendre << std::endl;
  }

}
