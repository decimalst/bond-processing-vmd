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
  for(int i=0; i<molecules_in->size(); i++){
    molecules->at(i).getMoleculeTrajectory()->set_number_of_frames(number_of_frames);
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
      angleholder[i] += fabs(((3/2)*pow(cos(molecules->at(i).getMoleculeTrajectory()
          ->get_trajectory_array()->at(j).get_angle_from(&zaxis)),2))-0.5);
    }
    angleholder[i] =angleholder[i] / (float)number_of_frames;
    average_over_ensemble += angleholder[i];
  }
  average_over_ensemble = average_over_ensemble / (float)number_of_molecules;
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
  int number_of_frames=get_number_of_frames();
  double x0xt_sum, y0yt_sum = 0.0, z0zt_sum = 0.0;
  //Second legendre sums
  double x0xt_squared_sum = 0.0, y0yt_squared_sum = 0.0, z0zt_squared_sum = 0.0;
  double x0y0xtyt_sum = 0.0, x0z0xtzt_sum = 0.0, y0z0ytzt_sum = 0.0;

  //Then, declare arrays for the output results vs t
  std::vector<double> first_x, first_y, first_z;
  std::vector<double> second_xsquare, second_ysquare, second_zsquare;
  std::vector<double> second_xy, second_xz, second_yz;
  std::vector<int> spacing_array;
  //Temp vectors to hold vectors we select
  vector_3d vector1, vector2;

  //String operations to make output filename
  std::string filename;
  filename = molecules->at(0).getName() + "_"
      + molecules->at(0).getMoleculeType() + "_legendre";

//Declare initial results arrays as empty
  for (int i = 0; i < number_of_frames; i++) {
    spacing_array.push_back(i);
    first_x.push_back(0.0);
    first_y.push_back(0.0);
    first_z.push_back(0.0);
    second_xsquare.push_back(0.0);
    second_ysquare.push_back(0.0);
    second_zsquare.push_back(0.0);
    second_xy.push_back(0.0);
    second_xz.push_back(0.0);
    second_yz.push_back(0.0);

  }

  //Step 1: For some molecular bond in simulation
  for (int i = 0; i < number_of_molecules; i++) {
    //Step 2: For all values of t-"spacing between frames" from t=1 to t=(n-1)
    //Let j represent t here
    for (int j = 1; j < number_of_frames; j++) {
      //Step 3: For all pairs of frames in simulation with spacing t.
      for (int k = 0; k < (number_of_frames - j); k++) {
        //Compute x(0)x(t), y(0)y(t), etc. for all involved pairs for C1, C2
        //this will be x[k]*x[k+j], y[k]y[k+j], etc.
        vector1 = molecules->at(i).getMoleculeTrajectory()->get_trajectory_array()->at(k);
        vector2 = molecules->at(i).getMoleculeTrajectory()->get_trajectory_array()->at(k+j);
        x0xt_sum += vector1.get_x() * vector2.get_x();
        y0yt_sum += vector1.get_y() * vector2.get_y();
        z0zt_sum += vector1.get_z() * vector2.get_z();
        x0xt_squared_sum += pow(vector1.get_x(), 2) * pow(vector2.get_x(), 2);
        y0yt_squared_sum += pow(vector1.get_y(), 2) * pow(vector2.get_y(), 2);
        z0zt_squared_sum += pow(vector1.get_z(), 2) * pow(vector2.get_z(), 2);
        x0y0xtyt_sum += vector1.get_x() * vector1.get_y() * vector2.get_x()
            * vector2.get_y();
        x0z0xtzt_sum += vector1.get_x() * vector1.get_z() * vector2.get_x()
            * vector2.get_z();
        y0z0ytzt_sum += vector1.get_y() * vector1.get_z() * vector2.get_y()
            * vector2.get_z();
        //Store in float sums
      }
      //Average float sums of products by (1/number_of_frames-spacing(j) )
//      std::cout << "j= " << j << " " << "x0xt_sum" << x0xt_sum << std::endl;
      x0xt_sum = x0xt_sum * (1.0 / (double)(number_of_frames-j));
//      std::cout << "j= " << j << " " << "x0xt_sum" << x0xt_sum << std::endl;
      y0yt_sum = y0yt_sum * (1.0 / (double)(number_of_frames-j));
      z0zt_sum = z0zt_sum * (1.0 / (double)(number_of_frames-j));
      x0xt_squared_sum = x0xt_squared_sum * (1.0 / (double)(number_of_frames-j));
      y0yt_squared_sum = y0yt_squared_sum * (1.0 / (double)(number_of_frames-j));
      z0zt_squared_sum = z0zt_squared_sum * (double)(1.0 / (number_of_frames-j));
      x0y0xtyt_sum = x0y0xtyt_sum * (1.0 / (double)(number_of_frames-j));
      x0z0xtzt_sum = x0z0xtzt_sum * (1.0 / (double)(number_of_frames-j));
      y0z0ytzt_sum = y0z0ytzt_sum * (1.0 / (double)(number_of_frames-j));

      //Store into results_array(t+=) as a result associated with the value of spacing
      first_x.at(j) = first_x.at(j)+ (x0xt_sum);
      first_y.at(j) += (y0yt_sum);
      first_z.at(j) += (z0zt_sum);
      second_xsquare.at(j) += (x0xt_squared_sum);
      second_ysquare.at(j) += (y0yt_squared_sum);
      second_zsquare.at(j) += (z0zt_squared_sum);
      second_xy.at(j) += (x0y0xtyt_sum);
      second_xz.at(j) += (x0z0xtzt_sum);
      second_yz.at(j) += (y0z0ytzt_sum);

      //Set float sum values back to zero.
      x0xt_sum = 0.0;
      y0yt_sum = 0.0;
      z0zt_sum = 0.0;
      x0xt_squared_sum = 0.0;
      y0yt_squared_sum = 0.0;
      z0zt_squared_sum = 0.0;
      x0y0xtyt_sum = 0.0;
      x0z0xtzt_sum = 0.0;
      y0z0ytzt_sum = 0.0;
    }
  }
  //Now that we have the ensemble sum stored into the arrays, we need to average it by ensemble
  //So multiply by 1/number_of_molecules each value we got for some time sep
  for (int i = 0; i < number_of_frames; i++) {
    first_x.at(i) = first_x.at(i) * (1.0 / ((double)number_of_molecules));
    first_y.at(i) = first_y.at(i) * (1.0 / ((double)number_of_molecules));
    first_z.at(i) = first_z.at(i) * (1.0 / ((double)number_of_molecules));
    second_xsquare.at(i) = second_xsquare.at(i) * (1.0 / ((double)number_of_molecules));
    second_ysquare.at(i) = second_ysquare.at(i) * (1.0 /((double)number_of_molecules));
    second_zsquare.at(i) = second_zsquare.at(i) * (1.0/((double)number_of_molecules));
    second_xy.at(i) = second_xy.at(i) * (1.0 / ((double)number_of_molecules));
    second_xz.at(i) = second_xz.at(i) * (1.0 / ((double)number_of_molecules));
    second_yz.at(i) = second_yz.at(i) * (1.0 / ((double)number_of_molecules));
  }

  //Now write all our values to a file
  std::ofstream write_to_file(filename);
  write_to_file << "t" << " " << "<x(0)x(t)>" << " "
                    << "<y(0)y(t)>" << " " << "<z(0)z(t)>" << " "
                    << "<x^2(0)x^2(t)>" << " " << "<y^2(0)y^2(t)>" << " "
                    << "<z^2(0)z^2(t)>" << " " << "<x(0)y(0)x(t)y(t)>" << " "
                    << "<x(0)z(0)x(t)z(t)>" << " " << "<y(0)z(0)y(t)z(t)>" << " "
                    << "C1" << " " << "C2"
                    << std::endl;
  float first_legendre, second_legendre;
  for (int i = 0; i < number_of_frames; i++) {
    first_legendre= first_x.at(i) + first_y.at(i)+first_z.at(i);
    second_legendre=1.5*(second_xsquare.at(i) + second_ysquare.at(i) + second_zsquare.at(i) + (2.0*second_xy.at(i))+ (2.0*second_xz.at(i))+(2.0*second_yz.at(i)))-0.5;
    write_to_file << spacing_array.at(i) << " " << first_x.at(i) << " "
                  << first_y.at(i) << " " << first_z.at(i) << " "
                  << second_xsquare.at(i) << " " << second_ysquare.at(i) << " "
                  << second_zsquare.at(i) << " " << second_xy.at(i) << " "
                  << second_xz.at(i) << " " << second_yz.at(i) << " "
                  << first_legendre << " " << second_legendre
                  << std::endl;
  }
}
