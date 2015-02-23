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
  for (unsigned int i = 0; i < molecules_in->size(); i++) {
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
  float average_over_sim = 0.0;
//Declare vector of size equal to the number of bonds in the simulation
//In this vector we will store the angle for each bond with z axis
//Once the vector has the angle of the entire simulation added, we
//average using the number of frames in the simulation.
  std::vector<float> angleholder(number_of_molecules);
  std::vector<float> scdholder(number_of_molecules);
  std::cout << angleholder.size() << " " << molecules->size() << std::endl;

  for (int j = 0; j < number_of_frames; j++) {
    for (int i = 0; i < number_of_molecules; i++) {
      scdholder[i] += (1.5 * pow(molecules->at(i).getMoleculeTrajectory()
                                ->get_trajectory_array()->at(j).get_dot_product(&zaxis), 2)) - 0.5;

    }
  }
  for (int i=0; i <number_of_molecules; i++){
    average_over_sim+=scdholder[i] / number_of_frames;
  }
  average_over_sim = average_over_sim / (float)number_of_molecules;
  std::cout << fabs(average_over_sim) << std::endl;

}
void simulation::perl_formatter() {
//Formats data for perl script
  //Calculates first and second legendre polynomials
  int number_of_bonds = number_of_molecules;
  ;
  int bonds_per_carbon = get_molecules()->at(0).getNumberOfBonds();
  int number_of_frames = get_number_of_frames();
  vector_3d vector1, vector2;

  std::string filename;
  filename = molecules->at(0).getName() + "_"
      + molecules->at(0).getMoleculeType() + "_perlformatted";
  std::ofstream write_to_file(filename);
  write_to_file.precision(10);
  std::cout.precision(10);
  for (int frame = 0; frame < number_of_frames; frame++) {
    write_to_file << "frame " << frame << std::endl;
    for (int bond = 0; bond < number_of_bonds; bond += bonds_per_carbon) {
      //std::cout << bond % bonds_per_carbon << std::endl;
      if (bonds_per_carbon == 1) {
        //then it's right bond
        vector1 = molecules->at(bond).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);
        write_to_file << " " << vector1.get_x() << " " << vector1.get_y() << " "
                      << vector1.get_z() << std::endl;
      }
      if (bonds_per_carbon == 2) {
        //then it's a left bond
        vector1 = molecules->at(bond).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);
        write_to_file << bond / (bonds_per_carbon) << " " << vector1.get_x()
                      << " " << vector1.get_y() << " " << vector1.get_z();
        vector1 = molecules->at(bond + 1).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);
        write_to_file << " " << vector1.get_x() << " " << vector1.get_y() << " "
                      << vector1.get_z() << std::endl;

      }
      if (bonds_per_carbon == 3) {
        //then it's a left bond
        vector1 = molecules->at(bond).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);

        vector1 = molecules->at(bond + 1).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);

        vector1 = molecules->at(bond + 2).getMoleculeTrajectory()
            ->get_trajectory_array()->at(frame);
      }
    }

  }

}
void simulation::calculate_second_legendre() {
//Calculates second legendre polynomial and writes it to a text file
}
void simulation::calculate_first_and_second_legendre(std::string outfile) {
//Calculates first and second legendre polynomials
  int number_of_bonds = number_of_molecules;
  int bonds_per_carbon = get_molecules()->at(0).getNumberOfBonds();
  int number_of_frames = get_number_of_frames();
  double x1, y1, z1, xx1, yy1, zz1, xy1, xz1, yz1;
  double x2, y2, z2, xx2, yy2, zz2, xy2, xz2, yz2;
  double x_all, y_all, z_all, xx_all, yy_all, zz_all, xy_all, xz_all, yz_all;
  double p1_sum = 0.0, p2_sum = 0.0;
  vector_3d vector1, vector2;

  std::string filename;
  filename = outfile + "_legendre";

  std::ofstream write_to_file(filename);
  write_to_file << "# t" << " " << "<x(0)x(t)>" << " " << "<y(0)y(t)>" << " "
                << "<z(0)z(t)>" << " " << "<x^2(0)x^2(t)>" << " "
                << "<y^2(0)y^2(t)>" << " " << "<z^2(0)z^2(t)>" << " "
                << "<x(0)y(0)x(t)y(t)>" << " " << "<x(0)z(0)x(t)z(t)>" << " "
                << "<y(0)z(0)y(t)z(t)>" << " " << "C1" << " " << "C2"
                << std::endl;
  write_to_file.precision(10);
  std::cout.precision(10);
  for (int spacing = 0; spacing < number_of_frames; spacing++) {
    std::cout << "spacing= " << spacing << " out of " << number_of_frames
              << std::endl;
    x_all = 0.0;
    y_all = 0.0;
    z_all = 0.0;
    xx_all = 0.0;
    yy_all = 0.0;
    zz_all = 0.0;
    xy_all = 0.0;
    xz_all = 0.0;
    yz_all = 0.0;
    for (int bond = 0; bond < number_of_bonds; bond++) {
      x1 = 0.0;
      x2 = 0.0;
      y1 = 0.0;
      y2 = 0.0;
      z1 = 0.0;
      z2 = 0.0;
      xx1 = 0.0;
      yy1 = 0.0;
      zz1 = 0.0;
      xy1 = 0.0;
      xz1 = 0.0;
      yz1 = 0.0;
      xx2 = 0.0;
      yy2 = 0.0;
      zz2 = 0.0;
      xy2 = 0.0;
      xz2 = 0.0;
      yz2 = 0.0;
      for (int frame = 0; frame < (number_of_frames - spacing); frame++) {
        //std::cout << bond % bonds_per_carbon << std::endl;
        if (bond % bonds_per_carbon == 0) {
          //then it's right bond
          vector1 = molecules->at(bond).getMoleculeTrajectory()
              ->get_trajectory_array()->at(frame);
          vector2 = molecules->at(bond).getMoleculeTrajectory()
              ->get_trajectory_array()->at(frame + spacing);
          x1 += vector1.get_x() * vector2.get_x();
          //std::cout << "x1 "<< x1 << std::endl;
          y1 += vector1.get_y() * vector2.get_y();
          z1 += vector1.get_z() * vector2.get_z();
          xx1 += pow(vector1.get_x(), 2) * pow(vector2.get_x(), 2);
          yy1 += pow(vector1.get_y(), 2) * pow(vector2.get_y(), 2);
          zz1 += pow(vector1.get_z(), 2) * pow(vector2.get_z(), 2);
          xy1 += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          xz1 += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          yz1 += vector1.get_y() * vector1.get_z() * vector2.get_y()
              * vector2.get_z();
        }
        if (bond % bonds_per_carbon == 1) {
          //then it's a left bond
          vector1 = molecules->at(bond).getMoleculeTrajectory()
              ->get_trajectory_array()->at(frame);
          vector2 = molecules->at(bond).getMoleculeTrajectory()
              ->get_trajectory_array()->at(frame + spacing);
          x2 += vector1.get_x() * vector2.get_x();
          y2 += vector1.get_y() * vector2.get_y();
          z2 += vector1.get_z() * vector2.get_z();
          xx2 += pow(vector1.get_x(), 2) * pow(vector2.get_x(), 2);
          yy2 += pow(vector1.get_y(), 2) * pow(vector2.get_y(), 2);
          zz2 += pow(vector1.get_z(), 2) * pow(vector2.get_z(), 2);
          xy2 += vector1.get_x() * vector1.get_y() * vector2.get_x()
              * vector2.get_y();
          xz2 += vector1.get_x() * vector1.get_z() * vector2.get_x()
              * vector2.get_z();
          yz2 += vector1.get_y() * vector1.get_z() * vector2.get_y()
              * vector2.get_z();
        }
      }
      x1 = x1 * (1.0 / (number_of_frames - spacing));
      x2 = x2 * (1.0 / (number_of_frames - spacing));
      y1 = y1 * (1.0 / (number_of_frames - spacing));
      y2 = y2 * (1.0 / (number_of_frames - spacing));
      z1 = z1 * (1.0 / (number_of_frames - spacing));
      z2 = z2 * (1.0 / (number_of_frames - spacing));
      xx1 = xx1 * (1.0 / (number_of_frames - spacing));
      yy1 = yy1 * (1.0 / (number_of_frames - spacing));
      zz1 = zz1 * (1.0 / (number_of_frames - spacing));
      xy1 = xy1 * (1.0 / (number_of_frames - spacing));
      xz1 = xz1 * (1.0 / (number_of_frames - spacing));
      yz1 = yz1 * (1.0 / (number_of_frames - spacing));
      xx2 = xx2 * (1.0 / (number_of_frames - spacing));
      yy2 = yy2 * (1.0 / (number_of_frames - spacing));
      zz2 = zz2 * (1.0 / (number_of_frames - spacing));
      xy2 = xy2 * (1.0 / (number_of_frames - spacing));
      xz2 = xz2 * (1.0 / (number_of_frames - spacing));
      yz2 = yz2 * (1.0 / (number_of_frames - spacing));
      p1_sum += (x1 + x2) + (y1 + y1) + (z1 + z2);
      p2_sum += (3
          * (xx1 + xx2 + yy1 + yy2 + zz1 + zz2 + 2 * (xy1 + xy2)
              + 2 * (xz1 + xz2) + 2 * (yz1 + yz2)) - 1) / 2;
      x_all += x1 + x2;
      y_all += y1 + y2;
      z_all += z1 + z2;
      xx_all += xx1 + xx2;
      yy_all += yy1 + yy2;
      zz_all += zz1 + zz2;
      xy_all += xy1 + xy2;
      xz_all += xz1 + xz2;
      yz_all += yz1 + yz2;
    }
    x_all = x_all * (1.0 / number_of_molecules);
    y_all = y_all * (1.0 / number_of_molecules);
    z_all = z_all * (1.0 / number_of_molecules);
    xx_all = xx_all * (1.0 / number_of_molecules);
    yy_all = yy_all * (1.0 / number_of_molecules);
    zz_all = zz_all * (1.0 / number_of_molecules);
    xy_all = xy_all * (1.0 / number_of_molecules);
    xz_all = xz_all * (1.0 / number_of_molecules);
    yz_all = yz_all * (1.0 / number_of_molecules);
    p1_sum = x_all + y_all + z_all;
    p2_sum = p2_sum * (1.0 / (number_of_molecules));
//    std::cout
//        << "p2_sum "
//        << p2_sum
//        << " vs "
//        << (3
//            * (xx_all + yy_all + zz_all + 2 * (xy_all) + 2 * (xz_all)
//                + 2 * (yz_all)) - 1) / 2;

    write_to_file << spacing << " " << x_all << " " << y_all << " " << z_all
                  << " " << xx_all << " " << yy_all << " " << zz_all << " "
                  << xy_all << " " << xz_all << " " << yz_all << " " << p1_sum
                  << " " << p2_sum << std::endl;
  }

}
