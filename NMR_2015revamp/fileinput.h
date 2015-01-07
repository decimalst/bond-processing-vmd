/*
 * fileinput.h
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */

#ifndef FILEINPUT_H_
#define FILEINPUT_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "molecule.h"
#include "trajectory.h"


const std::string xyz_extension = ".xyz";
class fileinput {
  std::string filename;
  std::string carbon_name;
  std::string lipid_name;
  bool isvalidfile;
  int number_of_lines;
  int lines_per_frame;
 public:
  void openfile(std::string);
  std::vector<molecule_bond>* loadmolecules();
  int find_frame_count();
  int find_atom_count();
  int get_line_number();
  fileinput(std::string, std::string, std::string);
  fileinput();
  void is_file_valid();
  int find_in_array(std::string, std::string*, int);

};

#endif /* FILEINPUT_H_ */
