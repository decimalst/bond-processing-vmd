/*
 * fileinput.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: byron
 */
#include "fileinput.h"

std::string cholesterolcarbons[] = { "C20", "C21", "C22", "C23", "C24", "C25",
    "C26", "C27" };
std::string dspccarbonssn1[] = { "C22", "C23", "C24", "C25", "C26", "C27",
    "C28", "C29", "C210", "C211", "C212", "C213", "C214", "C215", "C216",
    "C217", "C218" };
std::string* dspc_carbons_sn1_ptr = &dspccarbonssn1[0];
int dspccarbonssn1size = 17;
std::string dspccarbonssn2[] = { "C32", "C33", "C34", "C35", "C36", "C37",
    "C38", "C39", "C310", "C311", "C312", "C313", "C314", "C315", "C316",
    "C317", "C318" };
std::string* dspc_carbons_sn2_ptr = &dspccarbonssn2[0];
int dspccarbonssn2size = 17;
std::string popccarbons[] = { "C22", "C23", "C24", "C25", "C26", "C27", "C28",
    "C29", "C32", "C33", "C34", "C35", "C36", "C37", "C38", "C39", "C210",
    "C211", "C212", "C213", "C214", "C215", "C216", "C217", "C218", "C310",
    "C311", "C312", "C313", "C314", "C315", "C316" };

std::string dopccarbonssn1[] = { "C22", "C23", "C24", "C25", "C26", "C27",
    "C28", "C29", "C210", "C211", "C212", "C213", "C214", "C215", "C216",
    "C217", "C218" };
std::string* dopc_carbons_sn1_ptr = &dopccarbonssn1[0];
int dopccarbonssn1bonds[] =
    { 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 3 };
int dopccarbonssn1size = 17;

std::string dopccarbonssn2[] = { "C32", "C33", "C34", "C35", "C36", "C37",
    "C38", "C39", "C310", "C311", "C312", "C313", "C314", "C315", "C316",
    "C317", "C318" };
int dopccarbonssn2bonds[] = { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3 };
int dopccarbonssn2size = 17;
std::string* dopc_carbons_sn2_ptr = &dopccarbonssn2[0];

std::string dppccarbonssn1[] = { "C22", "C23", "C24", "C25", "C26", "C27",
    "C28", "C29", "C210", "C211", "C212", "C213", "C214", "C215", "C216" };
int dppccarbonssn1bonds[] =
    { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3 };
int dppccarbonssn1size = 15;
std::string* dppc_carbons_sn1_ptr = &dppccarbonssn1[0];

std::string dppccarbonssn2[] = { "C32", "C33", "C34", "C35", "C36", "C37",
    "C38", "C39", "C310", "C311", "C312", "C313", "C314", "C315", "C316" };
int dppccarbonssn2bonds[] =
    { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3 };
int dppccarbonssn2size = 15;
std::string* dppc_carbons_sn2_ptr = &dppccarbonssn2[0];

fileinput::fileinput(std::string file_name, std::string carbon_name_given,
                     std::string molecule_type) {
  filename = file_name;
  lipid_name = molecule_type;
  carbon_name = carbon_name_given;
  isvalidfile = false;
  number_of_lines = 0;
  lines_per_frame = 0;
}
fileinput::fileinput() {
  filename = ".null";
  lipid_name = "null";
  carbon_name = "null";
  isvalidfile = false;
  number_of_lines = 0;
  lines_per_frame = 0;
}
int fileinput::find_atom_count() {
  is_file_valid();
  if (!isvalidfile) {
    return 0;
  }
  int number_of_molecules = 0;
  std::ifstream file;
  file.open(
      filename.c_str());
  std::string templine;
  for (int i = 0; i <= 1; i++) {
    std::getline(
        file, templine);
    if (i == 0) {
      number_of_molecules = atoi(
          templine.c_str());
    }
  }
  lines_per_frame = number_of_molecules + 2;
  return number_of_molecules;
}

void fileinput::is_file_valid() {
  if (filename.find(
      xyz_extension, 0) != std::string::npos) {
    std::cout << "file is valid" << std::endl;
    isvalidfile = true;
  } else {
    std::cout << "file is not valid" << std::endl;
  }
}
int fileinput::get_line_number() {
  //Gets number of lines in a file
  is_file_valid();
  if (!isvalidfile) {
    return 0;
  }
  std::ifstream filetocount(
      filename);
  int i = 0;
  std::string line;
  while (getline(
      filetocount, line)) {
    i++;
  }
  return i;
}
int fileinput::find_in_array(std::string check_this, std::string* array_in,
                             int size) {
  int index = 0;
  for (int i = 0; i < size; i++) {
    if (array_in[i] == check_this) {
      index = i;
    }
    return index;
  }
}
std::vector<molecule_bond>* fileinput::loadmolecules() {
  int bonds_per_carbon = 0;
  int index_of_size = 0;
  int lines_in_file;
  //now, from the lipid_name and molecule_name, call find_in_array in order to find the number of bonds per
  //carbon for this specific carbon under consideration
  //pass a pointer to the start of one of the arrays listed at top of the header file
  if (lipid_name == "dspc") {
    index_of_size = find_in_array(
        carbon_name, dspc_carbons_sn1_ptr, dspccarbonssn1size);
    bonds_per_carbon = dopccarbonssn1bonds[index_of_size];
  }
  if (lipid_name == "dopc") {
    index_of_size = find_in_array(
        carbon_name, dopc_carbons_sn1_ptr, dopccarbonssn1size);
    bonds_per_carbon = dopccarbonssn1bonds[index_of_size];
  }
//TODO: Add sn1,sn2 chains for popc
//  if(lipid_name=="popc"){
//    index_of_size=find_in_array(carbon_name,popc_carbons_sn1_ptr,popccarbonssn1size);
//
//  }
  if (lipid_name == "dppc") {
    index_of_size = find_in_array(
        carbon_name, dppc_carbons_sn1_ptr, dppccarbonssn1size);
    bonds_per_carbon = dppccarbonssn1bonds[index_of_size];

  }
  if (lipid_name =="debug"){
    bonds_per_carbon=2;
  }
//Reads through .xyz file one line at a time, for each molecule, creates a trajectory for each C-H bond
//Control loop variables
  int current_frame = 0;
  int current_line_absolute = 0;
  int current_line_in_frame = 0;
  int per_line_count = 0;
  int bond_holder_index = 0;

  int atoms_per_frame = find_atom_count();
//Variables needed for file IO
  std::ifstream filetoload(
      filename);
  std::string templine;

//This holds the bonds that we load from the file
  std::vector<molecule_bond>* bondholder;
  bondholder = new std::vector<molecule_bond>;

//Holder variables used in inputting the float values from the file
  double temp_x, temp_y, temp_z;
//Holder variable for the string in the first column
  std::string identifier_holder;

//Holder variables for the bonds we load in, that are then put into the bondholder vector
  coords_3d temp_carbon, temp_hydrogen1, temp_hydrogen2, temp_hydrogen3;
  vector_3d temp_bond1, temp_bond2, temp_bond3;

  //Holder variables for molecule_bonds in the bondholder vector
  molecule_bond temp_mol_bond1, temp_mol_bond2, temp_mol_bond3;

  //Simple math defining how many bonds we have to load in given the number of atoms:
  //knowing we only have carbons and hydrogens and the number of bonds per carbon lets
  //us calculate how many bonds in each frame
  int number_of_bonds = (atoms_per_frame)
      * ((bonds_per_carbon) / (bonds_per_carbon + 1.0));
//  std::cout << "bonds_per_carbon" << bonds_per_carbon << std::endl;
//  std::cout << "number_of_bonds  " << number_of_bonds << std::endl;
//For each bond in the frame, instantiate it and add it to the bondholder vector
  for (int i = 0; i < number_of_bonds; i++) {
    bondholder->push_back(
        molecule_bond(
            carbon_name, lipid_name, bonds_per_carbon, new trajectory()));
  }
  std::cout << "bondhold size " << bondholder->size() << std::endl;
//For each line in the frame, get it, and store its value in templine
  while (std::getline(
      filetoload, templine)) {
    current_frame = floor(
        current_line_absolute / (atoms_per_frame + 2));
    current_line_in_frame = current_line_absolute % (atoms_per_frame + 2);
    current_line_absolute++;
//    std::cout <<"bond_holder_index "<< bond_holder_index << std::endl;
//    std::cout << "current_line_in_frame " << current_line_in_frame
//              << " current_frame " << current_frame
//              << " current_line_absolute  " << current_line_absolute
//              << std::endl;
    //   std::cout <<"current_frame"<< current_frame << std::endl;
//    std::cout << "per_line_count " << per_line_count << std::endl;
//    std::cout <<"current_line_absolute "<<current_line_absolute<<std::endl;
    //The first two lines of each frame will be the number of atoms, and the comment line
    //so we can disregard them
    if (current_line_in_frame > 1) {
      //Break the loaded line into a string stream for parsing
      std::stringstream iss(
          templine);
      iss >> identifier_holder >> temp_x >> temp_y >> temp_z;
      //Using the control variable "per_line_count" we know where to store the data
      if (per_line_count == 0 && per_line_count <= bonds_per_carbon) {
        temp_carbon = coords_3d(
            temp_x, temp_y, temp_z);
      }
      if (per_line_count == 1 && per_line_count <= bonds_per_carbon) {
        temp_hydrogen1 = coords_3d(
            temp_x, temp_y, temp_z);
      }
      if (per_line_count == 2 && per_line_count <= bonds_per_carbon) {
        temp_hydrogen2 = coords_3d(
            temp_x, temp_y, temp_z);
      }
      if (per_line_count == 3 && per_line_count <= bonds_per_carbon) {
        temp_hydrogen3 = coords_3d(
            temp_x, temp_y, temp_z);
      }
      //If we've loaded all the hydrogen coords associated with the first carbon
      //We can now calculate the unit vector from the coords and store it in the trajectory
      //associated with the molecule_bond
      if (per_line_count == bonds_per_carbon) {
        //Make a vector_3d of the carbon source to the two hydrogen atom targets
        //Add the vector_3ds to the trajectory of the molecular bond it belongs to
        //with the same index.
        if (bonds_per_carbon == 1) {
          temp_bond1 = vector_3d(
              temp_carbon, temp_hydrogen1);
          temp_mol_bond1 = bondholder->at(
              bond_holder_index);
          temp_mol_bond1.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond1);
          bond_holder_index++;
        }
        if (bonds_per_carbon == 2) {
          temp_bond1 = vector_3d(
              temp_carbon, temp_hydrogen1);
          temp_bond2 = vector_3d(
              temp_carbon, temp_hydrogen2);
          temp_mol_bond1 = bondholder->at(
              bond_holder_index);
          temp_mol_bond1.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond1);
          bond_holder_index++;
          temp_mol_bond2 = bondholder->at(
              bond_holder_index);
          temp_mol_bond2.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond2);
          bond_holder_index++;
        }
        if (bonds_per_carbon == 3) {
          temp_bond1 = vector_3d(
              temp_carbon, temp_hydrogen1);
          temp_bond2 = vector_3d(
              temp_carbon, temp_hydrogen2);
          temp_bond3 = vector_3d(
              temp_carbon, temp_hydrogen3);
          temp_bond1 = vector_3d(
              temp_carbon, temp_hydrogen1);
          temp_bond2 = vector_3d(
              temp_carbon, temp_hydrogen2);
          temp_mol_bond1 = bondholder->at(
              bond_holder_index);
          temp_mol_bond1.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond1);
          bond_holder_index++;
          temp_mol_bond2 = bondholder->at(
              bond_holder_index);
          temp_mol_bond2.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond2);
          bond_holder_index++;
          temp_mol_bond3 = bondholder->at(
              bond_holder_index);
          temp_mol_bond3.getMoleculeTrajectory()->get_trajectory_array()
              ->push_back(
              temp_bond3);
          bond_holder_index++;
        }
      }
      per_line_count++;
      if (per_line_count > bonds_per_carbon) {
        per_line_count = 0;
      }
      if (bond_holder_index == bondholder->size()) {
        bond_holder_index = 0;
      }

    }
    //Repeat per frame
  }
  return bondholder;
}

