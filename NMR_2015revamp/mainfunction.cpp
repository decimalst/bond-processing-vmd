//============================================================================
// Name        : NMR_2015revamp.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "vector_3d.h"
#include "coords_3d.h"
#include "fileinput.h"
#include "molecule.h"
#include "trajectory.h"
#include "simulation.h"

int coords_3d_test() {
  coords_3d* derp = new coords_3d(0.0, 0.0, 0.0);
  std::cout << derp->getX() << " " << derp->getY() << " " << derp->getZ()
            << std::endl;
  derp->setX(0.1);
  derp->setY(0.2);
  derp->setZ(0.3);
  std::cout << derp->getX() << " " << derp->getY() << " " << derp->getZ()
            << std::endl;
  delete derp;
  return 1;
}
int vector_3d_test() {
  vector_3d* test = new vector_3d();
  std::cout << test->get_x() << test->get_y() << test->get_z() << std::endl;
  delete test;
  return 1;
}
int fileinput_test() {
  fileinput* test = new fileinput("somefile.xyz", "test", "test");
  test->is_file_valid();
  fileinput* test2 = new fileinput("somefile.zyx", "test", "test");
  test2->is_file_valid();
  delete test;
  delete test2;
  fileinput* test3 = new fileinput("test.xyz", "test", "test");
  test3->is_file_valid();
  std::cout << test3->find_atom_count() << std::endl;
  std::cout << test3->get_line_number() << std::endl;
  delete test3;
  //Test loadmolecules() for "DOPC-C22.xyz"
  fileinput* test4 = new fileinput("DOPC-C22.xyz", "C22", "dopc");
  std::vector<molecule_bond> test_mol_bond;
  test_mol_bond = *test4->loadmolecules();
  std::cout << "test_mol_bond size " << test_mol_bond.size() << std::endl;
  molecule_bond testvectorcheck;
//  testvectorcheck = test_mol_bond.at(
//      1);
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          0).get_x()
//      << std::endl;
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          0).get_y()
//      << std::endl;
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          0).get_z()
//      << std::endl;
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          1).get_x()
//      << std::endl;
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          1).get_y()
//      << std::endl;
//  std::cout
//      << testvectorcheck.getMoleculeTrajectory()->get_trajectory_array()->at(
//          1).get_z()
//      << std::endl;
//  for (int i=0; i<152; i++){
//    std::cout << i << " "
//          << test_mol_bond.at(0).getMoleculeTrajectory()->get_trajectory_array()->at(i).get_x()
//          << " " << test_mol_bond.at(0).getMoleculeTrajectory()->get_trajectory_array()->at(i).get_y()
//          << " " << test_mol_bond.at(0).getMoleculeTrajectory()->get_trajectory_array()->at(i).get_z()
//          << std::endl;
//  }
  return 1;
}
int simulation_test() {
  fileinput* test_file = new fileinput("DOPC-C22.xyz", "C22", "dopc");
  std::vector<molecule_bond>* test_mol_bond;
  test_mol_bond = test_file->loadmolecules();
  simulation * newsim;
  newsim = new simulation(test_mol_bond);
  std::cout << "number of molecules " << newsim->get_number_of_molecules()
            << std::endl;
  std::cout << "number of frames " << newsim->get_number_of_frames()
            << std::endl;
  newsim->calculate_deuterium_order_parameter();
  newsim->calculate_first_and_second_legendre("newsim");
  return 1;
}
void debug_test() {
  fileinput* test_file = new fileinput("debug_test.xyz", "debug", "debug");
  std::vector<molecule_bond>* test_mol_bond;
  test_mol_bond = test_file->loadmolecules();
  simulation * newsim;
  newsim = new simulation(test_mol_bond);
  std::cout << "number of molecules " << newsim->get_number_of_molecules()
            << std::endl;
  std::cout << "number of frames " << newsim->get_number_of_frames()
            << std::endl;
  newsim->calculate_deuterium_order_parameter();
  newsim->calculate_first_and_second_legendre("test");
}

int main(int argc, char* argv[]) {
//  coords_3d_test();
//  vector_3d_test();
// fileinput_test();
// simulation_test();
//  debug_test();
  int choice;
  std::string inputchoice;
  bool isvalid = false;
  std::string filename;
  std::string carbonname;
  std::string lipidname;
  bool skip_prompt = false;
  if(argc > 1){
  if (std::string(argv[1])=="-np"){
    skip_prompt= true;
    filename = std::string(argv[2]);
    carbonname = std::string(argv[3]);
    lipidname = std::string(argv[4]);
    choice = atoi(argv[5]);
    std::cout << filename << " " << carbonname << " " << lipidname << std::endl;

  }
  }
  if (!skip_prompt) {
    std::cout << "Please input .xyz file name" << std::endl;

    std::cin >> filename;
    std::cout << "Please input carbon name(e.g. 'C22')" << std::endl;
    std::cin >> carbonname;
    std::cout << "Please input lipid name(e.g. 'dopc')" << std::endl;
    std::cin >> lipidname;
    std::cout << "You have selected " << lipidname << " carbon: " << carbonname
              << " from file " << filename << std::endl;
    std::cout
        << "Input 1 to calculate SCD parameter only, input 2 to calculate first and second Legendre polynomial only, input 3 to calculate both, input 4 to perl format the bonds"
        << std::endl;

    std::cin >> inputchoice;
    choice = stoi(inputchoice);
    if (choice == 1 || choice == 2 || choice == 3 || choice == 4) {
      isvalid = true;
    }
    while (!isvalid) {
      std::cout << "Invalid choice, try again" << std::endl;
      isvalid = false;
      std::cin >> inputchoice;
      choice = stoi(inputchoice);
      if (choice == 1 || choice == 2 || choice == 3) {
        isvalid = true;
      }
    }
  }
  fileinput* working = new fileinput(filename, carbonname, lipidname);
  working->is_file_valid();
  std::vector<molecule_bond>* mol_bond;
  mol_bond = working->loadmolecules();
  simulation worksim = simulation(mol_bond);
  if (choice == 1) {
    worksim.calculate_deuterium_order_parameter();
  }
  if (choice == 2) {
    worksim.calculate_first_and_second_legendre(filename);
  }
  if (choice == 3) {
    worksim.calculate_deuterium_order_parameter();
    worksim.calculate_first_and_second_legendre(filename);
  }
  if (choice == 4) {
    worksim.perl_formatter();
  }

}
