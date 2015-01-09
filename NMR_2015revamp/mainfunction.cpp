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
  coords_3d* derp = new coords_3d(
      0.0, 0.0, 0.0);
  std::cout << derp->getX() << " " << derp->getY() << " " << derp->getZ()
            << std::endl;
  derp->setX(
      0.1);
  derp->setY(
      0.2);
  derp->setZ(
      0.3);
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
  fileinput* test = new fileinput(
      "somefile.xyz", "test", "test");
  test->is_file_valid();
  fileinput* test2 = new fileinput(
      "somefile.zyx", "test", "test");
  test2->is_file_valid();
  delete test;
  delete test2;
  fileinput* test3 = new fileinput(
      "test.xyz", "test", "test");
  test3->is_file_valid();
  std::cout << test3->find_atom_count() << std::endl;
  std::cout << test3->get_line_number() << std::endl;
  delete test3;
  //Test loadmolecules() for "DOPC-C22.xyz"
  fileinput* test4 = new fileinput(
      "DOPC-C22.xyz", "C22", "dopc");
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
//          << test_mol_bond.at(i).getMoleculeTrajectory()->get_trajectory_array()->size()
//          << std::endl;
//  }
  return 1;
}
int simulation_test() {
  fileinput* test_file = new fileinput(
      "DOPC-C22.xyz", "C22", "dopc");
  std::vector<molecule_bond>* test_mol_bond;
  test_mol_bond = test_file->loadmolecules();
  simulation * newsim;
  newsim = new simulation(
      test_mol_bond);
  std::cout << "number of molecules " << newsim->get_number_of_molecules()
            << std::endl;
  std::cout << "number of frames " << newsim->get_number_of_frames()
            << std::endl;
  newsim->calculate_deuterium_order_parameter();
  newsim->calculate_first_and_second_legendre();
  return 1;
}
void debug_test(){
  fileinput* test_file = new fileinput(
       "debug_test.xyz", "debug", "debug");
   std::vector<molecule_bond>* test_mol_bond;
   test_mol_bond = test_file->loadmolecules();
   simulation * newsim;
   newsim = new simulation(
       test_mol_bond);
   std::cout << "number of molecules " << newsim->get_number_of_molecules()
             << std::endl;
   std::cout << "number of frames " << newsim->get_number_of_frames()
             << std::endl;
   newsim->calculate_deuterium_order_parameter();
   newsim->calculate_first_and_second_legendre();
}

int main(void) {
//  coords_3d_test();
//  vector_3d_test();
//  fileinput_test();
  simulation_test();
//  debug_test();
}
