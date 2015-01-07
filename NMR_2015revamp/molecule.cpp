/*
 * molecule.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */
#include "molecule.h"
  trajectory* molecule_bond::getMoleculeTrajectory() {
    return molecule_trajectory;
  }

  std::string molecule_bond::getMoleculeType() {
    return molecule_type;
  }

  std::string molecule_bond::getName() {
    return name;
  }

  int molecule_bond::getNumberOfBonds() {
    return number_of_bonds;
  }
  molecule_bond::molecule_bond(std::string to_name, std::string type, int numb_of_bonds, trajectory* trajectory){
    name=to_name;
    molecule_type=type;
    number_of_bonds=numb_of_bonds;
    molecule_trajectory=trajectory;
  }
  molecule_bond::molecule_bond(){
    name="null";
    molecule_type="null";
    number_of_bonds=0;
    molecule_trajectory=new trajectory();
  }


