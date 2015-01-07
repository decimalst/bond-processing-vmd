/*
 * molecule.h
 *
 *  Created on: Dec 22, 2014
 *      Author: byronl
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_
#include <string>
#include "trajectory.h"
class molecule_bond{
    std::string molecule_type,name;
    int number_of_bonds;
    trajectory* molecule_trajectory;
 public:
    molecule_bond(std::string,std::string,int,trajectory*);
    molecule_bond();
    trajectory* getMoleculeTrajectory();
    std::string getName();
    std::string getMoleculeType();
    int getNumberOfBonds();
};




#endif /* MOLECULE_H_ */
