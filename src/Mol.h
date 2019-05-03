/*
 * Mol.h
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#ifndef MOL_H_
#define MOL_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include <cmath>

using namespace std;

/*!
 *
 */
class Mol {
public:

// variaveis

	//! Number of atoms
	int N ;
	//! Number of residues
	int Nres;
	//! Number of atomtypes
	int Natomtypes;
	//! Charges of the atoms. Must be divided by 18.2223 to get electron charges
	vector<double>charges;
	//! Atomic masses for the system.
	vector<double>masses;
	//! Atom names for the system according to AMBER FF.
	vector<string>amberatoms;
	//! Atomic parameters
	vector<string> atomtypes_prm;
	//! Atomic names
	vector<string> atomnames;
	//! Atomic radii
	vector<double>radius;
	//! Welldepth parameters to each atom according to AMBER FF
	vector<double>welldepth;
	//! Atomic coordinates
	vector<vector<double> > xyz;
	//! Atomic coordinates after rotation/translation
	vector<vector<double> > new_xyz;
	//! Atomic coordinates from the last step accepted.
	vector<vector<double> > old_xyz;
	//! Atomic parameter for LJ computation
	vector<double>epsilons;
	//! Squared root of atoms epsilons
	vector<double> epsilons_sqrt;
	//! Atomic radii
	vector<double>radii;
	//! Mass value for an individidual atom
	double mass;
	//! Atomic charge
	double charge;
	//! C++ vector with the names of the residues
	vector<string> resnames;
	//! Pointers to the number of residues / number of atoms.
	vector<int> residue_pointer;
	//! Temporary string to parse prmtop file.
	string line;



//funcoes internas

	//! Class constructor
	//! @param prmtop PRMTOP file with molecule parameters;
	//! @param inpcrd INPCRD file with molecule coordinates;
	Mol(ifstream &prmtop, ifstream &inpcrd);

	//! Function to parse the number of atoms from the PRMTOP file
	//! @param prmtop PRMTOP file with molecule parameters;
	void get_N(ifstream &prmtop);

	//! Function to get the atomic names as defined in the PRMTOP file
	//! @param prmtop PRMTOP file with molecule parameters;
	//! @param N Number of atoms in the molecule;
	void get_atomnames(ifstream &prmtop, int N);

	//! Function to get the atomic charges
	//! @param prmtop PRMTOP file with molecule parameters;
	//! @param N Number of atoms in the molecule;
	void get_charges(ifstream &prmtop, int N);

	//! Function to get the atomic masses
	//! @param prmtop PRMTOP file with molecule parameters;
	//! @param N Number of atoms in the molecule;
	void get_masses(ifstream &prmtop, int N);

	//! Function to get the names of the residues in the PRMTOP file
	//! @param prmtop PRMTOP file with molecule parameters;
	void get_resnames(ifstream &prmtop);

	//! Function to parse the number of atoms in each residue;
	//! @param prmtop PRMTOP file with molecule parameters;
	void get_res_pointers(ifstream &prmtop);

	//! Function to print the atomic charges (unused)
	//! @param charges C++ vector (\a N-sized) with molecule atomic charges;
	//! @param N Number of atoms in the molecule;
	void print_charge(vector<double> charges, int N);

	//! Function to get the AMBER atom types in the PRMTOP file
	//! @param prmtop PRMTOP file with molecule parameters;
	//! @param N Number of atoms in the molecule;
	void get_amberatoms(ifstream &prmtop, int N);

	//! Function to parse the atomic parameter from AMBER vdw.param file
	void read_atomtypes_prm();

	//! Function to parse the epsilon parameter used for VDW computation
	//! @param atomtypes_prm VDW parameters parsed from \a vdw.param file;
	//! @param amberatoms C++ vector with atom names as defined in AMBER force fields;
	//! @param welldepth VDW values for welldepths according to AMBER force field.
	void get_epsilon(vector<string>atomtypes_prm, vector<string>amberatoms, vector<double>welldepth);

	//! Function to parse atomic radii used in VDW computation
	//! @param atomtypes_prm VDW parameters parsed from \a vdw.param file;
	//! @param amberatoms C++ vector with atom names as defined in AMBER force fields;
	//! @param radius VDW values for radii according to AMBER force field.
	void get_radius(vector<string>atomtypes_prm, vector<string>amberatoms, vector<double>radius);

	//! Function to get atomic coordinates from INPCRD file.
	//! @param inpcrd INPCRD file with molecule coordinates;
	//! @param N Number of atoms in the molecule;
	void get_xyz(ifstream &inpcrd, int N);

}; //EOClass


#endif /* MOL_H_ */
