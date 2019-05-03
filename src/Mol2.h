/*
 * Mol2.h
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#ifndef MOL2_H_
#define MOL2_H_


/*!
 * @brief This cless is designed to handle the information in a SYBYL MOL2 file.
 * @author Alessandro S. Nascimento
 * @date September, 2012.
 */

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<string.h>
#include<cstdio>
#include<cstdlib>
#include<cmath>
#include"Parser.h"

using namespace std;

/*!
 *
 */
class Mol2 {
public:


// Variables

	//! Number of atoms
	int N ;
	//! Number of residues
	int Nres;
	//! Number of atomtypes
	int Natomtypes;
	//! Number of bonds
	int Nbonds;
	//! Name for the molecule
	string molname;
	//! Charges of the atoms. Must be divided by 18.2223 to get electron charges
	vector<double>charges;
	//! Atomic masses for the system.
	vector<double>masses;
	//! Atom names for the system according to AMBER FF.
	vector<string>amberatoms;
	//! Atomic parameters
	vector<string>atomtypes_prm;
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
	//! C++ vector with the names of the residues
	vector<string> resnames;
	//! Pointers to the number of residues / number of atoms.
	vector<int> residue_pointer;
	//! Temporary string to parse prmtop file.
	string line;
	//! Temporary C string
	char str[80];
	//! String to get the environment variable
	char* elsa_dir_path;
	//! Keeps Vaa for RefMol/CompMol
	double self_obj_function;
	//! Atom names in SYBYL types
	vector<string> sybyl_atoms;
	//! Bonds defining molecular topology. Used to write MOL2 files
	vector<vector<string> >bonds;


	/*!
	 * Initializer. This class has, as arguments, a pointer to the class PARSER.
	 * The class uses some information given by the user there. The filename is also
	 * given as argument. This makes possible to use the same object to reference and
	 * comparing molecules.
	 */
	Mol2(Parser *Input, ifstream &mol2file); // not used anymore

	/*!
	 * Initializer. This method gets a Parser object and a C++ string with the MOL2 filename.
	 * This overloaded function is actually used throughout the code.
	 */
	Mol2(Parser *Input, string mol2);

	/*!
	 * This method is used to manually convert SYBYL atom types to AMBER (GAFF) atom
	 * types. The conversion does not need to be very accurate since the VDW parameters
	 * are the only parameters used here and only for alignment purposes.
	 */
	string convert2gaff(string atom);

	/*!
	 * This method parses the file "vdw.param" to get GAFF atomic VDW parameters. DEPRECATED.
	 */
	void read_atomtypes_prm();

	/*!
	 * This method initializes the GAFF atomic parameters without the need of a additional
	 * "vdw.param" file. This version is used in the current version of the program.
	 */
	void get_gaff_parameters();

	/*!
	 * This method reads the atomic types in the molecule and attributes GAFF vdw parameters
	 * to each atom.
	 * @param atomtypes_prm C++ vector with GAFF atom types;
	 * @param amberatom Atom type for atom "i"
	 * @param amberatom C++ vector with epsilons (welldepths) for GAFF atoms in the same order as in the
	 * atomtypes_prm vector.
	 */
	void get_epsilon(vector<string>atomtypes_prm, string amberatom, vector<double>welldepth);

	/*
	 * This method reads the atomic types in the molecule and attributes GAFF vdw parameters
	 * to each atom.
	 * @param atomtypes_prm C++ vector with GAFF atom types;
	 * @param amberatom Atom type for atom "i"
	 * @param radius C++ vector with atomic radii for GAFF atoms in the same order as in the
	 * atomtypes_prm vector.
	 */
	void get_radius(vector<string>atomtypes_prm, string amberatom, vector<double>radius);

	/*!
	 * This method attributes atomic masses to each of the atoms in the molecule. It is useful for
	 * the center of mass computation.
	 * @param atomname Atomic type.
	 */
	void get_masses(string atomname);

	//!
	int Parse_PDB(string molfile);

	//!
	double choose_charge(string resname);
};

#endif /* MOL2_H_ */
