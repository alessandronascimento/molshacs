/*
 * Parser.h
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#ifndef PARSER_H_
#define PARSER_H_

/*!
 * @brief This class handles the input information provided by the user or by the GUI.
 */

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstdlib>

#ifdef HAS_GUI
#include<QStringList>
#endif

using namespace std;

/*!
 * MolShaCS was designed to parse a text file with the (minimum) required information for a comparison computation.
 * The class PARSER handles this input. Here, all the required information is kept and used during the program execution.
 * There are two typical ways of using the class and/or using the program. First of all, the information provided in the
 * program GUI is processed and kept in this class. Alternatively, the user may provide a input file with the necessary information
 * for batch processing, for example. The file is processed and the necessary information is stored in this class.
 *
 * An important issue regarding the class (and the program) is that it was designed to abort the program execution if a incorrect
 * information is provided in the input file.
 */
class Parser {

public:
	//! Reference Molecule PRMTOP filename
	string refmol_prmtop;
	//! Reference Molecule Coordinates filename
	string refmol_inpcrd;
	//! Reference Molecule SYBYL MOL2 file
	string refmol_mol2;
	//! Mol2 file has AMBER (GAFF) atom types.
	bool mol2_aa;
	//! Parameter reading during parsing
	string param;
	//! Samping scheme during electrostactic potential computation
	double sampling;
	//! Angle step during rotational matching
	double step;
	//! Actually, half the box size, i.e., the distance from COM to the box edges in each direction.
	double box_size;
	//! File containing the comparing molecule file names.
	string molsfile;
	//! Whether the multifile comparison should be done or not.
	bool multi;
	//! Prefix for all the output files
	string output;
	//! Temporary string....
	string tmp;
	//! Number of minimization steps.
	int nsteps;
	//! Tolerance to gradient minimization.
	double tol;
	//! Delta for derivative computation in minimization
	double delta;
	//! Chooses the algorithm to spatial minimization
	string minimizer;
	//! Decides wether re-oriented coordinates will or not be written
	bool write_pdb;
	//! Scale of VDW term in minimization
	double vdw_scale;
	//! Scale of Electrostatic term in minimization
	double elec_scale;
	//! Timeout (in secods) for molecule minimization
	int timeout;
	//! Defines wether the algorithm should or should not align molecules first.
	bool align;
	//! Softcore term for VDW energy. Avoids numerical instability when r->0;
	double deltaij6;
	//! Softcore term for electrostatic energy. Avoids numerical instability when r->0;
	double deltaij_es6;
	//! Whether or not the program should use a threshold for writing atomic coordinates of the compared molecules.
	bool use_write_coord_threshold;
	//! Threshold for writing atomic coordinates of the compared molecules.
	float write_coord_threshold;
	//!
	bool macromol;
#ifdef HAS_GUI
	//! Qt String list with the list of comparing molecules
	QStringList comparing_molecules;
#endif



	//! Class constructor
	//! @param infile Text file with parameters to be parsed with program parameters
	Parser(char* infile);

	Parser(void);
	//! Function to get a individual parameter;
	//! @param infile Text file with parameters to be parsed with program parameters
	void set_parameters(char* infile);

	//! Function to compare the parameter token from input file and compare with known parameters.
	//! @param param Text string parsed in the \a infile with a program keyword
	//! @param input Text file with parameters to be parsed with program parameters
	void comparing(string param, ifstream &input);
};


#endif /* PARSER_H_ */
