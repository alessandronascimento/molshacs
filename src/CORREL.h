/*
 * CORREL.h
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#ifndef CORREL_H_
#define CORREL_H_

/**
 * @brief 	Class CORREL. This class handles some important function such as rotation and translation of a molecule. Here
 * a rotation is defined by a a rotation around the three Euler's angles alpha, beta and gamma.
 * @author Alessandro S. Nascimento
 * @date September, 2012
 */

#include<string>
#include<vector>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include "zlib.h"
#include "ElSA.h"

using namespace std;

/*!
 * Class CORREL. This class handles some important function such as rotation and translation of a molecule. Here
 * a rotation is defined by a a rotation around the three Euler's angles alpha, beta and gamma.
 * @author Alessandro S. Nascimento
 * @date September, 2012
 */
class CORREL {
public:

	//! Correlation index
	double correl;
	//! Electrostatic potential at a point
	double elec;
	//! Distance
	double r;
	//! X coordinate
	double x,
	//! Y coordinate
	y,
	//! Z Coordinate
	z;
	//! New coordinates after movement
	vector<vector<double> > new_xyz;
	//! Who knows?
	vector<double> temp;
	//!
	double sum_product;

// Methods
	/*!
	 * Default constructor. Does nothing.
	 */
	CORREL();

	//! Method for Rotation AND Translation
	//! @param xyz C++ vector of vector that stores molecule coordinates
	//! @param N Number of atoms in the molecule
	//! @param alpha Euler alpha angle
	//! @param beta Euler beta angle
	//! @param gamma Euler gamma angle
	//! @param transx X-shift in Angstrons
	//! @param transy Y-shift in Angstrons
	//! @param transz Z-shift in Angstrons
	vector<vector<double> > rototranslate(vector<vector<double> > xyz, int N, double alpha, double beta, double gamma, double transx, double transy, double transz);

	//! Method to compute the correlation index
	//! @param gridcoords C++ vector with the coordinates used to compute the electrostatic grid
	//! @param elec_potential C++ vector with the electrostatic grid
	//! @param N Number of atoms
	//! @param xyz C++ vector of a vector with atom coordinates
	//! @param charges C++ vector with atomic charges
	double compute_correl(vector<vector<double> > gridcoords, vector<double> elec_potential, int N, vector<vector<double> >xyz, vector<double> charges);

	//! Coordinates translation
	//! @param xyz C++ vector of a vector with atom coordinates
	//! @param N Number of atoms
	//! @param tx X-shift
	//! @param ty Y-shift
	//! @param tz Z-shift
	vector<vector<double> > translate(vector<vector<double> > xyz, int N, double tx, double ty, double tz);

	//! Method to write a PDB file of the transformed coordinates
	//! @param xyz C++ vector of a vector with atom coordinates
	//! @param N Number of atoms
	//! @param atomnames C++ vector of Strings with the names of the atoms as defined in PRMTOP file.
	//! @param name String with name of for PDB file output
	void write_pdb(vector<vector<double> >xyz, int N, vector<string> atomnames, string name);

	/*! Scalar product of two electrostatic potentials used to compute the Hodgkin index
	 * (Hodgkin, E.E., Richards, W.G. Int. J. Quant. Chem. Quant. Bio. Symp. (1987) 14:105-110.
	 * Cited in Blomberg N. et al. Proteins (1999) 37:379-387.
	 */
	//! @param p1 C++ vector with electrostatic grid of the first molecule
	//! @param p2 C++ vector with electrostatic grid of the second molecule
	double compute_correl_wade(vector<double> p1, vector<double> p2);

	//! Method to compute the electrostatic correlation as defined by Wade et al.
	//! @param elec_potential Electrostatic potential pre-computed for the reference molecule;
	//! @param p2 C++ vector with electrostatic grid of the second molecule
	//! @param sum_of_elec_potential Sum of the elements of the \a elec_potential vector. Useful to speed up computations.
	double compute_correl_wade(vector<double> elec_potential, vector<double> p2, double sum_of_elec_potential);
};


#endif /* CORREL_H_ */
