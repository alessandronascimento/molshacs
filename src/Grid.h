/*
 * Grid.h
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#ifndef GRID_H_
#define GRID_H_

#include<iostream>
#include<fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include"Mol.h"
#include"Mol2.h"
#include"Parser.h"

using namespace std;

/*!
 *
 */
class Grid {
public:

// variaveis

	//! Center of the computation box
	vector<double>center;
	//! Center of Mass
	vector<double>com;
	//! Lower limit for computation in X direction;
	double min_x,
	//! Lower limit for computation in Y direction;
	min_y,
	//! Lower limit for computation in Z direction;
	min_z,
	//! Upper limit for computation in X direction;
	max_x,
	//! Upper limit for computation in Y direction;
	max_y,
	//! Upper limit for computation in Z direction;
	max_z;
	//! Electrostactic correlation computed according to the method proposed by Rebecca Wade;
	double elec;
	//! X component of the Center of Mass;
	double comx;
	//! Y component of the Center of Mass;
	double comy;
	//! Z component of the Center of Mass;
	double comz;
	//! Sum of the atomic masses in the molecule;
	double sumofmasses;
	//! Distance \a r among two atoms;
	double r;
	//! C++ vector to store the electrostatic grid.
	vector<double> grid;
	//! Pointer to molecule PRMTOP class.
	Mol *Cmol;
	//! Pointer to Parser class.
	Parser *Input;


// Funcoes internas

	//! Class constructor. No actual function.
	Grid();

	//! Computes the distance among two atoms.
	//! @param x1 X coordinate of the first atom
	//! @param y1 Y coordinate of the first atom
	//! @param z1 Z coordinate of the first atom
	//! @param x2 X coordinate of the second atom
	//! @param y2 Y coordinate of the second atom
	//! @param Z2 Z coordinate of the second atom
	double dist(double x1, double x2, double y1, double y2, double z1, double z2);

	//! Computes the square of the distance among two atoms.
	//! @param x1 X coordinate of the first atom
	//! @param y1 Y coordinate of the first atom
	//! @param z1 Z coordinate of the first atom
	//! @param x2 X coordinate of the second atom
	//! @param y2 Y coordinate of the second atom
	//! @param Z2 Z coordinate of the second atom
	double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);

	//! Method t compute the center of mass of a molecule
	//! @param xyz XYZ C++ vector of a vector with molecule coordinates
	//! @param masses C++ vector with atomic masses
	vector<double>compute_com(vector<vector<double> >xyz, vector<double>masses);

	//! Method to compute the electrostactic grid
	//! @param xyz XYZ C++ vector of a vector with molecule coordinates
	//! @param charges C++ vector with atomic charges
	//! @param min_x Lower limit in X direction for electrostactic computation
	//! @param min_y Lower limit in Y direction for electrostactic computation
	//! @param min_z Lower limit in Z direction for electrostactic computation
	//! @param max_x Upper limit in X direction for electrostactic computation
	//! @param max_y Upper limit in Y direction for electrostactic computation
	//! @param max_z Upper limit in Z direction for electrostactic computation
	//! @param sampling defines the spacing among grid points in each direction
	//! @param N Number of atoms in the molecule
	vector<double> compute_grid(vector<vector<double> > xyz, vector<double> charges, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, double sampling, int N);

	/*!
	 * Method for grid computation. Not used anymore.
	 */
	vector<double> compute_grid(vector<vector<double> > xyz, vector<double> charges, double sampling, int N);

	/*!
	 * Method for grid computation. Not used anymore.
	 */
	vector<double> compute_grid(Mol* Cmol, Parser* Input);

	/*!
	 * Method for grid computation. Not used anymore.
	 */
	vector<double> compute_grid(Mol2* Cmol, Parser* Input);

	/*!
	 * Method for grid computation. Not used anymore.
	 */
	vector<double> compute_grid(vector<vector<double> > xyz, Mol* Cmol, Parser* Input);

	/*!
	 * Method for grid computation. Not used anymore.
	 */
	vector<double> compute_grid(vector<vector<double> > xyz, Mol2* Cmol, Parser* Input);

	//! Method to write a RCSB PDB file with the box that defines the grid
	//! @param min_x Lower limit in X direction for electrostactic computation
	//! @param min_y Lower limit in Y direction for electrostactic computation
	//! @param min_z Lower limit in Z direction for electrostactic computation
	//! @param max_x Upper limit in X direction for electrostactic computation
	//! @param max_y Upper limit in Y direction for electrostactic computation
	//! @param max_z Upper limit in Z direction for electrostactic computation
	void write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z);

	//! Method to compute the VDW repulsive grid
	//! @param Cmol pointer to PRMTOP class
	//! @param Input pointer to Parser class
	vector<double> compute_vdw_repulsive_grid(Mol *Cmol, Parser *Input);

	/*!
	 * Method for repulsive grid computation. Not used anymore.
	 */
	vector<double> compute_vdw_repulsive_grid(Mol2 *Cmol, Parser *Input);

	/*!
	 * Method for repulsive grid computation. Not used anymore.
	 */
	vector<double> compute_vdw_repulsive_grid(Mol *Cmol, Parser *Input, vector<vector<double> > xyz);

	/*!
	 * Method for repulsive grid computation. Not used anymore.
	 */
	vector<double> compute_vdw_repulsive_grid(Mol2 *Cmol, Parser *Input, vector<vector<double> > xyz);


};


#endif /* GRID_H_ */
