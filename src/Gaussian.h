/*
 * Gaussian.h
 *
 *  Created on: 18/11/2011
 *      Author: Nascimento
 */

#ifndef GAUSSIAN_H_
#define GAUSSIAN_H_

/**
 * @brief The Gaussian class is used to compute the Gaussian-based decriptors for molecular overlay and similarity computation.
 * @author Alessandro S. Nascimento
 * @date September, 2012.
 */

#include<stdio.h>
#include<iostream>
#include<cmath>
#include <ctime>
#include <vector>
#include"Mol2.h"
#include"ElSA.h"
#include"Parser.h"

using namespace std;

/*!
 * This class has the important function compute_shape_and_charge, where the descriptors are computed for
 * overlay or similarity analysis.
 */
class Gaussian {
public:

	//! Parameter for Gaussian amplitude
	double pi;
	//! Parameter for Gaussian amplitude
	double pj;
	//! Distance among two atoms
	double dij2;
	//!
	double alpha_i;
	//!
	double alpha_j;
	//!
	clock_t time0, time1;

	/*!
	 *
	 */
	Gaussian(void);
	/*!
	 *
	 */
	Gaussian(Mol2* RefMol, Mol2* CompMol, double* Vab);
	/**
	 *
	 */
	Gaussian(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, double* Vab);
	/*!
	 *
	 */
	double compute_si(Mol2* RefMol, Mol2* CompMol);
	/*!
	 *
	 */
	double compute_si(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_si_pos_charges(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_si_neg_charges(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/**
	 *
	 */
	double compute_shape_density(Mol2* RefMol, Mol2* CompMol);
	/*!
	 *
	 */
	double compute_shape_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_shape_density(Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_positive_charge_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);

	/*!
	 *
	 */
	double compute_positive_charge_density(Mol2* CompMol, vector<vector<double> > xyz);

	/*!
	 *
	 */
	double compute_negative_charge_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_negative_charge_density(Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_shape_and_charge_density(Parser *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz);
	/*!
	 *
	 */
	double compute_shape_and_charge_density(Parser *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, vector<vector<double> > cxyz);
	/*!
	 *
	 */
	double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);
	/*!
	 *
	 */
	virtual ~Gaussian();

};

#endif /* GAUSSIAN_H_ */
