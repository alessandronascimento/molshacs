/*
 * Minimizer2.h
 *
 *  Created on: 03/03/2011
 *      Author: Nascimento
 */

#ifndef MINIMIZER2_H_
#define MINIMIZER2_H_

#include<gsl/gsl_multimin.h>
#include"gsl/gsl_errno.h"
#include"gsl/gsl_math.h"
#include"gsl/gsl_vector.h"
#include<string.h>
#include<stdio.h>
#include<vector>
#include<cmath>
#include"Mol2.h"
#include"Parser.h"
#include"Grid.h"
#include "Gaussian.h"
#include "nlopt.h"
#include "nlopt.hpp"
#include<time.h>


/*!
 * Class Minimizer2. This is the second version of the Minimizer Class.
 * This Class has buit-in functions to handle the minimization of a objective function
 * f that has as parameters three rotation angles (Euler angles) and translation in cartesian axis
 * x, y and z. The function f(a,b,g,x,y,z) is minimized using different algorithms provided here
 * by both GNU Scientific Library (GSL) and (LGPL) NLOPT (Steven G. Johnson, The NLopt nonlinear-optimization package,
 * http://ab-initio.mit.edu/nlopt).
 */

class Minimizer2 {

public:

	//!
	static double alpha_i;
	//!
	static double alpha_j;
	//! Squared distance among atoms i and j
	static double dij2;
	//! Amplitude of the Gaussian for atom i
	static double pi;
	//! Amplitude of the Gaussian for atom j
	static double pj;
	//! Pointer to the MOL Class
	static Mol2 *Cmol;
	//! Pointer to the GRID Class
	static Grid *Cgrid;
	//! Pointer to the INPUT Class
	static Parser *Input;
	//! Pointer to the reference molecule object (type Mol2)
	static Mol2* RefMol;
	//! Pointer to the comparison molecule object
	static Mol2* CompMol;
	//! x, y and z coordinates
	static double x, y, z;
	//! Rototranslated coordinates.
	static vector<vector<double> >new_xyz;
	//! Objective function value
	static double t1;
	//! Minimum of the objective function reached after minimization.
	static double f_minimum;

	//! Class Constructor.
	//! @param _Cmol MOL Class pointer passed as a parameter
	//! @param _Input PARSER Class pointer passed as a parameter
	//! @param _CGrid GRID Class pointer passed as parameter
	Minimizer2(Mol2 *_RefMol, Mol2* _CompMol, Parser *_Input, Grid *_Cgrid);

	//! This function computes the distance among two atoms.
	/*!
	 * \f$(d_12)\f$ = \f$\sqrt{(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$
	 * @param x1 X coordinate of the first atom
	 * @param x2 X coordinate of the second atom
	 * @param y1 Y coordinate of the first atom
	 * @param y2 Y coordinate of the second atom
	 * @param z1 Z coordinate of the first atom
	 * @param z2 Z coordinate of the second atom
	 * \return The distance.
	 */
	static double dist(double x1, double x2, double y1, double y2, double z1, double z2);

	//! This function computes the square of the distance among two atoms.
	/*!
	 * \f$(d_12)\f$ = \f${(x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2}\f$
	 * @param x1 X coordinate of the first atom
	 * @param x2 X coordinate of the second atom
	 * @param y1 Y coordinate of the first atom
	 * @param y2 Y coordinate of the second atom
	 * @param z1 Z coordinate of the first atom
	 * @param z2 Z coordinate of the second atom
	 * \return The square of the distance.
	 */
	static double dist_squared(double x1, double x2, double y1, double y2, double z1, double z2);

	/*! Computes the center of mass of a molecular system
	 *  @param xyz C++ vector with molecular coordinates
	 *  @param Cmol Pointer to a Mol2 Object with Molecular information (massses, charges, etc)
	 *  @return A C++ vector with three elements (x, y and z) for the center of mass.
	 */
	static vector<double> compute_com(vector<vector<double> > xyz, Mol2 *Cmol);

	//! Computes new coordinates after a translation and a rotation.
	/*!
	 * @param Cmol MOL Class pointer with molecule paramters
	 * @param xyz Molecule coordinates
	 * @param alpha Euler alpha angle
	 * @param beta Euler beta angle
	 * @param gamma Euler gamma angle
	 * @param transx Shift in X direction
	 * @param transy Shift in Y direction
	 * @param transz Shift in Z direction
	 * @return A new vector of vector with new coordinates.
	 */
	static vector<vector<double> > rototranslate(Mol2 *Cmol, vector<vector<double> >xyz, double alpha, double beta, double gamma, double transx, double transy, double transz);


/*! This function computes the objective function for GSL optimizers. But is not
 * being used for now.
 * @param v GSL vector with the six variables(transx, transy, transz, alpha, beta and gamma)
 * @param params Coefficients for computation (void)
 * @return The objective function value
 */
	static double function_gaussian_shape(const gsl_vector *v, void *params);


	//! Minimizes the objective function using the method of moving asymptotes as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_mma();

	//! Minimizes the objective function using Improved Stochastic Ranking Evolution Strategy as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_isres();

	//! Minimizes the objective function using Subplex as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_subplex();

	//! Minimizes the objective function using Simplex as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_simplex();
	//! Minimizes the objective function using COBYLA as implemented in NLOPT
	//! @return Returns nothing.
	double minimize_nlopt_cobyla();

	//! Minimizes the objective function using BFGS as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_bfgs2();

	//! Minimizes the objective function using STOGO as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_stogo();

	//! Minimizes the objective function using DIRECT-L as implemented in NLOPT
	//! @return the minimum value of the objective function reached
	double minimize_nlopt_direct();

	/*! @brief Minimizes the objective function using local AUGLAG as implemented in NLOPT
	 * without using derivatives (LN)
	 * @return the minimum value of the objective function reached
	 */
	double minimize_nlopt_ln_auglag();

	/*! Objective function for Gaussian Shape and Charge computation
	 * @param x Rotation/Translation variables;
	 * @param grad Gradients
	 * @param data Coefficients (not used)
	 * @return The objective function values
	 */
	static double nlopt_func_gauss_density(const std::vector<double> &x, std::vector<double> &grad, void *data);
};

#endif /* MINIMIZER2_H_ */
