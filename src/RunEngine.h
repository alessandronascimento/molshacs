/*
 * RunEngine.h
 *
 *  Created on: 19/07/2011
 *      Author: Nascimento
 */

#ifndef RUNENGINE_H_
#define RUNENGINE_H_

/*!
 * @brief This class is the main engine of the code. Here, the basic functions are called for computation. The method is always called.
 * Can be called from GUI or after an input file is provided.
 */

#include "Parser.h"
#include <ctime>
#include<math.h>
#include<stdio.h>
#include<vector>
#include "Mol2.h"
#include "Grid.h"
#include "CORREL.h"
#include "Minimizer2.h"
#include "Writer.h"
#include "Parser.h"
#include "Gaussian.h"
#include "zlib.h"

#ifdef HAS_GUI

#include "GUI/QtWriter.h"
#include <QProgressBar>
#include <QPlainTextEdit>

#endif


using namespace std;
/*!
 * The RunEngine is the class where all the computation is done. It is always called from main.cpp.
 */
class RunEngine {
public:

	//! Distance
	double r;
	//! Electrostatic potential at a point
	double elect;
	//! Alpha angle
	double a,
	//! Beta angle
	b,
	//! Gamma angle
	g;
	//! Translation step in X direction
	double transx,
	//! Translation step in Y direction
	transy,
	//! Translation step in Z direction
	transz;
	//! String to parse the molecules (when running in the multimode) from the multifile.
	string indmol;
	//! Counter for the molecules in the multimode file.
	int count;
	//! Character array to set up a output name for the PDB file given as output
	char outname[50];
	//! Variable to get the timing of the execution.
	clock_t time0;
	//! Variable to get the timing of the execution.
	clock_t time1;
	//! Variable to get the timing of the execution.
	double total_time;
	//! Minimum value reached after minimization procedure.
	double f_min;
	//! Char array to print running information in the screen.
	char info[82];
	//! Correlation coefficient computed to each molecule after spatial matching.
	double cc;
	//!
	double Vab;

// FUNCTIONS

	RunEngine();

	/*!
	 * Function distance. Computes the distance between two atoms.
	 * @param x1 X-coordinate of the first atom.
	 * @param x2 X-coordinate of the second atom.
	 * @param y1 Y-coordinate of the first atom.
	 * @param y2 Y-coordinate of the second atom.
	 * @param z1 Z-coordinate of the first atom.
	 * @param z2 Z-coordinate of the second atom.
	 * @return the distance (double) between the atoms.
	 */
	double distance (double x1, double x2, double y1, double y2, double z1, double z2);

	/*!
	 * Function run. This function actually runs the program. It's where all the computation
	 * takes place.
	 * @param Input An Parser object with the parameters taken from user.
	 */
	void run(Parser Input);

	/*!
	 * Function ranker. This function writes a vector of integers with the indexes
	 * of the SI vector sorted decreasingly.
	 * @param CCs C++ vector with molecule SIs;
	 * @return A vector of integers with SIs indexes sorted.
	 */
	vector<int> ranker(vector<float> CCs);

#ifdef HAS_GUI
	/*!
	 * Function run (overloaded). This function actually runs the program. It's where all
	 * the computation takes place. This overloaded function was written to use with the
	 * Qt GUI. For this reason, takes some others parameters.
	 * @param Input A Parser object with user-entered parameters
	 * @param Ed A pointer to the QPlainTextEdit object, used to upgrade the log in the GUI.
	 * @param progressbar A pointer to he QProgressBar object used here to update the progress bar during
	 * program execution.
	 *
	 */
	void run(Parser Input, QPlainTextEdit* Ed, QProgressBar* progressbar);
#endif

};

#endif /* RUNENGINE_H_ */
