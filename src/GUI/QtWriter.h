/*
 * QtWriter.h
 *
 *  Created on: 22/07/2011
 *      Author: Nascimento
 */

#ifndef QTWRITER_H_
#define QTWRITER_H_

#include "../Parser.h"
#include "../Mol.h"
#include "../Mol2.h"
#include <QPlainTextEdit>
#include <QString>
#include<QApplication>
#include"zlib.h"
#include<stdio.h>

using namespace std;

class QtWriter {
public:
	//! logfile is defined as a C file to generate a written report (log) of the program run.
	FILE *logfile;
	//! Output file with aligned molecules in format SYBYL MOL2.
	gzFile outmol2;
	//! Pointer to a QPlainTextEdit class. This object is used to show the log in the screen (GUI)
	QPlainTextEdit* Editor;
	//! Pointer to Mol2 Object
	Mol *Cmol;

	//! Class constructor. Opens the logfile with the "output" name defined in the input file.
	//! @param Input Pointer to the Parser class.
	QtWriter(Parser *Input, QPlainTextEdit* Ed);

	//! This function writes down the parameters given by the user as input.
	//! @param Pointer to the Parser class
	void write_params(Parser *Input);
	//! This function writes a "welcome" message
	void write_welcome(void);
	//! This message shows a char array in the screen and writes it to the logfile.
	//! @param  info char array (82 char long) with a message to be written/shown.
	void write_to_log(char info[82]);
	//! This overloaded function writes a blank line in the logfile and in the screen.
	void write_to_log(void);
	//! Class destructor. Closes the logfile.

	/*! Method to write a PDB file (not used anymore).
	 * @param Cmol Pointer to a Mol2 Object with molecule parameters
	 * @param xyz C++ vector with molecule coordinates
	 * @param outname string for file output
	 */
	void write_pdb(Mol *Cmol, vector<vector<double> >xyz, string outname);

	/*!
	 * Class to write a PDB file (overloaded)
	 * @param Cmol Pointer to a Mol2 Object with molecule parameters
	 * @param xyz C++ vector with molecule coordinates
	 * @param outname string for file output
	 */
	void write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, string outname);

	/*!
	 * Class used to write Mol2 files
	 * @param Cmol Pointer to a Mol2 object with molecule parameters
	 * @param xyz C++ vector with molecule atomic coordinates
	 * @param cc Correlation coefficient written in mol2 file
	 */
	void writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double cc);

	/*!
	 * Deconstructor. Closes opened files.
	 */
	~QtWriter(void);
};

#endif /* QTWRITER_H_ */
