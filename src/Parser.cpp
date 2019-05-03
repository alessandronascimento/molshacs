/*
 * Parser.cpp
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#include "Parser.h"

using namespace std;

void Parser::comparing (string param, ifstream &input) {
	if(param == "refmol_mol2"){
		input >> Parser::refmol_mol2;
	}
	else if (param == "mol2_aa"){
		input >> tmp;
		if (tmp == "yes" or tmp=="YES" or tmp == "Yes"){
			Parser::mol2_aa = true;
		}
	}
	else if (param == "sampling") {
		input >> Parser::sampling;
	}
	else if (param == "step") {
		input >> Parser::step;
	}
	else if (param == "box_size") {
		input >> Parser::box_size;
		Parser::box_size = Parser::box_size/2.00;
	}
	else if (param == "multimode"){ // not necessary
		input >> tmp;
		if (tmp == "yes" or tmp=="YES" or tmp == "Yes"){
			Parser::multi = true;
		}
	}
	else if (param == "multimol"){
		input >> Parser::molsfile;
		this->multi = true;
	}
	else if (param == "output_prefix"){
		input >> Parser::output;
	}
	else if (param == "ntries"){
		input >> Parser::nsteps;
	}
	else if (param == "tol"){
		input >> Parser::tol;
	}
	else if (param == "delta"){
		input >> Parser::delta;
	}
	else if (param == "minimizer"){
		input >> Parser::minimizer;
	}
	else if (param == "vdw_scale"){
		input >> Parser::vdw_scale;
	}
	else if (param == "elec_scale"){
		input >> Parser::elec_scale;
	}
	else if (param == "align_molecules"){
		input >> tmp;
		if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
			Parser::align = true;
		}
		else {
			Parser::align = false;
		}
	}
	else if (param == "timeout"){
		input >> Parser::timeout;
	}
	else if (param == "deltaij6"){
		input >> Parser::deltaij6;
	}
	else if (param == "deltaij_es6"){
		input >> Parser::deltaij_es6;
	}
	else if (param == "write_coordinates"){
		input >> tmp;
		if (tmp == "yes" or tmp == "YES" or tmp == "Yes"){
			Parser::write_pdb = true;
		}
		else {
			Parser::write_pdb = false;
		}
	}

	else if (param == "write_coord_threshold"){
		Parser::use_write_coord_threshold = true;
		input >> Parser::write_coord_threshold;
	}
	else if (param == "macromol"){
		input >> tmp;
		if (tmp == "yes" or tmp == "YES" or tmp =="Yes"){
			this->macromol = true;
		}
		else {
			this->macromol = false;
		}
	}
	else {
		cout << "Unknown parameter: " << param << endl;
		exit(1);
	}
}

void Parser::set_parameters(char* infile){
	ifstream input(infile);
	while (!input.eof()){
		input >> param;
		Parser::comparing (param, input);
	}
	input.close();
}

Parser::Parser(char* infile){
	Parser::mol2_aa = false;
	Parser::multi = false;
	Parser::align = true;
	Parser::write_pdb = false;
	Parser::use_write_coord_threshold = false;
	Parser::step = 1.0;
	Parser::box_size = 12.5;
	Parser::nsteps = 50;
	Parser::tol = 0.1;
	Parser::delta = 1E-5;
	Parser::minimizer = "nlopt_mma";
	Parser::elec_scale = 1.0;
	Parser::vdw_scale = 1.0;
	Parser::output = "spca_out";
	Parser::timeout = 120;
	Parser::deltaij6 = 122.978496247489; //2,23^6
	Parser::deltaij_es6 = 28.722900390625;
	this->macromol = false;
	Parser::set_parameters(infile);
}

Parser::Parser(void){
	Parser::mol2_aa = false;
	Parser::multi = false;
	Parser::align = true;
	Parser::write_pdb = false;
	Parser::use_write_coord_threshold = false;
	Parser::step = 1.0;
	Parser::box_size = 12.5;
	Parser::nsteps = 50;
	Parser::tol = 0.1;
	Parser::delta = 1E-5;
	Parser::minimizer = "nlopt_mma";
	Parser::elec_scale = 1.0;
	Parser::vdw_scale = 1.0;
	Parser::output = "spca_out";
	Parser::timeout = 120;
	Parser::deltaij6 = 122.978496247489; //2,23^6
	Parser::deltaij_es6 = 28.722900390625;
	this->macromol = false;
}
