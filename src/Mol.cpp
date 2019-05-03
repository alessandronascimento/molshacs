/*
 * Mol.cpp
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#include "Mol.h"

Mol::Mol(ifstream &prmtop, ifstream &inpcrd) {
	if( (prmtop.is_open()) and (inpcrd.is_open()) ){
	this->get_N(prmtop);
	this->get_atomnames(prmtop, Mol::N);
	this->get_charges(prmtop, Mol::N);
	this->get_masses(prmtop, Mol::N);
	this->get_resnames(prmtop);
	this->get_res_pointers(prmtop);
	this->get_amberatoms(prmtop, Mol::N);
	this->read_atomtypes_prm();
	this->get_xyz(inpcrd, Mol::N);
	this->get_epsilon(Mol::atomtypes_prm,Mol::amberatoms,Mol::welldepth);
	this->get_radius(Mol::atomtypes_prm, Mol::amberatoms, Mol::radius);
	}
	else {
		printf("Could not open prmtop/inpcrd file.\n");
		printf("Please, check!\n");
		printf("Exiting...\n");
		exit(1);
	}
}

void Mol::get_N(ifstream &prmtop) {
	if (prmtop.is_open()){
		getline (prmtop, line);
		for (int i=1; i<=5; i++) {
			getline (prmtop, line); }
		prmtop >> Mol::N >> Mol::Natomtypes;
		for (int i=1; i<=9; i++){
			prmtop >> line;
		}
		prmtop >> Mol::Nres;
		getline (prmtop, line);
	}
	else {
		cout << "Couldn't open prmtop file;" << endl;
		exit(1);
	}
}

void Mol::get_atomnames(ifstream &prmtop, int N){
	string name, name2;
	getline (prmtop, line);
	while (line.size() < 6 or line.substr(6,9) != "ATOM_NAME")  {
		getline(prmtop, line);
	}
	getline(prmtop, line);
	int i=0;
	while (i<N){
		prmtop >> name;
		if (name.size() > 4){
			while (name.size()>4) {
				name2=name.substr(0,4);
				this->atomnames.push_back(name2);
				i++;
				name=name.substr(4);
			}
			this->atomnames.push_back(name);
			i++;
		}
		else {
			this->atomnames.push_back(name);
			i++;
		}
	}
}

void Mol::get_charges(ifstream &prmtop, int N) {
	double charge;
	getline (prmtop, line);
  	while (line.size() < 6 or line.substr(6,6) != "CHARGE")  {
  		getline (prmtop, line); }
  	getline (prmtop, line);		// FORMAT
  	for (int i=0; i<N; i++) {
    		prmtop >> charge;
    		this->charges.push_back(charge); }
}


void Mol::get_masses(ifstream &prmtop, int N) {
	double mass;
	getline (prmtop, line);
  	while (line.size() < 6 or line.substr(6,4) != "MASS")  {
  		getline (prmtop, line); }
  	getline (prmtop, line);		// FORMAT
  	for (int i=1; i<=N; i++) {
    		prmtop >> mass;
    		this->masses.push_back(mass); }
  	getline (prmtop, line);
  	getline (prmtop, line);
}

void Mol::get_resnames(ifstream &prmtop){
	getline (prmtop, line);
	string resname;
	while (line.size() < 6 or line.substr(6,13) != "RESIDUE_LABEL")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
	for (int i=1; i<= this->Nres; i++){
		prmtop >> resname;
		this->resnames.push_back(resname);
	}
	getline (prmtop, line);
}

void Mol::get_res_pointers(ifstream &prmtop){
	getline(prmtop, line);
	int respointer;
	while (line.size() < 6 or line.substr(6,15)!= "RESIDUE_POINTER")  {
		getline (prmtop, line);
	}
	getline (prmtop, line);		// FORMAT
	for (int i=0; i < this->Nres; i++){
		prmtop >> respointer;
		this->residue_pointer.push_back(respointer);
	}
	getline(prmtop, line);
}

void Mol::print_charge(vector<double>charges, int N) {
	cout << "charge: " << charges[N] << endl; }

void Mol::get_amberatoms(ifstream &prmtop, int N) {
	string atom;
	getline (prmtop, line);
	while (! prmtop.eof()) {
		getline (prmtop, line);
		if (line.size() >= 15) {
			if (line.substr(6,15) == "AMBER_ATOM_TYPE")  {
				getline (prmtop, line);		// FORMAT
				for (int i=1; i<=N; i++) {
					prmtop >> atom;
					this->amberatoms.push_back(atom); }
			}
		}
	}
}

void Mol::read_atomtypes_prm() {
	string atom;
	double rad, well;
	ifstream vdwprm("vdw.param");
	if (vdwprm.is_open()){
		getline(vdwprm, line);
		while (!vdwprm.eof()) {
			vdwprm >> atom >> rad >> well;
			this->atomtypes_prm.push_back(atom);
			this->radius.push_back(rad);
			this->welldepth.push_back(well); } }
	else {
		printf("Could not open vdw.param file. Please check.");
		exit(1);
	}
}

void Mol::get_epsilon(vector<string>atomtypes_prm, vector<string>amberatoms, vector<double>welldepth) {
	for (unsigned i=0; i<amberatoms.size(); i++) {
		for (unsigned j=0; j<atomtypes_prm.size(); j++) {
			if (amberatoms[i]== atomtypes_prm[j]) {
				this->epsilons.push_back(welldepth[j]);
				this->epsilons_sqrt.push_back(sqrt(welldepth[j]));
			}
		}
	}
}

void Mol::get_radius(vector<string>atomtypes_prm, vector<string>amberatoms, vector<double>radius) {
	for (unsigned i=0; i<amberatoms.size(); i++) {
		for (unsigned j=0; j<atomtypes_prm.size(); j++) {
			if (amberatoms[i]== atomtypes_prm[j]) {
				this->radii.push_back(radius[j]);
			}
		}
	}
}

void Mol::get_xyz(ifstream &inpcrd, int N) {
	double x, y, z;
	vector<double> xyz_tmp;
	getline (inpcrd, line);
	getline (inpcrd, line);
	for (int n=1; n<=N; n++) {
		xyz_tmp.clear();
		inpcrd >> x >> y >> z;
		xyz_tmp.push_back(x);
		xyz_tmp.push_back(y);
		xyz_tmp.push_back(z);
		this->xyz.push_back(xyz_tmp);
	}
}
