/*
 * Mol2.cpp
 *
 *  Created on: 10/10/2011
 *      Author: Nascimento
 */

#include "Mol2.h"

using namespace std;

Mol2::Mol2(Parser *Input, ifstream &mol2file) {

	if (mol2file.is_open()){
		getline(mol2file, line); 		//@<TRIPOS>MOLECULE
		getline(mol2file, line);		// Molecule name

		mol2file >> this->N;
		mol2file >> this->line; // number of bonds
		mol2file >> this->Nres;

		while (line != "@<TRIPOS>ATOM"){
			getline(mol2file, line);
		}

		double tx, ty, tz;
		vector<double> txyz;
		int tres;
		double tcharge;
		int count=0;
		this->read_atomtypes_prm();

		for (int i=0; i<N; i++){
			mol2file >> line; // atom number (i+1);

			mol2file >> line;
			this->atomnames.push_back(line);

			mol2file >> tx >> ty >> tz;
			txyz.push_back(tx);
			txyz.push_back(ty);
			txyz.push_back(tz);
			this->xyz.push_back(txyz);
			txyz.clear();

			mol2file >> line; //atom type (SYBYL)
			if (Input->mol2_aa){
				this->amberatoms.push_back(line);
/*
 * If atomtypes are given as gaff types, we *still* need to get VDW parameter
 * in the vdw.param file.
 * For SYBYL atom types there is no need to read the file, since the parameters
 * are parsed in the convert2gaff method.
 */

				this->get_epsilon(this->atomtypes_prm,this->amberatoms[i],this->welldepth);
				this->get_radius(this->atomtypes_prm, this->amberatoms[i], this->radius);
			}
			else {
				this->amberatoms.push_back(this->convert2gaff(line));
			}

			this->get_masses(amberatoms[i]);

			mol2file >> tres; // residue number;
			if (tres > count){
				this->residue_pointer.push_back(i+1);
				count = tres;
				mol2file >> line;
				this->resnames.push_back(line);
			}
			else {
				mol2file >> line; // residue name
			}

			mol2file >> tcharge;
			this->charges.push_back(tcharge);
		}
	}
	else {
		printf("Mol2 file could not be opened! Please check!\n");
		exit(1);
	}
}

Mol2::Mol2(Parser *Input, string mol2) {
	if (!Input->macromol){
		FILE* mol2file;
		mol2file=fopen(mol2.c_str(), "r");
		int tint;
		float tx, ty, tz;
		vector<double> txyz;
		int tres;
		float tcharge;
		int count=0;
		char tatomtype[10];
		char resname[10];
		if (Input->mol2_aa){
//			this->read_atomtypes_prm();
			this->get_gaff_parameters();
		}

		if (mol2file !=NULL){
			fgets(str, 80, mol2file); 	//@<TRIPOS>MOLECULE
			//		fgets(str, 80, mol2file); 	//Molecule name
			fscanf(mol2file, "%s", str);
			this->molname = str;
			fscanf(mol2file, "%d %d %d %d %d", &this->N, &this->Nbonds, &this->Nres, &tint, &tint);

#ifdef DEBUG
			printf("Number of atoms: %d\n", this->N);
#endif
			str[0] = '#';
			while (str[0] != '@'){
				fgets(str, 80, mol2file);
			}
			for (int i=0; i<this->N; i++){
				fscanf(mol2file, "%d %s %f %f %f %s %d %s %f\n", &tint, str, &tx, &ty, &tz, tatomtype, &tres, resname, &tcharge);

				txyz.push_back(tx);
				txyz.push_back(ty);
				txyz.push_back(tz);
				this->xyz.push_back(txyz);
				txyz.clear();

				this->charges.push_back(tcharge);
				this->atomnames.push_back(str);

				if (Input->mol2_aa){
					this->amberatoms.push_back(tatomtype);
					this->get_epsilon(this->atomtypes_prm,this->amberatoms[i],this->welldepth);
					this->get_radius(this->atomtypes_prm, this->amberatoms[i], this->radius);
				}
				else{
					this->amberatoms.push_back(this->convert2gaff(string(tatomtype)));
					this->sybyl_atoms.push_back(string(tatomtype));
				}

				this->get_masses(amberatoms[i]);

				if (tres > count){
					this->residue_pointer.push_back(i+1);
					count = tres;
					this->resnames.push_back(string(resname));
				}
			}
			fscanf(mol2file, "%s\n", str);
			if (str[0] != '@'){
				while (str[0] != '@'){
					fgets(str, 80, mol2file);
				}
			}

			vector<string> bond;
			char s1[6], s2[6], s3[5];
			for (int i=0; i<this->Nbonds; i++){
				fscanf(mol2file, "%d%s%s%s\n", &tint, s1, s2, s3);
//				cout << tint << " " << string(s1) << " " << string(s2) << " " << string(s3) << endl;
				bond.push_back(string(s1));
				bond.push_back(string(s2));
				bond.push_back(string(s3));
				this->bonds.push_back(bond);
//				printf("%6d%6.6s%6.6s%5.5s\n", i+1, bonds[i][0].c_str(),bonds[i][1].c_str(), bonds[i][2].c_str());
				bond.clear();
			}
		}

		else {
			printf("Mol2 file %s could not be opened! Please check!\n", mol2.c_str());
			exit(1);
		}
		fclose(mol2file);
	}



	else {
		this->Parse_PDB(mol2);
	}
}

int Mol2::Parse_PDB(string molfile){
	FILE* pdbin;
	pdbin = fopen(molfile.c_str(), "r");
	char keyword[6];
	int atomnumber;
	char atomname[4];
	char resname[4];
	char chain[2];
	int resnumber;
	float x, y, z, occ, bf;
	char atomtype[1];
	char look4chain[1];
	look4chain[0] = 'A';
	if (pdbin != NULL){
		str[0]='#';
		while (str[0] != 'A' or str[1] !='T' or str[2] !='O' or str[3] != 'M'){
			fgets(str, 80, pdbin);
		}
		/* Important: We lost the coordinates of the first atom. Usually the N for the
		 * residue 1
		 */
		vector<double> xyz_temp;
		int respointer =1;
		while (!feof(pdbin)){
			fscanf(pdbin, "%s %d %s %s %s %d %f %f %f %f %f %s", keyword, &atomnumber, atomname, resname, chain, &resnumber,
					&x, &y, &z, &occ, &bf, atomtype);
			string resn = string(resname);
			if (resn == "GLU" or resn == "ASP" or resn == "LYS" or resn == "ARG"){
				if (atomname[0] == 'C' and atomname[1] == 'A' and chain[0] == look4chain[0]){
					//				printf("%d %s %d %s %s %f %f %f\n", atomnumber, atomname, resnumber, chain, resname, x, y, z);
					xyz_temp.push_back(x);
					xyz_temp.push_back(y);
					xyz_temp.push_back(z);
					this->xyz.push_back(xyz_temp);
					xyz_temp.clear();
					this->radii.push_back(1.9080);
					this->epsilons.push_back(0.0860);
					this->charges.push_back(this->choose_charge(string(resname)));
					this->masses.push_back(12.01);
					this->residue_pointer.push_back(respointer);
					respointer++;
					this->resnames.push_back(string(resname));
					this->atomnames.push_back(string(atomname));

				}
			}
		}
		this->N = int(this->xyz.size());
	}
	return 0;
}

double Mol2::choose_charge(string resname){
	if (resname== "GLU" or resname== "ASP"){
			return -1.00;
	}
	else if (resname== "ARG" or resname== "LYS"){
		return 1.00;
	}
	else {
		return 0.00;
	}
}

string Mol2::convert2gaff(string atom){
	string gaff_atom;
	if (atom == "C.3"){
		gaff_atom = "c3";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.1094);
		this->epsilons_sqrt.push_back(sqrt(0.1094));
	}

	else if (atom =="C.2"){
		gaff_atom = "c2";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="C.1"){
		gaff_atom = "c1";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));

	}

	else if (atom =="C.ar"){
		gaff_atom = "ca";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="C.cat"){
		gaff_atom = "c";
		this->radii.push_back(1.9080);
		this->epsilons.push_back(0.0860);
		this->epsilons_sqrt.push_back(sqrt(0.0860));
	}

	else if (atom =="N.3"){
		gaff_atom = "n3";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.2"){
		gaff_atom = "n2";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.1"){
		gaff_atom = "n1";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.ar"){
		gaff_atom = "nh";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.am"){
		gaff_atom = "n";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.4"){
		gaff_atom = "n4";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="N.pl3"){
		gaff_atom = "na";
		this->radii.push_back(1.8240);
		this->epsilons.push_back(0.1700);
		this->epsilons_sqrt.push_back(sqrt(0.1700));
	}

	else if (atom =="O.3"){
		gaff_atom = "oh";
		this->radii.push_back(1.7210);
		this->epsilons.push_back(0.2104);
		this->epsilons_sqrt.push_back(sqrt(0.2104));
	}

	else if (atom =="O.2"){
		gaff_atom = "o";
		this->radii.push_back(1.6612);
		this->epsilons.push_back(0.2100);
		this->epsilons_sqrt.push_back(sqrt(0.2100));
	}

	else if (atom =="O.co2"){
		gaff_atom = "o";
		this->radii.push_back(1.6612);
		this->epsilons.push_back(0.2100);
		this->epsilons_sqrt.push_back(sqrt(0.2100));
	}

	else if (atom =="O.spc" or atom == "O.t3p"){
		gaff_atom = "ow";
		this->radii.push_back(1.7683);
		this->epsilons.push_back(0.1520);
		this->epsilons_sqrt.push_back(sqrt(0.1520));
	}

	else if (atom =="S.3"){
		gaff_atom = "sh"; //not sure... sh or ss
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.2"){
		gaff_atom = "s2";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.O" or atom == "S.o"){
		gaff_atom = "s4";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="S.O2" or atom == "S.o2"){
		gaff_atom = "s6";
		this->radii.push_back(2.0000);
		this->epsilons.push_back(0.2500);
		this->epsilons_sqrt.push_back(sqrt(0.2500));
	}

	else if (atom =="P.3"){
		gaff_atom = "p3";
		this->radii.push_back(2.1000);
		this->epsilons.push_back(0.2000);
		this->epsilons_sqrt.push_back(sqrt(0.2000));
	}

	else if (atom =="F"){
		gaff_atom = "f";
		this->radii.push_back(1.75);
		this->epsilons.push_back(0.0610);
		this->epsilons_sqrt.push_back(sqrt(0.0610));
	}

	else if (atom =="H"){
		gaff_atom = "hc";
		this->radii.push_back(0.600);
		this->epsilons.push_back(0.0157);
		this->epsilons_sqrt.push_back(sqrt(0.0157));
	}

	else if (atom =="H.spc" or atom=="H.t3p"){
		gaff_atom = "hw";
		this->radii.push_back(0.0000);
		this->epsilons.push_back(0.0000);
		this->epsilons_sqrt.push_back(0.0000);
	}

	else if (atom =="Cl"){
		gaff_atom = "cl";
		this->radii.push_back(1.948);
		this->epsilons.push_back(0.2650);
		this->epsilons_sqrt.push_back(sqrt(0.2650));
	}

	else if (atom =="Br"){
		gaff_atom = "br";
		this->radii.push_back(2.22);
		this->epsilons.push_back(0.3200);
		this->epsilons_sqrt.push_back(sqrt(0.3200));
	}

	else if (atom =="I"){
		gaff_atom = "i";
		this->radii.push_back(2.35);
		this->epsilons.push_back(0.4000);
		this->epsilons_sqrt.push_back(sqrt(0.4000));
	}

	else{
		printf("Atom type %s not found among GAFF parameters.\nPlease check Mol2.h source file.\n", atom.c_str());
		exit(1);
	}

	return(gaff_atom);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Mol2::read_atomtypes_prm() {
	char atom[3];
	float rad, well;
	char vdw[100];

	elsa_dir_path = getenv("MOLSHACS_DIR");
	if (elsa_dir_path != NULL){
		strcpy(vdw, elsa_dir_path);
		strcat(vdw, "/vdw.param");
		FILE* vdwprm;
		vdwprm = fopen(vdw, "r");
		if (vdwprm != NULL){
			fgets(str, 80, vdwprm);
			while(!feof(vdwprm)){
				fscanf(vdwprm, "%s %f %f", atom, &rad, &well);
				this->atomtypes_prm.push_back(string(atom));
				this->radius.push_back(rad);
				this->welldepth.push_back(well);
			}
			fclose(vdwprm);
		}
	}

//		ifstream vdwprm(vdw);
//		if (vdwprm.is_open()){
//			getline(vdwprm, line);
//			while (!vdwprm.eof()) {
//				vdwprm >> atom >> rad >> well;
//				this->atomtypes_prm.push_back(atom);
//				this->radius.push_back(rad);
//				this->welldepth.push_back(well);
//			}
//			vdwprm.close();
//		}

	else {																		// "MOLSHACS_DIR is not defined
//		cout << "Environment variable MOLSHACS_DIR is not defined." << endl;
//		cout << "Looking for vdw.param file in the local folder..." << endl;
/*		ifstream vdwprm("vdw.param");
		if (vdwprm.is_open()){
			getline(vdwprm, line);
			while (!vdwprm.eof()) {
				vdwprm >> atom >> rad >> well;
				this->atomtypes_prm.push_back(atom);
				this->radius.push_back(rad);
				this->welldepth.push_back(well);
			}
			vdwprm.close();
		}
*/
		FILE* vdwprm;
		vdwprm = fopen("vdw.param", "r");
		if (vdwprm != NULL){
			fgets(str, 80, vdwprm);
			while(!feof(vdwprm)){
				fscanf(vdwprm, "%s %f %f", atom, &rad, &well);
				this->atomtypes_prm.push_back(string(atom));
				this->radius.push_back(rad);
				this->welldepth.push_back(well);
			}
			fclose(vdwprm);
		}
		else {
			printf("Could not open vdw.param file %s. Please check.\n", vdw);
			exit(1);
		}
	}
}

void Mol2::get_epsilon(vector<string>atomtypes_prm, string amberatom, vector<double>welldepth) {
	for (unsigned i=0; i<atomtypes_prm.size(); i++) {
		if (atomtypes_prm[i] == amberatom) {
			this->epsilons.push_back(welldepth[i]);
			this->epsilons_sqrt.push_back(sqrt(welldepth[i]));
		}
	}
}

void Mol2::get_radius(vector<string>atomtypes_prm, string amberatom, vector<double>radius) {
	bool test = false;
	for (unsigned i=0; i<atomtypes_prm.size(); i++) {
		if (atomtypes_prm[i] == amberatom) {

#ifdef DEBUG
			printf("Atom %s == %s. Radii: %.5f\n", amberatom.c_str(), atomtypes_prm[i].c_str(), radius[i]);
#endif

			this->radii.push_back(radius[i]);
			test = true;
		}
	}
	if (! test){
		printf("Radius for atom %s was not found! Please check....\n", amberatom.c_str());
		exit(1);
	}
}

void Mol2::get_masses(string atomname){
	if(atomname.substr(0,1) == "h"){
		this->masses.push_back(1.008);
	}
	else if (atomname.substr(0,1) == "c" and atomname.substr(1,1) != "l"){ // carbon but not clorine
		this->masses.push_back(12.01);
	}
	else if (atomname.substr(0,1) == "n"){
		this->masses.push_back(14.01);
	}
	else if (atomname.substr(0,1) == "o"){
		this->masses.push_back(16.00);
	}
	else if (atomname.substr(0,1) == "p"){
		this->masses.push_back(30.97);
	}
	else if (atomname.substr(0,1) == "s"){
		this->masses.push_back(32.06);
	}
	else if (atomname.substr(0,1) == "f"){
		this->masses.push_back(19.00);
	}
	else if (atomname.substr(0,2) == "cl"){
		this->masses.push_back(32.06);
	}
	else if (atomname.substr(0,2) == "br"){
		this->masses.push_back(79.90);
	}
	else if (atomname.substr(0,1) == "i"){
		this->masses.push_back(126.9);
	}
	else {
		printf("Could not find atomic mass for atom %s\n. Please check Mol2.cpp file.\n", atomname.c_str());
		exit(1);
	}
}

void Mol2::get_gaff_parameters(){
	this->atomtypes_prm.clear();
	this->radius.clear();
	this->welldepth.clear();

// Reading GAFF atomtypes;
	this->atomtypes_prm.push_back("h1");		//0
	this->atomtypes_prm.push_back("h2");		//1
	this->atomtypes_prm.push_back("h3");		//2
	this->atomtypes_prm.push_back("h4");		//3
	this->atomtypes_prm.push_back("h5");		//4
	this->atomtypes_prm.push_back("ha");		//5
	this->atomtypes_prm.push_back("hc");		//6
	this->atomtypes_prm.push_back("hn");		//7
	this->atomtypes_prm.push_back("ho");		//8
	this->atomtypes_prm.push_back("hp");		//9
	this->atomtypes_prm.push_back("hs");		//10
	this->atomtypes_prm.push_back("hw");		//11
	this->atomtypes_prm.push_back("hx");		//12
	this->atomtypes_prm.push_back("o");			//13
	this->atomtypes_prm.push_back("oh");		//14
	this->atomtypes_prm.push_back("os");		//15
	this->atomtypes_prm.push_back("ow");		//16
	this->atomtypes_prm.push_back("c");			//17
	this->atomtypes_prm.push_back("cz");		//18
	this->atomtypes_prm.push_back("c1");		//19
	this->atomtypes_prm.push_back("c2");		//20
	this->atomtypes_prm.push_back("c3");		//21
	this->atomtypes_prm.push_back("ca");		//22
	this->atomtypes_prm.push_back("cc");		//23
	this->atomtypes_prm.push_back("cd");		//24
	this->atomtypes_prm.push_back("ce");		//25
	this->atomtypes_prm.push_back("cf");		//26
	this->atomtypes_prm.push_back("cg");		//27
	this->atomtypes_prm.push_back("ch");		//28
	this->atomtypes_prm.push_back("cp");		//29
	this->atomtypes_prm.push_back("cq");		//30
	this->atomtypes_prm.push_back("cu");		//31
	this->atomtypes_prm.push_back("cv");		//32
	this->atomtypes_prm.push_back("cx");		//33
	this->atomtypes_prm.push_back("cy");		//34
	this->atomtypes_prm.push_back("n");			//35
	this->atomtypes_prm.push_back("n1");		//36
	this->atomtypes_prm.push_back("n2");		//37
	this->atomtypes_prm.push_back("n3");		//38
	this->atomtypes_prm.push_back("n4");		//39
	this->atomtypes_prm.push_back("na");		//40
	this->atomtypes_prm.push_back("nb");		//41
	this->atomtypes_prm.push_back("nc");		//42
	this->atomtypes_prm.push_back("nd");		//43
	this->atomtypes_prm.push_back("ne");		//44
	this->atomtypes_prm.push_back("nf");		//45
	this->atomtypes_prm.push_back("nh");		//46
	this->atomtypes_prm.push_back("no");		//47
	this->atomtypes_prm.push_back("s");			//48
	this->atomtypes_prm.push_back("s2");		//49
	this->atomtypes_prm.push_back("s4");		//50
	this->atomtypes_prm.push_back("s6");		//51
	this->atomtypes_prm.push_back("sx");		//52
	this->atomtypes_prm.push_back("sy");		//53
	this->atomtypes_prm.push_back("sh");		//54
	this->atomtypes_prm.push_back("ss");		//55
	this->atomtypes_prm.push_back("p2");		//56
	this->atomtypes_prm.push_back("p3");		//57
	this->atomtypes_prm.push_back("p4");		//58
	this->atomtypes_prm.push_back("p5");		//59
	this->atomtypes_prm.push_back("pb");		//60
	this->atomtypes_prm.push_back("px");		//61
	this->atomtypes_prm.push_back("py");		//62
	this->atomtypes_prm.push_back("f");			//63
	this->atomtypes_prm.push_back("cl");		//64
	this->atomtypes_prm.push_back("br");		//65
	this->atomtypes_prm.push_back("i");			//66
	this->atomtypes_prm.push_back("DU");		//67

// Reading GAFF radii. Same order as above
	this->radius.push_back(1.3870);
	this->radius.push_back(1.2870);
	this->radius.push_back(1.1870);
	this->radius.push_back(1.4090);
	this->radius.push_back(1.3590);
	this->radius.push_back(1.4590);
	this->radius.push_back(1.4870);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(0.6000);
	this->radius.push_back(1.1000);
	this->radius.push_back(1.6612);
	this->radius.push_back(1.7210);
	this->radius.push_back(1.6837);
	this->radius.push_back(1.7683);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.9080);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(1.8240);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.0000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(2.1000);
	this->radius.push_back(1.75);
	this->radius.push_back(1.948);
	this->radius.push_back(2.22);
	this->radius.push_back(2.35);
	this->radius.push_back(0.0000);

// Reading GAFF welldepth. Same order as above
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0150);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0000);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.0000);
	this->welldepth.push_back(0.0157);
	this->welldepth.push_back(0.2100);
	this->welldepth.push_back(0.2104);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1520);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1094);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.0860);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.1700);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2500);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.2000);
	this->welldepth.push_back(0.0610);
	this->welldepth.push_back(0.2650);
	this->welldepth.push_back(0.3200);
	this->welldepth.push_back(0.4000);
	this->welldepth.push_back(0.0000);

	bool fail = true;
	if ((this->atomtypes_prm.size() == this->radius.size()) and (this->atomtypes_prm.size() == this->welldepth.size())){
		fail = false;
	}
	if (fail){
		printf("Size of vectors with GAFF parameters does not match. Please check Mol2.cpp file!\n");
		exit(1);
	}
}
