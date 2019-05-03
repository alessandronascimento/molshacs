/*
 * Grid.cpp
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#include "Grid.h"


Grid::Grid() {
}
double Grid::dist(double x1, double x2, double y1, double y2, double z1, double z2) {
	return (sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) ); }

double Grid::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ((((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

vector<double> Grid::compute_com(vector<vector<double> >xyz, vector<double>masses) {
	com.clear();
	sumofmasses=0.0;
	comx=0.0;
	comy=0.0;
	comz=0.0;
	for (unsigned i=0; i<masses.size(); i++) {
		comx= comx+(xyz[i][0]*masses[i]);
		comy= comy+(xyz[i][1]*masses[i]);
		comz= comz+(xyz[i][2]*masses[i]);
		sumofmasses=sumofmasses+masses[i]; }
	comx=comx/sumofmasses;
	comy=comy/sumofmasses;
	comz=comz/sumofmasses;
	com.push_back(comx);
	com.push_back(comy);
	com.push_back(comz);
	return(com);
}

vector<double> Grid::compute_grid(vector<vector<double> > xyz, vector<double> charges, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, double sampling, int N){
	grid.clear();
	for (double a=min_x; a<=max_x; a= a+sampling){
		for (double b=min_y; b<=max_y; b= b+sampling){
			for (double c=min_z; c<=max_z; c= c+sampling){
				elec=0.00;
				for (int i=0; i<N; i++){
					double r2 = Grid::dist_squared(a, xyz[i][0], b, xyz[i][1], c, xyz[i][2]);
					double rij6 = r2*r2*r2;
					elec+= (charges[i]) / (pow((rij6+28.722900390625), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

vector<double> Grid::compute_grid(vector<vector<double> > xyz, vector<double> charges, double sampling, int N){
	grid.clear();
	for (double a=this->min_x; a<=this->max_x; a= a+sampling){
		for (double b=this->min_y; b<=this->max_y; b= b+sampling){
			for (double c=this->min_z; c<=this->max_z; c= c+sampling){
				elec=0.00;
				for (int i=0; i<N; i++){
					double r2 = Grid::dist_squared(a, xyz[i][0], b, xyz[i][1], c, xyz[i][2]);
					double rij6 = r2*r2*r2;
					elec+= (charges[i]) / (pow((rij6+28.722900390625), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

vector<double> Grid::compute_grid(Mol* Cmol, Parser* Input){
	grid.clear();
	for (double a=this->min_x; a<=this->max_x; a= a+Input->sampling){
		for (double b=this->min_y; b<=this->max_y; b= b+Input->sampling){
			for (double c=this->min_z; c<=this->max_z; c= c+Input->sampling){
				elec=0.00;
				for (int i=0; i<Cmol->N; i++){
					double r2 = Grid::dist_squared(a, Cmol->xyz[i][0], b, Cmol->xyz[i][1], c, Cmol->xyz[i][2]);
					double rij6 = r2*r2*r2;

					elec+= (Cmol->charges[i]) / (pow((rij6+Input->deltaij_es6), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

vector<double> Grid::compute_grid(Mol2* Cmol, Parser* Input){
	grid.clear();
	for (double a=this->min_x; a<=this->max_x; a= a+Input->sampling){
		for (double b=this->min_y; b<=this->max_y; b= b+Input->sampling){
			for (double c=this->min_z; c<=this->max_z; c= c+Input->sampling){
				elec=0.00;
				for (int i=0; i<Cmol->N; i++){
					double r2 = Grid::dist_squared(a, Cmol->xyz[i][0], b, Cmol->xyz[i][1], c, Cmol->xyz[i][2]);
					double rij6 = r2*r2*r2;
					elec+= (Cmol->charges[i]) / (pow((rij6+Input->deltaij_es6), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

vector<double> Grid::compute_grid(vector<vector<double> > xyz, Mol* Cmol, Parser* Input){
	grid.clear();
	for (double a=this->min_x; a<=this->max_x; a= a+Input->sampling){
		for (double b=this->min_y; b<=this->max_y; b= b+Input->sampling){
			for (double c=this->min_z; c<=this->max_z; c= c+Input->sampling){
				elec=0.00;
				for (int i=0; i<Cmol->N; i++){
					double r2 = Grid::dist_squared(a, xyz[i][0], b, xyz[i][1], c, xyz[i][2]);
					double rij6 = r2*r2*r2;
					elec+= (Cmol->charges[i]) / (pow((rij6+Input->deltaij_es6), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

vector<double> Grid::compute_grid(vector<vector<double> > xyz, Mol2* Cmol, Parser* Input){
	grid.clear();
	for (double a=this->min_x; a<=this->max_x; a= a+Input->sampling){
		for (double b=this->min_y; b<=this->max_y; b= b+Input->sampling){
			for (double c=this->min_z; c<=this->max_z; c= c+Input->sampling){
				elec=0.00;
				for (int i=0; i<Cmol->N; i++){
					double r2 = Grid::dist_squared(a, xyz[i][0], b, xyz[i][1], c, xyz[i][2]);
					double rij6 = r2*r2*r2;
					elec+= (Cmol->charges[i]) / (pow((rij6+Input->deltaij_es6), (1.0/3.0)));
				}
				grid.push_back(elec);
			}
		}
	}
	return(grid);
}

void Grid::write_box(vector<double>center, double min_x, double min_y, double min_z, double max_x, double max_y, double max_z){

	FILE *box;
	box = fopen("box.pdb", "w");

	fprintf (box, "REMARK    CENTER OF THE BOX  %2.3f  %2.3f  %2.3f\n", center[0], center[1], center[2]);
	fprintf (box, "ATOM      1  DUA BOX     1     % 2.3f % 2.3f % 2.3f\n", min_x, min_y, min_z);
	fprintf (box, "ATOM      2  DUB BOX     1     % 2.3f % 2.3f % 2.3f\n", max_x, min_y, min_z);
	fprintf (box, "ATOM      3  DUC BOX     1     % 2.3f % 2.3f % 2.3f\n", max_x, min_y, max_z);
	fprintf (box, "ATOM      4  DUD BOX     1     % 2.3f % 2.3f % 2.3f\n", min_x, min_y, max_z);
	fprintf (box, "ATOM      5  DUE BOX     1     % 2.3f % 2.3f % 2.3f\n", min_x, max_y, min_z);
	fprintf (box, "ATOM      6  DUF BOX     1     % 2.3f % 2.3f % 2.3f\n", max_x, max_y, min_z);
	fprintf (box, "ATOM      7  DUG BOX     1     % 2.3f % 2.3f % 2.3f\n", max_x, max_y, max_z);
	fprintf (box, "ATOM      8  DUH BOX     1     % 2.3f % 2.3f % 2.3f\n", min_x, max_y, max_z);
	fprintf (box, "CONECT    1    2    4    5\n");
	fprintf (box, "CONECT    2    1    3    6\n");
	fprintf (box, "CONECT    3    2    4    7\n");
	fprintf (box, "CONECT    4    1    3    8\n");
	fprintf (box, "CONECT    5    1    6    8\n");
	fprintf (box, "CONECT    6    2    5    7\n");
	fprintf (box, "CONECT    7    3    6    8\n");
	fprintf (box, "CONECT    8    4    5    7\n");
	fprintf (box, "END\n");

	fclose(box);
}

vector<double> Grid::compute_vdw_repulsive_grid(Mol *Cmol, Parser *Input){
	vector<double> repulsive_potential;
	double repul;
	double r2, r6;
	for (double x=Grid::min_x; x<=Grid::max_x; x+=Input->sampling){
		for (double y=Grid::min_y; y<=Grid::max_y; y+=Input->sampling){
			for (double z=Grid::min_z; z<=Grid::max_z; z+=Input->sampling){
				repul=0.0;
				for (int i=0; i<Cmol->N; i++){
					r2=Grid::dist_squared(Cmol->xyz[i][0], x,Cmol->xyz[i][1], y, Cmol->xyz[i][2], z);
					r6=r2*r2*r2;
					repul+= ((Cmol->epsilons_sqrt[i])*(pow((Cmol->radii[i]), 12)))/(pow((r6+Input->deltaij6),2));
				}
				repulsive_potential.push_back(repul);
			}
		}
	}
	return(repulsive_potential);
}

vector<double> Grid::compute_vdw_repulsive_grid(Mol2 *Cmol, Parser *Input){
	vector<double> repulsive_potential;
	double repul;
	double r2, r6;
	for (double x=Grid::min_x; x<=Grid::max_x; x+=Input->sampling){
		for (double y=Grid::min_y; y<=Grid::max_y; y+=Input->sampling){
			for (double z=Grid::min_z; z<=Grid::max_z; z+=Input->sampling){
				repul=0.0;
				for (int i=0; i<Cmol->N; i++){
					r2=Grid::dist_squared(Cmol->xyz[i][0], x,Cmol->xyz[i][1], y, Cmol->xyz[i][2], z);
					r6=r2*r2*r2;
					repul+= ((Cmol->epsilons_sqrt[i])*(pow((Cmol->radii[i]), 12)))/(pow((r6+Input->deltaij6),2));
				}
				repulsive_potential.push_back(repul);
			}
		}
	}
	return(repulsive_potential);
}

vector<double> Grid::compute_vdw_repulsive_grid(Mol *Cmol, Parser *Input, vector<vector<double> > xyz){
	vector<double> repulsive_potential;
	double repul;
	double r2, r6;
	for (double x=Grid::min_x; x<=Grid::max_x; x+=Input->sampling){
		for (double y=Grid::min_y; y<=Grid::max_y; y+=Input->sampling){
			for (double z=Grid::min_z; z<=Grid::max_z; z+=Input->sampling){
				repul=0.0;
				for (int i=0; i<Cmol->N; i++){
					r2=Grid::dist_squared(xyz[i][0], x,xyz[i][1], y, xyz[i][2], z);
					r6=r2*r2*r2;
					repul+= ((Cmol->epsilons_sqrt[i])*(pow((Cmol->radii[i]), 12)))/(pow((r6+Input->deltaij6),2));
				}
				repulsive_potential.push_back(repul);
			}
		}
	}
	return(repulsive_potential);
}

vector<double> Grid::compute_vdw_repulsive_grid(Mol2 *Cmol, Parser *Input, vector<vector<double> > xyz){
	vector<double> repulsive_potential;
	double repul;
	double r2, r6;
	for (double x=Grid::min_x; x<=Grid::max_x; x+=Input->sampling){
		for (double y=Grid::min_y; y<=Grid::max_y; y+=Input->sampling){
			for (double z=Grid::min_z; z<=Grid::max_z; z+=Input->sampling){
				repul=0.0;
				for (int i=0; i<Cmol->N; i++){
					r2=Grid::dist_squared(xyz[i][0], x,xyz[i][1], y, xyz[i][2], z);
					r6=r2*r2*r2;
					repul+= ((Cmol->epsilons_sqrt[i])*(pow((Cmol->radii[i]), 12)))/(pow((r6+Input->deltaij6),2));
				}
				repulsive_potential.push_back(repul);
			}
		}
	}
	return(repulsive_potential);
}

