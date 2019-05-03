/*
 * Gaussian.cpp
 *
 *  Created on: 18/11/2011
 *      Author: Nascimento
 */

#include "Gaussian.h"

Gaussian::Gaussian(Mol2* RefMol, Mol2* CompMol, double* Vab) {
	time0=clock();
	pi=2.0*sqrt(2);
	pj=2.0*sqrt(2);
	printf("Refmol N: %d    CompMol N: %d\n", RefMol->N, CompMol->N);
	*Vab = this->compute_si(RefMol, CompMol);
	time1=clock();

#ifdef DEBUG
	printf("Gaussian computation took %7.3f seconds.\n", ((double(time1)-double(time0))/CLOCKS_PER_SEC));
#endif

}

Gaussian::Gaussian(void){
	pi=2.0*sqrt(2);
	pj=2.0*sqrt(2);
}

Gaussian::Gaussian(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, double* Vab) {
	time0=clock();
	pi=2.0*sqrt(2);
	pj=2.0*sqrt(2);
	printf("Refmol N: %d    CompMol N: %d\n", RefMol->N, CompMol->N);
	*Vab = this->compute_si(RefMol, CompMol, xyz);
	time1=clock();

#ifdef DEBUG
	printf("Gaussian computation took %7.3f seconds.\n", ((double(time1)-double(time0))/CLOCKS_PER_SEC));
#endif

}

double Gaussian::compute_shape_density(Mol2* RefMol, Mol2* CompMol){
	double V=0.0;
	for (int i=0; i<RefMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow( ((3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(RefMol->xyz[i][0], CompMol->xyz[j][0],RefMol->xyz[i][1], CompMol->xyz[j][1],RefMol->xyz[i][2], CompMol->xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_shape_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<RefMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_shape_density(Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<CompMol->N; i++){
		alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(CompMol->radii[i], 3))), (2.0/3.0)));
		for (int j=0; j<CompMol->N; j++){
			alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			dij2 = dist_squared(xyz[i][0], xyz[j][0],xyz[i][1], xyz[j][1],xyz[i][2], xyz[j][2]);
			V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
		}
	}
	return(V);
}

double Gaussian::compute_si(Mol2* RefMol, Mol2* CompMol){
	double si = ((2*this->compute_shape_density(RefMol, CompMol))/(this->compute_shape_density(RefMol, RefMol) + this->compute_shape_density(CompMol, CompMol)));
	printf("Shape density correlation: %.3f\n", si);
	return (si);
}

double Gaussian::compute_si(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double t1, t2, t3, si;
	t1= this->compute_shape_density(RefMol, CompMol, xyz);
	t2 = this->compute_shape_density(RefMol, RefMol);
	t3 = this->compute_shape_density(CompMol, xyz);
	si = 2*t1/(t2+t3);

#ifdef DEBUG
	printf("Shape density correlation: %.3f\n", si);
#endif

	return (si);
}

double Gaussian::compute_si_pos_charges(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double t1, t2, t3, si;
	t1= this->compute_positive_charge_density(RefMol, CompMol, xyz);
	t2 = this->compute_positive_charge_density(RefMol, RefMol, RefMol->xyz);
	t3 = this->compute_positive_charge_density(CompMol, xyz);
	si = 2*t1/(t2+t3);
	return (si);
}

double Gaussian::compute_si_neg_charges(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double t1, t2, t3, si;
	t1 = this->compute_negative_charge_density(RefMol, CompMol, xyz);
	t2 = this->compute_negative_charge_density(RefMol, RefMol, RefMol->xyz);
	t3 = this->compute_negative_charge_density(CompMol, xyz);
	si = 2*t1/(t2+t3);
	return (si);
}


double Gaussian::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2){
	return ((((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

double Gaussian::compute_positive_charge_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<RefMol->N; i++){
		if (RefMol->charges[i]> 0.0){
			alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
		}
		else {
			alpha_i = 0.0;
		}
		for (int j=0; j<CompMol->N; j++){
			if (CompMol->charges[j] > 0.00){
				alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
			}
			else {
				alpha_j = 0.0;
			}

			dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);

			if (alpha_i!= 0.0 and alpha_j != 0.0){
				V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
			}
		}
	}
	return(V);
}

double Gaussian::compute_positive_charge_density(Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<CompMol->N; i++){
		if (CompMol->charges[i]> 0.0){
			alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[i]), 3))), (2.0/3.0)));
		}
		else {
			alpha_i = 0.0;
		}
		for (int j=0; j<CompMol->N; j++){
			if (CompMol->charges[j] > 0.00){
				alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
			}
			else {
				alpha_j = 0.0;
			}

			dij2 = dist_squared(xyz[i][0], xyz[j][0],xyz[i][1], xyz[j][1],xyz[i][2], xyz[j][2]);

			if (alpha_i!= 0.0 and alpha_j != 0.0){
				V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
			}
		}
	}
	return(V);
}

double Gaussian::compute_negative_charge_density(Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<RefMol->N; i++){
		if (RefMol->charges[i] < 0.0){
			alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(abs(1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
		}
		else {
			alpha_i = 0.0;
		}
		for (int j=0; j<CompMol->N; j++){
			if (CompMol->charges[j] < 0.0){
				alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(abs(1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
			}
			else {
				alpha_j= 0.0;
			}

			dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);

			if (alpha_i !=0.0 and alpha_j != 0.0){
				V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
			}
		}
	}
	return(V);
}

double Gaussian::compute_negative_charge_density(Mol2* CompMol, vector<vector<double> > xyz){
	double V=0.0;
	for (int i=0; i<CompMol->N; i++){
		if (CompMol->charges[i] < 0.0){
			alpha_i = PI*( pow( ( (3*pi)/(4*PI*pow(abs(1.00 + CompMol->charges[i]), 3))), (2.0/3.0)));
		}
		else {
			alpha_i = 0.0;
		}
		for (int j=0; j<CompMol->N; j++){
			if (CompMol->charges[j] < 0.0){
				alpha_j = PI*pow(  ( (3*pi)/(4*PI*pow(abs(1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
			}
			else {
				alpha_j= 0.0;
			}

			dij2 = dist_squared(xyz[i][0], xyz[j][0],xyz[i][1], xyz[j][1],xyz[i][2], xyz[j][2]);

			if (alpha_i!= 0.00 and alpha_j != 0.0){
				V+= pi*pj*(exp((-alpha_i*alpha_j*dij2)/(alpha_i+alpha_j)))*(pow((PI/(alpha_i+alpha_j)), (3.0/2.0)));
			}
		}
	}
	return(V);
}

double Gaussian::compute_shape_and_charge_density(Parser *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz){

	double Vshape=0.0, Vpos=0.0, Vneg=0.0;

	double alphai_shape=0.0, alphai_pos=0.0, alphai_neg=0.0;

	double alphaj_shape=0.0, alphaj_pos=0.0, alphaj_neg=0.0;

	for (int i=0; i<RefMol->N; i++){

		if (RefMol->radii[i] != 0.0){
			alphai_shape = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		}
		else {
			alphai_shape=0.0;
		}


		if (RefMol->charges[i] > 0.0){
			alphai_pos = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
			alphai_neg=0.0;
		}

		else if (RefMol->charges[i] < 0.0) {
			alphai_neg = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + abs(RefMol->charges[i])), 3))), (2.0/3.0)));
			alphai_pos=0.0;
		}

		for (int j=0; j<CompMol->N; j++){
			if (CompMol->radii[j] != 0.0000){
				alphaj_shape = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			}
			else {
				alphaj_shape=0.0000;
			}


			if (CompMol->charges[j] > 0.00){
				alphaj_pos = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
				alphaj_neg=0.0;
			}

			else if (CompMol->charges[j] < 0.00){
				alphaj_neg = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + abs(CompMol->charges[j])), 3))), (2.0/3.0));
				alphaj_pos = 0.0;
			}

			dij2 = dist_squared(RefMol->xyz[i][0], xyz[j][0],RefMol->xyz[i][1], xyz[j][1],RefMol->xyz[i][2], xyz[j][2]);

			if (alphai_shape != 0.0 and alphaj_shape != 0.0 ){
				Vshape += pi*pj*(exp((-alphai_shape*alphaj_shape*dij2)/(alphai_shape+alphaj_shape)))*(pow((PI/(alphai_shape+alphaj_shape)), (3.0/2.0)));
			}

			if (alphai_pos != 0.0 and alphaj_pos != 0.0 ){
				Vpos += pi*pj*(exp((-alphai_pos*alphaj_pos*dij2)/(alphai_pos+alphaj_pos)))*(pow((PI/(alphai_pos+alphaj_pos)), (3.0/2.0)));
			}

			if (alphai_neg != 0.0 and alphaj_neg != 0.0 ){
				Vneg += pi*pj*(exp((-alphai_neg*alphaj_neg*dij2)/(alphai_neg+alphaj_neg)))*(pow((PI/(alphai_neg+alphaj_neg)), (3.0/2.0)));
			}
		}
	}

#ifdef DEBUG
	printf("Vshape: %.3f    Vpos: %.3f    Vneg: %.3f    DIFF:%.3f\n", Input->vdw_scale*Vshape, Input->elec_scale*Vpos, Input->elec_scale*Vneg, (Input->vdw_scale*Vshape)-(Input->elec_scale*Vpos)-(Input->elec_scale*Vneg));
#endif

	return((Input->vdw_scale*Vshape)+(Input->elec_scale*Vpos)+(Input->elec_scale*Vneg));
}

double Gaussian::compute_shape_and_charge_density(Parser *Input, Mol2* RefMol, Mol2* CompMol, vector<vector<double> > xyz, vector<vector<double> > cxyz){

	double Vshape=0.0, Vpos=0.0, Vneg=0.0;

	double alphai_shape=0.0, alphai_pos=0.0, alphai_neg=0.0;

	double alphaj_shape=0.0, alphaj_pos=0.0, alphaj_neg=0.0;

	for (int i=0; i<RefMol->N; i++){

		if (RefMol->radii[i] != 0.0){
			alphai_shape = PI*( pow( ( (3*pi)/(4*PI*pow(RefMol->radii[i], 3))), (2.0/3.0)));
		}
		else {
			alphai_shape=0.0;
		}


		if (RefMol->charges[i] > 0.0){
			alphai_pos = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + RefMol->charges[i]), 3))), (2.0/3.0)));
			alphai_neg=0.0;
		}

		else if (RefMol->charges[i] < 0.0) {
			alphai_neg = PI*( pow( ( (3*pi)/(4*PI*pow((1.00 + abs(RefMol->charges[i])), 3))), (2.0/3.0)));
			alphai_pos=0.0;
		}

		for (int j=0; j<CompMol->N; j++){
			if (CompMol->radii[j] != 0.0000){
				alphaj_shape = PI*pow(  ( (3*pi)/(4*PI*pow(CompMol->radii[j], 3))), (2.0/3.0));
			}
			else {
				alphaj_shape=0.0000;
			}


			if (CompMol->charges[j] > 0.00){
				alphaj_pos = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + CompMol->charges[j]), 3))), (2.0/3.0));
				alphaj_neg=0.0;
			}

			else if (CompMol->charges[j] < 0.00){
				alphaj_neg = PI*pow(  ( (3*pi)/(4*PI*pow((1.00 + abs(CompMol->charges[j])), 3))), (2.0/3.0));
				alphaj_pos = 0.0;
			}

			dij2 = dist_squared(xyz[i][0], cxyz[j][0],xyz[i][1], cxyz[j][1],xyz[i][2], cxyz[j][2]);

			if (alphai_shape != 0.0 and alphaj_shape != 0.0 ){
				Vshape += pi*pj*(exp((-alphai_shape*alphaj_shape*dij2)/(alphai_shape+alphaj_shape)))*(pow((PI/(alphai_shape+alphaj_shape)), (3.0/2.0)));
			}

			if (alphai_pos != 0.0 and alphaj_pos != 0.0 ){
				Vpos += pi*pj*(exp((-alphai_pos*alphaj_pos*dij2)/(alphai_pos+alphaj_pos)))*(pow((PI/(alphai_pos+alphaj_pos)), (3.0/2.0)));
			}

			if (alphai_neg != 0.0 and alphaj_neg != 0.0 ){
				Vneg += pi*pj*(exp((-alphai_neg*alphaj_neg*dij2)/(alphai_neg+alphaj_neg)))*(pow((PI/(alphai_neg+alphaj_neg)), (3.0/2.0)));
			}
		}
	}

#ifdef DEBUG
	printf("Vshape: %.3f    Vpos: %.3f    Vneg: %.3f    DIFF:%.3f\n", Input->vdw_scale*Vshape, Input->elec_scale*Vpos, Input->elec_scale*Vneg, (Input->vdw_scale*Vshape)-(Input->elec_scale*Vpos)-(Input->elec_scale*Vneg));
#endif

	return((Input->vdw_scale*Vshape)+(Input->elec_scale*Vpos)+(Input->elec_scale*Vneg));
}

Gaussian::~Gaussian() {
}
