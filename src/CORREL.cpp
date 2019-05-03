/*
 * CORREL.cpp
 *
 *  Created on: 20/08/2010
 *      Author: Nascimento
 */

#include "CORREL.h"

double dist(double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) ); }

CORREL::CORREL(){

}

vector<vector<double> >CORREL::rototranslate(vector<vector<double> > xyz, int N, double alpha, double beta, double gamma, double transx, double transy, double transz){
	new_xyz.clear();
	for(int i=0; i < N ; i++){
		x=xyz[i][0];
		y=xyz[i][1];
		z=xyz[i][2];
		temp.push_back((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx);
		temp.push_back((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy);
		temp.push_back((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz);
		new_xyz.push_back(temp);
		temp.clear();
	}
	return(new_xyz);
}

double CORREL::compute_correl(vector<vector<double> > gridcoords, vector<double> elec_potential, int N, vector<vector<double> >xyz, vector<double> charges){
	correl = 0.00;
	for (unsigned i=0; i<gridcoords.size(); i++){
		elec=0.00;
		for (int j=0; j < N; j++){
			r = dist(gridcoords[i][0], xyz[j][0], gridcoords[i][1], xyz[j][1], gridcoords[i][2], xyz[j][2]);
			elec+=((charges[j]/18.2223)/r);
		}
		correl+=(elec*elec_potential[i]);
	}
	return(correl);
}

double CORREL::compute_correl_wade(vector<double> p1, vector<double> p2){
	sum_product=0.00;
	if (p1.size() != p2.size()){
		printf("Vector P1 size: %10d    Vector P2 size: %10d\n", p1.size(), p2.size());
		exit(1);
	}
	else {
		for (unsigned i=0; i<p1.size(); i++){
			sum_product+=p1[i]*p2[i];
		}
		return(sum_product);
	}
}

double CORREL::compute_correl_wade(vector<double> elec_potential, vector<double> p2, double sum_of_elec_potential){
	double t1 = CORREL::compute_correl_wade(elec_potential, p2);
	double t2 = CORREL::compute_correl_wade(p2, p2);
	return ((2*t1)/(t2+sum_of_elec_potential));
}


vector<vector<double> > CORREL::translate(vector<vector<double> > xyz, int N, double tx, double ty, double tz){
	for (int i=0; i<N; i++){
		xyz[i][0]=xyz[i][0]+tx;
		xyz[i][1]=xyz[i][1]+ty;
		xyz[i][2]=xyz[i][2]+tz;
	}
	return(xyz);
}

void CORREL::write_pdb(vector<vector<double> >xyz, int N, vector<string> atomnames, string name){
	gzFile outpdb;
	outpdb = gzopen(name.c_str(), "w");
	gzprintf(outpdb, "MODEL\n");
	for (int i=0; i<N; i++){
		gzprintf(outpdb, "ATOM    %3d %4s LIG     1    % 8.3f% 7.3f% 7.3f\n", i, atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2]);
	}
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
}


