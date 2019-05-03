/*
 * Minimizer2.cpp
 *
 *  Created on: 03/03/2011
 *      Author: Nascimento
 */

#include "Minimizer2.h"

Minimizer2::Minimizer2(Mol2 *_RefMol, Mol2* _CompMol, Parser *_Input, Grid *_Cgrid) {
	Minimizer2::Cmol = _CompMol;
	Minimizer2::RefMol = _RefMol;
	Minimizer2::CompMol = _CompMol;
	Minimizer2::Cgrid = _Cgrid;
	Minimizer2::Input = _Input;
	Minimizer2::new_xyz = Cmol->xyz;
	Minimizer2::pi = 2.0*sqrt(2);
	Minimizer2::pj = 2.0*sqrt(2);
}


double Minimizer2::dist(double x1, double x2, double y1, double y2, double z1, double z2) {
	return (sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

double Minimizer2::dist_squared(double x1, double x2, double y1, double y2, double z1, double z2){
	return ((((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) );
}

vector<double> Minimizer2::compute_com(vector<vector<double> > xyz, Mol2 *Cmol){
	vector<double> com;
	double total_mass =0.00;
	double comx=0.00, comy=0.00, comz=0.00;
	for (int i=0; i<Cmol->N; i++){
		comx += xyz[i][0]*Cmol->masses[i];
		comy += xyz[i][1]*Cmol->masses[i];
		comz += xyz[i][2]*Cmol->masses[i];
		total_mass += Cmol->masses[i];
	}
	com.push_back(comx/total_mass);
	com.push_back(comy/total_mass);
	com.push_back(comz/total_mass);
	return(com);
}

vector<vector<double> > Minimizer2::rototranslate(Mol2 *Cmol, vector<vector<double> >xyz, double alpha, double beta, double gamma, double transx, double transy, double transz){
	vector<vector<double> >nxyz;
	vector<double> temp;
	vector<double> com = compute_com(xyz, Cmol);
	for(int i=0; i < Cmol->N ; i++){
		x=xyz[i][0]-com[0];
		y=xyz[i][1]-com[1];
		z=xyz[i][2]-com[2];
		temp.push_back((((x)*(((cos(alpha*PI/180))*(cos(gamma*PI/180)))-((sin(alpha*PI/180))*(cos(beta*PI/180))*sin(gamma*PI/180)))) + ((y)*(((-cos(alpha*PI/180))*(sin(gamma*PI/180)))-(sin(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180))))+ ((z)*(sin(beta*PI/180)*sin(alpha*PI/180))))+transx+com[0]);
		temp.push_back((((x)*((sin(alpha*PI/180)*cos(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*sin(gamma*PI/180)))) + ((y)*((-sin(alpha*PI/180)*sin(gamma*PI/180))+(cos(alpha*PI/180)*cos(beta*PI/180)*cos(gamma*PI/180)))) + ((z)*(-sin(beta*PI/180)*cos(alpha*PI/180))))+transy+com[1]);
		temp.push_back((((x)*(sin(beta*PI/180)*sin(gamma*PI/180))) + ((y)*sin(beta*PI/180)*cos(gamma*PI/180)) + ((z)*cos(beta*PI/180)))+transz+com[2]);
		nxyz.push_back(temp);
		temp.clear();
	}
	return(nxyz);
}

double Minimizer2::function_gaussian_shape(const gsl_vector *v, void *params){

	double f;
	double a, b, g, transx, transy, transz;
	vector<double> vr;

	a = gsl_vector_get(v,0);
	b = gsl_vector_get(v,1);
	g = gsl_vector_get(v,2);
	transx = gsl_vector_get(v,3);
	transy = gsl_vector_get(v,4);
	transz= gsl_vector_get(v,5);

	new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, a, b, g, transx, transy, transz);

// The Gaussian functions

	Gaussian* Gauss = new Gaussian;
	f = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
	delete Gauss;

	return(-f);
}


double Minimizer2::nlopt_func_gauss_density(const std::vector<double> &x, std::vector<double> &grad, void *data){

	double t2;
	Gaussian* Gauss = new Gaussian;

	new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1], x[2], x[3], x[4], x[5]);
	t1 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);

	if (!grad.empty()){
		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0] + Minimizer2::Input->delta , x[1], x[2], x[3], x[4], x[5]);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[0] = (t2-t1)/Minimizer2::Input->delta;

		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1] + Minimizer2::Input->delta , x[2], x[3], x[4], x[5]);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[1] = (t2-t1)/Minimizer2::Input->delta;

		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1], x[2]+ Minimizer2::Input->delta, x[3], x[4], x[5]);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[2] = (t2-t1)/Minimizer2::Input->delta;

		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1], x[2], x[3]+ Minimizer2::Input->delta, x[4], x[5]);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[3] = (t2-t1)/Minimizer2::Input->delta;

		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1], x[2], x[3], x[4]+ Minimizer2::Input->delta, x[5]);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[4] = (t2-t1)/Minimizer2::Input->delta;

		new_xyz = Minimizer2::rototranslate(Minimizer2::Cmol, Cmol->xyz, x[0], x[1], x[2], x[3], x[4], x[5]+ Minimizer2::Input->delta);
		t2 = Gauss->compute_shape_and_charge_density(Input, RefMol, CompMol, new_xyz);
		grad[5] = (t2-t1)/Minimizer2::Input->delta;
	}
	return (t1);
}

double Minimizer2::minimize_nlopt_ln_auglag(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::LN_AUGLAG,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_mma(){
	srand(time(NULL));
	nlopt::opt *opt = new nlopt::opt(nlopt::LD_MMA,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_isres(){
	srand(time(NULL));

	nlopt::opt* opt = new nlopt::opt(nlopt::GN_ISRES,6);
	const nlopt::opt local_opt(nlopt::LN_NELDERMEAD, 6);
	opt->set_local_optimizer(local_opt);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));
	f_minimum=0.00;

	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_subplex(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::LN_SBPLX,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_simplex(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::LN_NELDERMEAD,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_cobyla(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::LN_COBYLA,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);

	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));
	f_minimum=0.00;

	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_bfgs2(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::LD_LBFGS,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;

	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_direct(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::GN_DIRECT_L,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}

double Minimizer2::minimize_nlopt_stogo(){
	srand(time(NULL));
	nlopt::opt* opt = new nlopt::opt(nlopt::GD_STOGO,6);
	vector<double> lb(6);
	lb[0] = -180.0;
	lb[1] = -90.0;
	lb[2] = -180.0;
	lb[3] = Minimizer2::Cgrid->min_x;
	lb[4] = Minimizer2::Cgrid->min_y;
	lb[5] = Minimizer2::Cgrid->min_z;
	vector<double> ub(6);
	ub[0]= 180.0;
	ub[1] = 90.0;
	ub[2] = 180.0;
	ub[3] = Minimizer2::Cgrid->max_x;
	ub[4] = Minimizer2::Cgrid->max_y;
	ub[5] = Minimizer2::Cgrid->max_z;

	opt->set_lower_bounds(lb);
	opt->set_upper_bounds(ub);
	opt->set_max_objective(Minimizer2::nlopt_func_gauss_density, NULL);
	opt->set_xtol_rel(Minimizer2::Input->tol);
	opt->set_maxtime(Minimizer2::Input->timeout);

	vector<double> x(6);
	double rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[0] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[1] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[2] = (-10.0 + (1.0 *(rnumber*20.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[3] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[4] = (-1.0 + (1.0 *(rnumber*2.0)));
	rnumber = double(rand())/(double(RAND_MAX)+1.0);
	x[5] = (-1.0 + (1.0 *(rnumber*2.0)));

	f_minimum=0.00;
	nlopt::result* result = new nlopt::result(opt->optimize(x,f_minimum));
	delete result;
	delete opt;
	return (f_minimum);
}
