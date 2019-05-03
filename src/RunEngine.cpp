/*
 * RunEngine.cpp
 *
 *  Created on: 19/07/2011
 *      Author: Nascimento
 */

#include "RunEngine.h"

RunEngine::RunEngine(){

}

double RunEngine::distance (double x1, double x2, double y1, double y2, double z1, double z2) {
	return ( sqrt(((x2-x1)*(x2-x1))+((y2-y1)*(y2-y1))+((z2-z1)*(z2-z1))) ); }

#ifdef HAS_GUI

void RunEngine::run(Parser Input, QPlainTextEdit* Ed, QProgressBar* progressbar){

	time0=clock();

	QtWriter* Writer = new QtWriter(&Input, Ed);
	Writer->Editor->updatesEnabled();
	Writer->write_welcome();
	Writer->write_params(&Input);
	Writer->Editor->update();
	CORREL* Correl = new CORREL;

	Mol2* RefMol = new Mol2(&Input, Input.refmol_mol2);
	Grid* RefGrid = new Grid;


	RefGrid->center=RefGrid->compute_com(RefMol->xyz, RefMol->masses);
	RefMol->xyz = Correl->translate(RefMol->xyz, RefMol->N, (-RefGrid->center[0]), (-RefGrid->center[1]), (-RefGrid->center[2]));
	Writer->write_pdb(RefMol, RefMol->xyz, "refmol");
	RefGrid->center=RefGrid->compute_com(RefMol->xyz, RefMol->masses);

	RefGrid->min_x= (RefGrid->center[0]-Input.box_size);
	RefGrid->max_x= (RefGrid->center[0]+Input.box_size);
	RefGrid->min_y= (RefGrid->center[1]-Input.box_size);
	RefGrid->max_y= (RefGrid->center[1]+Input.box_size);
	RefGrid->min_z= (RefGrid->center[2]-Input.box_size);
	RefGrid->max_z= (RefGrid->center[2]+Input.box_size);
	sprintf(info, "Limits of the computation box: %.2f %.2f %.2f %.2f %.2f %.2f",RefGrid->min_x, RefGrid->max_x, RefGrid->min_y, RefGrid->max_y, RefGrid->min_z, RefGrid->max_z);
	Writer->write_to_log(info);

#ifdef DEBUG
	RefGrid->write_box(RefGrid->center, RefGrid->min_x, RefGrid->min_y, RefGrid->min_z, RefGrid->max_x, RefGrid->max_y, RefGrid->max_z);
#endif

/*
 * Ligand Part
 *
 * Ligand part starts here. I want to define the center of mass (COM) for ligand
 * and use it to superpose with the reference molecule by matching both COMs. After that,
 * the ligands are fully searched for a perfect matching by rotation through Euler angles.
*/

	sprintf(info, "%7s %20s %20s %6s %6s %6s %6s", "#", "Mol 1", "Mol 2", "SIElec", "SIVDW", "SI","f_min");
	Writer->write_to_log(info);

	Gaussian* Gauss = new Gaussian;

	FILE *ccfile;
	ccfile = fopen((Input.output+".cc.dat").c_str(), "w");

	RefMol->self_obj_function = Gauss->compute_shape_and_charge_density(&Input, RefMol, RefMol, RefMol->xyz);

	for (signed i=0; i<Input.comparing_molecules.size(); i++){
		indmol = Input.comparing_molecules.at(i).toStdString();
		total_time=0.00;

#ifdef DEBUG
		printf("Opening %s file...\n", indmol.c_str());
#endif


		Mol2 *CompMol = new Mol2(&Input, indmol);
		Grid *CompGrid = new Grid;


		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);
		CompMol->xyz = Correl->translate(CompMol->xyz, CompMol->N, -RefGrid->center[0], -RefGrid->center[1], -RefGrid->center[2]);
		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);

		transx = CompGrid->center[0]-RefGrid->center[0];
		transy = CompGrid->center[1]-RefGrid->center[1];
		transz = CompGrid->center[2]-RefGrid->center[2];

		CompMol->xyz=Correl->translate(CompMol->xyz, CompMol->N, -transx, -transy, -transz);
		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);

/*!
* Minimization starts here by calling a new Minimizer class "Min"
* Through parameter passage (coordinates, charges, N, etc) the class
* is set up and the minimizer actually starts with the method minimize.
* The minimum difference function value reached during the minimization
* is returned to the "f_min" variable and printed below
*/


#ifdef DEBUG
		printf("Initializing the minimizer...\n");
#endif

		Minimizer2 *Min = new Minimizer2(RefMol, CompMol, &Input, RefGrid);

#ifdef DEBUG
		printf("Minimizing...\n");
#endif

		if (Input.align == true ){

			if (Input.minimizer == "nlopt_ln_auglag") {
				f_min = Min->minimize_nlopt_ln_auglag();
			}
			else if(Input.minimizer == "nlopt_mma"){
				f_min = Min->minimize_nlopt_mma();
			}

			else if(Input.minimizer == "nlopt_isres"){
				f_min = Min->minimize_nlopt_isres();
			}
			else if (Input.minimizer == "nlopt_subplex"){
				f_min = Min->minimize_nlopt_subplex();
			}
			else if (Input.minimizer == "nlopt_simplex"){
				f_min = Min->minimize_nlopt_simplex();
			}
			else if (Input.minimizer == "nlopt_lbfgs2"){
				f_min = Min->minimize_nlopt_bfgs2();
			}
			else if (Input.minimizer == "nlopt_cobyla"){
				f_min = Min->minimize_nlopt_cobyla();
			}
			else {
				printf("Minimization method %s not implemented!\n", Input.minimizer.c_str());
				printf("Available minimizers:\n");
				printf("\t \"nlopt_simplex\"     ->  Nelder-Mead Simplex method as implemented in NLOPT\n");
				printf("\t \"nlopt_subplex\"     ->  Subplex - a variant of Nelder-Mead simplex that uses Nelder-Mead on a sequence of subspaces as implemented in NLOPT\n\n");
				printf("\t \"nlopt_mma\"         ->  Method of Moving Asymptotes as implemented in NLOPT\n");
				printf("\t \"nlopt_lbfgs2\"      ->  Broyden-Fletcher-Goldfarb-Shanno method as implemented in NLOPT\n");
				printf("\t \"nlopt_cobyla\"  ->  NLOPT's COBYLA minimization of gaussian sum of atomic volumes and charges\n");
				printf("\t \"nlopt_ln_auglag\"      ->  Augmented Lagrangian method as implemented in NLOPT\n");
				printf("\t \"nlopt_isres\"       ->  Improved Stochastic Ranking Evolution Strategy as implemented in NLOPT\n");
				printf("Please check!\n");
				printf("Exiting...\n");
				exit(1);
			}
		}

#ifdef DEBUG
		printf("Computing Correlation coefficient...\n");
#endif

		Vab = Gauss->compute_shape_and_charge_density(&Input, RefMol, CompMol, Min->new_xyz);
		CompMol->self_obj_function = Gauss->compute_shape_and_charge_density(&Input, CompMol, CompMol, Min->new_xyz, Min->new_xyz);

		cc = 2*Vab/(RefMol->self_obj_function+CompMol->self_obj_function);

/*
* End of the minimization procedure.
*
* Some output is written below.
*/

		sprintf(info, "%7d %63.63s %6.3f %6.3e", int(i+1), indmol.c_str(), cc, f_min);
		Writer->write_to_log(info);

		fprintf(ccfile, "%d ; %6.3e ; %.5f\n", int(i+1), f_min, cc);

		sprintf(outname, "%s.%d", Input.output.c_str(), i+1);

#ifdef DEBUG
		printf("Writing coordinates...\n");
#endif

		if (Input.write_pdb ==true){
			if (Input.use_write_coord_threshold){
				if ( cc >= Input.write_coord_threshold){
					Writer->writeMol2(CompMol, Min->new_xyz, cc);
				}
			}
			else {
				Writer->writeMol2(CompMol, Min->new_xyz, cc);
			}
		}

		delete CompMol;
		delete CompGrid;
		delete Min;

	progressbar->setValue(int((i+1)*100.0/Input.comparing_molecules.size()));

	}

	fclose(ccfile);

	Writer->write_to_log();
	Writer->write_to_log();

	time1 = clock();

	total_time = (double(time1)-double(time0))/CLOCKS_PER_SEC;

	sprintf(info, "Time elapsed : %10.3f seconds (%10.2f minutes)", total_time, (total_time/60.0));
	Writer->write_to_log(info);
	Writer->write_to_log();

	sprintf(info, "%s", "Calculation Finished!");
	Writer->write_to_log(info);
	Writer->write_to_log();

	sprintf(info, "%s", "Quitting gracefully...");
	Writer->write_to_log(info);

	sprintf(info, "%s", "*********************************************************************************");
	Writer->write_to_log(info);
	Writer->write_to_log();

	delete Writer;
	delete Correl;
	delete RefMol;
	delete RefGrid;
//	Writer.~QtWriter();

}

#endif

void RunEngine::run(Parser Input){

	time0=clock();

	Grid* RefGrid = new Grid;;
	CORREL* Correl = new CORREL;

	Writer* BWriter = new Writer(&Input);
	BWriter->write_welcome();
	BWriter->write_params(&Input);

	Mol2* RefMol = new Mol2(&Input, Input.refmol_mol2);

	RefGrid->center=RefGrid->compute_com(RefMol->xyz, RefMol->masses);
	RefMol->xyz = Correl->translate(RefMol->xyz, RefMol->N, (-RefGrid->center[0]), (-RefGrid->center[1]), (-RefGrid->center[2]));
	BWriter->writeMol2(RefMol, RefMol->xyz, "refmol");
	RefGrid->center=RefGrid->compute_com(RefMol->xyz, RefMol->masses);

	RefGrid->min_x= (RefGrid->center[0]-Input.box_size);
	RefGrid->max_x= (RefGrid->center[0]+Input.box_size);
	RefGrid->min_y= (RefGrid->center[1]-Input.box_size);
	RefGrid->max_y= (RefGrid->center[1]+Input.box_size);
	RefGrid->min_z= (RefGrid->center[2]-Input.box_size);
	RefGrid->max_z= (RefGrid->center[2]+Input.box_size);
	sprintf(info, "Limits of the computation box: %.2f %.2f %.2f %.2f %.2f %.2f",RefGrid->min_x, RefGrid->max_x, RefGrid->min_y, RefGrid->max_y, RefGrid->min_z, RefGrid->max_z);
	BWriter->write_to_log(info);

#ifdef DEBUG
	RefGrid->write_box(RefGrid.center, RefGrid.min_x, RefGrid.min_y, RefGrid.min_z, RefGrid.max_x, RefGrid.max_y, RefGrid.max_z);
#endif


	/*
	 * Ligand Part
	 *
	 * Ligand part starts here. I want to define the center of mass (COM) for ligand
	 * and use it to superpose with the reference molecule by matching both COMs. After that,
	 * the ligands are fully searched for a perfect matching by rotation through Euler angles.
	 */

	ifstream multifile(Input.molsfile.c_str());
	sprintf(info, "%-7s %-56s %-6s %-6s", "#", "Mol", "SI","f_min");
	BWriter->write_to_log(info);
	multifile >> indmol;
	int count=0;
	total_time=0.00;

	Gaussian* Gauss = new Gaussian;

	FILE *ccfile;
	ccfile = fopen((Input.output+".cc.dat").c_str(), "w");

	BWriter->write_is_running();

	RefMol->self_obj_function = Gauss->compute_shape_and_charge_density(&Input, RefMol, RefMol, RefMol->xyz);

	while (indmol != "EOF" && !multifile.eof()){
		count++;

		Mol2 *CompMol = new Mol2(&Input, indmol);
		Grid *CompGrid = new Grid;

		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);
		CompMol->xyz = Correl->translate(CompMol->xyz, CompMol->N, -RefGrid->center[0], -RefGrid->center[1], -RefGrid->center[2]);
		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);

		transx = CompGrid->center[0]-RefGrid->center[0];
		transy = CompGrid->center[1]-RefGrid->center[1];
		transz = CompGrid->center[2]-RefGrid->center[2];

		CompMol->xyz=Correl->translate(CompMol->xyz, CompMol->N, -transx, -transy, -transz);
		CompGrid->center=CompGrid->compute_com(CompMol->xyz, CompMol->masses);


		/*!
		 * Minimization starts here by calling a new Minimizer class "Min"
		 * Through parameter passage (coordinates, charges, N, etc) the class
		 * is set up and the minimizer actually starts with the method minimize.
		 * The minimum difference function value reached during the minimization
		 * is returned to the "f_min" variable and printed below
		 */


		Minimizer2 *Min = new Minimizer2(RefMol, CompMol, &Input, RefGrid);

		if (Input.align == true ){

			if (Input.minimizer == "nlopt_ln_auglag") {
				f_min = Min->minimize_nlopt_ln_auglag();
			}
			else if(Input.minimizer == "nlopt_mma"){
				f_min = Min->minimize_nlopt_mma();
			}
			else if(Input.minimizer == "nlopt_isres"){
				f_min = Min->minimize_nlopt_isres();
			}
			else if (Input.minimizer == "nlopt_subplex"){
				f_min = Min->minimize_nlopt_subplex();
			}
			else if (Input.minimizer == "nlopt_simplex"){
				f_min = Min->minimize_nlopt_simplex();
			}
			else if (Input.minimizer == "nlopt_lbfgs2"){
				f_min = Min->minimize_nlopt_bfgs2();
			}
			else if(Input.minimizer == "nlopt_stogo"){
				f_min = Min->minimize_nlopt_stogo();
			}
			else if(Input.minimizer == "nlopt_direct-l"){
				f_min = Min->minimize_nlopt_direct();
			}
			else if(Input.minimizer == "nlopt_cobyla"){
				f_min = Min->minimize_nlopt_cobyla();
			}
			else {
				printf("Minimization method %s not implemented!\n", Input.minimizer.c_str());
				printf("Available minimizers:\n");
				printf("\t \"nlopt_simplex\"     ->  Nelder-Mead Simplex method as implemented in NLOPT\n");
				printf("\t \"nlopt_subplex\"     ->  Subplex - a variant of Nelder-Mead simplex that uses Nelder-Mead on a sequence of subspaces as implemented in NLOPT\n\n");
				printf("\t \"nlopt_mma\"         ->  Method of Moving Asymptotes as implemented in NLOPT\n");
				printf("\t \"nlopt_lbfgs2\"      ->  Broyden-Fletcher-Goldfarb-Shanno method as implemented in NLOPT\n");
				printf("\t \"nlopt_cobyla\"      ->  NLOPT's COBYLA minimization\n");
				printf("\t \"nlopt_ln_auglag\"   ->  Augmented Lagrangian method as implemented in NLOPT\n");
				printf("\t \"nlopt_isres\"       ->  Improved Stochastic Ranking Evolution Strategy as implemented in NLOPT\n");
				printf("\t \"nlopt_stogo\"       ->  Global optimization algorithm that systematically divides the search space, as implemented in NLOPT\n");
				printf("\t \"nlopt_direct-l\"    ->  DIviding RECTangles algorithm for global optimization, as implemented in NLOPT\n");
				printf("Please check!\n");
				printf("Exiting...\n");
				exit(1);
			}
		}


		Vab = Gauss->compute_shape_and_charge_density(&Input, RefMol, CompMol, Min->new_xyz);

		CompMol->self_obj_function = Gauss->compute_shape_and_charge_density(&Input, CompMol, CompMol, Min->new_xyz, Min->new_xyz);

		cc = 2*Vab/(RefMol->self_obj_function+CompMol->self_obj_function);

		sprintf(info, "%7d %-56.56s %6.3f %6.3e", count, indmol.c_str(), cc, f_min);
		BWriter->write_to_log(info);

		// Writing CC file...
		fprintf(ccfile, "%d ; %6.3e ; %.5f\n", count, f_min, cc);

		//written
		sprintf(outname, "%s.%d", Input.output.c_str(), count);

		if (Input.write_pdb == true){
			if (Input.use_write_coord_threshold){
				if ( cc >= Input.write_coord_threshold){
					BWriter->writeMol2(CompMol, Min->new_xyz, cc);
				}
			}
			else {
				BWriter->writeMol2(CompMol, Min->new_xyz, cc);
			}
		}

		delete CompMol;
		delete CompGrid;
		delete Min;

		multifile >> indmol;
	}

	multifile.close();

	BWriter->write_to_log();
	BWriter->write_to_log();

	time1=clock();
	total_time=(double(time1)-double(time0))/CLOCKS_PER_SEC;

	sprintf(info, "Time elapsed for minimization: %10.3f seconds (%10.2f minutes)", total_time, (total_time/60.0));
	BWriter->write_to_log(info);
	BWriter->write_to_log();

	fclose(ccfile);


	BWriter->write_to_log();

	sprintf(info, "%s", "Calculation Finished!");
	BWriter->write_to_log(info);
	BWriter->write_to_log();

	sprintf(info, "%s", "Quitting gracefully...");
	BWriter->write_to_log(info);

	BWriter->close();
	delete BWriter;
	delete RefMol;
	delete RefGrid;
	delete Correl;
}

vector<int> RunEngine::ranker(vector<float> CCs){
	vector<int> sorted_indexes;
	int tmp;
	float dtmp;

// Filling the int vector with its own indexes
	for (unsigned i=0; i<CCs.size(); i++){
		sorted_indexes.push_back(i);
	}

//Ranking the int vector

	for (unsigned i=0; i<CCs.size()-1; i++){
		for (unsigned j=i+1; j<CCs.size(); j++){
			if (CCs[i] < CCs[j]){
				tmp = sorted_indexes[i];
				sorted_indexes[i] = sorted_indexes[j];
				sorted_indexes[j] = tmp;

				dtmp = CCs[i];
				CCs[i] = CCs[j];
				CCs[j] = dtmp;
			}
		}
	}

//End of ranking

	return (sorted_indexes);
}
