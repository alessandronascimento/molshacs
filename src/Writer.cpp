/*
 * Writer.cpp
 *
 *  Created on: 03/03/2011
 *      Author: Nascimento
 */

#include "Writer.h"

Writer::Writer(Parser *Input) {
	logfile= fopen((Input->output+".log").c_str(), "w");
	outmol2 = gzopen((Input->output+".mol2.gz").c_str(), "w");
}

void Writer::write_params(Parser *Input){
	printf("* %-82s *\n", "Parsed parameters:");
	fprintf(logfile, "* %-82s *\n", "Parsed parameters:");

	printf("* %82s *\n", "");
	fprintf(logfile,"* %82s *\n", "");

	printf("*%-20s %-63s*\n", "refmol_mol2", Input->refmol_mol2.c_str());
	fprintf(logfile,"*%-20s %-63s*\n", "refmol_mol2", Input->refmol_mol2.c_str());

	printf("*%-20s %-20d %42s*\n", "mol2_aa", Input->mol2_aa, "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "mol2_aa", Input->mol2_aa, "");

	printf("*%-20s %-20d %42s*\n", "mol2 amber atoms?", Input->mol2_aa, "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "mol2 amber atoms", Input->mol2_aa, "");

	printf("*%-20s %-20d %42s*\n", "align_molecules", Input->align, "");
	fprintf(logfile,"*%-20s %-20d %42s*\n", "align_molecules", Input->align, "");

	printf("*%-20s %-20.3e %42s*\n", "search step", Input->step , "");
	fprintf(logfile,"*%-20s %-20.ef %42s*\n", "search step", Input->step , "");

	printf("*%-20s %-20.2e %42s*\n", "minimization tol", Input->tol , "");
	fprintf(logfile, "*%-20s %-20.2e %42s*\n", "minimization tol", Input->tol , "");

	printf("*%-20s %-20d %42s*\n", "minimization steps", Input->nsteps , "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "minimization steps", Input->nsteps , "");

	printf("*%-20s %-20.4e %42s*\n", "minimization delta", Input->delta, "");
	fprintf(logfile, "*%-20s %-20.4e %42s*\n", "minimization delta", Input->delta, "");

	printf("*%-20s %-63s*\n", "minimizer", Input->minimizer.c_str());
	fprintf(logfile,"*%-20s %-63s*\n", "minimizer", Input->minimizer.c_str());

	printf("*%-20s %-20d %42s*\n", "multi", Input->multi, "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "multi", Input->multi, "");

	printf("*%-20s %-63s*\n", "molsfile", Input->molsfile.c_str());
	fprintf(logfile, "*%-20s %-63s*\n", "molsfile", Input->molsfile.c_str());

	printf("*%-20s %-63s*\n", "output", Input->output.c_str());
	fprintf(logfile, "*%-20s %-63s*\n", "output", Input->output.c_str());

	printf("*%-20s %-20d %42s*\n", "write coordinates", Input->write_pdb, "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "write coordinates", Input->write_pdb, "");

	printf("*%-20s %-20.2f %42s*\n", "electrostatic scale", Input->elec_scale, "");
	fprintf(logfile, "*%-20s %-20.2f %42s*\n", "electrostatic scale", Input->elec_scale, "");

	printf("*%-20s %-20.2f %42s*\n", "VDW scale", Input->vdw_scale, "");
	fprintf(logfile, "*%-20s %-20.2f %42s*\n", "VDW scale", Input->vdw_scale, "");

	printf("*%-20s %-20d %42s*\n", "timeout", Input->timeout, "");
	fprintf(logfile, "*%-20s %-20d %42s*\n", "timeout", Input->timeout, "");

	if (Input->use_write_coord_threshold){
		printf("*%-20s %-20f %42s*\n", "write_threshold", Input->write_coord_threshold, "");
		fprintf(logfile, "*%-20s %-20f %42s*\n", "write_threshold", Input->write_coord_threshold, "");
	}

	printf("* %82s *\n", "");
	fprintf(logfile, "* %82s *\n", "");
	printf("* %82s *\n", "");
	printf("**************************************************************************************\n");
	fprintf(logfile,"**************************************************************************************\n");
	printf("* %82s *\n", "");
}

void Writer::write_welcome(void){
	printf("**************************************************************************************\n");
	printf("*                                                                                    *\n");
	printf("*                                   MolShaCS v. 1.0                                  *\n");
	printf("*                        Molecular Shape and Charge Similarity                       *\n");
	printf("*                                                                                    *\n");
	printf("* Written by Alessandro S. Nascimento-Jun/2010 -  alessandro.nascimento@ufabc.edu.br *\n");
	printf("*                                                                                    *\n");
	printf("*                   Universidade Federal do ABC - Brazil                             *\n");
	printf("*                                                                                    *\n");
	printf("*               http://nascimento.ufabc.edu.br/wordpress/?page_id=14                 *\n");
	printf("*                                                                                    *\n");
	printf("**************************************************************************************\n");
	printf("* %82s *\n", ""); 						// I am using a 85 character long field to print all messages!

	fprintf(logfile, "**************************************************************************************\n");
	fprintf(logfile, "*                                 MolShaCS v. 1.0                                    *\n");
	fprintf(logfile, "*                        Molecular Shape and CHarge Similarity                       *\n");
	fprintf(logfile, "*                                                                                    *\n");
	fprintf(logfile, "* Written by Alessandro S. Nascimento-Jun/2010 -  alessandro.nascimento@ufabc.edu.br *\n");
	fprintf(logfile, "*                                                                                    *\n");
	fprintf(logfile, "*                   Universidade Federal do ABC - Brazil                             *\n");
	fprintf(logfile, "*                                                                                    *\n");
	fprintf(logfile, "*               http://nascimento.ufabc.edu.br/wordpress/?page_id=14                 *\n");
	fprintf(logfile, "*                                                                                    *\n");
	fprintf(logfile, "**************************************************************************************\n");
	fprintf(logfile, "* %82s *\n", ""); 						// I am using a 85 character long field to print all messages!

}

void Writer::write_to_log(char info[82]){
	fprintf(logfile, "* %-82s *\n", info);
}

void Writer::write_is_running(void){
	printf("* %-82s *\n", "Running.... please wait");
}

void Writer::write_to_log(void){
	fprintf(logfile, "* %-82s *\n", "");
}

void Writer::write_pdb(Mol *Cmol, vector<vector<double> >xyz, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
	gzprintf(outpdb, "MDL\n");
	gzprintf(outpdb, "REMARK\n");
	int i=0;
	int resn=0;

	while (resn < int(Cmol->residue_pointer.size()-1)){
		while(i < Cmol->residue_pointer[resn+1]-1){
			gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000    0.00  1\n", i+1, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
			i++;
		}
		resn++;
	}

	while (i < Cmol->N){
		gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000    0.00  1\n", i+1, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2]);
		i++;
	}
	gzprintf(outpdb, "TER\n");
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
}

void Writer::write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
	gzprintf(outpdb, "MDL\n");
	gzprintf(outpdb, "REMARK %s\n", Cmol->molname.c_str());
	int i=0;
	int resn=0;

	while (resn < int(Cmol->residue_pointer.size()-1)){
		while(i < Cmol->residue_pointer[resn+1]-1){
			gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000   % 4.2f  1\n", i+1, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2], Cmol->charges[i]);
			i++;
		}
		resn++;
	}

	while (i < Cmol->N){
		gzprintf(outpdb, "ATOM   %4d%4s  %3.3s  %4d    % 8.3f % 7.3f % 7.3f   0.000   % 4.2f  1\n", i+1, Cmol->atomnames[i].c_str(), Cmol->resnames[resn].c_str(), resn+1, xyz[i][0], xyz[i][1], xyz[i][2], Cmol->charges[i]);
		i++;
	}
	gzprintf(outpdb, "TER\n");
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
}

void Writer::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, string outname){
	gzFile outmol2;
	outmol2 = gzopen((outname+".mol2.gz").c_str(), "a");
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
				gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
				gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
			gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
			gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}
	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
	gzclose(outmol2);
}

void Writer::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double cc){
	gzprintf(outmol2, "@<TRIPOS>MOLECULE\n");
	gzprintf(outmol2, "%s\n", Cmol->molname.c_str());
	gzprintf(outmol2, "%d %d %d\n", Cmol->N, Cmol->Nbonds, Cmol->Nres);
	gzprintf(outmol2, "SMALL\n");
	gzprintf(outmol2, "USER_CHARGES\n");
	gzprintf(outmol2, "Similarity_Index: %7.3f\n", cc);
	gzprintf(outmol2, "@<TRIPOS>ATOM\n");
	int i=0;
	unsigned resn=0;

	while(resn < Cmol->residue_pointer.size()-1){
		while(i < Cmol->residue_pointer[resn+1]-1){
			if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
				gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			else {
				gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
			}
			i++;
		}
		resn++;
	}
	while(i < Cmol->N){
		if (int(Cmol->sybyl_atoms.size()) == Cmol->N){
			gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->sybyl_atoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		else{
			gzprintf(outmol2, "%7d %-3.3s      %9.4f %9.4f %9.4f %-5.5s %5d %5s %8.4f\n", i+1, Cmol->atomnames[i].c_str(), xyz[i][0], xyz[i][1], xyz[i][2], Cmol->amberatoms[i].c_str(), resn+1, Cmol->resnames[resn].c_str(), Cmol->charges[i]);
		}
		i++;
	}

	gzprintf(outmol2, "@<TRIPOS>BOND\n");

	for (unsigned j=0; j<Cmol->bonds.size(); j++){
		gzprintf(outmol2, "%6d%6.6s%6.6s%5.5s\n", j+1, Cmol->bonds[j][0].c_str(), Cmol->bonds[j][1].c_str(),Cmol->bonds[j][2].c_str());
	}
}


Writer::~Writer(void){
}

void Writer::close(void){
	fprintf(logfile,"**************************************************************************************\n");
	fclose(logfile);
	gzclose(outmol2);
}

