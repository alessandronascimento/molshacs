/*
 * QtWriter.cpp
 *
 *  Created on: 22/07/2011
 *      Author: Nascimento
 */

#include "QtWriter.h"

QtWriter::QtWriter(Parser *Input, QPlainTextEdit* Ed) {
	Editor = Ed;
	logfile= fopen((Input->output+".log").c_str(), "w");
	outmol2 = gzopen((Input->output + ".mol2.gz").c_str(), "w");
}

void QtWriter::write_params(Parser *Input){
	fprintf(logfile,"**************************************************************************************\n");
	fprintf(logfile, "* %82s *\n", "Parsed parameters:");

	Editor->appendPlainText(QString("%1\n").arg("Parsed parameters:"));
	fprintf(logfile,"* %82s *\n", "");

	Editor->appendPlainText(QString("%1: %2").arg("refmol_mol2").arg(Input->refmol_mol2.c_str()));
	fprintf(logfile,"*%-20s %-63.63s*\n", "refmol_mol2", Input->refmol_mol2.c_str());

	Editor->appendPlainText(QString("%1: %2 %3").arg("Amber atomtypes?").arg(Input->mol2_aa).arg(""));
	fprintf(logfile, "*%-20s %-20d %42s*\n", "Amber atomtypes?", Input->mol2_aa, "");

//	Editor->appendPlainText(QString("%1: %2").arg("sampling").arg(Input->sampling));
//	fprintf(logfile, "*%-20s %-20.4f %42s*\n", "sampling", Input->sampling, "");

//	Editor->appendPlainText(QString("%1: %2").arg("box_size").arg(Input->box_size*2));
//	fprintf(logfile,"*%-20s %-20.2f %42s*\n", "box_size", Input->box_size*2, "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("align_molecules").arg(Input->align).arg(""));
	fprintf(logfile,"*%-20s %-d %42s*\n", "align_molecules", Input->align, "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("search step").arg(Input->step).arg(""));
	fprintf(logfile,"*%-20s %-20.3f %42s*\n", "search step", Input->step , "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("minimization tol").arg(Input->tol).arg(""));
	fprintf(logfile, "*%-20s %-20.2f %42s*\n", "minimization tol", Input->tol , "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("minimization steps").arg(Input->nsteps).arg(""));
	fprintf(logfile, "*%-20s %-20d %42s*\n", "minimization steps", Input->nsteps , "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("minimization delta").arg(Input->delta).arg(""));
	fprintf(logfile, "*%-20s %-20.4f %42s*\n", "minimization delta", Input->delta, "");

	Editor->appendPlainText(QString("%1: %2").arg("minimizer").arg(Input->minimizer.c_str()));
	fprintf(logfile,"*%-20s %-63s*\n", "minimizer", Input->minimizer.c_str());

	Editor->appendPlainText(QString("%1: %2 %3").arg("multi").arg(Input->multi).arg(""));
	fprintf(logfile, "*%-20s %-20d %42s*\n", "multi", Input->multi, "");

	Editor->appendPlainText(QString("%1: %2").arg("molsfile").arg(Input->molsfile.c_str()));
	fprintf(logfile, "*%-20s %-63s*\n", "molsfile", Input->molsfile.c_str());

	Editor->appendPlainText(QString("%1: %2").arg("output").arg(Input->output.c_str()));
	fprintf(logfile, "*%-20s %-63s*\n", "output", Input->output.c_str());

	Editor->appendPlainText(QString("%1: %2 %3").arg("write coordinates").arg(Input->write_pdb).arg(""));
	fprintf(logfile, "*%-20s %-20d %42s*\n", "write coordinates", Input->write_pdb, "");

	if (Input->use_write_coord_threshold){
		Editor->appendPlainText(QString("%1: %2 %3").arg("Writing threshold: ").arg(Input->write_coord_threshold).arg(""));
		fprintf(logfile, "*%-20s %-20.2f %42s*\n", "Writing threshold:", Input->write_coord_threshold, "");
	}

	Editor->appendPlainText(QString("%1: %2 %3").arg("electrostatic scale").arg(Input->elec_scale).arg(""));
	fprintf(logfile, "*%-20s %-20.2f %42s*\n", "electrostatic scale", Input->elec_scale, "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("VDW scale").arg(Input->vdw_scale).arg(""));
	fprintf(logfile, "*%-20s %-20.2f %42s*\n", "VDW scale", Input->vdw_scale, "");

	Editor->appendPlainText(QString("%1: %2 %3").arg("timeout").arg(Input->timeout).arg(""));
	fprintf(logfile, "*%-20s %-20d %42s*\n", "timeout", Input->timeout, "");

	if (Input->use_write_coord_threshold){
		fprintf(logfile, "*%-20s %-20f %42s*\n", "write_threshold", Input->write_coord_threshold, "");
		Editor->appendPlainText(QString("%1: %2 %3").arg("write_threshold").arg(Input->write_coord_threshold).arg(""));
	}

	Editor->appendPlainText(QString("%1").arg(""));
	fprintf(logfile, " %82s \n", "");

	Editor->appendPlainText(QString("%1").arg(""));

	Editor->appendPlainText("**************************************************************************************");
	fprintf(logfile,"**************************************************************************************\n");

	Editor->appendPlainText(QString("%1").arg(""));
	QApplication::processEvents();
}

void QtWriter::write_welcome(void){
	Editor->clear();
	Editor->appendPlainText("**************************************************************************************");
	Editor->appendPlainText("                                                                                       ");
	Editor->appendPlainText("                                   MolShaCS v. 1.0                                     ");
	Editor->appendPlainText("                      Molecular Shape and Charge Similarity                            ");
	Editor->appendPlainText("                                                                                      ");
	Editor->appendPlainText("  Written by Alessandro S. Nascimento-Jun/2010 -  alessandro.nascimento@ufabc.edu.br  ");
	Editor->appendPlainText("                                                                                      ");
	Editor->appendPlainText("                        Universidade Federal do ABC - Brazil                          ");
	Editor->appendPlainText("                                                                                      ");
	Editor->appendPlainText("                           http://nascimento.ufabc.edu.br/                            ");
	Editor->appendPlainText("                                                                                      ");
	Editor->appendPlainText("**************************************************************************************");
	Editor->appendPlainText("");

	fprintf(logfile, "**************************************************************************************\n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "                                   MolShaCS v. 1.0                                    \n");
	fprintf(logfile, "                       Molecular Shape and Charge Similarity                          \n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "  Written by Alessandro S. Nascimento-Jun/2010 -  alessandro.nascimento@ufabc.edu.br  \n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "                        Universidade Federal do ABC - Brazil                          \n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "                 http://nascimento.ufabc.edu.br/wordpress/?page_id=14                 \n");
	fprintf(logfile, "                                                                                      \n");
	fprintf(logfile, "**************************************************************************************\n");

	QApplication::processEvents();
}

void QtWriter::write_to_log(char info[82]){
	string str = info;
	Editor->appendPlainText(QString("%1").arg(str.c_str()));
	Editor->update();
	fprintf(logfile, "* %-82s *\n", info);
	QApplication::processEvents();
}

void QtWriter::write_to_log(void){
	Editor->appendPlainText("");
	Editor->update();
	fprintf(logfile, "* %-82s *\n", "");
	QApplication::processEvents();
}

void QtWriter::write_pdb(Mol *Cmol, vector<vector<double> >xyz, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "a");
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
	QApplication::processEvents();
}

void QtWriter::write_pdb(Mol2 *Cmol, vector<vector<double> >xyz, string outname){
	gzFile outpdb;
	outpdb = gzopen((outname+".pdb.gz").c_str(), "w");
	gzprintf(outpdb, "MDL\n");
	gzprintf(outpdb, "REMARK %s\n", Cmol->molname.c_str());
	int i=0;
	int resn=0;

	while (resn < int(Cmol->residue_pointer.size()-1)){
		while(i < Cmol->residue_pointer[resn+1]-1){
			gzprintf(outpdb, "%6s%5d %-4s%1s%3.3s%1s%4d%5s%8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM  ", i+1, Cmol->atomnames[i].c_str(), " ",Cmol->resnames[resn].c_str()," ",resn+1," ",xyz[i][0], xyz[i][1], xyz[i][2], 1.0, Cmol->charges[i]);
			i++;
		}
		resn++;
	}

	while (i < Cmol->N){
		gzprintf(outpdb, "%6s%5d %-4s%1s%3.3s%1s%4d%5s%8.3f%8.3f%8.3f%6.2f%6.2f\n", "ATOM  ", i+1, Cmol->atomnames[i].c_str(), " ",Cmol->resnames[resn].c_str()," ",resn+1," ",xyz[i][0], xyz[i][1], xyz[i][2], 1.0, Cmol->charges[i]);
		i++;
	}
	gzprintf(outpdb, "TER\n");
	gzprintf(outpdb, "ENDMDL\n");
	gzclose(outpdb);
	QApplication::processEvents();
}

void QtWriter::writeMol2(Mol2* Cmol, vector<vector<double> >xyz, double cc){
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

QtWriter::~QtWriter(void){
	fprintf(logfile,"**************************************************************************************\n");
	QApplication::processEvents();
	fclose(logfile);
	gzclose(outmol2);
}
