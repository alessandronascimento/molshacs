#include "Widget.h"

using namespace std;

Widget::Widget(Parser *Input, QPlainTextEdit* Ed)
{
	Editor = Ed;
	INP = Input;
	mainLayout = new QVBoxLayout(this);
	editLayout = new QGridLayout;
	buttonLayout = new QHBoxLayout;
	messageLayout = new QHBoxLayout;
	progress_layout = new QVBoxLayout;
	mainLayout->addLayout(editLayout);
	mainLayout->addStretch();
	mainLayout->addLayout(buttonLayout);

	rmol2 = new QPushButton(tr("Reference molecule file (MOL2 format): "));
	choose_ref_mol2 = new QLineEdit(tr("Choose File..."));
	refmol_mol2_ql = new QLabel;

	mol2_aa_lab = new QLabel(tr("MOL2 Amber Atoms?"));
	mol2_aa = new QComboBox;
	mol2_aa->addItem(tr("True"));
	mol2_aa->addItem(tr("False"));
	mol2_aa->setCurrentIndex(1);

	alignmol_label = new QLabel(tr("Align Molecules ? "));
	alignmol_combo = new QComboBox;
	alignmol_combo->addItem(tr("True"));
	alignmol_combo->addItem(tr("False"));

	step_label = new QLabel(tr("Searching Step: "));
	step_spin = new QDoubleSpinBox;
	step_spin->setMaximum(180.0);
	step_spin->setMinimum(0.1);
	step_spin->setSingleStep(0.1);
	step_spin->setDecimals(1);
	step_spin->setValue(1.0);

	tol_label = new QLabel(tr("Searching Tolerance: "));
	tol_spin = new QDoubleSpinBox;
	tol_spin->setMaximum(5.0);
	tol_spin->setMinimum(0.00001);
	tol_spin->setSingleStep(0.00001);
	tol_spin->setDecimals(5);
	tol_spin->setValue(0.01);

	ntries_label = new QLabel(tr("Searching Cycles: "));
	ntries_spin = new QSpinBox;
	ntries_spin->setMaximum(1000);
	ntries_spin->setMinimum(1);
	ntries_spin->setValue(50);

	delta_label = new QLabel(tr("Gradient Delta: "));
	delta_spin = new QDoubleSpinBox;
	delta_spin->setMaximum(1.0);
	delta_spin->setMinimum(1E-10);
	delta_spin->setDecimals(10);
	delta_spin->setSingleStep(1E-5);
	delta_spin->setValue(1E-10);

	minimizer_label = new QLabel(tr("Minimizer: "));
	minimizer_combo = new QComboBox;
	minimizer_combo->addItem(tr("AUGLAG")); //0
	minimizer_combo->addItem(tr("MMA")); //1
	minimizer_combo->addItem(tr("ISRES")); //2
	minimizer_combo->addItem(tr("SUBPLEX")); //3
	minimizer_combo->addItem(tr("SIMPLEX")); //4
	minimizer_combo->addItem(tr("LBFGS2")); //5
	minimizer_combo->addItem(tr("STOGO")); //6
	minimizer_combo->addItem(tr("DIRECT-L")); //7
	minimizer_combo->addItem(tr("COBYLA")); //8
	minimizer_combo->setCurrentIndex(0);

	multifile_label = new QPushButton(tr("Comparing Molecule files (MOL2 format): "));
	choose_multifile = new QLineEdit(tr("Click to choose the files..."));

	output_label = new QLabel(tr("Output prefix: "));
	output_linedit = new QLineEdit(tr("MolShaCS"));

	writecoor_label = new QLabel(tr("Write coordinates ? "));
	writecoor_combo = new QComboBox;
	writecoor_combo->addItem("True");
	writecoor_combo->addItem("False");

	use_writing_threshold_label = new QLabel(tr("Use writing threshold?"));
	use_writing_threshold = new QComboBox;
	use_writing_threshold->addItem("True"); // 0
	use_writing_threshold->addItem("False"); // 1
	use_writing_threshold->setCurrentIndex(1);

	writting_threshold_label = new QLabel(tr("Writing threshold: "));
	writting_threshold = new QDoubleSpinBox;
	writting_threshold->setMinimum(-1.0);
	writting_threshold->setMaximum(1.0);
	writting_threshold->setDecimals(2);
	writting_threshold->setValue(0.85);

	timeout_label = new QLabel("Timeout: ");
	timeout_spin = new QSpinBox;
	timeout_spin->setMaximum(2147483647);
	timeout_spin->setMinimum(1);
	timeout_spin->setValue(120);

	elecscale_label = new QLabel("Electrostatic Scale: ");
	vdwscale_label = new QLabel("VDW Scale: ");
	elecscale_spin = new QDoubleSpinBox;
	vdwscale_spin = new QDoubleSpinBox;
	elecscale_spin->setMaximum(1.0);
	elecscale_spin->setMinimum(0.0);
	elecscale_spin->setDecimals(2);
	elecscale_spin->setValue(1.0);
	vdwscale_spin->setMaximum(1.0);
	vdwscale_spin->setMinimum(0.0);
	vdwscale_spin->setDecimals(2);
	vdwscale_spin->setValue(1.0);

	progressbar = new QProgressBar;
	progressbar->setRange(0, 100);
	progressbar->setValue(0);


	editLayout->addWidget(rmol2, 0, 0);
	editLayout->addWidget(choose_ref_mol2, 0, 1);
	editLayout->addWidget(mol2_aa_lab, 1, 0);
	editLayout->addWidget(mol2_aa, 1, 1);
	editLayout->addWidget(alignmol_label, 2, 0);
	editLayout->addWidget(alignmol_combo, 2, 1);
//	editLayout->addWidget(step_label, 3,0);
//	editLayout->addWidget(step_spin, 3, 1);
	editLayout->addWidget(tol_label, 4,0);
	editLayout->addWidget(tol_spin, 4, 1);
//	editLayout->addWidget(ntries_label, 5,0);
//	editLayout->addWidget(ntries_spin, 5, 1);
//	editLayout->addWidget(delta_label, 6,0);
//	editLayout->addWidget(delta_spin, 6, 1);
	editLayout->addWidget(minimizer_label, 7, 0);
	editLayout->addWidget(minimizer_combo, 7, 1);
	editLayout->addWidget(multifile_label, 8, 0);
	editLayout->addWidget(choose_multifile, 8, 1);
	editLayout->addWidget(output_label, 9,0);
	editLayout->addWidget(output_linedit, 9, 1);
	editLayout->addWidget(writecoor_label, 10,0);
	editLayout->addWidget(writecoor_combo, 10,1);
	editLayout->addWidget(use_writing_threshold_label, 11, 0);
	editLayout->addWidget(use_writing_threshold, 11, 1);
	editLayout->addWidget(writting_threshold_label, 12, 0);
	editLayout->addWidget(writting_threshold, 12, 1);
	editLayout->addWidget(timeout_label, 13,0);
	editLayout->addWidget(timeout_spin, 13,1);
	editLayout->addWidget(elecscale_label, 14,0);
	editLayout->addWidget(vdwscale_label, 15,0);
	editLayout->addWidget(elecscale_spin, 14,1);
	editLayout->addWidget(vdwscale_spin, 15,1);

	progress_label = new QLabel(tr("<h2>Progress</h2>"));
	progress_label->setAlignment(Qt::AlignHCenter);

	progress_layout->addWidget(progress_label, Qt::AlignCenter);
	progress_layout->addWidget(progressbar);

	mainLayout->addStretch();
	mainLayout->addLayout(progress_layout);
	mainLayout->addStretch();

	Runbutton = new QPushButton(tr("Run"));
	writebutton = new QPushButton(tr("Write INP File"));

	buttonLayout->addStretch();
	buttonLayout->addWidget(writebutton, Qt::AlignHCenter);
	buttonLayout->addWidget(Runbutton, Qt::AlignHCenter);

	Runbutton->setDefault(true);
	setWindowTitle(tr("MolShaCS by Alessandro Nascimento"));

	connect(rmol2, SIGNAL(clicked()), this, SLOT(choose_ref_mol2_file()));
//	connect(multifile_label, SIGNAL(clicked()), this, SLOT(choose_dat_file()));
	connect(multifile_label, SIGNAL(clicked()), this, SLOT(choose_comparing_mol2()));
	connect(Runbutton, SIGNAL(clicked()), this, SLOT(run()));
	connect(writebutton, SIGNAL(clicked()), this, SLOT(write_inp_file()));
}

void Widget::choose_ref_mol2_file(void){
	QString mol2 = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("MOL2 Files (*.mol2)"));
	choose_ref_mol2->setText(mol2.toUtf8());
	editLayout->addWidget(refmol_mol2_ql, 0, 1);
}

void Widget::choose_dat_file(void){
	QString datfile = QFileDialog::getOpenFileName(this, tr("Choose a File"), "", tr("DAT Files (*.dat)"));
	choose_multifile->setText(datfile.toUtf8());
	editLayout->addWidget(choose_multifile, 10, 1);
}

void Widget::choose_comparing_mol2(void){
	QStringList ql = QFileDialog::getOpenFileNames(this, tr("Choose MOL2 Files"), "", tr("MOL2 Files (*.mol2)"));
	INP->comparing_molecules = ql;
	QString tmp = QString("%1 Molecules read!").arg(INP->comparing_molecules.size());
	choose_multifile->setText(tmp);
	editLayout->addWidget(choose_multifile, 10, 1);
}

void Widget::write_inp_file(void){
	QFile output;
	output.setFileName("MolShaCS.inp");
	if(!output.open(QIODevice::WriteOnly|QIODevice::Text)) {
		printf("Failed to open file molshacs.inp for writing");
	}

	QString line = ("refmol_mol2 " + choose_ref_mol2->text() + "\n");
	output.write(line.toUtf8());
	line = 	choose_ref_mol2->text();
	INP->refmol_mol2 = line.toStdString();

	QString combo_ans;
	if (mol2_aa->currentIndex() == 0){
		combo_ans = "yes";
		INP->mol2_aa = true;
	}
	else {
		combo_ans = "no";
		INP->mol2_aa = false;
	}
	line = ("mol2_aa " + combo_ans + "\n");
	output.write(line.toUtf8());

	if(alignmol_combo->currentIndex() == 0){
		combo_ans="yes";
		INP->align = true;
	}
	else {
		combo_ans="no";
		INP->align = false;
	}
	line = ("align_molecules " + combo_ans + "\n");
	output.write(line.toUtf8());

	line = ("step " + QString::number(step_spin->value()) + "\n");
	INP->step = step_spin->value();
	output.write(line.toUtf8());

	line = ("tol " + QString::number(tol_spin->value()) + "\n");
	INP->tol = tol_spin->value();
	output.write(line.toUtf8());

	line = ("ntries " + QString::number(ntries_spin->value()) + "\n");
	INP->nsteps = ntries_spin->value();
	output.write(line.toUtf8());

	line = ("delta " + QString::number(delta_spin->value()) + "\n");
	INP->delta = delta_spin->value();
	output.write(line.toUtf8());

	switch(minimizer_combo->currentIndex()){
	case 0:
		combo_ans = "ln_auglag";
		INP->minimizer = "ln_auglag";
		break;
	case 1:
		combo_ans = "nlopt_mma";
		INP->minimizer = "nlopt_mma";
		break;
	case 2:
		combo_ans = "nlopt_isres";
		INP->minimizer = "nlopt_isres";
		break;
	case 3:
		combo_ans = "nlopt_subplex";
		INP->minimizer = "nlopt_subplex";
		break;
	case 4:
		combo_ans = "nlopt_simplex";
		INP->minimizer = "nlopt_simplex";
		break;
	case 5:
		combo_ans = "nlopt_lbfgs2";
		INP->minimizer = "nlopt_lbfgs2";
		break;
	case 6:
		combo_ans = "nlopt_stogo";
		INP->minimizer = "nlopt_stogo";
		break;
	case 7:
		combo_ans = "nlopt_direct-l";
		INP->minimizer = "nlopt_direct-l";
		break;
	case 8:
		combo_ans = "nlopt_cobyla";
		INP->minimizer = "nlopt_cobyla";
		break;
	}
	line = ("minimizer " + combo_ans + "\n");
	output.write(line.toUtf8());

	line = ("multimode yes\n");
	INP->multi = true;
	output.write(line.toUtf8());


	QFile multimol;
	multimol.setFileName("multimol.dat");
	if(!multimol.open(QIODevice::WriteOnly|QIODevice::Text)) {
		printf("Failed to open file multimol.dat for writing");
	}
	for (signed i=0; i<INP->comparing_molecules.size(); i++){
		multimol.write(((INP->comparing_molecules.at(i)).toUtf8()) + "\n");
	}
	multimol.write("EOF"); // Just to make sure!! ;)
	multimol.close();
	line = ("multimol multimol.dat \n");
	output.write(line.toUtf8());

	line = ("output_prefix " + output_linedit->text() + "\n");
	INP->output = output_linedit->text().toStdString();
	output.write(line.toUtf8());

	if(writecoor_combo->currentIndex() == 0){
		combo_ans="yes";
		INP->write_pdb = true;
	}
	else {
		combo_ans="no";
		INP->write_pdb = false;
	}
	line = ("write_coordinates " + combo_ans + "\n");
	output.write(line.toUtf8());

	if (use_writing_threshold->currentIndex() == 0){
		INP->use_write_coord_threshold = true;
		INP->write_coord_threshold = float(writting_threshold->value());
	}
	line = ("write_coord_threshold " + QString::number(writting_threshold->value()) + " \n");
	output.write(line.toUtf8());

	line = ("timeout " + QString::number(timeout_spin->value()) + "\n");
		INP->timeout = timeout_spin->value();
	output.write(line.toUtf8());

	line = ("elec_scale " + QString::number(elecscale_spin->value()) + "\n");
	INP->elec_scale = elecscale_spin->value();
	output.write(line.toUtf8());

	line = ("vdw_scale " + QString::number(vdwscale_spin->value()) + "\n");
	INP->vdw_scale = vdwscale_spin->value();
	output.write(line.toUtf8());

	output.close();

	done = new QLabel(tr("Parameter file MolShaCS.inp written successfully!"));
	editLayout->addWidget(done, 18, 0);
}


void Widget::set_parameters(void){
	INP->refmol_mol2 = choose_ref_mol2->text().toStdString();

	if (mol2_aa->currentIndex() == 0){
		INP->mol2_aa = true;
	}

	else {
		INP->mol2_aa = false;
	}

	if(alignmol_combo->currentIndex() == 0){
		INP->align = true;
	}
	else {
		INP->align = false;
	}
	INP->step = step_spin->value();
	INP->tol = tol_spin->value();
	INP->nsteps = ntries_spin->value();
	INP->delta = delta_spin->value();
	switch(minimizer_combo->currentIndex()){
	case 0:
		INP->minimizer = "nlopt_ln_auglag";
		break;
	case 1:
		INP->minimizer = "nlopt_mma";
		break;
	case 2:
		INP->minimizer = "nlopt_isres";
		break;
	case 3:
		INP->minimizer = "nlopt_subplex";
		break;
	case 4:
		INP->minimizer = "nlopt_simplex";
		break;
	case 5:
		INP->minimizer = "nlopt_lbfgs2";
		break;
	case 6:
		INP->minimizer = "nlopt_stogo";
		break;
	case 7:
		INP->minimizer = "nlopt_direct-l";
		break;
	case 8:
		INP->minimizer = "nlopt_cobyla";
		break;
	}
	INP->multi = true;
	INP->molsfile = choose_multifile->text().toStdString();
	INP->output = output_linedit->text().toStdString();
	if(writecoor_combo->currentIndex() == 0){
		INP->write_pdb = true;
		if (use_writing_threshold->currentIndex() == 0){
			INP->use_write_coord_threshold = true;
			INP->write_coord_threshold = float(writting_threshold->value());
		}
	}

	else {
		INP->write_pdb = false;
	}
	INP->timeout = timeout_spin->value();
	INP->elec_scale = elecscale_spin->value();
	INP->vdw_scale = vdwscale_spin->value();
}

Widget::~Widget()
{
}

void Widget::run(){
	set_parameters();
	MolShaCS = new RunEngine;
	MolShaCS->run(*INP, Editor, progressbar);
	delete MolShaCS;
}
