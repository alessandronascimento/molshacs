#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <iostream>
#include <QFileDialog>
#include <QPlainTextEdit>
#include <QEventLoop>
#include <QStringList>
#include <QProgressBar>
#include "../Parser.h"
#include "../RunEngine.h"
#include <cstdio>

using namespace std;

class Widget : public QWidget
{
    Q_OBJECT

public:
	Widget(Parser *Input, QPlainTextEdit* Ed);
    ~Widget();

    QPlainTextEdit* Editor;
    Parser* INP;

    QGridLayout* editLayout;
    QVBoxLayout* mainLayout;
    QHBoxLayout* buttonLayout;
    QHBoxLayout* messageLayout;

    QVBoxLayout* progress_layout;
    QProgressBar *progressbar;
    QLabel* progress_label;

    QLabel* done;

    QPushButton* rmol2;
    QComboBox* mol2_aa;
    QLabel* mol2_aa_lab;

    QPushButton* Runbutton;
    QPushButton* writebutton;

    QLabel* refmol_mol2_ql;
    QLabel* refmol_inpcrd;

    QLineEdit* choose_ref_mol2;

    QLabel* alignmol_label;
    QComboBox* alignmol_combo;

    QDoubleSpinBox* step_spin;
    QLabel* step_label;

    QDoubleSpinBox* tol_spin;
    QLabel* tol_label;

    QSpinBox* ntries_spin;
    QLabel* ntries_label;

    QDoubleSpinBox* delta_spin;
    QLabel* delta_label;

    QLabel* minimizer_label;
    QComboBox* minimizer_combo;

    QPushButton* multifile_label;
    QLineEdit* choose_multifile;

    QLabel* output_label;
    QLineEdit* output_linedit;

    QLabel* writecoor_label;
    QComboBox* writecoor_combo;

    QLabel* use_writing_threshold_label;
    QComboBox* use_writing_threshold;

    QLabel* writting_threshold_label;
    QDoubleSpinBox* writting_threshold;

    QLabel* timeout_label;
    QSpinBox* timeout_spin;

    QLabel* elecscale_label;
    QLabel* vdwscale_label;
    QDoubleSpinBox* elecscale_spin;
    QDoubleSpinBox* vdwscale_spin;

    RunEngine* MolShaCS;

    void set_parameters(void);

public slots:
    void write_inp_file(void);

private slots:
    void choose_ref_mol2_file(void);
    void choose_dat_file(void);
    void choose_comparing_mol2(void);
    void run(void);

};

#endif // WIDGET_H
