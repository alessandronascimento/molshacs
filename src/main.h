/*
 * main.h
 *
 *  Created on: 19/08/2010
 *      Author: Nascimento
 */

#ifndef MAIN_H_
#define MAIN_H_
#ifdef HAS_GUI
#include <QApplication>
#include <QCoreApplication>
#include <QSplashScreen>
#include "../src/GUI/GUI.h"
#include <QTimer>
#endif
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include "Parser.h"
#include "RunEngine.h"
#include "ElSA.h"

using namespace std;

double Minimizer2::x, Minimizer2::y, Minimizer2::z;
vector<vector<double> > Minimizer2::new_xyz;
double Minimizer2::t1;
double Minimizer2::f_minimum;
Parser* Minimizer2::Input;
Mol2* Minimizer2::Cmol;
Mol2* Minimizer2::RefMol;
Mol2* Minimizer2::CompMol;
Grid* Minimizer2::Cgrid;

double Minimizer2::dij2;
double Minimizer2::alpha_i;
double Minimizer2::alpha_j;
double Minimizer2::pi;
double Minimizer2::pj;

#endif /* MAIN_H_ */
