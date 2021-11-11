/*
 * main.cpp
 *  Created on: 19/08/2010
 *  Author: Nascimento
 */
/*
 **********************************************************************************************
                MolShaCS - Molecular Shape and Charge Similarity v. 1.0

This software is dedicated to structural analysis of potential interaction sites given a
dataset of superposed ligands.

                 ReWritten by Alessandro S. Nascimento - june/2010
 ***********************************************************************************************
 */

#include "../src/main.h"

using namespace std;

int main (int argc, char* argv[]) {

	if (argc > 1) {
		Parser Input(argv[1]);
		RunEngine Engine;
		Engine.run(Input);
		return(0);
	}
	else {

#ifdef HAS_GUI

        QApplication app(argc, argv);
        QCoreApplication::setOrganizationName(ORGANIZATION);
//		app.setStyle("gtk");
		app.setApplicationName(NAME);
		app.setApplicationVersion(VERSION);


		QSplashScreen *splash = new QSplashScreen;
		splash->setPixmap(QPixmap(":/elsa.png"));
		splash->show();

		splash->showMessage(app.organizationName());

		GUI w(argc, argv);
		QIcon windowIcon(":/elsa.png");
//		w.setWindowIcon(windowIcon);

		QTimer::singleShot(2000, splash, SLOT(close()));
		QTimer::singleShot(2000, &w, SLOT(show()));

		return app.exec();

#endif
	}
}
