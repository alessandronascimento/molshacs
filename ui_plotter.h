/********************************************************************************
** Form generated from reading UI file 'plotter.ui'
**
** Created: Tue Jul 19 20:52:47 2011
**      by: Qt User Interface Compiler version 4.7.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PLOTTER_H
#define UI_PLOTTER_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_plotterClass
{
public:

    void setupUi(QWidget *plotterClass)
    {
        if (plotterClass->objectName().isEmpty())
            plotterClass->setObjectName(QString::fromUtf8("plotterClass"));
        plotterClass->resize(400, 300);

        retranslateUi(plotterClass);

        QMetaObject::connectSlotsByName(plotterClass);
    } // setupUi

    void retranslateUi(QWidget *plotterClass)
    {
        plotterClass->setWindowTitle(QApplication::translate("plotterClass", "plotter", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class plotterClass: public Ui_plotterClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PLOTTER_H
