/********************************************************************************
** Form generated from reading UI file 'spca_gui.ui'
**
** Created: Tue Jul 19 20:52:47 2011
**      by: Qt User Interface Compiler version 4.7.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SPCA_GUI_H
#define UI_SPCA_GUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SPCA_GUIClass
{
public:
    QWidget *centralwidget;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *SPCA_GUIClass)
    {
        if (SPCA_GUIClass->objectName().isEmpty())
            SPCA_GUIClass->setObjectName(QString::fromUtf8("SPCA_GUIClass"));
        SPCA_GUIClass->resize(800, 600);
        centralwidget = new QWidget(SPCA_GUIClass);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        SPCA_GUIClass->setCentralWidget(centralwidget);
        menubar = new QMenuBar(SPCA_GUIClass);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 800, 21));
        SPCA_GUIClass->setMenuBar(menubar);
        statusbar = new QStatusBar(SPCA_GUIClass);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        SPCA_GUIClass->setStatusBar(statusbar);

        retranslateUi(SPCA_GUIClass);

        QMetaObject::connectSlotsByName(SPCA_GUIClass);
    } // setupUi

    void retranslateUi(QMainWindow *SPCA_GUIClass)
    {
        SPCA_GUIClass->setWindowTitle(QApplication::translate("SPCA_GUIClass", "MainWindow", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class SPCA_GUIClass: public Ui_SPCA_GUIClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SPCA_GUI_H
