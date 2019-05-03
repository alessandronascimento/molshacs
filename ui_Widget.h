/********************************************************************************
** Form generated from reading UI file 'Widget.ui'
**
** Created: Tue Jul 19 20:52:47 2011
**      by: Qt User Interface Compiler version 4.7.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_WIDGET_H
#define UI_WIDGET_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHeaderView>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_WidgetClass
{
public:

    void setupUi(QWidget *WidgetClass)
    {
        if (WidgetClass->objectName().isEmpty())
            WidgetClass->setObjectName(QString::fromUtf8("WidgetClass"));
        WidgetClass->resize(400, 300);

        retranslateUi(WidgetClass);

        QMetaObject::connectSlotsByName(WidgetClass);
    } // setupUi

    void retranslateUi(QWidget *WidgetClass)
    {
        WidgetClass->setWindowTitle(QApplication::translate("WidgetClass", "Widget", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class WidgetClass: public Ui_WidgetClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_H
