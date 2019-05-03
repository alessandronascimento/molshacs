#ifndef PLOTTER_H
#define PLOTTER_H

#include <vector>
#include <Qwidget>
#include <QFile>
#include <QPainter>
#include <QFileDialog>
#include<QPaintEvent>

using namespace std;

class plotter : public QWidget
{
    Q_OBJECT

public:
    plotter(QWidget *parent = 0);
    ~plotter();
    vector<double> x, fmin, si;
    enum { Margin = 50 };

private:
    void paintEvent(QPaintEvent* event);
};

#endif // PLOTTER_H
