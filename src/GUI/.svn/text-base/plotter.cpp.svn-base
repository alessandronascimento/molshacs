#include "plotter.h"

plotter::plotter(QWidget *parent)
: QWidget(parent)
{
	setBackgroundRole(QPalette::Dark);
	setAutoFillBackground(true);
	setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	setFocusPolicy(Qt::StrongFocus);
	QString filename = QFileDialog::getOpenFileName(this, tr("Choose a Results File"),
			"", tr("DAT Files (*.dat)"));
	QFile graph_file(filename);
	int count = 0;
	if(graph_file.open(QIODevice::ReadOnly|QIODevice::Text)){
		while (!graph_file.atEnd()){
			QString line = graph_file.readLine();
			QStringList fields = line.split(';');
			if (fields.size() == 3){
				x.push_back(fields.takeFirst().toDouble());
				fmin.push_back(fields.takeFirst().toDouble());
				si.push_back(fields.takeFirst().toDouble());
				count++;
			}
		}
	}
}

plotter::~plotter()
{
}

void plotter::paintEvent(QPaintEvent *){
    if (x.size() > 0) {
        QPolygonF data_si(x.size());
        QPainter painter(this);
        painter.setRenderHint(QPainter::Antialiasing, true);
        QRect rect(Margin, Margin, width()-2 * Margin, height()-2 * Margin);

        int xmin=1, xmax=x.size();
        double ymin=-1.0, ymax=1.0; 			// util para desenhar o grid.

        for (unsigned i=0; i<x.size(); i++){
            double dx = x[i] - xmin;
            double dy = si[i] - ymin;
            double xd = rect.left() + (dx * (rect.width()-1)/(xmax-xmin));
            double yd = rect.bottom() - (dy * (rect.height()-1)/(ymax-ymin));
           data_si[i]= QPointF(xd, yd);
        }

        int numXticks = int(xmax/2);
        int numYticks = 10;

        for (int i=0; i<=numXticks; i++){
            int x = rect.left()  + ( i* rect.width() -1 )/numXticks;
            double label = xmin + (i  * (xmax-xmin)/numXticks);
            painter.setPen(palette().dark().color().light());
            painter.drawLine(x, rect.top(), x, rect.bottom());
            painter.setPen(palette().light().color());
            painter.drawLine(x, rect.bottom(), x, rect.bottom()+5);
            painter.drawText(x-50, rect.bottom()+5, 100, 15,
                             Qt::AlignHCenter | Qt::AlignTop,
                             QString::number(label));

        }

        for (int j=0; j<=numYticks; j++){
            int y = rect.bottom() - (j*rect.height()-1)/numYticks;
            double label = (ymin + (j* ((ymax-ymin)/numYticks)));
            painter.setPen(palette().dark().color().light());
            painter.drawLine(rect.left(), y, rect.right(), y);
            painter.setPen(palette().light().color());
            painter.drawLine(rect.left()-5, y, rect.left(), y);
            painter.drawText(rect.left()-Margin, y-10, Margin-5, 20,
                             Qt::AlignRight | Qt::AlignVCenter,
                             QString::number(label));
        }

        painter.setPen(Qt::red);
        painter.drawText(25, rect.bottom()+30, 100, 15, Qt::AlignHCenter | Qt::AlignTop, "SI");

        painter.setPen(Qt::white);

        painter.drawRect(rect.adjusted(0, 0, -1, -1));
        painter.setClipRect(rect.adjusted(+1, +1, -1, -1));

        painter.setPen(Qt::red);
        painter.drawPolyline(data_si);

    }
}
