#include <QCoreApplication>
#include <QApplication>
#include <iostream>
#include <armadillo/include/armadillo>
#include "examples.hpp"
#include <string>
#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QChar>
#include <QTime>
#include <QDate>
#include <QVector>


#include <iostream>
#include <conio.h>
#include "Load_Files/include_cpp/libxl.h"

#include "kalmanfilter.hpp"
#include "qcustomplot.h"





int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);
    QApplication a(argc, argv);
    arma::arma_rng::set_seed_random();

    double Tm = 0.1;
    BasicMatrix A = { {1, Tm},
                       {0,  1} };

    BasicMatrix C = { {1, 0},
                      {0, 1} };

    BasicMatrix P = { {1, 0},
                      {0, 1} };

    BasicMatrix Q = { {0, 0},
                      {0, 0} };

    BasicMatrix R = { {1, 0},
                      {0, 1} };


//    BasicVector XK1 = arma::linspace(0, 100);
    BasicVector XK1 = BasicVector(101,1); XK1[0] = 0;
    BasicVector Time = BasicVector(101,1); Time[0] = 0;
    for (arma::uword i = 1 ; i <= 100 ; i++ ) {
        XK1[i] = XK1[i-1] + 1;
        Time[i] = Time[i-1] + 0.1;
    }

    BasicVector XK2 = XK1;  XK2.fill(10);

    //BasicVector RANDON = arma::randu<arma::vec>(XK1.n_rows);
    BasicVector RANDON = arma::randn<arma::vec>(XK1.n_rows);

    BasicVector ZK1 = XK1 + RANDON;
    BasicVector ZK2 = XK2 + RANDON;

//    std::cout << XK1 << "\n" << XK2 << std::endl;
//    //std::cout << RANDON << std::endl;
//    std::cout << ZK1 << "\n" << ZK2 << std::endl;
    //BasicMatrix XK; XK.insert_cols(XK1);

    BasicMatrix XK; XK.insert_cols(0,XK1); XK.insert_cols(1, XK2);
    BasicMatrix ZK; ZK.insert_cols(0,ZK1); ZK.insert_cols(1, ZK2);


    std::cout << XK << std::endl;
    std::cout << ZK << std::endl;
//    std::cout << Time << std::endl;


    KalmanFilter filter({0, 5}, Q, R);
    KalmanParameters params = filter.Configuration();
    params.set_A(A).set_C(C);
    filter.SetConfiguration(params);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK.n_rows,1);
    BasicVector XK2_corregida = BasicVector(ZK.n_rows,1);


    for (arma::uword i = 0 ; i <= ZK.n_rows -1 ; i++) {

        BasicVector input = ZK.row(i).t();
        BasicVector output = filter.update(input);
        XK1_corregida[i] = output.at(0);
        XK2_corregida[i] = output.at(1);

//        std::cout << input  << " - " << ZK.row(i).t() << std::endl;
//        std::cout << input << " \n-------- " << output << std::endl;

    }


    // Converting BasicVector to QVector

    QVector<double> VTime;
    QVector<double> VZK1, VXK1, VXKC1;
    QVector<double> VZK2, VXK2, VXKC2;


    for (arma::uword i = 0 ; i <= ZK.n_rows -1 ; i++) {

        VTime.append( Time(i) );

        VZK1.append( ZK(i, 0) );
        VXK1.append( XK(i, 0) );

        VZK2.append( ZK(i, 1) );
        VXK2.append( XK(i, 1) );

        VXKC1.append( XK1_corregida(i) );
        VXKC2.append( XK2_corregida(i) );

    }

//    qDebug() << VZK2.empty();
//    qDebug() << VXK2.empty();
//    qDebug() << VXKC2.empty();

/*  PLOTTING  */

    QCustomPlot plot;
    plot.addGraph();
    plot.graph(0)->setData(VTime, VZK1);
    plot.graph(0)->setPen(QColor(Qt::blue));

    plot.addGraph();
    plot.graph(1)->setData(VTime, VXK1);
    plot.graph(1)->setPen(QColor(Qt::black));

    plot.addGraph();
    plot.graph(2)->setData(VTime, VXKC1);
    plot.graph(2)->setPen(QColor(Qt::red));
    // give the axes some labels:
    plot.xAxis->setLabel("Time");
    plot.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot.xAxis->setRange(-1, 12);
    plot.yAxis->setRange(-5, 120);
    plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot.replot();
    plot.setMinimumSize(QSize(500, 500));
    plot.show();


    QCustomPlot plot2;
    plot2.addGraph();
    plot2.graph(0)->setData(VTime, VZK2);
    plot2.graph(0)->setPen(QColor(Qt::blue));

    plot2.addGraph();
    plot2.graph(1)->setData(VTime, VXK2);
    plot2.graph(1)->setPen(QColor(Qt::black));

    plot2.addGraph();
    plot2.graph(2)->setData(VTime, VXKC2);
    plot2.graph(2)->setPen(QColor(Qt::red));
    // give the axes some labels:
    plot2.xAxis->setLabel("Time");
    plot2.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot2.xAxis->setRange(-0.5, 10.5);
    plot2.yAxis->setRange(5, 20);
    plot2.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot2.replot();
    plot2.setMinimumSize(QSize(500, 500));
    plot2.show();

    return a.exec();
}
