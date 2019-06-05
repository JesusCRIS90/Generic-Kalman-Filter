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

#include <QFile>
#include <QByteArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QDebug>

/**/
double ECM_2(QVector<double>& input1, QVector<double>& input2){
    BasicVector in1(input1.toStdVector());
    BasicVector in2(input2.toStdVector());
    BasicVector result( in1 - in2 );
    BasicVector r = result%result;          // Schur Product
//    std::cout << r << std::endl;
    return ( qSqrt(arma::mean(r)) );
}
/**/

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);
    QApplication a(argc, argv);

// -------------------- EXAMPLE 1 -------------------------------------------
/*
    arma::arma_rng::set_seed_random();

    double Tm = 0.1;
    BasicMatrix A = { {1, Tm},
                       {0,  1} };

    BasicMatrix C = { {1, 0} };

    BasicMatrix P = { {1, 0},
                      {0, 1} };

    BasicMatrix Q = { {0, 0},
                      {0, 0} };

//    BasicMatrix R = { {1, 0},
//                      {0, 1} };

    BasicMatrix R = { 1 };

    BasicVector XK1 = BasicVector(101,1); XK1[0] = 0;
    BasicVector Time = BasicVector(101,1); Time[0] = 0;

    for (arma::uword i = 1 ; i <= 100 ; i++ ) {
        XK1[i] = XK1[i-1] + 1;
        Time[i] = Time[i-1] + 0.1;
    }

    BasicVector XK2 = XK1;  XK2.fill(10);

    BasicVector RANDON = arma::randn<arma::vec>(XK1.n_rows);

    BasicVector ZK1 = XK1 + RANDON;
    BasicVector ZK2 = XK2 + RANDON;

    BasicMatrix XK; XK.insert_cols(0,XK1); XK.insert_cols(1, XK2);
    BasicMatrix ZK; ZK.insert_cols(0,ZK1); ZK.insert_cols(1, ZK2);

    BasicVector X0 = {0, 5};
    BasicVector Z0 = {5};
    BasicVector U0 = {1};

    KalmanFilter filter(X0, Z0, Q, R);
    KalmanParameters params = filter.Configuration();
    params.set_A(A).set_C(C);
    filter.SetConfiguration(params);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK.n_rows,1);
    BasicVector XK2_corregida = BasicVector(ZK.n_rows,1);


    for (arma::uword i = 0 ; i <= ZK.n_rows -1 ; i++) {

//        BasicVector input = ZK.row(i).t();
        BasicVector input { ZK1.at(i) };
        BasicVector output = filter.update(input);
        XK1_corregida[i] = output.at(0);
        XK2_corregida[i] = output.at(1);
    }


//   Converting BasicVector to QVector

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


//  PLOTTING

    QCustomPlot plot;

    plot.legend->setVisible(true);
    QFont legendFont;
    legendFont.setPointSize(9);
    plot.legend->setFont(legendFont);
    plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);


//    plot.addGraph(plot.yAxis, plot.xAxis);
    plot.graph(0)->setData(VTime, VZK1);
    plot.graph(0)->setPen(QColor(Qt::blue));
    plot.graph(0)->setName("Medida Sensor");

    plot.addGraph();
    plot.graph(1)->setData(VTime, VXK1);
    plot.graph(1)->setPen(QColor(Qt::black));
    plot.graph(1)->setName("Medida Real");

    plot.addGraph();
    plot.graph(2)->setData(VTime, VXKC1);
    plot.graph(2)->setPen(QColor(Qt::red));
    plot.graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot.graph(2)->setName("Medida Filtro Kalman");

    // give the axes some labels:
    plot.xAxis->setLabel("Time");
    plot.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot.xAxis->setRange(-1, 12);
    plot.yAxis->setRange(-5, 120);
    plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot.replot();
    plot.setMinimumSize(QSize(500, 500));
    plot.setWindowTitle("Distancia recorrida del Tren");
    plot.show();


    qDebug() << ECM_(VXK1, VZK1);
    qDebug() << ECM_(VXK1, VXKC1);

    QCustomPlot plot2;

    plot2.legend->setVisible(true);
    plot2.legend->setFont(legendFont);
    plot2.legend->setBrush(QBrush(QColor(255,255,255,230)));
    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    plot2.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

//    plot2.addGraph();
//    plot2.graph(0)->setData(VTime, VZK2);
//    plot2.graph(0)->setPen(QColor(Qt::blue));
//    plot2.graph(0)->setName("Medida Sensor");

    plot2.addGraph();
    plot2.graph(0)->setData(VTime, VXK2);
    plot2.graph(0)->setPen(QColor(Qt::black));
    plot2.graph(0)->setName("Medida Real");

    plot2.addGraph();
    plot2.graph(1)->setData(VTime, VXKC2);
    plot2.graph(1)->setPen(QColor(Qt::red));
    plot2.graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot2.graph(1)->setName("Medida Filtro Kalman");

    // give the axes some labels:
    plot2.xAxis->setLabel("Time");
    plot2.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot2.xAxis->setRange(-0.5, 10.5);
    plot2.yAxis->setRange(5, 20);
    plot2.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot2.replot();
    plot2.setMinimumSize(QSize(500, 500));
    plot2.setWindowTitle("Velocidad del Tren");
    plot2.show();

    qDebug() << "--------------------------";
    qDebug() << ECM_(VXK2, VZK2);
    qDebug() << ECM_(VXK2, VXKC2);


*/

// --------------------------------------------------------------------------


// ------------------- EXAMPLE 2 --------------------------------------------
/*

    arma::arma_rng::set_seed_random();

    double Tm = 0.1;
    BasicMatrix A = { {1, Tm},
                       {0,  1} };

    BasicMatrix B = { {1, 0},
                      {0, 1} };

//    BasicMatrix B = { {0, 0} };

//    B = B.t();

    BasicMatrix C = { {1, 0} };

    BasicMatrix P = { {1, 0},
                      {0, 1} };

    BasicMatrix Q = { {0, 0},
                      {0, 0} };

//    BasicMatrix R = { {1, 0},
//                      {0, 1} };

    BasicMatrix R = { 1 };

    BasicVector XK1 = BasicVector(101,1); XK1[0] = 0;
    BasicVector Time = BasicVector(101,1); Time[0] = 0;

    for (arma::uword i = 1 ; i <= 100 ; i++ ) {
        XK1[i] = XK1[i-1] + 1;
        Time[i] = Time[i-1] + 0.1;
    }

    BasicVector XK2 = XK1;  XK2.fill(10);

    BasicVector RANDON = arma::randn<arma::vec>(XK1.n_rows);

    BasicVector ZK1 = XK1 + RANDON;
    BasicVector ZK2 = XK2 + RANDON;

    BasicMatrix XK; XK.insert_cols(0,XK1); XK.insert_cols(1, XK2);
    BasicMatrix ZK; ZK.insert_cols(0,ZK1); ZK.insert_cols(1, ZK2);

    BasicVector X0 = {0, 5};
    BasicVector Z0 = { ZK1[0] };
    BasicVector U0 = {0, 0};

    KalmanFilter filter(A, B, C, Q, R, X0, Z0, P, {}, U0);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK.n_rows,1);
    BasicVector XK2_corregida = BasicVector(ZK.n_rows,1);


    for (arma::uword i = 0 ; i <= ZK.n_rows -1 ; i++) {

//        BasicVector output = filter.GetState();
//        BasicVector xnext = A*output + B*U0;
//        BasicVector Zknext = C*xnext + RANDON[i];
//        filter.update(Zknext);
//
        //BasicVector input = ZK.row(i).t();
        BasicVector input { ZK1.at(i) };
        filter.update(input);

        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
        XK2_corregida[i] = output.at(1);
    }


//   Converting BasicVector to QVector

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


//  PLOTTING

    QCustomPlot plot;

    plot.legend->setVisible(true);
    QFont legendFont;
    legendFont.setPointSize(9);
    plot.legend->setFont(legendFont);
    plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);


//    plot.addGraph(plot.yAxis, plot.xAxis);
    plot.addGraph();
    plot.graph(0)->setData(VTime, VZK1);
    plot.graph(0)->setPen(QColor(Qt::blue));
    plot.graph(0)->setName("Medida Sensor");

    plot.addGraph();
    plot.graph(1)->setData(VTime, VXK1);
    plot.graph(1)->setPen(QColor(Qt::black));
    plot.graph(1)->setName("Medida Real");

    plot.addGraph();
    plot.graph(2)->setData(VTime, VXKC1);
    plot.graph(2)->setPen(QColor(Qt::red));
    plot.graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot.graph(2)->setName("Medida Filtro Kalman");

    // give the axes some labels:
    plot.xAxis->setLabel("Time");
    plot.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot.xAxis->setRange(-1, 12);
    plot.yAxis->setRange(-5, 120);
    plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot.replot();
    plot.setMinimumSize(QSize(500, 500));
    plot.setWindowTitle("Distancia recorrida del Tren");
    plot.show();


    qDebug() << ECM_2(VXK1, VZK1);
    qDebug() << ECM_2(VXK1, VXKC1);

    QCustomPlot plot2;

    plot2.legend->setVisible(true);
    plot2.legend->setFont(legendFont);
    plot2.legend->setBrush(QBrush(QColor(255,255,255,230)));
    // by default, the legend is in the inset layout of the main axis rect. So this is how we access it to change legend placement:
    plot2.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

//    plot2.addGraph();
//    plot2.graph(0)->setData(VTime, VZK2);
//    plot2.graph(0)->setPen(QColor(Qt::blue));
//    plot2.graph(0)->setName("Medida Sensor");

    plot2.addGraph();
    plot2.graph(0)->setData(VTime, VXK2);
    plot2.graph(0)->setPen(QColor(Qt::black));
    plot2.graph(0)->setName("Medida Real");

    plot2.addGraph();
    plot2.graph(1)->setData(VTime, VXKC2);
    plot2.graph(1)->setPen(QColor(Qt::red));
    plot2.graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot2.graph(1)->setName("Medida Filtro Kalman");

    // give the axes some labels:
    plot2.xAxis->setLabel("Time");
    plot2.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot2.xAxis->setRange(-0.5, 10.5);
    plot2.yAxis->setRange(5, 20);
    plot2.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    plot2.replot();
    plot2.setMinimumSize(QSize(500, 500));
    plot2.setWindowTitle("Velocidad del Tren");
    plot2.show();

    qDebug() << "--------------------------";
    qDebug() << ECM_2(VXK2, VZK2);
    qDebug() << ECM_2(VXK2, VXKC2);




//*/

// --------------------------------------------------------------------------


// ------------------- EXAMPLE 3 --------------------------------------------
/*
    arma::arma_rng::set_seed_random();

 //   double Tm = 0.1;
    BasicMatrix A = { 1 };

    BasicMatrix B = { 0 };

    BasicMatrix C = { 1 };

    BasicMatrix P = { 100 };

    BasicMatrix Q = { 0.00001 };

    BasicMatrix R = { 0.01 };

    int indice = 300;

    BasicVector XK1 = BasicVector(indice, 1); XK1[0] = -0.3773;
    BasicVector Time = BasicVector(indice,1); Time[0] = 0;

    for (arma::uword i = 1 ; i <= indice-1 ; i++ ) {
        XK1[i] = -0.3773;
        Time[i] = Time[i-1] + 1;
    }

    BasicVector RANDON = arma::randn<arma::vec>(XK1.n_rows);
    BasicVector ZK1 = XK1 + RANDON/10;

    BasicVector X0 = {0};
    BasicVector Z0 = {0};
    BasicVector U0 = {10, 10, 100};

    std::cout << "U0 Rows:" << U0.n_rows << "\t U0 Cols:" << U0.n_cols << std::endl;
    //KalmanFilter filter(A, B, C, Q, R, X0, Z0, P);
    KalmanFilter filter(X0, Z0, Q, R, {}, P);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK1.n_rows,1);      // Estado Corregido

    for (arma::uword i = 0 ; i <= ZK1.n_rows -1 ; i++) {

        BasicVector input { ZK1.at(i) };
        filter.update(input);
        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
    }

    std::cout << "FINISHING" << std::endl;

    //   Converting BasicVector to QVector

        QVector<double> VTime;
        QVector<double> V_Real, VZ_RUIDO, VXKC1;

        for (arma::uword i = 0 ; i <= ZK1.n_rows -1 ; i++) {

            VTime.append( Time(i) );

            V_Real.append( XK1(i) );

            VZ_RUIDO.append( ZK1(i) );

            VXKC1.append( XK1_corregida(i) );

        }

        std::cout << "FINISHING 2" << std::endl;

//  PLOTTING

            QCustomPlot plot;

            plot.legend->setVisible(true);
            QFont legendFont;
            legendFont.setPointSize(9);
            plot.legend->setFont(legendFont);
            plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
            plot.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);


        //    plot.addGraph(plot.yAxis, plot.xAxis);
            plot.addGraph();
            plot.graph(0)->setData(VTime, VZ_RUIDO);
            plot.graph(0)->setPen(QColor(Qt::blue));
            plot.graph(0)->setName("Medida Sensor");

            plot.addGraph();
            plot.graph(1)->setData(VTime, V_Real);
            plot.graph(1)->setPen(QColor(Qt::black));
            plot.graph(1)->setName("Medida Real");

            plot.addGraph();
            plot.graph(2)->setData(VTime, VXKC1);
            plot.graph(2)->setPen(QColor(Qt::red));
            plot.graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
            plot.graph(2)->setName("Estimacion Filtro Kalman");

            // give the axes some labels:
            plot.xAxis->setLabel("Time");
            plot.yAxis->setLabel("Value");
            // set axes ranges, so we see all data:
            plot.xAxis->setRange(0, 320);
            plot.yAxis->setRange(0, -1);
            plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
            plot.replot();
            plot.setMinimumSize(QSize(500, 500));
            plot.setWindowTitle("Resultados");
            plot.show();


            qDebug() << ECM_2(V_Real, VZ_RUIDO);
            qDebug() << ECM_2(V_Real, VXKC1);

//*/

// --------------------------------------------------------------------------



// ---------------------- EJEMPLO CON SIMULACION REAL -----------------------
/**/

    arma::arma_rng::set_seed_random();
    double m = 25;
    double k = 24;
    double b = 8;
    arma::uword N_Simulation = 200;


    BasicMatrix A = { {0, 1},
                      {-0.96,  -0.32} };

    BasicMatrix B = { 0, 0.04 };
    B = B.t();

    BasicMatrix C = { {1, 0} };

    BasicMatrix P = { {1, 0},
                      {0, 1} };

    BasicMatrix Q = { {0, 0},
                      {0, 0} };

//    BasicMatrix R = { {1, 0},
//                      {0, 1} };

    BasicMatrix R = { 1 };


    BasicVector U = BasicVector(N_Simulation+ 1, 1);

    for (arma::uword i = 1 ; i <= N_Simulation ; i++ ) {

        if( (i > N_Simulation*0.25) && (i < N_Simulation*0.75) )
            U[i] = 8;
        else {
            U[i] = -5;
        }

    }

    BasicVector X0 = {1, 1};
    BasicVector Z0 = { 0 };
    BasicVector U0 = { 0 };

    BasicVector X_NEXT;
    BasicVector X_CURRENT = X0;
    BasicVector U_CURRENT = { U[0] };
    BasicVector Z_CURRENT;

    // FOR REPRESENTATION OF THE SYSTEM
    BasicVector X1 = BasicVector(N_Simulation + 1, 1);
    BasicVector X2 = BasicVector(N_Simulation + 1, 1);
    BasicVector Z1 = BasicVector(N_Simulation+ 1, 1);



    // SYSTEM SIMULATION

    for (arma::uword idx = 0 ; idx <= N_Simulation ; idx++) {

        // X(k+1) = A*X(k) + B*U(k)
        X_NEXT = A*X_CURRENT + B*U_CURRENT;
        //X_NEXT = A*X_CURRENT + B*U0;

        // Y(k+1) = C*X(k+1)
        Z_CURRENT = C*X_NEXT;

        X_CURRENT = X_NEXT;
        U_CURRENT = U[idx];

//        std::cout << "ESTADO/n:" << X_CURRENT << std::endl;
//        std::cout << "MEDIDA:" <<Z_CURRENT << std::endl;
//        std::cout << "ACTUACION:" <<U_CURRENT << std::endl;

        // MAKING REGISTERS
        X1.at(idx) = X_CURRENT[0];
        X2.at(idx) = X_CURRENT[1];
        Z1.at(idx) = Z_CURRENT.at(0);
    }

    BasicVector RANDON = arma::randn<arma::vec>(Z1.n_rows);
    Z1 += 1*RANDON;

    // KALMAN FILTER

    X0 = {-2, 2};
    KalmanFilter filter(A, B, C, Q, R, X0, Z0, P, {}, U0);
    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(Z1.n_rows,1);
    BasicVector XK2_corregida = BasicVector(Z1.n_rows,1);


    for (arma::uword i = 0 ; i <= Z1.n_rows -1 ; i++) {
//
        //BasicVector input = ZK.row(i).t();
        BasicVector input { Z1.at(i) };
        std::cout << "MEDIDA:" << input << std::endl;
        filter.update(input, U_CURRENT);
        //filter.update(input, U0);

        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
        XK2_corregida[i] = output.at(1);

        U_CURRENT = U[i];

    }

    // REPRESENTATION

    BasicVector Time = BasicVector(N_Simulation+ 1, 1); Time[0] = 0;
    for (arma::uword i = 1 ; i <= N_Simulation ; i++ ) {
        Time[i] = Time[i-1] + 0.05;
    }


    QVector<double> VTime;
    QVector<double> VX1, VX2, VZ, VU;
    QVector<double> V_Kalman_X1, V_Kalman_X2;


    for (arma::uword i = 0 ; i <= Time.n_rows -1 ; i++) {

        VTime.append( Time(i) );
        VX1.append( X1(i) );
        VX2.append( X2(i) );
        VZ.append( Z1(i) );
        VU.append( U(i) );

        V_Kalman_X1.append( XK1_corregida(i) );
        V_Kalman_X2.append( XK2_corregida(i) );
    }

    QCustomPlot plot1;
    QFont legendFont;
    plot1.legend->setVisible(true);
    legendFont.setPointSize(9);
    plot1.legend->setFont(legendFont);
    plot1.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot1.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot1.addGraph();
    plot1.graph(0)->setData(VTime, VX1);
    plot1.graph(0)->setPen(QColor(Qt::blue));
    plot1.graph(0)->setName("Estado X1 Sistema");

    plot1.addGraph();
    plot1.graph(1)->setData(VTime, V_Kalman_X1);
    plot1.graph(1)->setPen(QColor(Qt::red));
    plot1.graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot1.graph(1)->setName("Estado X1 Kalman");

    plot1.addGraph();
    plot1.graph(2)->setData(VTime, VU);
    plot1.graph(2)->setPen(QColor(Qt::green));
    plot1.graph(2)->setName("Actuacion U");


    plot1.xAxis->setLabel("Time");
    plot1.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot1.xAxis->setRange(-1, 11);
    plot1.yAxis->setRange(-1, 15);
    plot1.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot1.replot();
    plot1.setMinimumSize(QSize(500, 500));
    plot1.setWindowTitle("ESTADO X1");
    plot1.show();

    std::cout << "X1_ECM:" << ECM_2(VX1, V_Kalman_X1) << std::endl;

    QCustomPlot plot2;
    plot2.legend->setVisible(true);
    legendFont.setPointSize(9);
    plot2.legend->setFont(legendFont);
    plot2.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot2.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot2.addGraph();
    plot2.graph(0)->setData(VTime, VX2);
    plot2.graph(0)->setPen(QColor(Qt::blue));
    plot2.graph(0)->setName("Estado X2");

    plot2.addGraph();
    plot2.graph(1)->setData(VTime, V_Kalman_X2);
    plot2.graph(1)->setPen(QColor(Qt::red));
    plot2.graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot2.graph(1)->setName("Estado X2 Kalman");

    plot2.addGraph();
    plot2.graph(2)->setData(VTime, VU);
    plot2.graph(2)->setPen(QColor(Qt::green));
    plot2.graph(2)->setName("Actuacion U");


    plot2.xAxis->setLabel("Time");
    plot2.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot2.xAxis->setRange(-1, 11);
    plot2.yAxis->setRange(-1, 5);
    plot2.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot2.replot();
    plot2.setMinimumSize(QSize(500, 500));
    plot2.setWindowTitle("MEDIDA X");
    plot2.show();


    std::cout << "X2_ECM:" << ECM_2(VX2, V_Kalman_X2) << std::endl;

//    std::cout << A << "\n" << B << "\n" << C << std::endl;
//    std::cout << "Maximo:" << arma::max(X2) << std::endl;
//    std::cout << "Size:" << Time.n_rows << std::endl;

//*/


// --------------------------------------------------------------------------
    return a.exec();
}
