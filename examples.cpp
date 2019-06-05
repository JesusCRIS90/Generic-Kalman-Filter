#include "examples.hpp"
#include <armadillo/include/armadillo>
#include <iostream>
#include <qcustomplot.h>

QString QString2String(std::string str){
   return QString::fromStdString(str);
}


arma::mat EmptyInizialization(arma::mat received = {})
{
    return received;
}

double ECM_(QVector<double>& input1, QVector<double>& input2){
    BasicVector in1(input1.toStdVector());
    BasicVector in2(input2.toStdVector());
    BasicVector result( in1 - in2 );
    BasicVector r = result%result;          // Schur Product
//    std::cout << r << std::endl;
    return ( qSqrt(arma::mean(r)) );
}

void example1()
{

    // Initialize the random generator
            arma::arma_rng::set_seed_random();

        // Create a 4x4 random matrix and print it on the screen
            arma::Mat<double> A = arma::randu(4, 4);
            std::cout << "A:\n" << A << "\n";

        // Multiply A with his transpose:
            std::cout << "A * A.t() =\n";
            std::cout << A * A.t() << "\n";

        // Access/Modify rows and columns from the array:
             A.row(0) = A.row(1) + A.row(3);
             A.col(3).zeros();
             std::cout << "add rows 1 and 3, store result in row 0, also fill 4th column with zeros:\n";
             std::cout << "A:\n" << A << "\n";

        // Create a new diagonal matrix using the main diagonal of A:
             arma::Mat<double>B = arma::diagmat(A);
             std::cout << "B:\n" << B << "\n";

        // Save matrices A and B:
             A.save("A_mat.txt", arma::arma_ascii);
             B.save("B_mat.txt", arma::arma_ascii);

}

void example2(){

    std::string Data = "1 2 3; 4 5 6; 7 8 9";

    arma::Mat<double> C(Data);

//    std::cout << "A:\n" << C << "\n";
//    std::cout << arma::det(C) << std::endl;


    //arma::mat A = arma::randu<mat>(5,5);
    arma::mat B = C.t()*C;  // generate a symmetric matrix

    std::cout << B << std::endl;
    std::cout << arma::inv(C) << std::endl;

    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, B);

    std::cout << eigval << std::endl;
    std::cout << eigvec << std::endl;
}

void example3(){
    libxl::Book *book = xlCreateBook();

    if(book){
        if(book->load(L"C:\\Users\\JESUS\\Desktop\\BORRAME\\Project_Qt_Armadillo\\QT_Armadillo\\Files\\DemandaReal.xls"))
        {
          libxl::Sheet* sheet = book->getSheet(0);
          if(sheet){

              int last = sheet->lastRow() - 1;          // Last = sheet->lastRow() - 1;
              const wchar_t* s = sheet->readStr(last,5);
              /* Convert from wchar to wstring */
              std::wstring ws(s);
              std::string test(ws.begin(), ws.end());

              if(s) {std::wcout << s << std::endl;}

              double val = sheet->readNum(last,4);

              std::cout << val << std::endl;
              std::cout << test << std::endl;
              qDebug() << "All Rigth";

              QString Date_Time = QString2String(test).replace(QLatin1String("T"),QLatin1String(" ")).replace(QLatin1String("+"),QLatin1String(" "));
              //Date_Time = Date_Time.replace(QLatin1String("T"),QLatin1String(" ")).replace(QLatin1String("+"),QLatin1String(" "));
              //Date_Time = Date_Time.replace(QLatin1String("+"),QLatin1String(" "));

              qDebug() << Date_Time;

//              QStringList list = Date_Time.split(' ');

//              qDebug() << list;

              //QStringList Date = list[0].split('-');
              //QStringList Time = list[1].split(':');

              QStringList Date = Date_Time.split(' ')[0].split('-');
              QStringList Time = Date_Time.split(' ')[1].split(':');


              qDebug() << Date << "   " << Time;
//              qDebug() << Date[0].toInt() << "  " << Date[1].toInt() << "  " << Date[2].toInt();
//              qDebug() << Time[0].toInt() << "  " << Time[1].toInt() << "  " << Time[2].toInt();
              QDateTime startDate(QDate(Date[0].toInt(), Date[1].toInt(), Date[2].toInt()),
                                  QTime( Time[0].toInt(), Time[1].toInt(), Time[2].toInt() ) );

              qDebug() << startDate << "   " << startDate.toSecsSinceEpoch();


               QList<date_time> FinalList;
               date_time temp;
               temp.consumo = val;
               temp.TimeStamp = startDate.toSecsSinceEpoch();

               FinalList.append(temp);


          }
        }
        else{
            qDebug() << "Error Reading the File";
        }
    }

    else {
        qFatal("Error creating a Book");
    }
}

void example4()
{
    //    arma::arma_rng::set_seed_random();

        arma::mat p = EmptyInizialization( arma::randu(4, 4) );
        arma::mat* ptr = new arma::mat(arma::randu(4, 4));
        arma::mat p2 = EmptyInizialization();

    //    p.is_empty() ? std::cout << "\nYES" : std::cout << "\nNO";
    //    p2.is_empty() ? std::cout << "\nYES" : std::cout << "\nNO";

    //    KalmanParameters param;
    //    BasicMatrix A = arma::randu(4,4);
    //    param.set_A(A).set_B(A).set_C().set_D().set_Q().set_R().set_Pk();

        arma::mat copy(*ptr);

        std::cout << *ptr << "\n\n" << copy << std::endl;

        std::cout << p2.is_empty() << "\t" << p2.n_cols << "\t" << p2.n_rows << std::endl;
}

void example5(){

//    BasicMatrix Q = { {1, 0, 0},
//                      {0, 1, 0},
//                      {0, 0, 1} };

////    BasicMatrix Q = {1};

////    BasicMatrix R = {1};

//    BasicMatrix R = { {1, 0, 0},
//                      {0, 1, 0},
//                      {0, 0, 1} };

//    BasicVector Xk = {1, 2, 3};
//    BasicVector Uk = {1, 0, 1};
//    //std::cout << Q << "\n\n" << R << "\n\n" << Xk << std::endl;



//    KalmanParameters params(Xk, Q, R, Uk);
////    std::cout << params << std::endl;

//    KalmanFilter filter(Xk, Q, R);

//    std::cout << filter << std::endl;

//    KalmanParameters newParams = filter.Configuration();

//    newParams.copyParameters(params);

//    filter.SetConfiguration(newParams);

//    std::cout << "\n\n\n" << filter << std::endl;

//    BasicVector input = {2, 3, 4};
//    BasicVector output = filter.update(input);

//    std::cout << input << std::endl;
//    std::cout << output << std::endl;

////    filter.SetConfiguration(params);

////    std::cout << "\n\n\n" << filter << std::endl;

////    // Trying Copy Constructor
////    KalmanParameters copy(params);
////    std::cout << "\n\n\n ------------------------------------- \n " << copy << std::endl;

////    std::unique_ptr<KalmanParameters> pt = std::make_unique<KalmanParameters>(Xk, Q, R);
////    KalmanParameters* pt = new KalmanParameters(Xk, Q, R, Uk);
////    std::cout << *pt << std::endl;

////    KalmanParameters copy(pt);
////    std::cout << "\n\n\n ------------------------------------- \n " << copy << std::endl;
}

// Kalman Filter. Train Example (SIMPLE)
void example6(){

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
        filter.update(input);
        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
        XK2_corregida[i] = output.at(1);
    }


/*   Converting BasicVector to QVector   */

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


/*  PLOTTING  */

    QCustomPlot plot;

    plot.legend->setVisible(true);
    QFont legendFont;
    legendFont.setPointSize(9);
    plot.legend->setFont(legendFont);
    plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);


    plot.addGraph(/*plot.yAxis, plot.xAxis*/);
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


    //std::getchar();

}

// Kalman Filter. Train Example (Complex)
void example7()
{
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
    BasicVector Z0 = {5};
    BasicVector U0 = {0, 0};

    KalmanFilter filter(A, B, C, Q, R, X0, Z0, P, {}, U0);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK.n_rows,1);
    BasicVector XK2_corregida = BasicVector(ZK.n_rows,1);


    for (arma::uword i = 0 ; i <= ZK.n_rows -1 ; i++) {

//        BasicVector input = ZK.row(i).t();
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
}

// Kalman Filter. Sensor Example (Simple)
void example8()
{
    arma::arma_rng::set_seed_random();

    double Tm = 0.1;
    BasicMatrix A = { 1 };

    BasicMatrix B = { 1 };

    BasicMatrix C = { 1 };

    BasicMatrix P = { 1 };

    BasicMatrix Q = { 0.001 };

    BasicMatrix R = { 1 };

    BasicVector XK1 = BasicVector(101, 1); XK1[0] = -0.3773;
    BasicVector Time = BasicVector(101,1); Time[0] = 0;

    for (arma::uword i = 1 ; i <= 100 ; i++ ) {
        XK1[i] = -0.3773;
        Time[i] = Time[i-1] + 1;
    }

    BasicVector RANDON = arma::randn<arma::vec>(XK1.n_rows);
    BasicVector ZK1 = XK1 + RANDON;

    BasicVector X0 = {0};
    BasicVector Z0 = {5};
    BasicVector U0 = {0, 0};

    KalmanFilter filter(A, B, C, Q, R, X0, Z0, P);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(ZK1.n_rows,1);      // Estado Corregido

    for (arma::uword i = 0 ; i <= ZK1.n_rows -1 ; i++) {

        BasicVector input { ZK1.at(i) };
        filter.update(input);
        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
    }

    //   Converting BasicVector to QVector

        QVector<double> VTime;
        QVector<double> V_Real, VZ_RUIDO, VXKC1;

        for (arma::uword i = 0 ; i <= ZK1.n_rows -1 ; i++) {

            VTime.append( Time(i) );

            V_Real.append( XK1(i, 0) );

            VZ_RUIDO.append( ZK1(i, 1) );

            VXKC1.append( XK1_corregida(i) );

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
            plot.xAxis->setRange(-1, 12);
            plot.yAxis->setRange(-5, 120);
            plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
            plot.replot();
            plot.setMinimumSize(QSize(500, 500));
            plot.setWindowTitle("Distancia recorrida del Tren");
            plot.show();


            qDebug() << ECM_(V_Real, VZ_RUIDO);
            qDebug() << ECM_(V_Real, VXKC1);
}


void exampleJSon()
{
    // PROBANDO JSON

        QFile file("C:/Users/JESUS/Desktop/BORRAME/Project_Qt_Armadillo/QT_Armadillo/JsonPrueba.json");
        file.open(QIODevice::ReadOnly | QIODevice::Text);
        QByteArray jsonData = file.readAll();
        file.close();

        QJsonDocument document = QJsonDocument::fromJson(jsonData);
        QJsonObject object = document.object();

        qDebug() << object;
        qDebug() << object["A"];
        qDebug() << object["A"].toObject().take("NoRows").toInt();
        qDebug() << object["A"].toObject().take("Row1").toArray();

        QJsonArray list = object["A"].toObject().take("Row1").toArray();

        foreach (const QJsonValue & v, list){
                    qDebug() << v.toDouble();
            }
}


void exampleCar()
{
    arma::arma_rng::set_seed_random();

    BasicVector Time = BasicVector(101,1); Time[0] = 0;
    BasicVector Uk = BasicVector(101,1); Uk[0] = 0;

    double dt = 0.1;
    double SensorStd = 0.03;
    double SpeedStd = 0.2;

//    BasicVector x = ( 0.05*Time )%( arma::cos(2*Time) );
    BasicVector x = BasicVector(101, 1); x[0] = 0.05*Time[0]*std::cos(2*Time[0]);


    BasicVector RANDON = arma::randn<arma::vec>(Time.n_rows);

    for (arma::uword i = 1 ; i <= 100 ; i++ ) {
        Time[i] = Time[i-1] + dt;
        x[i] = 0.05*Time[i]*std::cos(2*Time[i]);
    }

    BasicVector u = arma::diff(x) + SpeedStd*arma::randn<arma::vec>(Time.n_rows - 1);

    //    std::cout << u.n_rows;
    for (arma::uword i = 0 ; i <= u.n_rows ; i++ ) {
            Uk[i+1] = u[i];
    }

    Uk += RANDON/4;

    BasicVector Zk = x - SensorStd*RANDON;


// ------------------- Kalman Filter

    BasicMatrix A = { 0.1 };

    BasicMatrix C = { 0.1 };

    BasicMatrix B = { dt };

    BasicMatrix P = { 0 };

    BasicMatrix Q = { SpeedStd*SpeedStd };

    BasicMatrix R = { SensorStd*SensorStd };

    BasicVector X0 = {0};
    BasicVector Z0 = {0};
    BasicVector U0 = { Uk[0] };

    KalmanFilter filter(A, B, C, Q, R, X0, Z0, P, {}, U0);

    std::cout << filter << std::endl;

    BasicVector XK1_corregida = BasicVector(Time.n_rows,1);


    for (arma::uword i = 0 ; i <= Time.n_rows -1 ; i++) {

        BasicVector input { Zk.at(i) };
        BasicVector input_uk { Uk.at(i) };
        filter.update(input, input_uk );
        BasicVector output = filter.GetState();
        XK1_corregida[i] = output.at(0);
//        XK2_corregida[i] = output.at(1);

    }




// Types Conversion


    QVector<double> VTime;
    QVector<double> VX, VUk, VZK, VX_Corregido;


    for (arma::uword i = 0 ; i <= Time.n_rows -1 ; i++) {

        VTime.append( Time(i) );

        VX.append( x(i) );
        VUk.append( Uk(i) );
        VZK.append( Zk(i) );
        VX_Corregido.append( XK1_corregida(i) );

    }

//    std::cout << x << std::endl;
//    std::cout << Time << std::endl;


// ----------- PLOTTING ---------------

    QCustomPlot plot;

    plot.legend->setVisible(true);
    QFont legendFont;
    legendFont.setPointSize(9);
    plot.legend->setFont(legendFont);
    plot.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot.addGraph();
    plot.graph(0)->setData(VTime, VX);
    plot.graph(0)->setPen(QColor(Qt::blue));
    plot.graph(0)->setName("Vector XK Estimado");

//    plot.addGraph();
//    plot.graph(1)->setData(VTime, VXK1);
//    plot.graph(1)->setPen(QColor(Qt::black));
//    plot.graph(1)->setName("Medida Real");

//    plot.addGraph();
//    plot.graph(2)->setData(VTime, VXKC1);
//    plot.graph(2)->setPen(QColor(Qt::red));
//    plot.graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
//    plot.graph(2)->setName("Medida Filtro Kalman");

    // give the axes some labels:
    plot.xAxis->setLabel("Time");
    plot.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot.xAxis->setRange(-1, 12);
    plot.yAxis->setRange(-1, 1);
    plot.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot.replot();
    plot.setMinimumSize(QSize(500, 500));
    plot.setWindowTitle("Vector XK");
    plot.show();



    QCustomPlot plot2;

    plot2.legend->setVisible(true);
    legendFont.setPointSize(9);
    plot2.legend->setFont(legendFont);
    plot2.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot2.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot2.addGraph();
    plot2.graph(0)->setData(VTime, VUk);
    plot2.graph(0)->setPen(QColor(Qt::blue));
    plot2.graph(0)->setName("Vector Uk");

    plot2.xAxis->setLabel("Time");
    plot2.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot2.xAxis->setRange(-1, 12);
    plot2.yAxis->setRange(-1, 1);
    plot2.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot2.replot();
    plot2.setMinimumSize(QSize(500, 500));
    plot2.setWindowTitle("Vector Uk");
    plot2.show();


    QCustomPlot plot3;

    plot3.legend->setVisible(true);
    legendFont.setPointSize(9);
    plot3.legend->setFont(legendFont);
    plot3.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot3.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot3.addGraph();
    plot3.graph(0)->setData(VTime, VZK);
    plot3.graph(0)->setPen(QColor(Qt::blue));
    plot3.graph(0)->setName("Vector Uk");

    plot3.xAxis->setLabel("Time");
    plot3.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot3.xAxis->setRange(-1, 12);
    plot3.yAxis->setRange(-1, 1);
    plot3.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot3.replot();
    plot3.setMinimumSize(QSize(500, 500));
    plot3.setWindowTitle("Vector ZK");
    plot3.show();



    QCustomPlot plot4;

    plot4.legend->setVisible(true);
    legendFont.setPointSize(9);
    plot4.legend->setFont(legendFont);
    plot4.legend->setBrush(QBrush(QColor(255,255,255,230)));
    plot4.axisRect()->insetLayout()->setInsetAlignment(0, Qt::AlignBottom|Qt::AlignRight);

    plot4.addGraph();
    plot4.graph(0)->setData(VTime, VX);
    plot4.graph(0)->setPen(QColor(Qt::blue));
    plot4.graph(0)->setName("Vector VX");

    plot4.addGraph();
    plot4.graph(1)->setData(VTime, VX_Corregido);
    plot4.graph(1)->setPen(QColor(Qt::red));
    plot4.graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssStar, Qt::red, Qt::white, 7));
    plot4.graph(1)->setName("Vector VX_Corregida");

    plot4.xAxis->setLabel("Time");
    plot4.yAxis->setLabel("Value");
    // set axes ranges, so we see all data:
    plot4.xAxis->setRange(-1, 12);
    plot4.yAxis->setRange(-1, 1);
    plot4.setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    plot4.replot();
    plot4.setMinimumSize(QSize(500, 500));
    plot4.setWindowTitle("Vector ZK");
    plot4.show();

}
