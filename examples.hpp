#ifndef EXAMPLES_HPP
#define EXAMPLES_HPP

#include <QFile>
#include <QTextStream>
#include <QDebug>
#include <QChar>
#include <QTime>
#include <QDate>


#include <string>
#include <iostream>
#include <conio.h>
#include "Load_Files/include_cpp/libxl.h"
#include "kalmanfilter.hpp"


struct date_time
{
    double consumo;
    qint64 TimeStamp;
};


void example1();
void example2();
void example3();
void example4();
void example5();
void example6();
void example7();
void example8();
void exampleJSon();
void exampleCar();









#endif // EXAMPLES_HPP
