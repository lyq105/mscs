#include "mainwindow.h"
#include <QApplication>
#include <stdio.h>
#include <iostream>
int main(int argc, char *argv[])
{
//    std::streambuf *psbuf, *backup;
//    std::ofstream filestr;
//    filestr.open ("test.txt");

//    backup = std::cout.rdbuf();     // back up cout's streambuf

//    psbuf = filestr.rdbuf();        // get file's streambuf
//    std::cout.rdbuf(psbuf);         // assign streambuf to cout
    QApplication a(argc, argv);
 //   QTextCodec::setCodecForTr(QTextCodec::codecForName("utf8"));
    MainWindow w;
    w.show();
    return a.exec();
}
