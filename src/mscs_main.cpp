#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
 //   QTextCodec::setCodecForTr(QTextCodec::codecForName("utf8"));
    MainWindow w;
    w.show();

    return a.exec();
}
