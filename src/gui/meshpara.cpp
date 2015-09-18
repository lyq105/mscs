#include "meshpara.h"
#include "ui_meshpara.h"

MeshPara::MeshPara(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::MeshPara)
{
    ui->setupUi(this);
}

MeshPara::~MeshPara()
{
    delete ui;
}
