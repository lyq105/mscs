#include "setcelldialog.h"
#include "ui_setcelldialog.h"

SetCellDialog::SetCellDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetCellDialog)
{
    ui->setupUi(this);
}

SetCellDialog::~SetCellDialog()
{
    delete ui;
}
