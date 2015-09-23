#include "setcelldialog.h"
#include "ui_setcelldialog.h"

SetCellDialog::SetCellDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetCellDialog)
{
    ui->setupUi(this);
    set_default_value();
}

SetCellDialog::~SetCellDialog()
{
    delete ui;
}


void SetCellDialog::set_default_value()
{
    ui->lE_P1x->setText("0");
    ui->lE_P1y->setText("0");
    ui->lE_P1z->setText("0");
    ui->lE_P2x->setText("1");
    ui->lE_P2y->setText("1");
    ui->lE_P2z->setText("1");
    ui->dSB_volRation->setValue(0.42);

    ui->dSB_aMax_coin->setValue(0.1);
    ui->dSB_aMin_coin->setValue(0.1);
    ui->dSB_bMin_coin->setValue(0.1);
    ui->dSB_cValue_coin->setValue(0.1);

    ui->dSB_aMin_ellipsoid->setValue(0.05);
    ui->dSB_aMax_ellipsoid->setValue(0.05);
    ui->dSB_bMin_ellipsoid->setValue(0.05);

    ui->dSB_aMin_fiber->setValue(0.05);
    ui->dSB_aMax_fiber->setValue(0.05);
    ui->dSB_bValue_fiber->setValue(0.05);

}

void SetCellDialog::get_cell_info()
{
    ui->lE_P1x->setText("0");
    ui->lE_P1y->setText("0");
    ui->lE_P1z->setText("0");
    ui->lE_P2x->setText("1");
    ui->lE_P2y->setText("1");
    ui->lE_P2z->setText("1");
    ui->dSB_volRation->setValue(0.42);

    ui->dSB_aMax_coin->setValue(0.1);
    ui->dSB_aMin_coin->setValue(0.1);
    ui->dSB_bMin_coin->setValue(0.1);
    ui->dSB_cValue_coin->setValue(0.1);

    ui->dSB_aMin_ellipsoid->setValue(0.05);
    ui->dSB_aMax_ellipsoid->setValue(0.05);
    ui->dSB_bMin_ellipsoid->setValue(0.05);

    ui->dSB_aMin_fiber->setValue(0.05);
    ui->dSB_aMax_fiber->setValue(0.05);
    ui->dSB_bValue_fiber->setValue(0.05);

}
