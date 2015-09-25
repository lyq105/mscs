#include "setcelldialog.h"
#include "ui_setcelldialog.h"

SetCellDialog::SetCellDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetCellDialog)
{
    ui->setupUi(this);
    ui->tabWidget->setCurrentWidget(ui->tab_globalInfo);
    set_default_value();
}

SetCellDialog::~SetCellDialog()
{
    delete ui;
}

/// 设置表单初始值
void SetCellDialog::set_default_value()
{
    ui->lE_P1x->setText("0");
    ui->lE_P1y->setText("0");
    ui->lE_P1z->setText("0");
    ui->lE_P2x->setText("1");
    ui->lE_P2y->setText("1");
    ui->lE_P2z->setText("1");
    ui->cb_coating->setChecked(0);
    ui->cb_ellipsoid->setChecked(1);
    ui->dsb_ellipsoid_ratio->setValue(1.0);
    ui->dsb_volume_ratio->setValue(0.4);
    ui->dsb_g_meshsize->setValue(0.02);

    ui->rb_d_uniform->setChecked(1);
    ui->dsb_amax->setValue(0.02);
    ui->dsb_amin->setValue(0.01);
    ui->dsb_bmax->setValue(0.02);
    ui->dsb_bmin->setValue(0.01);

    ui->rb_cd_uiniform->setChecked(1);

    ui->dsb_ellipsoid_amax->setValue(0.02);
    ui->dsb_ellipsoid_amin->setValue(0.01);
    ui->dsb_ellipsoid_bmin->setValue(0.005);
    ui->dsb_ellipsoid_cmin->setValue(0.002);
}

void SetCellDialog::get_cell_info()
{
    ui->lE_P1x->setText("0");
    ui->lE_P1y->setText("0");
    ui->lE_P1z->setText("0");
    ui->lE_P2x->setText("1");
    ui->lE_P2y->setText("1");
    ui->lE_P2z->setText("1");
//    ui->dSB_volRation->setValue(0.42);

//    ui->dSB_aMax_coin->setValue(0.1);
//    ui->dSB_aMin_coin->setValue(0.1);
//    ui->dSB_bMin_coin->setValue(0.1);
//    ui->dSB_cValue_coin->setValue(0.1);

//    ui->dSB_aMin_ellipsoid->setValue(0.05);
//    ui->dSB_aMax_ellipsoid->setValue(0.05);
//    ui->dSB_bMin_ellipsoid->setValue(0.05);

//    ui->dSB_aMin_fiber->setValue(0.05);
//    ui->dSB_aMax_fiber->setValue(0.05);
//    ui->dSB_bValue_fiber->setValue(0.05);

}
