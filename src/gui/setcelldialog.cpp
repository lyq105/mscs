#include "setcelldialog.h"
#include "ui_setcelldialog.h"
#include "sotsinterface.h"

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
    ui->dsb_volume_ratio->setValue(0.1);
    ui->dsb_g_meshsize->setValue(0.02);

    ui->rb_d_uniform->setChecked(1);
    on_rb_d_uniform_clicked();
    ui->dsb_amax->setValue(1.57);
    ui->dsb_amin->setValue(0.0);
    ui->dsb_bmax->setValue(3.14);
    ui->dsb_bmin->setValue(0.0);

    ui->rb_cd_uiniform->setChecked(1);
    on_rb_cd_uiniform_clicked();
    ui->dsb_ellipsoid_amax->setValue(0.08);
    ui->dsb_ellipsoid_amin->setValue(0.08);
    ui->dsb_ellipsoid_bmin->setValue(0.04);
    ui->dsb_ellipsoid_cmin->setValue(0.04);

}

void SetCellDialog::save_cell_info(std::string cell_file)
{
    std::ofstream ofile(cell_file.c_str());
    /// 单胞边界
    ofile << ui->lE_P1x->text().toStdString() << "  "
          << ui->lE_P2x->text().toStdString() << "  "
          << ui->lE_P1y->text().toStdString() << "  "
          << ui->lE_P2y->text().toStdString() << "  "
          << ui->lE_P1z->text().toStdString() << "  "
          << ui->lE_P2z->text().toStdString() << std::endl;

    unsigned int distribution_sign1,distribution_sign2;

    if(ui->rb_d_uniform->isChecked()) distribution_sign1 = 0;
    if(ui->rb_d_normal->isChecked()) distribution_sign1 = 1;
    ofile << distribution_sign1 << std::endl;

    unsigned int shape1=0,shape2=0,shape3=0;
    if(ui->cb_fibre->isChecked()) shape1 =1;
    if(ui->cb_ellipsoid->isChecked()) shape2 =1;
    if(ui->cb_money->isChecked()) shape3 =1;
    ofile << shape1 <<" " << shape2 <<" "<< shape3 <<std::endl;

    ofile << ui->dsb_ellipsoid_amin->text().toStdString() << " "
          << ui->dsb_ellipsoid_amax->text().toStdString() << " "
          << ui->dsb_ellipsoid_bmin->text().toStdString() << " "
          << ui->dsb_ellipsoid_cmin->text().toStdString() << std::endl;


    if(ui->rb_d_uniform->isChecked())
    {
        ofile << ui->dsb_amin->text().toStdString() << " "
              << ui->dsb_amax->text().toStdString() << " "
              << ui->dsb_bmin->text().toStdString() << " "
              << ui->dsb_bmax->text().toStdString() << std::endl;
    }

    ofile << shape1 <<" " << shape2 <<" "<< shape3 <<" ";


    if(ui->rb_d_normal->isChecked())
    {
        ofile << ui->dsb_sigma1->text().toStdString() << " "
              << ui->dsb_mu1->text().toStdString() << " "
              << ui->dsb_sigma2->text().toStdString() << " "
              << ui->dsb_mu2->text().toStdString() << std::endl;
    }
    ofile << distribution_sign2 << std::endl;
    if(ui->rb_cd_uiniform->isChecked())  // 中心坐标均匀分布
    {
        distribution_sign2 = 0;
        ofile << distribution_sign2 << std::endl;
    }
    if(ui->rb_cd_power->isChecked())  // 中心坐标指数分布
    {
        distribution_sign2 = 1;
        ofile << distribution_sign2 << std::endl;
        ofile << ui->dsb_cd_power_1->text().toStdString() << " "
              << ui->dsb_cd_power_1->text().toStdString() << " "
              << ui->dsb_cd_power_1->text().toStdString() << std::endl;
    }
    if(ui->rb_cd_normal->isChecked())  // 中心坐标正态分布
    {
        distribution_sign2 = 2;
        ofile << distribution_sign2 << std::endl;
        ofile << ui->dsb_cd_normal_sigma1->text().toStdString() << " "
              << ui->dsb_cd_normal_mu1->text().toStdString() << " "
              << ui->dsb_cd_normal_sigma2->text().toStdString() << " "
              << ui->dsb_cd_normal_mu2->text().toStdString() << " "
              << ui->dsb_cd_normal_sigma2->text().toStdString() << " "
              << ui->dsb_cd_normal_mu2->text().toStdString() << std::endl;
    }
    ofile << ui->dsb_volume_ratio->text().toStdString() <<std::endl;
    ofile.close();
}

void SetCellDialog::get_cell_info(PMCell_Info* cell_info)
{
    cell_info->amax[0] = ui->dsb_ellipsoid_amax->value();
    cell_info->amax[1] = ui->dsb_amin->value();
    cell_info->amax[2] = ui->dsb_ellipsoid_amax->value();
    cell_info->amax[0] = ui->dsb_ellipsoid_amax->value();
    cell_info->amax[0] = ui->dsb_ellipsoid_amax->value();
    cell_info->amax[0] = ui->dsb_ellipsoid_amax->value();
    //   cell_info->amax[0] =
    ui->lE_P1x->setText("0");
    ui->lE_P1y->setText("0");
    ui->lE_P1z->setText("0");
    ui->lE_P2x->setText("1");
    ui->lE_P2y->setText("1");
    ui->lE_P2z->setText("1");
 }

void SetCellDialog::on_rb_cd_uiniform_clicked()
{
    ui->dsb_cd_normal_mu1->setEnabled(0);
    ui->dsb_cd_normal_mu2->setEnabled(0);
    ui->dsb_cd_normal_mu3->setEnabled(0);
    ui->dsb_cd_normal_sigma1->setEnabled(0);
    ui->dsb_cd_normal_sigma2->setEnabled(0);
    ui->dsb_cd_normal_sigma3->setEnabled(0);

    ui->dsb_cd_power_1->setEnabled(0);
    ui->dsb_cd_power_2->setEnabled(0);
    ui->dsb_cd_power_3->setEnabled(0);

}

void SetCellDialog::on_rb_cd_power_clicked()
{
    ui->dsb_cd_power_1->setEnabled(1);
    ui->dsb_cd_power_2->setEnabled(1);
    ui->dsb_cd_power_3->setEnabled(1);

    ui->dsb_cd_normal_mu1->setEnabled(0);
    ui->dsb_cd_normal_mu2->setEnabled(0);
    ui->dsb_cd_normal_mu3->setEnabled(0);
    ui->dsb_cd_normal_sigma1->setEnabled(0);
    ui->dsb_cd_normal_sigma2->setEnabled(0);
    ui->dsb_cd_normal_sigma3->setEnabled(0);
}

void SetCellDialog::on_rb_cd_normal_clicked()
{
    ui->dsb_cd_normal_mu1->setEnabled(1);
    ui->dsb_cd_normal_mu2->setEnabled(1);
    ui->dsb_cd_normal_mu3->setEnabled(1);
    ui->dsb_cd_normal_sigma1->setEnabled(1);
    ui->dsb_cd_normal_sigma2->setEnabled(1);
    ui->dsb_cd_normal_sigma3->setEnabled(1);

    ui->dsb_cd_power_1->setEnabled(0);
    ui->dsb_cd_power_2->setEnabled(0);
    ui->dsb_cd_power_3->setEnabled(0);
}

void SetCellDialog::on_rb_d_uniform_clicked()
{
    ui->dsb_mu1->setEnabled(0);
    ui->dsb_mu2->setEnabled(0);
    ui->dsb_sigma1->setEnabled(0);
    ui->dsb_sigma2->setEnabled(0);
    ui->dsb_amax->setEnabled(1);
    ui->dsb_amin->setEnabled(1);
    ui->dsb_bmax->setEnabled(1);
    ui->dsb_bmin->setEnabled(1);
}

void SetCellDialog::on_rb_d_normal_clicked()
{
    ui->dsb_mu1->setEnabled(1);
    ui->dsb_mu2->setEnabled(1);
    ui->dsb_sigma1->setEnabled(1);
    ui->dsb_sigma2->setEnabled(1);
    ui->dsb_amax->setEnabled(0);
    ui->dsb_amin->setEnabled(0);
    ui->dsb_bmax->setEnabled(0);
    ui->dsb_bmin->setEnabled(0);
}
