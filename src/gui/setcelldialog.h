#ifndef SETCELLDIALOG_H
#define SETCELLDIALOG_H

#include <QDialog>

namespace Ui {
class SetCellDialog;
}
class PMCell_Info;
class SetCellDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SetCellDialog(QWidget *parent = 0);
    ~SetCellDialog();
    /// 读取单胞几何参数
    void get_cell_info(PMCell_Info*);

private slots:
    void on_rb_cd_uiniform_clicked();
    void on_rb_cd_power_clicked();
    void on_rb_cd_normal_clicked();
    void on_rb_d_uniform_clicked();
    void on_rb_d_normal_clicked();

private:
    void set_default_value();

    Ui::SetCellDialog *ui;
};

#endif // SETCELLDIALOG_H
