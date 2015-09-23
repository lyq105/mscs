#ifndef SETCELLDIALOG_H
#define SETCELLDIALOG_H

#include <QDialog>

namespace Ui {
class SetCellDialog;
}

class SetCellDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SetCellDialog(QWidget *parent = 0);
    ~SetCellDialog();

private:
    void set_default_value();
    void get_cell_info();
    Ui::SetCellDialog *ui;
};

#endif // SETCELLDIALOG_H
