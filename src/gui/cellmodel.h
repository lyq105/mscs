#ifndef CELLMODEL_H
#define CELLMODEL_H

#include "ui_cellmodel.h"

class CellModel : public QMainWindow
{
    Q_OBJECT

public:
    explicit CellModel(QWidget *parent = 0);

protected:
    void changeEvent(QEvent *e);

private:
    Ui::CellModel ui;
};

#endif // CELLMODEL_H
