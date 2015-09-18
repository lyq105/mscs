#include "cellmodel.h"

CellModel::CellModel(QWidget *parent) :
    QMainWindow(parent)
{
    ui.setupUi(this);
}

void CellModel::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui.retranslateUi(this);
        break;
    default:
        break;
    }
}
