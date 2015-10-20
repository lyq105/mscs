#ifndef CELLSOLVEROPTION_H
#define CELLSOLVEROPTION_H

#include <QDialog>

namespace Ui {
class CellSolverOption;
}

class CellSolverOption : public QDialog
{
    Q_OBJECT

public:
    explicit CellSolverOption(QWidget *parent = 0);
    ~CellSolverOption();
    Ui::CellSolverOption *ui;


private:
};

#endif // CELLSOLVEROPTION_H
