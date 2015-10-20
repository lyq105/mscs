#include "cellsolveroption.h"
#include "ui_cellsolveroption.h"

CellSolverOption::CellSolverOption(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::CellSolverOption)
{
    ui->setupUi(this);
}

CellSolverOption::~CellSolverOption()
{
    delete ui;
}
