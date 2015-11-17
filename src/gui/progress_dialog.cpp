#include "progress_dialog.h"
#include "ui_progress_dialog.h"
#include <QTimer>

Progress_Dialog::Progress_Dialog(QString qst, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Progress_Dialog)
{
    ui->setupUi(this);
    ui->progressBar->setMaximum(100);
    ui->progressBar->setTextVisible(false);
    ui->label->setText(qst);
    ui->label->setAlignment(Qt::AlignCenter);
    timer = new QTimer(this);
    connect(timer,SIGNAL(timeout()),this,SLOT(update_value()));
    timer->start(10);
    var = 0;
}

Progress_Dialog::~Progress_Dialog()
{
    delete ui;
}

void Progress_Dialog::update_value()
{
    var += 1;
    ui->progressBar->setValue(var);
    var = var%100;
}

