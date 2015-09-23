#ifndef MATERIEL_H
#define MATERIEL_H

#include "ui_materiel.h"

class materiel : public QDialog
{
    Q_OBJECT

public:
    explicit materiel(QWidget *parent = 0);

protected:
    void changeEvent(QEvent *e);

private slots:
    void onChanged(int);

private:
    Ui::materiel ui;
};

#endif // MATERIEL_H
