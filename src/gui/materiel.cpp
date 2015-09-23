#include "materiel.h"

materiel::materiel(QWidget *parent) :
    QDialog(parent)
{
    ui.setupUi(this);
 //   ui.comboBox->
    ui.tab->setDisabled(1);
    ui.tab_2->setDisabled(0);
    ui.tab_3->setDisabled(1);
    ui.tabWidget->setCurrentIndex(0);
    connect(ui.comboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(onChanged(int)));
 //   ui.tableWidget->item(0,1)->text();

}

void materiel::changeEvent(QEvent *e)
{
    QDialog::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui.retranslateUi(this);
        break;
    default:
        break;
    }
}

void materiel::onChanged(int index)
{
    if (index == 0)
    {
    ui.tab->setDisabled(1);
    ui.tab_2->setDisabled(0);
    ui.tab_3->setDisabled(1);
    ui.tabWidget->setCurrentIndex(0);
    }
    if (index == 1)
    {
    ui.tab->setDisabled(0);
    ui.tab_2->setDisabled(1);
    ui.tab_3->setDisabled(1);
    ui.tabWidget->setCurrentIndex(1);
    }
    if (index == 2)
    {
    ui.tab->setDisabled(1);
    ui.tab_2->setDisabled(1);
    ui.tab_3->setDisabled(0);
    ui.tabWidget->setCurrentIndex(2);
    }
}
