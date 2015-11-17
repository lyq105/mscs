#ifndef PROGRESS_DIALOG_H
#define PROGRESS_DIALOG_H

#include <QDialog>

namespace Ui {
class Progress_Dialog;
}
class QTimer;

class Progress_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Progress_Dialog(QString qst,QWidget *parent = 0);
    ~Progress_Dialog();
private slots:
    void update_value();

private:
    Ui::Progress_Dialog *ui;
    QTimer* timer;
    unsigned  var;
};

#endif // PROGRESS_DIALOG_H
