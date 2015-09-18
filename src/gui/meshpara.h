#ifndef MESHPARA_H
#define MESHPARA_H

#include <QDialog>

namespace Ui {
class MeshPara;
}

class MeshPara : public QDialog
{
    Q_OBJECT

public:
    explicit MeshPara(QWidget *parent = 0);
    ~MeshPara();

private:
    Ui::MeshPara *ui;
};

#endif // MESHPARA_H
