#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}
class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkAppendPolyData;


class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void creatSlots();
    void creatMenu();
    void creatToolbar();
    void openSTL(QString stl_filename);
    void openVtu(QString vtu_filename);
    void read_neutral_format(QString filename );
    void meshSTL(QString stl_filename);
private slots:
    void on_action_Mscs_triggered();

    void on_action_3_triggered();

    void on_Openfile_triggered();

    void on_action_show_massage_box_triggered();

    void on_action_show_massage_box_triggered(bool checked);

private:
    Ui::MainWindow *ui;
    QVTKWidget* qvtkWidget;
    vtkRenderer* renderer;
};

#endif // MAINWINDOW_H
