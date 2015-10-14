#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "sotsinterface.h"

namespace Ui {
class MainWindow;
}
class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkAppendPolyData;
class vtkOrientationMarkerWidget;

class QTextEdit;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();


private slots:
    void about();
    void show_message_box();
    void show_message_box(bool checked);
//    void show_cell();
//    void show_mesh();

//    void open_input();
//    void solve_cell_problem();
//    void solve_homo_problem();
//    void load_cell_mesh();
//    void add_material();
    void import_mesh();
    void import_geo();
    void new_project();

    /// 直接导入计算文件计算
    void import_inp();

    /// 新建材料
    void new_material();

    /// 导入单胞网格
    void import_cell_mesh();
    void import_cell_geo();
    void set_cell();

    /// view controls
    void show_pcmcell();
    void show_axes(bool flags);
    void show_XY();
    void show_XZ();
    void show_YZ();


private:
    void createSlots();
    void createStatusBar();
    void createToolbar();
    void createVTKview();
    void createMessageBox();


    void openSTL(QString stl_filename);
    void read_neutral_format(QString filename );
    void meshSTL(QString stl_filename);
    vtkUnstructuredGrid* import_bdf(QString bdf_filename);
    //    void openVtu(QString vtu_filename);
    void plot_cell();


    SOTSinterface sots;
    Ui::MainWindow *ui;
    QTextEdit* message_content;
    QVTKWidget* qvtkWidget;
    vtkRenderer* renderer;
    vtkOrientationMarkerWidget* widget;

    vtkActor* cell_mesh_actor;
    vtkActor* cell_geo_actor;
    vtkActor* cell_actor;
};

#endif // MAINWINDOW_H
