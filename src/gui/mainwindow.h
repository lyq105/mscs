#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTextStream>
#include "sotsinterface.h"

namespace Ui {
class MainWindow;
}

class VTKviewer;
class QTextEdit;
class QDebugStream;



class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    SOTSinterface sots;


private slots:
    void about();
    void show_message_box();
    void show_message_box(bool checked);
    void show_sidebar(bool checked);
//    void show_cell();
//    void show_mesh();

//    void open_input();
//    void solve_cell_problem();
//    void solve_homo_problem();
//    void load_cell_mesh();
//    void add_material();
    void set_cell_solver();
    void analysis_type();
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

    void show_matrix();
    void show_reinforcement();

private:
    void createSlots();
    void createStatusBar();
    void createToolbar();
    void createVTKview();
    void createMessageBox();

    Ui::MainWindow *ui;
    QTextEdit* message_content;
    VTKviewer* viewer;
    QDebugStream* qcout;

 //   QTextStream cout;
};

#endif // MAINWINDOW_H
