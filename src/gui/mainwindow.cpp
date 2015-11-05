#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "aboutdialog.h"
#include <QFileDialog>
#include <QLabel>
#include "QVTKWidget.h"
#include <QFile>
#include <QFileInfo>
#include <meshpara.h>
#include <QTextEdit>
#include <QMessageBox>
#include "analysistype.h"
#include "materiel.h"
#include <sstream>
#include <QScrollBar>
#include <QtCore>
#include <QFuture>
#include "setcelldialog.h"
#include "vtkviewer.h"
#include "mylogger.h"
#include "qstream.h"
#include "cellsolveroption.h"


/// 构造和析构函数
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowIcon(QIcon(":/images/elisa_128px_540496_easyicon.net.png"));

    ui->tool_box->hide();
    //    ui->sidebar->setWidget();
    createMessageBox();
    createStatusBar();
    createVTKview();
    createSlots();
}

MainWindow::~MainWindow()
{
    delete message_content;
    delete ui;
}

/// 设置输出窗口
void MainWindow::createMessageBox()
{
    message_content = new QTextEdit;
    qcout = new QDebugStream(std::cout,message_content);
    message_content->setReadOnly(1);
    ui->message_box->hide();
    ui->message_box->setWidget(message_content);
    //    message_content->setText("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nkai shi le xie xie");
    //    message_content->append(QString::fromLocal8Bit("<font color=red>错误</font>"));
    QScrollBar *pScroll = message_content->verticalScrollBar();
    pScroll->setSliderPosition(pScroll->maximum());
}

/// 设置vtk显示窗口
void MainWindow::createVTKview()
{
    viewer = new VTKviewer;
    setCentralWidget(viewer->getVTKWidget());

}


/// 设置状态栏
void MainWindow::createStatusBar()
{
    QLabel* _pQLabel = new QLabel;
    _pQLabel->show();
    ui->statusBar->addPermanentWidget(_pQLabel);
    ui->statusBar->showMessage("Initializing...",1000);
}

/// 设置信号和槽
void MainWindow::createSlots()
{
    // about
    connect(ui->action_Mscs, SIGNAL(triggered()), this, SLOT(about()));

    // view menu
    connect(ui->action_show_cell, SIGNAL(triggered()), this, SLOT(show_pcmcell()));
    connect(ui->action_show_massage_box, SIGNAL(triggered(bool)), this, SLOT(show_message_box(bool)));
    connect(ui->action_show_axes, SIGNAL(triggered(bool)), this, SLOT(show_axes(bool)));
    ui->action_show_axes->setChecked(1);
    connect(ui->action_show_XY,SIGNAL(triggered()), this, SLOT(show_XY()));
    connect(ui->action_show_XZ,SIGNAL(triggered()), this, SLOT(show_XZ()));
    connect(ui->action_show_YZ,SIGNAL(triggered()), this, SLOT(show_YZ()));
    connect(ui->action_show_sidebar,SIGNAL(triggered(bool)), this, SLOT(show_sidebar(bool)));
    ui->action_show_sidebar->setChecked(1);
    connect(ui->action_show_reinforcement,SIGNAL(triggered()), this, SLOT(show_reinforcement()));
    connect(ui->action_show_matrix,SIGNAL(triggered()), this, SLOT(show_matrix()));
    ui->action_show_matrix->setEnabled(false);
    ui->action_show_reinforcement->setEnabled(false);

    // model menu
    connect(ui->action_anatype, SIGNAL(triggered()), this, SLOT(analysis_type()));

    connect(ui->action_import_mesh, SIGNAL(triggered()), this, SLOT(import_mesh()));
    connect(ui->action_import_geo, SIGNAL(triggered()), this, SLOT(import_geo()));

    connect(ui->action_import_cell_mesh, SIGNAL(triggered()), this, SLOT(import_cell_mesh()));
    connect(ui->action_import_cell_geo, SIGNAL(triggered()), this, SLOT(import_cell_geo()));

    connect(ui->action_new_mat, SIGNAL(triggered()), this, SLOT(new_material()));
    connect(ui->action_view_mat, SIGNAL(triggered()), this, SLOT(import_geo()));
    connect(ui->action_gen_cell, SIGNAL(triggered()), this, SLOT(set_cell()));

    // solve menu
    connect(ui->action_cell_solver_option, SIGNAL(triggered()), this, SLOT(set_cell_solver()));
    // file menu
    connect(ui->action_new_project,SIGNAL(triggered()),this,SLOT(new_project()));
    connect(ui->action_openfile,SIGNAL(triggered()),this,SLOT(import_inp()));

}

// slots ====================================================================================

void MainWindow::set_cell_solver()
{
    CellSolverOption csoption;
    if(csoption.exec() == QDialog::Accepted)
    {
        std::cout << "Cell solver has been set.\n";
        std::cout << "单胞求解器已设置为" <<std::endl;
        message_content->append(QString::fromLocal8Bit("<font color=red>单胞求解器已设置为</font>"));
        return;
    }
}

void MainWindow::analysis_type()
{
    AnalysisType anatype(this);
    if(anatype.exec() == QDialog::Accepted)
    {
        sots.set_analysis_type(anatype.get_analysis_type().toStdString());
        //       sots.set_prj_folder(anatype.get_prj_folder().toStdString());
        sots.set_prj_name(anatype.get_prj_name().toStdString());
    }
    //    QString logfile = anatype.get_prj_folder() + QString("/") + anatype.get_prj_name() +".log";
}


void MainWindow::show_matrix()
{
    viewer->show_matrix();
}
void MainWindow::show_reinforcement()
{
    viewer->show_reinforcement();
}

/// 设置显示侧边栏

void MainWindow::show_sidebar(bool checked)
{
    //if (checked) ui->sidebar->show();
    //else ui->sidebar->hide();
    //viewer->show_cell(std::string("/home/yzh/dataOutEllipDat.dat"));
}



/// 显示XY平面

void MainWindow::show_XY()
{
    viewer->show_XY();
}

/// 显示XZ平面
void MainWindow::show_XZ()
{
    viewer->show_XZ();
}

/// 显示YZ平面
void MainWindow::show_YZ()
{
    viewer->show_YZ();
}


/// 显示坐标轴
void MainWindow::show_axes(bool flags)
{
    if(flags){
        viewer->axes_on();
    }
    else{
        viewer->axes_off();
    }
}

/// 显示颗粒增强单胞
void MainWindow::show_pcmcell()
{
    if (sots.celldata == NULL)return;
    viewer->show_ug_scalar(sots.celldata,std::string("subdomainid"));
}

/// 设置单胞
void MainWindow::set_cell()
{
    SetCellDialog setcell(this);

    if(setcell.exec() == QDialog::Accepted)
    {
        QString path;
        QDir dir;
        path=dir.currentPath();
        QString cellfile(path+"/cellinfo.txt");
        setcell.save_cell_info(cellfile.toStdString().c_str());
        QFileInfo temp(cellfile);
        QString outfile= temp.absolutePath()+"/out.dat";
        QString datafile= temp.absolutePath()+"/data.dat";
        sots.set_pmcell(cellfile.toStdString(),outfile.toStdString(),datafile.toStdString());
        QFuture<int> sumf = QtConcurrent::run(&sots,&SOTSinterface::build_pmcell);
        sumf.waitForFinished();
        //sots.build_pmcell(cellfile.toStdString(),outfile.toStdString(),datafile.toStdString());
        viewer->show_cell(outfile.toStdString());
        //        QMessageBox::warning(0,"PATH",path,QMessageBox::Yes);//查看路径
        //        sots.mesh_pmcell(outfile.toStdString(),0.02);

        //        //sumf.waitForFinished();
        //        viewer->show_ug_mesh(sots.celldata);
        //        viewer->write_ug(sots.celldata,cellfile.toStdString()+".vtk");
    }
}


/// 导入单胞几何
void MainWindow::import_cell_geo()
{
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;MSCS Input Files(*.inp)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;MSCS Input Files(*.inp)");
    QString infile = QFileDialog::getOpenFileName(
                this, QString::fromLocal8Bit("直接导入计算文件"),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);
    if(infile.isNull()) return;
    QFileInfo temp(infile);
    QString outfile= temp.absolutePath()+"/out.dat";
    QString datafile= temp.absolutePath()+"/data.dat";
    //sots.solve(infile.toStdString(),outfile.toStdString(),datafile.toStdString());
    //sots.build_pmcell(infile.toStdString(),datafile.toStdString());
    sots.build_pmcell(infile.toStdString(),outfile.toStdString(),datafile.toStdString());
    viewer->show_cell(outfile.toStdString());
}

/// 直接导入单胞网格
void MainWindow::import_cell_mesh()
{
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;BDF Files(*.bdf)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;BDF Files(*.bdf)");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString::fromLocal8Bit("导入单胞网格文件"),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);
    //   cout << sots.celldata<< endl;
    if(!sFileName.isNull())
    {
        sots.celldata = viewer->load_bdf(sFileName.toStdString());
        //       cout << sots.celldata << endl;
        //show_pcmcell();
        viewer->show_ug_mesh(sots.celldata);
        ui->action_show_matrix->setEnabled(true);
        ui->action_show_reinforcement->setEnabled(true);
    }
}

/// 创建材料
void MainWindow::new_material()
{
    materiel mat(this);
    mat.exec();
}

extern int build_cell(MainWindow* qw,std::string infile,std::string datafile )
{
    qw->sots.build_pmcell(infile,datafile);
    return 0;
}
/// 直接导入计算文件
void MainWindow::import_inp()
{
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;MSCS Input Files(*.inp)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;MSCS Input Files(*.inp)");
    QString infile = QFileDialog::getOpenFileName(
                this, QString::fromLocal8Bit("直接导入计算文件"),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);
    if(infile.isNull()) return;
    QFileInfo temp(infile);
    QString outfile= temp.absolutePath()+"/out.dat";
    QString datafile= temp.absolutePath()+"/data.dat";
    //sots.solve(infile.toStdString(),outfile.toStdString(),datafile.toStdString());
    //sots.build_pmcell(infile.toStdString(),datafile.toStdString());
    //sots.build_pmcell(infile.toStdString(),outfile.toStdString(),datafile.toStdString());
    //viewer->show_cell(outfile.toStdString());

    sots.solve_cell(infile.toStdString(),datafile.toStdString());
    //   sots.mesh_pmcell(outfile.toStdString(),0.02);
    //QFuture<int> sumf =QtConcurrent::run(build_cell,this,infile.toStdString(),datafile.toStdString());
    //sumf.waitForFinished();
    //    viewer->show_ug_mesh(sots.celldata);
    viewer->write_ug(sots.celldata,infile.toStdString()+".vtk");
    //QMessageBox::warning(0,"PATH",path,QMessageBox::Yes);
    infile = "均匀化求解完毕！\n 单胞解写入文件 "+infile +".vtk";
    QMessageBox::about(NULL, "  ", QString::fromLocal8Bit(infile.toAscii()));
    viewer->show_ug_scalar(sots.celldata,std::string("subdomainid"));
    // ui->action_show_matrix->setEnabled(true);
    // ui->action_show_reinforcement->setEnabled(true);

    return;
}

/// 新建工程
void MainWindow::new_project()
{

    QFileDialog* openFilePath = new QFileDialog( this, QString::fromLocal8Bit("请选择一个文件夹"), "");
    openFilePath->setFileMode( QFileDialog::DirectoryOnly);
    if ( openFilePath->exec() == QDialog::Accepted )
    {
        sots.set_prj_folder(openFilePath->selectedFiles()[0].toStdString());
    }
    std::cout << "Choose " << openFilePath->selectedFiles()[0].toStdString()
              << "  as work directory" << endl;
    std::cout << "当前工作目录是 " << openFilePath->selectedFiles()[0].toStdString()
              << endl;
    delete openFilePath;
    return;
}

/// 关于对话框
void MainWindow::about()
{
    AboutDialog aboutdialog(this);
    aboutdialog.exec();
}

void MainWindow::show_message_box()
{
    ui->message_box->show();
}

void MainWindow::show_message_box(bool checked)
{
    if (checked)
    {
        ui->message_box->show();
    }
    else
        ui->message_box->hide();
}

/// 读入网格
void MainWindow::import_mesh()
{
    //QString homeName = getHomePath();
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;MSH Files(*.mesh);;BDF File(*.bdf)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;MSH Files(*.mesh);;BDF File(*.bdf)");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString(tr("Open Mesh")).toLocal8Bit(),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);

    if (!sFileName.isNull())
    {
        std::cout << sFileName.toStdString() << std::endl;
        viewer->show_ug_mesh(viewer->load_netgen_mesh(sFileName.toStdString()));
    }
}


void MainWindow::import_geo()
{
    //QString homeName = getHomePath();
    QString sDefaultName = tr("");
    QString selectedFilter = tr("STL Files(*.stl)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;STL Files(*.stl)");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString(tr("Open Geometry")).toLocal8Bit(),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);

    if (!sFileName.isNull())
    {
        std::cout << sFileName.toStdString() << std::endl;
        viewer->show_stl(sFileName.toStdString());
    }
}




