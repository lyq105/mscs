#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "aboutdialog.h"
#include <QFileDialog>
#include <QLabel>
#include "QVTKWidget.h"
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkProperty.h>
#include <vtkPropPicker.h>
#include <vtkCallbackCommand.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCleanPolyData.h>
#include <vtkAxesActor.h>
#include <vtkSmartPointer.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkTransform.h>
#include <meshpara.h>
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkUnstructuredGridVolumeRayCastMapper.h>
#include <QTextEdit>
#include "analysistype.h"
#include "materiel.h"
#include <sstream>
#include <QScrollBar>
#include <vtkCellData.h>
#include "setcelldialog.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowIcon(QIcon(":/images/elisa_128px_540496_easyicon.net.png"));


    ui->menu_mesh->setEnabled(0);

    /// 设置输出窗口
    message_content = new QTextEdit;
    message_content->setReadOnly(1);
    ui->message_box->hide();
    ui->message_box->setWidget(message_content);
    QScrollBar *pScroll = message_content->verticalScrollBar();
    pScroll->setSliderPosition(pScroll->maximum());
    message_content->setText("kai shi le xie xie");
    message_content->append(QString::fromLocal8Bit("<font color=red>错误</font>"));

    createStatusBar();
    createSlots();
    createVTKview();
}

MainWindow::~MainWindow()
{
    //    delete renderer;
    delete message_content;
    delete qvtkWidget;
    delete ui;
}

void MainWindow::createVTKview()
{
    qvtkWidget = new QVTKWidget(this);
    setCentralWidget(qvtkWidget);
    renderer = vtkRenderer::New();
    // Setup VTK window
    vtkSmartPointer<vtkRenderWindow> renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
    renderer->ResetCamera();
    //qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();
    // qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
    qvtkWidget->SetRenderWindow(renderWindow);
    renderer->SetBackground(169/255., 169./255, 169./255);
    renderer->SetBackground(1, 1, 1);

    //renderer->SetBackground2(1,1,1);
    //renderer->SetGradientBackground(1);
    vtkSmartPointer<vtkAxesActor> axes =
            vtkSmartPointer<vtkAxesActor>::New();
    //renderer->AddActor(axes);
    //renderer->RemoveAllViewProps();

    vtkSmartPointer<vtkOrientationMarkerWidget> widget =
            vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget->SetOutlineColor( 0, 0, 0 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( qvtkWidget->GetInteractor() );
    widget->SetViewport( 0.0, 0.0, 1, 1 );
    widget->SetEnabled( 1 );
    widget->InteractiveOn();
    qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();

//    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
//    ren->SetViewport(0,0,0.2,0.2);
//    ren->AddActor(axes);
//    ren->SetBackground(1,1,1);
//    qvtkWidget->GetRenderWindow()->AddRenderer(ren);

    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
    qvtkWidget->GetRenderWindow()->Render();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
    //    QCoreApplication::processEvents();
}


////================================================================================================
void MainWindow::createStatusBar()
{
    QLabel* _pQLabel = new QLabel;
    _pQLabel->show();
    ui->statusBar->addPermanentWidget(_pQLabel);
    ui->statusBar->showMessage("Initializing...",1000);
}

void MainWindow::createSlots()
{
    // about
    connect(ui->action_Mscs, SIGNAL(triggered()), this, SLOT(about()));

    // view menu
    connect(ui->action_show_massage_box, SIGNAL(triggered()), this, SLOT(show_message_box()));
    connect(ui->action_show_massage_box, SIGNAL(triggered(bool)), this, SLOT(show_message_box(bool)));
    // model menu
    connect(ui->action_import_mesh, SIGNAL(triggered()), this, SLOT(import_mesh()));
    connect(ui->action_import_geo, SIGNAL(triggered()), this, SLOT(import_geo()));

    connect(ui->action_import_cell_mesh, SIGNAL(triggered()), this, SLOT(import_cell_mesh()));
    connect(ui->action_import_cell_geo, SIGNAL(triggered()), this, SLOT(import_cell_geo()));

    connect(ui->action_new_mat, SIGNAL(triggered()), this, SLOT(new_material()));
    connect(ui->action_view_mat, SIGNAL(triggered()), this, SLOT(import_geo()));
    connect(ui->action_gen_cell, SIGNAL(triggered()), this, SLOT(set_cell()));


    // file menu
    connect(ui->action_new_project,SIGNAL(triggered()),this,SLOT(new_project()));
    connect(ui->action_openfile,SIGNAL(triggered()),this,SLOT(import_inp()));

}


// slots
void MainWindow::set_cell()
{
    SetCellDialog setcell(this);
    setcell.exec();
}



void MainWindow::import_cell_geo()
{

}

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
    if(!sFileName.isNull())import_bdf(sFileName);
}
void MainWindow::new_material()
{
    materiel mat(this);
    mat.exec();
}

void MainWindow::import_inp()
{
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;MSH Files(*.inp)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr("");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString::fromLocal8Bit("直接导入计算文件"),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);
    return;
}

void MainWindow::new_project()
{
    AnalysisType anatype(this);
    if(anatype.exec() == QDialog::Accepted)
    {
        sots.set_analysis_type(anatype.get_analysis_type().toStdString());
        sots.set_prj_folder(anatype.get_prj_folder().toStdString());
        sots.set_prj_name(anatype.get_prj_name().toStdString());
    }
    //freopen()
}
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
        ui->message_box->show();
    else
        ui->message_box->hide();
}

//void MainWindow::on_action_3_triggered()
//{
//    MeshPara meshpara(this);
//    meshpara.exec();
//}

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
        read_neutral_format(sFileName);
    }
}
void MainWindow::import_geo()
{
    //QString homeName = getHomePath();
    QString sDefaultName = tr("");
    QString selectedFilter = tr(";;STL Files(*.stl)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr(";;STL Files(*.stl)");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString(tr("Open Geometry")).toLocal8Bit(),
                sDefaultName,filter,&selectedFilter,
                QFileDialog::DontUseNativeDialog);

    if (!sFileName.isNull())
    {
        std::cout << sFileName.toStdString() << std::endl;
        openSTL(sFileName);
    }
}
void MainWindow::openSTL(QString stl_filename)
{
    std::string inputFilename = stl_filename.toStdString();

    vtkSmartPointer<vtkSTLReader> reader =
            vtkSmartPointer<vtkSTLReader>::New();
    reader->SetFileName(inputFilename.c_str());
    reader->Update();

    // Visualize
    vtkSmartPointer<vtkPolyDataMapper> mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    //  actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    // actor->GetProperty()->EdgeVisibilityOn();
    renderer->RemoveAllViewProps();
    renderer->AddActor(actor);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
    //   return 1;

}

void MainWindow::import_bdf(QString bdf_filename)
{

    vtkIdType number_of_points = 0, number_of_tetra = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIntArray> intValue = vtkSmartPointer<vtkIntArray>::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("subdomainid");
    //     intValue->InsertNextValue(5);
    //    vtkSmartPointer<vtkDoubleArray> mat =
    //        vtkSmartPointer<vtkDoubleArray>::New();

    //      // Create the data to store (here we just use (0,0,0))
    //      double locationValue[3] = {0,0,0};

    //      location->SetNumberOfComponents(3);
    //      location->SetName("MyDoubleArray");
    //      location->InsertNextTuple(locationValue);
    // The data is added to FIELD data (rather than POINT data as usual)
    //      polydata->GetFieldData()->AddArray(location);

    std::ifstream infile(bdf_filename.toStdString().c_str());
    std::string s,l;
    while (getline(infile,s))
    {
        //跳过注释行
        if(s.substr(0,1)=="$")getline(infile,s);
        std::istringstream ss(s);
        //       std::cout << s << std::endl;
        ss >> l;
        if (l == "GRID")number_of_points ++;
        if (l == "CTETRA" )number_of_tetra++;
    }
    points->SetNumberOfPoints(number_of_points);
    std::cout << number_of_points << std::endl;
    std::cout << number_of_tetra << std::endl;
    infile.close();
    {
        std::ifstream infile(bdf_filename.toStdString().c_str());
        std::string s,l;
        while (getline(infile,s))
        {
            //跳过注释行
            if(s.substr(0,1)=="$")getline(infile,s);
            std::istringstream ss(s);
            //std::cout << s << std::endl;
            ss >> l;
            if (l == "GRID")
            {
                int n;
                double x, y, z;
                ss >> n >> x >> y >> z;
                points->SetPoint(n-1, x, y, z);
            }
            if (l == "CTETRA" )
            {
                vtkIdType a, b, c, d;
                vtkIdType n,subdomainId;
                ss >> n >> subdomainId >>  a >> b >> c >> d;
                if(subdomainId == 1) continue;
                //                std::cout << n << subdomainId << a << b << c <<d << std::endl;
                vtkSmartPointer<vtkTetra> tetra =
                        vtkSmartPointer<vtkTetra>::New();
                intValue->InsertNextValue(subdomainId);

                tetra->GetPointIds()->SetId(0, a-1);
                tetra->GetPointIds()->SetId(1, b-1);
                tetra->GetPointIds()->SetId(2, c-1);
                tetra->GetPointIds()->SetId(3, d-1);
                //                std::cout << n << std::endl;
                cellArray->InsertNextCell(tetra);
            }
        }

    }
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);
    unstructuredGrid->GetCellData()->AddArray(intValue);

    // Write file
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName((bdf_filename+".vtu").toStdString().c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();

    // Read and display file for verification that it was written correclty
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
            vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName((bdf_filename+".vtu").toStdString().c_str());
    reader->Update();

    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //       actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();
    renderer->RemoveAllViewProps();
    renderer->AddActor(actor);

    //      qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();

}

void MainWindow::meshSTL(QString stl_filename)
{
    std::string inputFilename = stl_filename.toStdString();

    //   GmshInitialize(NULL, 0);
    //    GmshSetOption("Mesh", "Algorithm3D", 0.1);
    //    GModel *m = new GModel();
    //    m->readSTL(inputFilename.c_str());
    //GmshMergeFile("../../tutorial/t5.geo"); // will also set the bbox
    //    m->mesh(3);
    //    for(GModel::riter it = m->firstRegion(); it != m->lastRegion(); ++it){
    //      GRegion *r = *it;
    //      printf("volume %d contains %d elements:\n", r->tag(), r->getNumMeshElements());
    //      for(unsigned int i = 0; i < r->getNumMeshElements(); i++)
    //        printf(" %d", r->getMeshElement(i)->getNum());
    //      printf("\n");
    //    }
    //    m->writeVTK((inputFilename+".vtu").c_str());
    // //   m->writeUNV("test.unv");
    //    delete m;
    //    GmshFinalize();

    vtkSmartPointer<vtkGenericDataObjectReader> reader =
            vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader->SetFileName((inputFilename+".vtu").c_str());
    reader->Update();

    // Visualize
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();
    renderer->RemoveAllViewProps();
    renderer->AddActor(actor);


    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
    //   return 1;

}

//void MainWindow::openVtu(QString vtu_filename)
//{
//read all the data from the file
//      vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
//        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
//      reader->SetFileName(filename.c_str());
//      reader->Update();

//      //Create a mapper and actor
//      vtkSmartPointer<vtkDataSetMapper> mapper =
//        vtkSmartPointer<vtkDataSetMapper>::New();
//      mapper->SetInputConnection(reader->GetOutputPort());

//      vtkSmartPointer<vtkActor> actor =
//        vtkSmartPointer<vtkActor>::New();
//      actor->SetMapper(mapper);

//      //Create a renderer, render window, and interactor
//      vtkSmartPointer<vtkRenderer> renderer =
//        vtkSmartPointer<vtkRenderer>::New();
//      vtkSmartPointer<vtkRenderWindow> renderWindow =
//        vtkSmartPointer<vtkRenderWindow>::New();
//      renderWindow->AddRenderer(renderer);
//      vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//        vtkSmartPointer<vtkRenderWindowInteractor>::New();
//      renderWindowInteractor->SetRenderWindow(renderWindow);

//      //Add the actor to the scene
//      renderer->AddActor(actor);
//      renderer->SetBackground(.3, .6, .3); // Background color green

//      //Render and interact
//      renderWindow->Render();
//      renderWindowInteractor->Start();
//}

void MainWindow::read_neutral_format(QString filename )
{
    std::ifstream infile(filename.toStdString().c_str());
    vtkIdType number_of_points, number_of_tetra;
    infile >> number_of_points;

    vtkSmartPointer<vtkPoints> points
            = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(number_of_points);
    for (vtkIdType i = 0; i < number_of_points; i++)
    {
        double x, y, z;
        infile >> x >> y >> z;
        points->SetPoint(i, x, y, z);
    }

    infile >> number_of_tetra;
    vtkSmartPointer<vtkCellArray> cellArray
            = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType i = 0; i < number_of_tetra; i++)
    {
        vtkIdType a, b, c, d;
        vtkIdType subdomainId;
        infile >> subdomainId >> a >> b >> c >> d;
        //       std::cout << a << b << c << d << std::endl;
        vtkSmartPointer<vtkTetra> tetra =
                vtkSmartPointer<vtkTetra>::New();

        tetra->GetPointIds()->SetId(0, a - 1);
        tetra->GetPointIds()->SetId(1, b - 1);
        tetra->GetPointIds()->SetId(2, c - 1);
        tetra->GetPointIds()->SetId(3, d - 1);
        cellArray->InsertNextCell(tetra);
    }
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);

    // Write file
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName((filename+".vtu").toStdString().c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();

    // Read and display file for verification that it was written correclty
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
            vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName((filename+".vtu").toStdString().c_str());
    reader->Update();

    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputConnection(reader->GetOutputPort());

    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //       actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();
    renderer->RemoveAllViewProps();
    renderer->AddActor(actor);

    //      qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();
    //    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
}


