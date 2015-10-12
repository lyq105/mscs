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
#include <vtkBorderWidget.h>

#include <vtkSmartPointer.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkTextMapper.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkParametricFunctionSource.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkParametricEllipsoid.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyData.h>
#include <vtkStripper.h>
#include <vtkFeatureEdges.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkSphereSource.h>
#include <vtkBox.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkCubeSource.h>
#include <vtkLookupTable.h>
#include <vtkAssignAttribute.h>


/// 构造和析构函数
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowIcon(QIcon(":/images/elisa_128px_540496_easyicon.net.png"));
    ui->menu_mesh->setEnabled(0);
    ui->menu_mesh->hide();
    createMessageBox();
    createStatusBar();
    createVTKview();
    createSlots();
}

MainWindow::~MainWindow()
{
    //    delete renderer;
    delete message_content;
    delete qvtkWidget;
    delete ui;
}

/// 设置输出窗口
void MainWindow::createMessageBox()
{
    message_content = new QTextEdit;
    message_content->setReadOnly(1);
    ui->message_box->hide();
    ui->message_box->setWidget(message_content);
    message_content->setText("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nkai shi le xie xie");
    message_content->append(QString::fromLocal8Bit("<font color=red>错误</font>"));
    QScrollBar *pScroll = message_content->verticalScrollBar();
    pScroll->setSliderPosition(pScroll->maximum());
}

/// 设置vtk显示窗口
void MainWindow::createVTKview()
{
    qvtkWidget = new QVTKWidget(this);
    setCentralWidget(qvtkWidget);
    renderer = vtkRenderer::New();
    // Setup VTK window

    renderer->SetBackground(0.5, 0.5, 1);
    renderer->SetBackground2(1,1,1);
    renderer->SetGradientBackground(1);
    qvtkWidget->GetRenderWindow()->AddRenderer(renderer);

    vtkAxesActor* axes = vtkAxesActor::New();
    widget = vtkOrientationMarkerWidget::New();
    widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    widget->SetOrientationMarker( axes );
    widget->SetInteractor( qvtkWidget->GetInteractor() );
    widget->SetViewport( 0.0, 0.0, 0.3, 0.3 );
    widget->On();
    widget->InteractiveOn();

    renderer->ResetCamera();
    vtkSmartPointer<vtkCamera> camera =  renderer->GetActiveCamera();
    camera->SetPosition(1,1,1);
    qvtkWidget->GetRenderWindow()->Render();
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

// slots ====================================================================================

/// 显示XY平面

void MainWindow::show_XY()
{
    vtkSmartPointer<vtkCamera> camera =  renderer->GetActiveCamera();
    camera->SetPosition(0,0,1);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,1,0);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
}

/// 显示XZ平面
void MainWindow::show_XZ()
{
    vtkSmartPointer<vtkCamera> camera =  renderer->GetActiveCamera();
    camera->SetPosition(0,1,0);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,0,1);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
}

/// 显示YZ平面
void MainWindow::show_YZ()
{
    vtkSmartPointer<vtkCamera> camera =  renderer->GetActiveCamera();
    camera->SetPosition(1,0,0);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,0,1);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
}


/// 显示坐标轴
void MainWindow::show_axes(bool flags)
{
    if(flags){
        widget->On();
        qvtkWidget->GetRenderWindow()->Render();
    }
    else{
        widget->Off();
        qvtkWidget->GetRenderWindow()->Render();
    }
}

/// 显示颗粒增强单胞
void MainWindow::show_pcmcell()
{
    if (sots.celldata == NULL)return;


    vtkSmartPointer<vtkAssignAttribute> mat =
            vtkSmartPointer<vtkAssignAttribute>::New();
    mat->SetInput(sots.celldata);
    mat->Assign("SubdomainID", vtkDataSetAttributes::SCALARS,
                vtkAssignAttribute::CELL_DATA);
    // Create a lookup table to map cell data to colors
    vtkSmartPointer<vtkLookupTable> lut =
            vtkSmartPointer<vtkLookupTable>::New();
    //    int tableSize = std::max(resolution*resolution + 1, 10);
    lut->SetNumberOfTableValues(2);
    lut->Build();

    // Fill in a few known colors, the rest will be generated if needed
    lut->SetTableValue(0, 1, 0, 0, 1);
    lut->SetTableValue(1, 0, 0, 1, 1);
    //    lut->SetTableValue(2, 1.0000, 0.3882, 0.2784, 1); // Tomato
    //    lut->SetTableValue(3, 0.9608, 0.8706, 0.7020, 1); // Wheat
    //    lut->SetTableValue(4, 0.9020, 0.9020, 0.9804, 1); // Lavender
    //    lut->SetTableValue(5, 1.0000, 0.4900, 0.2500, 1); // Flesh
    //    lut->SetTableValue(6, 0.5300, 0.1500, 0.3400, 1); // Raspberry
    //    lut->SetTableValue(7, 0.9804, 0.5020, 0.4471, 1); // Salmon
    //    lut->SetTableValue(8, 0.7400, 0.9900, 0.7900, 1); // Mint
    //    lut->SetTableValue(9, 0.2000, 0.6300, 0.7900, 1); // Peacock
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInput(sots.celldata);
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(1, 2);
    mapper->ScalarVisibilityOn();
    //   cout << sots.celldata;
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();
    /*
    vtkSmartPointer<vtkParametricEllipsoid> ellip = vtkSmartPointer<vtkParametricEllipsoid>::New();
    ellip->SetXRadius(1.0);
    ellip->SetYRadius(2.0);
    ellip->SetZRadius(1.0);

    //vtkSmartPointer<vtkSphereSource> ellipSource =
    //vtkSmartPointer<vtkSphereSource>::New();
    //ellipSource->SetCenter(.5, 0.5, 0);
    //ellipSource->SetThetaResolution(600);
    //ellipSource->SetPhiResolution(600);
    //ellipSource->Update();

    vtkSmartPointer<vtkPolyData> polyData;

    vtkSmartPointer<vtkCubeSource> cubeSource =
            vtkSmartPointer<vtkCubeSource>::New();
    cubeSource->SetCenter(.5, 0.5, 0.5);
    cubeSource->Update();

    vtkSmartPointer<vtkBox> implicitCube =
            vtkSmartPointer<vtkBox>::New();
    implicitCube->SetBounds(cubeSource->GetOutput()->GetBounds());


    vtkSmartPointer<vtkParametricFunctionSource> ellipSource =
            vtkSmartPointer<vtkParametricFunctionSource>::New();
    ellipSource->SetParametricFunction(ellip);

    vtkSmartPointer<vtkClipPolyData> clipper =
            vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputConnection(ellipSource->GetOutputPort());
    clipper->SetClipFunction(implicitCube);
    //clipper->SetValue(0);
    clipper->InsideOutOn();
    clipper->Update();


    polyData = clipper->GetOutput();


    vtkSmartPointer<vtkFeatureEdges> boundaryEdges =
            vtkSmartPointer<vtkFeatureEdges>::New();
#if VTK_MAJOR_VERSION <= 5
    boundaryEdges->SetInput(polyData);
#else
    boundaryEdges->SetInputData(polyData);
#endif
    boundaryEdges->BoundaryEdgesOn();
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();

    vtkSmartPointer<vtkStripper> boundaryStrips =
            vtkSmartPointer<vtkStripper>::New();
    boundaryStrips->SetInputConnection(boundaryEdges->GetOutputPort());
    boundaryStrips->Update();

    // Change the polylines into polygons
    vtkSmartPointer<vtkPolyData> boundaryPoly =
            vtkSmartPointer<vtkPolyData>::New();
    boundaryPoly->SetPoints(boundaryStrips->GetOutput()->GetPoints());
    boundaryPoly->SetPolys(boundaryStrips->GetOutput()->GetLines());


    vtkSmartPointer<vtkPolyDataMapper> boundaryMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
    boundaryMapper->SetInput(boundaryPoly);
#else
    boundaryMapper->SetInputData(boundaryPoly);
#endif

    vtkSmartPointer<vtkActor> boundaryActor =
            vtkSmartPointer<vtkActor>::New();
    boundaryActor->SetMapper(boundaryMapper);
    boundaryActor->GetProperty()->SetColor(0.,1,1);
    boundaryActor->GetProperty()->SetEdgeColor(0,0,0);
    //boundaryActor->GetProperty()->EdgeVisibilityOn();
    //boundaryActor->GetProperty()->SetRepresentationToWireframe();


    vtkSmartPointer<vtkPolyDataMapper>  mapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    //mapper->SetInputConnection(ellipSource->GetOutputPort());
    mapper->SetInput(polyData);

    vtkSmartPointer<vtkActor>  actor =
            vtkSmartPointer<vtkActor>::New() ;
    actor->SetMapper(mapper);

    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    //actor->GetProperty()->SetEdgeColor(0,0,0);
    //actor->GetProperty()->EdgeVisibilityOn();
  */
    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
    //   renderer->AddActor(boundaryActor);
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
}

/// 设置单胞
void MainWindow::set_cell()
{
    SetCellDialog setcell(this);
    setcell.exec();
}


/// 导入单胞几何并剖分
void MainWindow::import_cell_geo()
{

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
        sots.celldata = import_bdf(sFileName);
 //       cout << sots.celldata << endl;
    }
}

/// 创建材料
void MainWindow::new_material()
{
    materiel mat(this);
    mat.exec();
}

/// 直接导入计算文件
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

/// 新建工程
void MainWindow::new_project()
{
    AnalysisType anatype(this);
    if(anatype.exec() == QDialog::Accepted)
    {
        sots.set_analysis_type(anatype.get_analysis_type().toStdString());
        sots.set_prj_folder(anatype.get_prj_folder().toStdString());
        sots.set_prj_name(anatype.get_prj_name().toStdString());
    }
    QString logfile = anatype.get_prj_folder() + QString("/") + anatype.get_prj_name() +".log";
    cout << logfile.toStdString().c_str() << endl;
    //    freopen(logfile.toStdString().c_str(),"w",stdout);
    cout << "这是一个分析文件" << endl;
    //freopen("CON", "w", stdout);
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
        read_neutral_format(sFileName);
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

vtkUnstructuredGrid* MainWindow::import_bdf(QString bdf_filename)
{

    vtkIdType number_of_points = 0, number_of_tetra = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkIntArray> intValue = vtkSmartPointer<vtkIntArray>::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("subdomainid");

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
    //   std::cout << number_of_points << std::endl;
    //   std::cout << number_of_tetra << std::endl;
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
                //                intValue->InsertNextValue(subdomainId);
                intValue->InsertNextValue(a);

                tetra->GetPointIds()->SetId(0, a-1);
                tetra->GetPointIds()->SetId(1, b-1);
                tetra->GetPointIds()->SetId(2, c-1);
                tetra->GetPointIds()->SetId(3, d-1);
                //                std::cout << n << std::endl;
                cellArray->InsertNextCell(tetra);
            }
        }

    }
    vtkUnstructuredGrid* unstructuredGrid = vtkUnstructuredGrid::New();
//    cout << unstructuredGrid << endl;
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
    return unstructuredGrid;
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
    //mapper->SetInputConnection(reader->GetOutputPort());
    mapper->SetInput(unstructuredGrid);
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
    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
}

void MainWindow::plot_cell()
{

}

