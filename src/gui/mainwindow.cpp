#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "aboutdialog.h"
#include <QFileDialog>
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
#include "Gmsh.h"
#include "GModel.h"
#include "MElement.h"
#include <vtkGenericDataObjectReader.h>
#include <vtkUnstructuredGridVolumeRayCastMapper.h>
#include <QTextEdit>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setWindowIcon(QIcon(":/images/elisa_128px_540496_easyicon.net.png"));
    qvtkWidget = new QVTKWidget(this);
    setCentralWidget(qvtkWidget);
    QTextEdit* message_content = new QTextEdit;
    message_content->setReadOnly(1);
    ui->message_box->hide();
    ui->message_box->setWidget(message_content);

    renderer = vtkRenderer::New();
    renderer->SetBackground(169/255., 169./255, 169./255);

    qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
    qvtkWidget->GetRenderWindow()->Render();
    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
    qvtkWidget->show();
}

MainWindow::~MainWindow()
{
    delete qvtkWidget;
    delete ui;

}

void MainWindow::on_action_Mscs_triggered()
{
    AboutDialog aboutdialog(this);
    aboutdialog.exec();
}

void MainWindow::on_action_3_triggered()
{
    MeshPara meshpara(this);
    meshpara.exec();
}

void MainWindow::on_Openfile_triggered()
{
    //QString homeName = getHomePath();
    QString sDefaultName = tr("/");
    QString selectedFilter = tr("All Files(*.*)");
    //    QDir defaultDir = QFSFileEngine::homePath();
    QString filter = tr("All Files(*.*);;Dat Files(*.dat)");
    filter += tr(";;MSH Files(*.msh);;BDF File(*.bdf)");
    filter += tr(";;STL Files(*.stl);;GMSH POS File(*.pos)");
   // filter += tr(";;IGES Files(*.igs *.iges)");
   // filter += tr(";;STEP Files(*.stp *.step)");
   // filter += tr(";;BREP Files(*.brep)");
    QString sFileName = QFileDialog::getOpenFileName(
                this, QString(tr("打开文件")).toLocal8Bit(),
                sDefaultName,
                filter,
                /*NULL*/&selectedFilter,
                QFileDialog::DontUseNativeDialog);
    std::cout << sFileName.toStdString() << std::endl;
   // openSTL(sFileName);
   read_neutral_format(sFileName);
   // meshSTL(sFileName);
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
    actor->GetProperty()->EdgeVisibilityOn();

    renderer->AddActor(actor);

    qvtkWidget->GetRenderWindow()->Render();
    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
 //   return 1;

}

void MainWindow::meshSTL(QString stl_filename)
{
    std::string inputFilename = stl_filename.toStdString();

    GmshInitialize(NULL, 0);
    GmshSetOption("Mesh", "Algorithm3D", 0.1);
    GModel *m = new GModel();
    m->readSTL(inputFilename.c_str());
    //GmshMergeFile("../../tutorial/t5.geo"); // will also set the bbox
    m->mesh(3);
    for(GModel::riter it = m->firstRegion(); it != m->lastRegion(); ++it){
      GRegion *r = *it;
      printf("volume %d contains %d elements:\n", r->tag(), r->getNumMeshElements());
      for(unsigned int i = 0; i < r->getNumMeshElements(); i++)
        printf(" %d", r->getMeshElement(i)->getNum());
      printf("\n");
    }
    m->writeVTK((inputFilename+".vtu").c_str());
 //   m->writeUNV("test.unv");
    delete m;
    GmshFinalize();

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
    actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();
    renderer->AddActor(actor);


    renderer->ResetCamera();
    qvtkWidget->GetRenderWindow()->Render();
    qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
 //   return 1;

}

void MainWindow::openVtu(QString vtu_filename)
{
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
}

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
       actor->GetProperty()->SetRepresentationToWireframe();
       actor->GetProperty()->SetColor(0,1,1);
       actor->GetProperty()->SetEdgeColor(0,0,0);
       actor->GetProperty()->EdgeVisibilityOn();
       renderer->AddActor(actor);

   //      qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
       renderer->ResetCamera();
       qvtkWidget->GetRenderWindow()->Render();
       qvtkWidget->GetRenderWindow()->GetInteractor()->Initialize();
       qvtkWidget->GetRenderWindow()->GetInteractor()->Start();
}

void MainWindow::on_action_show_massage_box_triggered()
{
    ui->message_box->show();
}

void MainWindow::on_action_show_massage_box_triggered(bool checked)
{
    if (checked)
        ui->message_box->show();
    else
        ui->message_box->hide();
}
