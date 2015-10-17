#include "vtkviewer.h"
#include <QVTKWidget.h>
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
#include <vtkCellData.h>
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

#include <sstream>



VTKviewer::VTKviewer()
{
    qvtkWidget = new QVTKWidget;
    renderer = vtkRenderer::New();
    // Setup VTK window

    renderer->SetBackground(0.5, 0.5, 1);
    renderer->SetBackground2(1,1,1);
    renderer->SetGradientBackground(1);
    qvtkWidget->GetRenderWindow()->AddRenderer(renderer);

    vtkAxesActor* axes = vtkAxesActor::New();
    axes_widget = vtkOrientationMarkerWidget::New();
    axes_widget->SetOutlineColor( 0.9300, 0.5700, 0.1300 );
    axes_widget->SetOrientationMarker( axes );
    axes_widget->SetInteractor( qvtkWidget->GetInteractor() );
    axes_widget->SetViewport( 0.0, 0.0, 0.3, 0.3 );
    axes_widget->On();
    axes_widget->InteractiveOn();

    renderer->ResetCamera();
    vtkSmartPointer<vtkCamera> camera
            =  renderer->GetActiveCamera();
    camera->SetPosition(1,1,1);
    qvtkWidget->GetRenderWindow()->Render();
}

void
VTKviewer::render()
{
    qvtkWidget->GetRenderWindow()->Render();
}


void
VTKviewer::show_XY()
{
    vtkSmartPointer<vtkCamera> camera
            =  renderer->GetActiveCamera();
    camera->SetPosition(0,0,1);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,1,0);
    renderer->ResetCamera();
    render();
}
void
VTKviewer::show_XZ()
{
    vtkSmartPointer<vtkCamera> camera
            =  renderer->GetActiveCamera();
    camera->SetPosition(0,1,0);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,0,1);
    renderer->ResetCamera();
    render();
}
void
VTKviewer::show_YZ()
{
    vtkSmartPointer<vtkCamera> camera
            =  renderer->GetActiveCamera();
    camera->SetPosition(1,0,0);
    camera->SetFocalPoint(0,0,0);
    camera->SetViewUp(0,0,1);
    renderer->ResetCamera();
    render();
}

void
VTKviewer::axes_on()
{
    axes_widget->On();
    render();
}
void
VTKviewer::axes_off()
{
    axes_widget->Off();
    render();
}

/// 载入数据
vtkPolyData*
VTKviewer::load_stl(std::string filename)
{
    return NULL;
}

vtkUnstructuredGrid*
VTKviewer::load_netgen_mesh(std::string filename)
{
    std::ifstream infile(filename.c_str());
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
    vtkUnstructuredGrid* unstructuredGrid
            = vtkUnstructuredGrid::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);

    return unstructuredGrid;
}


vtkUnstructuredGrid*
VTKviewer::load_bdf(std::string bdf_file)
{
    vtkIdType number_of_points = 0, number_of_tetra = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> matArray = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkCellArray> reinArray = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkIntArray> intValue = vtkSmartPointer<vtkIntArray>::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("subdomainid");

    std::ifstream infile(bdf_file.c_str());
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
        std::ifstream infile(bdf_file.c_str());
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
                //if(subdomainId == 2) continue;
                //                std::cout << n << subdomainId << a << b << c <<d << std::endl;
                vtkSmartPointer<vtkTetra> tetra =
                        vtkSmartPointer<vtkTetra>::New();
                intValue->InsertNextValue(subdomainId);
                //intValue->InsertNextValue(a);

                tetra->GetPointIds()->SetId(0, a-1);
                tetra->GetPointIds()->SetId(1, b-1);
                tetra->GetPointIds()->SetId(2, c-1);
                tetra->GetPointIds()->SetId(3, d-1);
                //                std::cout << n << std::endl;
                cellArray->InsertNextCell(tetra);
                if(subdomainId == 1) matArray->InsertNextCell(tetra);
                else reinArray->InsertNextCell(tetra);
            }
        }

    }
    vtkUnstructuredGrid* unstructuredGrid = vtkUnstructuredGrid::New();
    //    cout << unstructuredGrid << endl;
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);
    unstructuredGrid->GetCellData()->SetScalars(intValue);

    matrix = vtkUnstructuredGrid::New();
    matrix->SetPoints(points);
    matrix->SetCells(VTK_TETRA, matArray);

    reinforcement = vtkUnstructuredGrid::New();
    reinforcement->SetPoints(points);
    reinforcement->SetCells(VTK_TETRA, reinArray);

    return unstructuredGrid;
}


void
VTKviewer::write_ug(vtkUnstructuredGrid* ug,std::string ug_file)
{
    // Write file
    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
            vtkSmartPointer<vtkUnstructuredGridWriter>::New();
    writer->SetFileName(ug_file.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(ug);
#else
    writer->SetInputData(ug);
#endif
    writer->Write();
}

/// 显示
void
VTKviewer::show_ug_mesh(vtkUnstructuredGrid* ug)
{
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(ug);
#else
    mapper->SetInputData(ug);
#endif
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();

    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
    renderer->ResetCamera();
    render();
}

void
VTKviewer::show_ug_scalar(vtkUnstructuredGrid* ug,std::string scalar)
{
    //    vtkSmartPointer<vtkAssignAttribute> mat =
    //                vtkSmartPointer<vtkAssignAttribute>::New();
    //        mat->SetInput(sots.celldata);
    //       // mat->Assign("subdomainid", vtkDataSetAttributes::SCALARS,
    //       //             vtkAssignAttribute::CELL_DATA);
    //        // Create a lookup table to map cell data to colors
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
    double subrange[2];
    ug->GetCellData()->GetArray(scalar.c_str())->GetRange(subrange);
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    // mapper->SetLookupTable(lut);

#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(ug);
#else
    mapper->SetInputData(ug);
#endif
    mapper->SelectColorArray(scalar.c_str());
    mapper->SetScalarRange(subrange);

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();

    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
    renderer->ResetCamera();
    render();
}


void
VTKviewer::show_surface(vtkUnstructuredGrid*ug)
{

}
void VTKviewer::show_stl(std::string stl_file)
{
    std::string inputFilename = stl_file;

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
}


void VTKviewer::show_matrix()
{
    // ug->GetCellData()->GetArray("MaterialID");
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(matrix);
#else
    mapper->SetInputData(matrix);
#endif
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,1);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();

    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
    renderer->ResetCamera();
    render();
}

void VTKviewer::show_reinforcement()
{
    // ug->GetCellData()->GetArray("MaterialID");
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();

#if VTK_MAJOR_VERSION <= 5
    mapper->SetInput(reinforcement);
#else
    mapper->SetInputData(reinforcement);
#endif
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToSurface();
    //actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetColor(0,1,0);
    actor->GetProperty()->SetEdgeColor(0,0,0);
    actor->GetProperty()->EdgeVisibilityOn();

    renderer->RemoveAllViewProps();

    renderer->AddActor(actor);
    renderer->ResetCamera();
    render();
}

void VTKviewer::show_cell()
{

}

