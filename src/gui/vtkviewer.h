#ifndef VTKVIEWER_H
#define VTKVIEWER_H

//前置声明
class QVTKWidget;
class VTKActor;
class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkAppendPolyData;
class vtkOrientationMarkerWidget;

class VTKviewer
{
public:
    VTKviewer();
    ~VTKviewer();

    void render();
    QVTKWidget* getVTKWidget(){return qvtk;}

private:
    QVTKWidget* qvtk;
    vtkRenderer* renderer;
    vtkOrientationMarkerWidget* widget;

    vtkActor* cell_mesh_actor;
    vtkActor* cell_geo_actor;
    vtkActor* cell_actor;
};

#endif // VTKVIEWER_H
