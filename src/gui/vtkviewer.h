#ifndef VTKVIEWER_H
#define VTKVIEWER_H


#include <string>

//前置声明
class QVTKWidget;
class VTKActor;
class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkPolyData;
class vtkAppendPolyData;
class vtkUnstructuredGrid;
class vtkOrientationMarkerWidget;


class VTKviewer
{
public:
    VTKviewer();
    ~VTKviewer();

    void render();

    /// 坐标轴管理
    void show_XY();
    void show_XZ();
    void show_YZ();
    void axes_on();
    void axes_off();


    /// 载入数据
    vtkPolyData* load_stl(std::string filename);
    vtkUnstructuredGrid* load_netgen_mesh(std::string filename);
    vtkUnstructuredGrid* load_bdf(std::string bdf_file);


    /// 显示
    void show_ug_mesh(vtkUnstructuredGrid* ug);
    void show_surface(vtkUnstructuredGrid* ug);
    void show_stl(std::string stl_file);
    void write_ug(vtkUnstructuredGrid *ug, std::string ug_file);


    QVTKWidget* getVTKWidget(){return qvtkWidget;}

private:
    QVTKWidget* qvtkWidget;
    vtkRenderer* renderer;
    vtkOrientationMarkerWidget* axes_widget;

    vtkActor* cell_mesh_actor;
    vtkActor* cell_geo_actor;
    vtkActor* cell_actor;
};

#endif // VTKVIEWER_H
