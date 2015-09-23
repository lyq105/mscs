#ifndef SOTSINTERFACE_H
#define SOTSINTERFACE_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <string>


struct Cell_info
{
    std::string cell_type;
    double x0[3],x1[3];
    double value[100];
};

struct Material
{
    int number;
    std::string mat_type;
    double value[100];
};

// This class provide an interface for interchange data between gui
// and SOTS solver.

class SOTSinterface
{
public:
    SOTSinterface();
    ~SOTSinterface();
    void set_cell_mesh();
    void set_material(int mat_index);

    void set_prj_name(std::string name){ prj_name = name;}
    void set_prj_folder(std::string folder){ prj_folder = folder;}
    void set_analysis_type(std::string type){ analysis_type = type;}
    void set_cell_info(Cell_info& cell_info);


    vtkSmartPointer<vtkUnstructuredGrid> get_cell_data();
    void get_material(int mat_index);

    void write_cell_solution(std::string filename);
    void write_homo_solution(std::string filename);
private:
    vtkSmartPointer<vtkUnstructuredGrid> celldata;
    vtkSmartPointer<vtkUnstructuredGrid> homodata;
    std::string analysis_type;
    std::string prj_folder;
    std::string prj_name;
    Cell_info cell_info;
};


//class SOTSSolver
//{
//public:
//    SOTSSolver(){};
//    ~SOTSSolver(){};
//}
#endif // SOTSINTERFACE_H
