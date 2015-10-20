#ifndef SOTSINTERFACE_H
#define SOTSINTERFACE_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <string>
#include <QThread>

class Mesher;
class PMMesher;

struct PMCell_Info
{
    /// 单胞类型
    std::string cell_type;
    /// 单胞对角顶点
    double x0[3],x1[3];
    /// 颗粒类型
    int type[3];
    /// 是否添加界面层
    bool coating;
    /// 网剖分参数
    double global_meshsize;
    ///轴长参数设置
    double amax[3],amin[3],bmin[3],bmax[3];
    /// 倾角分布参数
    int d_sign;
    double d_value[4];
    /// 中心坐标分布参数
    int coord_sign;
    double coord_value[6];
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
    void set_pmcell_info(PMCell_Info& cell_info);


    vtkUnstructuredGrid* get_cell_data();
    PMCell_Info* get_pmcell_info(){return cell_info;}
    void get_material(int mat_index);

    void write_cell_solution(std::string filename);
    void write_homo_solution(std::string filename);
    /// 生成颗粒增强单胞
    int build_pmcell(std::string in_file, std::string data_file, int CNum=1);
    int build_pmcell(std::string ellips_para, std::string ellips_file,std::string data_file);
    int mesh_pmcell(std::string meshpara);
    std::string Get_Line(std::ifstream &infile)const;

    /// 直接求解单胞问题
    int solve(std::string in_file,std::string out_file,std::string data_file);
    int process_pmcell(const PMMesher *mesh);

    vtkUnstructuredGrid* celldata;
    vtkUnstructuredGrid* homodata;
private:
    std::string analysis_type;
    std::string prj_folder;
    std::string prj_name;
    PMCell_Info* cell_info;
};


//class SOTSSolver
//{
//public:
//    SOTSSolver(){};
//    ~SOTSSolver(){};
//}
#endif // SOTSINTERFACE_H
