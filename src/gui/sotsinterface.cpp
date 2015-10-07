#include "sotsinterface.h"

SOTSinterface::SOTSinterface()
{
}

SOTSinterface::~SOTSinterface()
{

}
void SOTSinterface::set_cell_mesh()
{
}
//void SOTSinterface::
void SOTSinterface::set_material(int mat_index)
{

}
void SOTSinterface::get_material(int mat_index)
{

}
void SOTSinterface::write_cell_solution(std::string filename)
{

}
void SOTSinterface::write_homo_solution(std::string filename)
{

}
void SOTSinterface::set_pmcell_info(PMCell_Info& mycell_info)
{
    //cell_info.coating = mycell_info.coating;

}
vtkSmartPointer<vtkUnstructuredGrid> SOTSinterface::get_cell_data()
{
    return celldata;
}
