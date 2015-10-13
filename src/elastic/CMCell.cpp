//===========================================================================
// CMCell.cpp
// 复合材料单胞类成员函数
// Member Functions in a Class of Cell of Composite Material
//
// 复合材料单胞的定义：包括基体和增强材料的信息和几何参数
// 复合材料单胞类是一个基类，之后将由它派生出许多子类
// 取长宽高都等于1的小立方体作为单元单胞(unit cell)，它是默认值
// 取杨氏模量和泊松比分别为200和0.3作为基体的默认属性
// 取杨氏模量和泊松比分别为200和0.3作为增强材料的默认属性
//===========================================================================
#include "CMCell.h"

//---------------------------------------------------------------------------
//构造函数，缺省为单胞长宽高都是1
CMCell::CMCell(double length, double width, double height)
{
	//长宽高
	clength = length ? length : 1;
	cwidth  = width  ? width  : 1;
	cheight = height ? height : 1;
	//起始点坐标
	origin_x=0;
	origin_y=0;
	origin_z=0;
	//生成单胞立方体的8个顶点的一组点向量，向量元素是每个点的三个坐标值
	//生成单胞立方体的12条边的一组边向量，向量元素是边的两个端点，由点向量描述
	generate();
}

//---------------------------------------------------------------------------
//设置单胞立方体的初始点信息，并生成单胞的点向量和边向量
void CMCell::set_origin(double x,double y,double z)
{
	origin_x=x;
	origin_y=y;
	origin_z=z;

	generate();
}

//---------------------------------------------------------------------------
//生成单胞立方体的8个顶点和12条边信息
void CMCell::generate()
{
	//生成单胞的8个点向量
	Point corner_points[8] = {{origin_x,origin_y,origin_z},
							 {origin_x+clength,origin_y,origin_z},
							 {origin_x+clength,origin_y+cwidth,origin_z},
							 {origin_x,origin_y+cwidth,origin_z},
							 {origin_x,origin_y,origin_z+cheight},
							 {origin_x+clength,origin_y,origin_z+cheight},
							 {origin_x+clength,origin_y+cwidth,origin_z+cheight},
							 {origin_x,origin_y+cwidth,origin_z+cheight}};

	points_vec.clear();
	for (int i=0; i<8; i++)
	{
		points_vec.push_back(corner_points[i]);
	}

	//生成单胞的12条边向量
	Line edge_lines[12] = {	Line(points_vec[0],points_vec[1]),
							Line(points_vec[1],points_vec[2]),
							Line(points_vec[2],points_vec[3]),
							Line(points_vec[3],points_vec[0]),
							Line(points_vec[0],points_vec[4]),
							Line(points_vec[1],points_vec[5]),
							Line(points_vec[2],points_vec[6]),
							Line(points_vec[3],points_vec[7]),
							Line(points_vec[4],points_vec[5]),
							Line(points_vec[5],points_vec[6]),
							Line(points_vec[6],points_vec[7]),
							Line(points_vec[7],points_vec[4])	};

	lines_vec.clear();
	for (int i=0; i<12; i++)
	{
		lines_vec.push_back(edge_lines[i]);
	}
}

//---------------------------------------------------------------------------
//返回单胞的一些信息（虚函数）
double CMCell::length() { return clength; }	//返回单胞的长
double CMCell::width()  { return cwidth;  }	//返回单胞的宽
double CMCell::height() { return cheight; }	//返回单胞的高

int CMCell::get_mat_num() { return material_num; }		//返回材料种类数

double CMCell::volume() { return clength*cwidth*cheight; }	//返回单胞的体积

//===========================================================================
