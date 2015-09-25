//===========================================================================
// CMCell.h
// 复合材料单胞类头文件
// A Class of Cell of Composite Material
//
// 复合材料单胞的定义：包括基体和增强材料的信息和几何参数
// 复合材料单胞类是一个基类，之后将由它派生出许多子类
//===========================================================================

#ifndef CMCELL_H
#define CMCELL_H

#include <vector>
#include "Geometry.h"
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
class CMCell
{
	public:
		//构造函数
		CMCell(double length=0, double width=0, double height=0); 
		//虚函数
		virtual double length();						//返回单胞的长
		virtual double width(); 						//返回单胞的宽
		virtual double height() ;						//返回单胞的高

		virtual int get_mat_num();						//返回材料种类数

		virtual double volume();						//返回单胞的体积

		void set_origin(double x,double y,double z);
		void print();									//输出单胞信息

		//单胞的几何参数
		double clength, cwidth, cheight;
		double origin_x,origin_y,origin_z;				//单胞起点坐标
		double Unitcell_V;

		vector<Point> points_vec;						//存储点的信息，首先是单胞立方体的8个顶点
		vector<Line>  lines_vec;							//存储所有的线段信息, 首先是单胞立方体的12条边线段信息  
		vector<Surface*> surfaces_vec;				//存储所有的面信息, 增强项的界面

	protected:                
		int material_num ;

	private:
		void generate();       //生成单胞立方体的8个顶点和12条边信息
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
