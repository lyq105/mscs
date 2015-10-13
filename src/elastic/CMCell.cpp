//===========================================================================
// CMCell.cpp
// ���ϲ��ϵ������Ա����
// Member Functions in a Class of Cell of Composite Material
//
// ���ϲ��ϵ����Ķ��壺�����������ǿ���ϵ���Ϣ�ͼ��β���
// ���ϲ��ϵ�������һ�����֮࣬�������������������
// ȡ����߶�����1��С��������Ϊ��Ԫ����(unit cell)������Ĭ��ֵ
// ȡ����ģ���Ͳ��ɱȷֱ�Ϊ200��0.3��Ϊ�����Ĭ������
// ȡ����ģ���Ͳ��ɱȷֱ�Ϊ200��0.3��Ϊ��ǿ���ϵ�Ĭ������
//===========================================================================
#include "CMCell.h"

//---------------------------------------------------------------------------
//���캯����ȱʡΪ��������߶���1
CMCell::CMCell(double length, double width, double height)
{
	//�����
	clength = length ? length : 1;
	cwidth  = width  ? width  : 1;
	cheight = height ? height : 1;
	//��ʼ������
	origin_x=0;
	origin_y=0;
	origin_z=0;
	//���ɵ����������8�������һ�������������Ԫ����ÿ�������������ֵ
	//���ɵ����������12���ߵ�һ�������������Ԫ���Ǳߵ������˵㣬�ɵ���������
	generate();
}

//---------------------------------------------------------------------------
//���õ���������ĳ�ʼ����Ϣ�������ɵ����ĵ������ͱ�����
void CMCell::set_origin(double x,double y,double z)
{
	origin_x=x;
	origin_y=y;
	origin_z=z;

	generate();
}

//---------------------------------------------------------------------------
//���ɵ����������8�������12������Ϣ
void CMCell::generate()
{
	//���ɵ�����8��������
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

	//���ɵ�����12��������
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
//���ص�����һЩ��Ϣ���麯����
double CMCell::length() { return clength; }	//���ص����ĳ�
double CMCell::width()  { return cwidth;  }	//���ص����Ŀ�
double CMCell::height() { return cheight; }	//���ص����ĸ�

int CMCell::get_mat_num() { return material_num; }		//���ز���������

double CMCell::volume() { return clength*cwidth*cheight; }	//���ص��������

//===========================================================================
