//===========================================================================
// CMCell.h
// ���ϲ��ϵ�����ͷ�ļ�
// A Class of Cell of Composite Material
//
// ���ϲ��ϵ����Ķ��壺�����������ǿ���ϵ���Ϣ�ͼ��β���
// ���ϲ��ϵ�������һ�����֮࣬�������������������
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
		//���캯��
		CMCell(double length=0, double width=0, double height=0); 
		//�麯��
        virtual double length();						//���ص����ĳ�
        virtual double width(); 						//���ص����Ŀ�
        virtual double height() ;						//���ص����ĸ�

        virtual int get_mat_num();						//���ز���������

        virtual double volume();						//���ص��������

        void set_origin(double x,double y,double z);
        void print();									//���������Ϣ

		//�����ļ��β���
		double clength, cwidth, cheight;
        double origin_x,origin_y,origin_z;				//�����������
		double Unitcell_V;

        vector<Point> points_vec;						//�洢�����Ϣ�������ǵ����������8������
        vector<Line>  lines_vec;							//�洢���е��߶���Ϣ, �����ǵ����������12�����߶���Ϣ  
        vector<Surface*> surfaces_vec;				//�洢���е�����Ϣ, ��ǿ��Ľ���

	protected:                
        int material_num ;

	private:
        void generate();       //���ɵ����������8�������12������Ϣ
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
