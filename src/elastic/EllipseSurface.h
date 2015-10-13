//===========================================================================
// EllipseSurface.h
// ��������ͷ�ļ�
// A class of Ellipse Surface
//===========================================================================
#ifndef ELLIPSESURFACE_H
#define ELLIPSESURFACE_H

#include "Geometry.h"
#include <iostream>

//---------------------------------------------------------------------------
class EllipseSurface : public Surface 
{
	public:
        EllipseSurface(double xx, double yy, double zz, double xy, double yz,
                       double zx, double x,  double y,  double z,  double c );
        EllipseSurface(double     x0, double     y0, double     z0,
                       double      a, double      b, double      c,
                       double alpha1, double alpha2, double alpha3,
                       double  beta1, double  beta2, double  beta3,
                       double  gama1, double  gama2, double  gama3 );
        //�жϸ������Ƿ������ڣ����ظ���������������򷽳̵�ֵ  
        int is_contain_usually( Point *point, int mode = 0 );	//һ�㷽��
        int is_contain( Point *point );

        int intersect(Point &point1, Point &point2, Point &rpoint);
        int intersect(Line &line, Point &rpoint);
        double dis_to(Point *point);     
        
        //����������꣨x,y,z)���������淽�̵�ֵ--F(x,y,z)
        double fvalue(double x, double y, double z);

        //to compute the point which the given point project to the curve on the normal direction;
		//����������ظ�������ͶӰ�������ϵ�ͶӰ��
        Point project_normal(Point *thepoint, int *error=NULL);
        //to compute the point which the given point project to the surface on the given direction;
		//����������ظ�������ͶӰ�����ϵ�ͶӰ��
        Point project(Point *thepoint, TDVector *vec, int *error=NULL, int mod = 0);

		//����һ�������������򣬶�������ÿ����������任
		int change_coor(Point *opoint, Point *point, double ratio);

        void print( int mod = 0 );
	private:
        double xx_c,yy_c,zz_c,xy_c,yz_c,zx_c,x_c,y_c,z_c,c_c;
        double x0, y0, z0;
        double a, b, c;

        double coe1,coe2,coe3,coe4,coe5,coe6,rhs;


        //�Ѹ�����15������ת��������һ�㷽�̵�10������
        void translate_para(double     x0, double     y0, double     z0,
                            double      a, double      b, double      c,
                            double alpha1, double alpha2, double alpha3,
                            double  beta1, double  beta2, double  beta3,
                            double  gama1, double  gama2, double  gama3 );

        //��һԪ���η��̣��޽ⷵ��0,���������x������ֵ��
        int sol_2order_equ(double a, double b, double c, double* x);

        //��һ���������Է����飨A*x=B),�������x��
        int sol_3order_lina_equs(double a[3][3], double b[3], double x[3]);

        //����3������ʽ��ֵ������������Է�����ʱʹ�ã�
        double det_3order(double a[3][3]);         

        //����Fx(x,y,z)
        double fxvalue(double x, double y, double z);
        //����Fy(x,y,z)
        double fyvalue(double x, double y, double z);
        //����Fz(x,y,z)
        double fzvalue(double x, double y, double z);
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
