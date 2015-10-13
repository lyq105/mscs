//===========================================================================
// EllipseSurface.h
// 椭球面类头文件
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
        //判断给定点是否在面内，返回给定点坐标带入椭球方程的值  
        int is_contain_usually( Point *point, int mode = 0 );	//一般方程
        int is_contain( Point *point );

        int intersect(Point &point1, Point &point2, Point &rpoint);
        int intersect(Line &line, Point &rpoint);
        double dis_to(Point *point);     
        
        //计算给定坐标（x,y,z)代入椭球面方程的值--F(x,y,z)
        double fvalue(double x, double y, double z);

        //to compute the point which the given point project to the curve on the normal direction;
		//计算给定点沿给定方向投影到曲线上的投影点
        Point project_normal(Point *thepoint, int *error=NULL);
        //to compute the point which the given point project to the surface on the given direction;
		//计算给定点沿给定方向投影到面上的投影点
        Point project(Point *thepoint, TDVector *vec, int *error=NULL, int mod = 0);

		//根据一定比例紧缩椭球，对椭球上每个点做坐标变换
		int change_coor(Point *opoint, Point *point, double ratio);

        void print( int mod = 0 );
	private:
        double xx_c,yy_c,zz_c,xy_c,yz_c,zx_c,x_c,y_c,z_c,c_c;
        double x0, y0, z0;
        double a, b, c;

        double coe1,coe2,coe3,coe4,coe5,coe6,rhs;


        //把给定的15个参数转换成椭球一般方程的10个参数
        void translate_para(double     x0, double     y0, double     z0,
                            double      a, double      b, double      c,
                            double alpha1, double alpha2, double alpha3,
                            double  beta1, double  beta2, double  beta3,
                            double  gama1, double  gama2, double  gama3 );

        //解一元二次方程，无解返回0,结果保存在x（两个值）
        int sol_2order_equ(double a, double b, double c, double* x);

        //解一个三阶线性方程组（A*x=B),结果放在x中
        int sol_3order_lina_equs(double a[3][3], double b[3], double x[3]);

        //计算3阶行列式的值（求解三阶线性方程组时使用）
        double det_3order(double a[3][3]);         

        //计算Fx(x,y,z)
        double fxvalue(double x, double y, double z);
        //计算Fy(x,y,z)
        double fyvalue(double x, double y, double z);
        //计算Fz(x,y,z)
        double fzvalue(double x, double y, double z);
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
