//===========================================================================
// Geometry.h
// 几何点、线、面、坐标系和几何算法等类头文件
// Classes of Point,line，surface and coordinates and so on
//===========================================================================

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <stdlib.h>
#include "Hns.h"
using namespace hns;

enum Curve_Type { line=1, arc , curve };
//---------------------------------------------------------------------------
//定义点
class Point
{
	public:
        Point operator+( Point &pt );
        Point operator+( double d );
        Point operator-( double d );
        Point operator-( Point &pt );
        Point operator*( double d );
        Point operator/( double d );
        double distance_to(Point &pt);

		//-------------------------------------------------------------------
		double x, y, z;
};

//---------------------------------------------------------------------------
//定义三维向量
class TDVector 
{
	public:
		TDVector( double xx=0, double yy=0, double zz=0 );		//构造函数
        TDVector( Point &p1, Point &p2 );

		double dot_product( TDVector *tdvec );					//求向量点积
		TDVector  cro_product( TDVector *tdvec );				//求向量叉积
		double angle_between(  TDVector *tdvec );				//两向量夹角（0～PI,弧度）

		TDVector& unitize();									//向量单位化

		TDVector operator+ ( TDVector tdv );					//向量运算
        TDVector operator+ ( double d );
		TDVector operator- ( TDVector tdv );
        TDVector operator- ( double d );
		TDVector operator* ( TDVector tdv );
		TDVector operator/ ( TDVector tdv );
		TDVector operator* ( double m );
		TDVector operator/ ( double d );

		bool operator== ( TDVector tdv ) ;

		TDVector vertical_vec();								//返回一个与本向量垂直的向量
     
		double length();										//返回向量的长度值

		//-------------------------------------------------------------------
		double x,y,z;

	private:
     //   bool operator!= ( TDVector tdv ) ;
};

//---------------------------------------------------------------------------
//定义坐标系
class Coors
{
	public:
        Coors();                                      //构造函数
        Coors(Point &p);
		Coors(TDVector &v1, TDVector &v2, TDVector &v3);
        Coors(Point &p, TDVector &v1, TDVector &v2, TDVector &v3);

		void print() ;								  //输出

		//-------------------------------------------------------------------
		Point ori;                                    //坐标原点
        TDVector vec_x,vec_y,vec_z;                   //三个坐标轴向量
 
};  

//---------------------------------------------------------------------------
//定义曲线
class Curve 
{
	public:
		//int point1_num, point2_num;
        Point point1, point2;
        Curve(Point p1, Point p2){ point1=p1; point2=p2; };
		//Curve(int pn1, int pn2){ point1_num=pn1; point2_num=pn2; };
        virtual double xx(double tt){ return 0.0; };
        virtual double yy(double tt){ return 0.0; };
        virtual double zz(double tt){ return 0.0; };
        virtual double xx_(double tt){ return 0.0; };  //the derivative
        virtual double yy_(double tt){ return 0.0; };
        virtual double zz_(double tt){ return 0.0; };
        virtual double length(){return 0.0; };
        virtual Curve_Type curve_type() const { return Curve_Type(0); };
        virtual double nearest_point(Point*)=0;   //曲线上距离给定点最近的点
};

//---------------------------------------------------------------------------
class Line : public Curve 
{
	public:
		//int point1,point2;
        Line (Point p1, Point p2);
        double xx(double tt)
		{
			return tt*A+point1.x;
        };
        double yy(double tt)
		{
			return tt*B+point1.y;
        };
        double zz(double tt)
		{
			return tt*C+point1.z;
        };
        double xx_(double tt)
		{
			return A;
        };
        double yy_(double tt)
		{
			return B;
        };
        double zz_(double tt)
		{
			return C;
        };
        double length()
		{
			return sqrt((point2.x-point1.x)*(point2.x-point1.x)+
						(point2.y-point1.y)*(point2.y-point1.y)+
						(point2.z-point1.z)*(point2.z-point1.z));
        };
        Curve_Type curve_type () const { return line; };
        double nearest_point(Point*);         //直线上距离给定点最近的点
protected:
        double A, B, C ;       //直线方向
};

//---------------------------------------------------------------------------
class Arc : public Curve
{
	public:
        //In this case, point1 means the start point, point2 means the end point
        //the point3 means center point.
        Point point3;
        Arc(Point p1, Point p2, Point p3) : Curve(p1,p2) { point3=p3 ; };
        Curve_Type curve_type () const { return arc; };
};

//---------------------------------------------------------------------------
//定义立方体
class Cubic 
{
	public:
        double x_min,x_max,y_min,y_max,z_min,z_max;
        Cubic(Point* p1, Point* p2)
		{
			x_min = (p1->x<p2->x) ? p1->x : p2->x ;
			x_max = (p1->x>p2->x) ? p1->x : p2->x ;
			y_min = (p1->y<p2->y) ? p1->y : p2->y ;
			y_max = (p1->y>p2->y) ? p1->y : p2->y ;
			z_min = (p1->z<p2->z) ? p1->z : p2->z ;
			z_max = (p1->z>p2->z) ? p1->z : p2->z ;
        };
        Cubic(double xx_min,double xx_max,double yy_min,double yy_max,double zz_min,double zz_max)
		{
			x_min = xx_min;
			x_max = xx_max;
			y_min = yy_min;
			y_max = yy_max;
			z_min = zz_min;
			z_max = zz_max;
        };
};

//---------------------------------------------------------------------------
//定义面
class Surface 
{
	public:
        virtual double length_u(){ return 0.0; }   //u向长度
        virtual double length_v(){ return 0.0; }
        virtual int is_contain(Point *point){ return -2;}
        virtual Point pp(double u, double v);
        virtual int boundary_p(double u, double v, Point* thepoint, Cubic* cb) ;
        virtual Point project_normal(Point* thepoint,int *error=NULL) ;
        virtual Point project(Point* thepoint, TDVector* vec, int *error=NULL, int mod=0) ;
        virtual int intersect(Point &point1, Point &point2, Point &rpoint) = 0 ; //求交点
        virtual int intersect(Line &line, Point &rpoint) = 0 ;
        virtual double dis_to(Point* thepoint){ return 0.0; }
        virtual double fvalue(double x, double y, double z){ return 0.0; }

		//根据一定比例紧缩椭球，对椭球上每个点做坐标变换
		virtual int change_coor(Point *opoint, Point *point, double ratio){ return 0; };

        virtual void print() const{} ;

        int material_num ;
        Curve *sscurve;
};

//---------------------------------------------------------------------------
//定义实体 
struct Solid 
{
	double x_min,x_max,y_min,y_max,z_min,z_max;
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
