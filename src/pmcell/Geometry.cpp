//===========================================================================
// Geometry.cpp
// 几何点、线、面、坐标系和几何算法等类成员函数
// Member Functions in classes of Point,line，surface and coordinates and so on
//===========================================================================
#include "Geometry.h"

//向量运算类TDVector
//---------------------------------------------------------------------------
//构造函数
TDVector::TDVector( double xx, double yy, double zz )
{
	x = xx;
	y = yy;
	z = zz;
}
//---------------------------------------------------------------------------
TDVector::TDVector( Point &p1, Point &p2 )
{
	x = p2.x - p1.x;
	y = p2.y - p1.y;
	z = p2.z - p1.z;
}
//---------------------------------------------------------------------------
//向量点乘
double TDVector::dot_product( TDVector *tdvec )
{
	return x * tdvec->x + y * tdvec->y + z * tdvec->z ;
}
//---------------------------------------------------------------------------
//向量叉乘
TDVector TDVector::cro_product( TDVector *tdvec )
{
	double xx = y * tdvec->z - z * tdvec->y ;
	double yy = z * tdvec->x - x * tdvec->z ;
	double zz = x * tdvec->y - y * tdvec->x ;

	return TDVector(xx,yy,zz);
}
//---------------------------------------------------------------------------
//两个向量之间的夹角
double TDVector::angle_between( TDVector *tdvec )
{
	double len1 = length();
	if( len1 == 0 )
	{
		return 0.0;
	}

	double len2 = tdvec->length();
	if( len2 == 0 ) 
	{
		return 0.0;
	}

	double acos_v = dot_product(tdvec)/(len1 * len2);

	if( acos_v > 1.0 ) 
	{
		acos_v = 1.0;
	}
	else if( acos_v < -1.0 ) 
	{
		acos_v = -1.0;
	}

	return acos(acos_v);
}
//---------------------------------------------------------------------------
//向量单位化
TDVector& TDVector::unitize()
{
	double len = length();
	x = x/len;
	y = y/len;
	z = z/len;

	return *this;
}
//---------------------------------------------------------------------------
//返回向量的长度值
double TDVector::length()
{
	return sqrt(x * x + y * y + z * z );
}

//---------------------------------------------------------------------------
//向量运算
TDVector TDVector::operator+ ( TDVector tdv )						//向量加法
{
	return TDVector( x + tdv.x, y + tdv.y, z + tdv.z );
}

TDVector TDVector::operator+ ( double d )							//向量加法(同加一个数)
{
	return TDVector( x + d, y + d, z + d );
}

TDVector TDVector::operator- ( TDVector tdv )						//向量减法		
{
	return TDVector( x - tdv.x, y - tdv.y, z - tdv.z );
}

TDVector TDVector::operator- ( double d )							//向量减法(同减一个数)
{
	return TDVector( x - d, y - d, z - d );
}

TDVector TDVector::operator* ( TDVector tdv )						//向量乘法
{
	return TDVector( x * tdv.x, y * tdv.y, z * tdv.z );
}

TDVector TDVector::operator/ ( TDVector tdv )						//向量除法
{
	return TDVector( x / tdv.x, y / tdv.y, z / tdv.z );
}

TDVector TDVector::operator* ( double m )							//向量乘法（放大m倍）
{
	return TDVector( x * m, y * m, z * m );
}

TDVector TDVector::operator/ ( double d )							//向量除法（缩小d倍）
{
	return TDVector( x / d, y / d, z / d );
}

bool TDVector::operator== ( TDVector tdv )				//向量相等
{
	if( fabs(tdv.x) + fabs(tdv.y) + fabs(tdv.z) == 0)
	{
		return fabs(x) + fabs(y) + fabs(z) == 0 ;
	}
	if( fabs(x) + fabs(y) + fabs(z) == 0 )
	{
		return fabs(tdv.x) + fabs(tdv.y) + fabs(tdv.z) == 0 ;
	}

	double t;
	if( x != 0 ) t = tdv.x / x;
	else if( y != 0 ) t = tdv.y /y;
	else if( z != 0 ) t = tdv.z /z;
	return tdv.x == t * x && tdv.y == t * y && tdv.z == t * z;
}

//---------------------------------------------------------------------------
//返回一个与本向量垂直的向量
TDVector TDVector::vertical_vec()
{
	if( x == 0 ) return TDVector( 1, 0, 0 );
	if( y == 0 ) return TDVector( 0, 1, 0 );
	if( z == 0 ) return TDVector( 0, 0, 1 ); 

	double xx = -x;
	double yy = -y;
	double zz = (x*x+y*y)/z;
	return TDVector(xx, yy, zz);
}

//*************************************************************************//
//点类
//---------------------------------------------------------------------------
Point Point::operator+( Point &pt )
{
	Point rp = { x + pt.x, y + pt.y, z + pt.z };
	return rp;
}
//---------------------------------------------------------------------------
Point Point::operator+( double d )
{
	Point rp = { x + d, y + d, z + d };
	return rp;
}
//---------------------------------------------------------------------------
Point Point::operator-( Point &pt )
{
	Point rp = { x - pt.x, y - pt.y, z - pt.z };
	return rp;
}
//---------------------------------------------------------------------------
Point Point::operator-( double d )
{
	Point rp = { x - d, y - d, z - d };
	return rp;
}
//---------------------------------------------------------------------------
Point Point::operator*( double d )
{
	Point rp = { x*d, y*d, z*d };
	return rp;
}
//---------------------------------------------------------------------------
Point Point::operator/( double d )
{
	Point rp = { x/d, y/d, z/d };
	return rp;
}
//---------------------------------------------------------------------------
double Point::distance_to( Point &pt )
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y)+(z-pt.z)*(z-pt.z);
	return sqrt(rv2);
}

//*************************************************************************//
//坐标系类
//---------------------------------------------------------------------------
//构造函数
Coors::Coors()
{
	ori.x = 0;
	ori.y = 0;
	ori.z = 0;
	vec_x = TDVector(1,0,0);
	vec_y = TDVector(0,1,0);
	vec_z = TDVector(0,0,1);
}
//---------------------------------------------------------------------------
//重载
Coors::Coors(Point &p)
{
	ori = p;
	vec_x = TDVector(1,0,0);
	vec_y = TDVector(0,1,0);
	vec_z = TDVector(0,0,1);
}

//---------------------------------------------------------------------------
//重载
Coors::Coors(TDVector &v1, TDVector &v2, TDVector &v3)
{
	ori.x = 0;
	ori.y = 0;
	ori.z = 0;
	vec_x = v1;
	vec_y = v2;
	vec_z = v3;
}

//---------------------------------------------------------------------------
//重载
Coors::Coors(Point &p, TDVector &v1, TDVector &v2, TDVector &v3)
{
	ori = p;
	vec_x = v1;
	vec_y = v2;
	vec_z = v3;
}

//---------------------------------------------------------------------------

//*************************************************************************//
//线类
//---------------------------------------------------------------------------
Line::Line(Point p1, Point p2) : Curve(p1,p2)
{
	A = p2.x - p1.x ;
	B = p2.y - p1.y ;
	C = p2.z - p1.z ;
}
//---------------------------------------------------------------------------
//点到直线的距离
double Line::nearest_point(Point* point)
{
	return (A*(point->x - point1.x) +
			B*(point->y - point1.y) +
			C*(point->z - point1.z))/
			(A*A + B*B + C*C ) ;
}

//*************************************************************************//
//面类
//---------------------------------------------------------------------------
Point Surface::pp(double u, double v)
{
	Point rp={0,0,0};
	return rp;
}

//---------------------------------------------------------------------------
Point Surface::project_normal(Point* thepoint, int *error)
{
	Point rp={0,0,0};
	return rp;
}
//---------------------------------------------------------------------------
Point Surface::project(Point* thepoint, TDVector* vec, int *error, int mod)
{
	Point rp={0,0,0};
	return rp;
}
//---------------------------------------------------------------------------
int Surface::boundary_p(double u, double v, Point* thepoint, Cubic* cb)
{
	return 0;
}

//===========================================================================
