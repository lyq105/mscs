//===========================================================================
// EllipseSurface.cpp
// 椭球面类成员函数
// Member Functions in a class of Ellipse Surface
//===========================================================================
#include "EllipseSurface.h"
#define ZERO 1.0e-12

//---------------------------------------------------------------------------
//构造函数，已知椭球面一般方程的情况(此构造函数不完整，小心使用）
EllipseSurface::EllipseSurface(	double xx, double yy, double zz, double xy, double yz,
												double zx, double x,  double y,  double z,  double c	)
{
	xx_c=xx;
	yy_c=yy;
	zz_c=zz;
	xy_c=xy;
	yz_c=yz;
	zx_c=zx;
	x_c =x;
	y_c =y;
	z_c =z;
	c_c =c;
}
//---------------------------------------------------------------------------
//构造函数
EllipseSurface::EllipseSurface(	double     x0, double     y0, double     z0,
												double      a, double      b, double      c,
												double alpha1, double alpha2, double alpha3,
												double  beta1, double  beta2, double  beta3,
												double  gama1, double  gama2, double  gama3	)
{
	this->x0 = x0;
	this->y0 = y0;
	this->z0 = z0;
	this->a  = a ;
	this->b  = b ;
	this->c  = c ;

	translate_para(	x0,y0,z0,a,b,c,alpha1,alpha2,alpha3,
							beta1,beta2,beta3,gama1,gama2,gama3	);
}
//---------------------------------------------------------------------------
//计算给定坐标（x,y,z）代入椭球面方程的值
double EllipseSurface::fvalue(double x, double y, double z)
{
	double xx = x-x0;
	double yy = y-y0;
	double zz = z-z0;
	double ff = coe1*xx*xx+coe2*yy*yy+coe3*zz*zz+2.0*coe4*xx*yy+2.0*coe5*yy*zz+
				2.0*coe6*zz*xx-rhs;
//	hout << "ff=" << ff << coe1 <<" "  << coe2 << " "  << coe3 << " "  << coe4 << " "  << coe5 << " "  << coe6 << " "  << rhs <<" "  << endl;
	return ff;
}
//---------------------------------------------------------------------------
//计算Fx(x,y,z)
double EllipseSurface::fxvalue(double x, double y, double z)
{
	double xx = x-x0;
	double yy = y-y0;
	double zz = z-z0;
	double fx = (coe1*xx+coe4*yy+coe6*zz)*2.0;
	return fx;
}
//---------------------------------------------------------------------------
//计算Fy(x,y,z)
double EllipseSurface::fyvalue(double x, double y, double z)
{
	double xx = x-x0;
	double yy = y-y0;
	double zz = z-z0;
	double fy = (coe4*xx+coe2*yy+coe5*zz)*2.0;
	return fy;
}
//---------------------------------------------------------------------------
//计算Fz(x,y,z)
double EllipseSurface::fzvalue(double x, double y, double z)
{
	double xx = x-x0;
	double yy = y-y0;
	double zz = z-z0;
	double fz = (coe6*xx+coe5*yy+coe3*zz)*2.0;
	return fz;
}
//---------------------------------------------------------------------------
//计算3阶行列式的值（求解三阶线性方程组时使用）
double EllipseSurface::det_3order(double a[3][3])
{
	double s1 = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]);
	double s2 = a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]);
	double s3 = a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
	return s1-s2+s3;
}
//---------------------------------------------------------------------------
//给定点到椭球面的距离
double EllipseSurface::dis_to(Point *point)
{
	Point rp = project_normal(point);
	return point->distance_to(rp);
}
//---------------------------------------------------------------------------
//点和椭球面的位置关系（使用椭球面一般方程）
int EllipseSurface::is_contain_usually( Point* thepoint, int mod )
{
	double xx=thepoint->x;
	double yy=thepoint->y;
	double zz=thepoint->z;
	double ff=xx_c*xx*xx+yy_c*yy*yy+zz_c*zz*zz+xy_c*xx*yy+yz_c*yy*zz+
			  zx_c*zz*xx+x_c*xx+y_c*yy+z_c*zz+c_c;
	if( mod == 1)
	{
		hout << " ff: " << ff << endl;
	}   
    if( ff > ZERO ) return 1;
	else if( ff < -ZERO ) return -1;
	else return 0;
}
//---------------------------------------------------------------------------
//点和椭球面的位置关系
int EllipseSurface::is_contain( Point* thepoint )
{
    double ff = fvalue(thepoint->x,thepoint->y,thepoint->z);
	if( ff > ZERO ) return 1;
	else if( ff < -ZERO ) return -1;
	else return 0;
}
//---------------------------------------------------------------------------
//寻找线段与椭球面的交点
int EllipseSurface::intersect(Point &point1, Point &point2, Point &rpoint)
{
	int flag1 = is_contain( &point1 );
	if ( flag1 == 0 )		 //point1 在面上
	{   
		rpoint = point1;
		return 1;
	}
	int flag2 = is_contain( &point2 );
	if ( flag2 == 0 )		//point2 在面上
	{   
		rpoint = point2;
		return 1;
	}

	if( flag1 * flag2 > 0 ) return 0;   //point1 point2 在面同侧

	//直接投影到椭球面上
	TDVector pvec(point1, point2);
	int error;
	rpoint = project(&point1,&pvec,&error,1);
	return error;                               
}
//---------------------------------------------------------------------------
//寻找线段与椭球面的交点
int EllipseSurface::intersect(Line &line, Point &rpoint)
{
	Point point1 = line.point1;
	Point point2 = line.point2;
	return intersect( point1, point2, rpoint ) ;
}

//---------------------------------------------------------------------------
//求给定点到椭球面的投影（沿法向投影）
Point EllipseSurface::project_normal( Point* thepoint, int *error )
{
	bool debuging = false;

	double xx = thepoint->x;
	double yy = thepoint->y;
	double zz = thepoint->z;
        
	int ic = is_contain(thepoint);
	//自身在椭球面上
	if( ic == 0 )
	{
		if( error != NULL ) *error = 1;
		return *thepoint;
	}

	//计算给定点与椭球中心点连线与椭球面的交点
	Point p0 = {x0, y0, z0};
	//检查p0和thepoint是不是一个点（有可能要投影的点就是球心点）
	TDVector ppv( p0, *thepoint );
	if( ppv.length() <= ZERO )
	{
		p0 = p0 + 1.0;
	}
	Point p1;
	if( ic < 0 )
	{
		//椭球内
		TDVector pvec( p0, *thepoint );
		p1 = project(thepoint, &pvec);
	}
	else
	{
		//椭球外
		TDVector pvec( *thepoint, p0 );
		p1 = project( thepoint, &pvec );
	}

	if( debuging )
	{
		hout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		hout << "The point: " << xx << " " << yy << " " << zz << endl;
	}

	//采用梯度法求解，BB方法
	//初值
	TDVector xv0(xx,yy,zz);
	//初始函数值（三个函数值）
	double f10 = fvalue(xv0.x, xv0.y, xv0.z);
	double f20 = 0.0;
	double f30 = 0.0;
	TDVector fv0(f10,f20,f30);
	//步长初值
	double alpha = 1.0;

	TDVector xv1(p1.x,p1.y,p1.z);

	double f11 = fvalue(xv1.x,xv1.y,xv1.z);
	double f21 = (xv1.x-xx)*fyvalue(xv1.x,xv1.y,xv1.z)-
				 (xv1.y-yy)*fxvalue(xv1.x,xv1.y,xv1.z);
	double f31 = (xv1.y-yy)*fzvalue(xv1.x,xv1.y,xv1.z)-
                 (xv1.z-zz)*fyvalue(xv1.x,xv1.y,xv1.z);
	TDVector fv1(f11,f21,f31);

	if( debuging )
	{
		hout << "xv0: " <<  "     " << xv0.x
								<<  "     " << xv0.y
								<<   "     "  << xv0.z << endl;
		hout << "fv0: " <<  "     " << fv0.x
								<<   "     " << fv0.y
								<<  "     " << fv0.z << endl;
		hout << "xv1: " <<  "     " << xv1.x
								<<  "     " << xv1.y
								<<  "     " << xv1.z << endl;
		hout << "fv1: " <<  "     " << fv1.x
								<<  "     " << fv1.y
								<<  "     " << fv1.z << endl;
	}

	//迭代开始
	//记录两次迭代产生的点的最小距离，防止发散
	//（当距离增大到100倍时停止，因为梯度法本身也不是单调的）
	double dis_min = 1E200;
	int is_failed = 0;
	int cl_num=0;
	TDVector dxv = xv1-xv0;
	while( (cl_num == 0 || dxv.length() > ZERO) && cl_num++ < 1500 )
	{
		if( debuging ) hout << "----------------------------------------------------" << endl;
		if( debuging ) hout << "cl_num : " << cl_num << endl;

		TDVector yk = fv1 - fv0;
		TDVector sk = xv1 - xv0;
		alpha = sk.dot_product(&yk)/(yk.x*yk.x+yk.y*yk.y+yk.z*yk.z);

		TDVector xv2 = fv1*alpha;
	    xv2 = xv1-xv2;
		dxv = xv2-xv1;

		double f12 = fvalue(xv2.x,xv2.y,xv2.z);
		double f22 = (xv2.x-xx)*fyvalue(xv2.x,xv2.y,xv2.z)-
					 (xv2.y-yy)*fxvalue(xv2.x,xv2.y,xv2.z);
		double f32 = (xv2.y-yy)*fzvalue(xv2.x,xv2.y,xv2.z)-
					 (xv2.z-zz)*fyvalue(xv2.x,xv2.y,xv2.z);
		TDVector fv2(f12,f22,f32);

		//更新变量
		xv0 = xv1;
		xv1 = xv2;
		fv0 = fv1;
		fv1 = fv2;

		dxv = xv1 - xv0;
		double dxv_len = dxv.length();
		if( dxv_len < dis_min ) dis_min = dxv_len;
		else if( dxv_len > dis_min*1000 )
		{
			//发散了，不幸啊
			if( debuging ) hout << "Failed. dis: " <<dxv_len <<" dis_min:" <<dis_min << endl << endl;
			is_failed = 1;
			break;
		}

		if( debuging )
		{
			hout << "alpha: " <<  "     " << alpha << endl;
			hout << "xv0: " <<  "     " << xv0.x
									<<  "     " << xv0.y
									<<  "     " << xv0.z << endl;
			hout << "fv0: " <<  "     " << fv0.x
									<<  "     " << fv0.y
									<<  "     " << fv0.z << endl;
			hout << "xv1: " <<  "     " << xv1.x
									<<  "     " << xv1.y
									<<  "     " << xv1.z << endl;
			hout << "fv1: " <<  "     " << fv1.x
									<<  "     " << fv1.y
									<<  "     " << fv1.z << endl;
			hout << "dis: " <<  "     " << dxv.length() << " dis_min: " << dis_min << endl;
		}
	}
	if( is_failed == 0 )
	{
		//迭代后的结果有时候不能满足精度要求，所以再沿点的方向投影一下
		Point rp = {xv1.x, xv1.y, xv1.z};
		TDVector p2svec(*thepoint,rp);
		Point rp1 = project(thepoint,&p2svec);

		if( debuging )
		{
			hout << "Before project: " <<  "        " << xx << " " <<  "        " << yy << " " <<  "        " << zz << endl;
			hout << "Middle point  : " <<  "        " << rp.x << " " <<  "        " << rp.y << " " <<  "        " << rp.z
				 << " fvalue: " << fvalue(rp.x,rp.y,rp.z) << endl;
			hout << "Return point  : " <<  "        " << rp1.x << " " <<  "        " << rp1.y << " " <<  "        " << rp1.z
				 << " fvalue: " << fvalue(rp1.x,rp1.y,rp1.z) << endl;

		}

		if( error != NULL ) *error = 1;
		return rp1;
	}
	else
	{                      
		if( error != NULL ) *error = 0;
		return *thepoint;
	 }
}
//---------------------------------------------------------------------------
//求给定点沿给定方向到椭球面的投影
Point EllipseSurface::project(Point* thepoint, TDVector* vec, int *error, int mod)
{
	//单位化投影向量
	TDVector lvec = *vec;
	if( lvec.length() < 1.0e-8 )
	{
		if( error != NULL ) *error = 0;
		return *thepoint;
	}
	lvec = lvec/lvec.length();

	Point rpoint=*thepoint;
        
	double A=lvec.x;
	double B=lvec.y;
	double C=lvec.z;
	double x1=thepoint->x;
	double y1=thepoint->y;
	double z1=thepoint->z;
	double dx=x1-x0;
	double dy=y1-y0;
	double dz=z1-z0;
	double a=2.0*coe4*A*B+2.0*coe5*B*C+2.0*coe6*C*A+coe1*A*A+coe2*B*B+coe3*C*C;
	double b=2.0*(coe1*dx*A+coe2*dy*B+coe3*dz*C+coe4*(dx*B+dy*A)+coe5*(dy*C+dz*B)+coe6*(dz*A+dx*C));
	double c=coe1*dx*dx+coe2*dy*dy+coe3*dz*dz+2.0*coe4*dx*dy+2.0*coe5*dy*dz+2.0*coe6*dz*dx-rhs;
	double t[2];

	if( sol_2order_equ(a,b,c,t) == 1 )
	{
		//还是用mod来控制吧，mod==0:  返回最近的，mod==1: 先按方向返回
		if( mod == 0 )
		{
			//看来还是要选取距离最近的点，否则光顺处理的时候搞不清方向
			//先这样了，回头再说
			Point point1={A*t[0]+x1,B*t[0]+y1,C*t[0]+z1};
			Point point2={A*t[1]+x1,B*t[1]+y1,C*t[1]+z1};
			double dis1 = thepoint->distance_to(point1);
			double dis2 = thepoint->distance_to(point2);

			if( dis1 < dis2 )
			{
				rpoint = point1;
			}
			else    
			{
				rpoint = point2;
			}
		}
		else if( mod == 1 )
		{
			if( t[0] < 0 )
			{
				Point point2={A*t[1]+x1,B*t[1]+y1,C*t[1]+z1};
				rpoint = point2;
			}
			else
			{
				Point point1={A*t[0]+x1,B*t[0]+y1,C*t[0]+z1};
				if( t[1] < 0 )
				{
					rpoint = point1;
				}
				else
				{
					Point point2={A*t[1]+x1,B*t[1]+y1,C*t[1]+z1};
					//取距离近的点
					if( thepoint->distance_to(point1)<thepoint->distance_to(point2))
					{
						rpoint = point1;
					}
					else    
					{
						rpoint = point2;
					}
				}
			}
		}
		if( error != NULL ) *error = 1;
	}
	else
	{
		if( error != NULL ) *error = 0;
	}

	return rpoint;
}
//---------------------------------------------------------------------------
//解一元二次方程，无解返回0,结果保存在x（两个值）
int EllipseSurface::sol_2order_equ(double a, double b, double c, double* x)
{
	double dd = b*b-4.0*a*c;
	if( dd < 0 ) return 0;
	if( a == 0 )
	{
		if( b == 0 ) return 0;
		x[0] = -c/b;
		x[1] = x[0];
		return 1;
	};
	x[0] = (-b+sqrt(dd))/(2.0*a);  
	x[1] = (-b-sqrt(dd))/(2.0*a);
	return 1;
}
//---------------------------------------------------------------------------
//解一个三阶线性方程组（A*x=B),结果放在x中，卡莱姆法则
int EllipseSurface::sol_3order_lina_equs(double a[3][3], double b[3], double x[3])
{
	double D = det_3order(a);
	if( D == 0 ) return 0;
	double a1[3][3];
	for( int i=0; i<3; i++ )
	{
		for( int j=0; j<3; j++ )
		{
			a1[i][j] = a[i][j];
		}
		a1[i][0]=b[i];
	}
	double D1 = det_3order(a1);
	x[0] = D1/D;

	for( int i=0; i<3; i++ )
	{
		for( int j=0; j<3; j++ )
		{
			a1[i][j] = a[i][j];
		}
		a1[i][1]=b[i];
	}
	double D2 = det_3order(a1);
	x[1] = D2/D;

	for( int i=0; i<3; i++ )
	{
		for( int j=0; j<3; j++ )
		{
			a1[i][j] = a[i][j];
		}
		a1[i][2]=b[i];
	}
	double D3 = det_3order(a1);
	x[2] = D3/D;

	return 1;
}
//---------------------------------------------------------------------------
//把给定的15个参数转换成椭球一般方程的10个参数
void EllipseSurface::translate_para(double     x0, double     y0, double     z0,
													double      a, double      b, double      c,
													double alpha1, double alpha2, double alpha3,
													double  beta1, double  beta2, double  beta3,
													double  gama1, double  gama2, double  gama3 )
{
	double A1 = alpha1,	A2 =	alpha2,	A3 =	alpha3;
	double B1 =	beta1,	B2 =	beta2,	B3 =	beta3;
	double C1 =	gama1,	C2 =	gama2,	C3 =	gama3;

	coe1 = a*a*b*b*A3*A3 + a*a*c*c*A2*A2 + b*b*c*c*A1*A1;
	coe2 = a*a*b*b*B3*B3 + a*a*c*c*B2*B2 + b*b*c*c*B1*B1;
	coe3 = a*a*b*b*C3*C3 + a*a*c*c*C2*C2 + b*b*c*c*C1*C1;
	coe4 = a*a*b*b*A3*B3 + a*a*c*c*A2*B2 + b*b*c*c*A1*B1;
	coe5 = a*a*b*b*C3*B3 + a*a*c*c*C2*B2 + b*b*c*c*C1*B1;
	coe6 = a*a*b*b*A3*C3 + a*a*c*c*A2*C2 + b*b*c*c*A1*C1;
	rhs  = a*a*b*b*c*c;

	//化简，使coe1=1;
	coe2 = coe2/coe1;
	coe3 = coe3/coe1;
	coe4 = coe4/coe1;
	coe5 = coe5/coe1;
	coe6 = coe6/coe1;
	rhs  = rhs /coe1;
	coe1 = 1.0;

	xx_c = coe1;
	yy_c = coe2;
	zz_c = coe3;
	xy_c = 2.0*coe4;
	yz_c = 2.0*coe5;
	zx_c = 2.0*coe6;
	x_c  = -2.0*coe1*x0-2.0*coe4*y0-2.0*coe6*z0;
	x_c  = -2.0*coe4*x0-2.0*coe2*y0-2.0*coe5*z0;
	x_c  = -2.0*coe6*x0-2.0*coe5*y0-2.0*coe3*z0;
	c_c  = coe1*x0*x0+coe2*y0*y0+coe3*z0*z0+2.0*coe4*x0*y0+2.0*coe5*y0*z0+
		   2.0*coe6*z0*x0-rhs;
} 
//---------------------------------------------------------------------------
void EllipseSurface::print( int mod )
{
	if( mod == 0 )
	{
		hout << "coes: " << coe1 << " " << coe2 << " " << coe3 << " "
			 << coe4 << " " << coe5 << " " << coe6 << " " << rhs << endl;
	}
	else
	{
		hout << "10参数:"	<< xx_c << " " << yy_c << " " << zz_c << endl
							<< xy_c << " " << yz_c << " " << zx_c << endl
							<< x_c  << " " << y_c  << " " << z_c  << endl
							<< c_c  << endl;
	}
}						

//---------------------------------------------------------------------------
//根据一定比例紧缩椭球，对椭球上每个点做坐标变换
int EllipseSurface::change_coor(Point *opoint, Point *point, double ratio)
{
	if (ratio>1.0)
	{
		hout << "      椭球收缩比例竟然比1.0还大，错误！" << endl;
		return 0;
	}

	point->x = x0 + (opoint->x - x0)*(1.0-ratio);
	point->y = y0 + (opoint->y - y0)*(1.0-ratio);
	point->z = z0 + (opoint->z - z0)*(1.0-ratio);

	return 1;
}

//===========================================================================
