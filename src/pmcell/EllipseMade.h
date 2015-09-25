//===========================================================================
// EllipseMade.h
// 椭球颗粒随机生成类头文件（排除法）
// A class of ellipse mading with random
//===========================================================================
#ifndef ELLIPSEMADE_H
#define ELLIPSEMADE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <vector>
#include "Hns.h"
using namespace hns;

class EllipseMade
{
public:
	//-----------------------------------------------------------------
	//成员函数
	EllipseMade() {};
	void ellip_generation(ifstream &infile,string output_file,string data_file, double mod=0);	//mod=0表示不画tecplot图


	// liyq add @ 2015-09-23  
	void ellip_generation(const string input,string output_file,string data_file, double mod=0);	//mod=0表示不画tecplot图
	void ellip_generation(double* input_value,string output_file,string data_file, double mod=0);	//mod=0表示不画tecplot图
private:
	//-----------------------------------------------------------------
	//数据成员
	struct elliparam		
		{	
			double x;
			double y;
			double z;
			double a;
			double b;
			double c;
			double alpha1;				// the all cos value in the following 
			double alpha2;
			double alpha3;
			double beta1;
			double beta2;
			double beta3;
			double gamma1;
			double gamma2;
			double gamma3;
		};
	struct	space_valu		
		{	
			double x0;
			double y0;
			double z0;
			double	x;
			double	y;
			double	z;
		};

	vector<struct elliparam> ellip;
	struct space_valu space;
	double ellip_ratio;
	//-----------------------------------------------------------------
	//成员函数
	string Get_Line(ifstream &)const;																//读取输入文件一行，跳过注释行（以%开始）
	double ellipse_volume(const struct elliparam &,const double&);					//计算椭球在单胞中的体积
	void draw_tecplot(const string&,const double&)const;								//
	void output_EllipseData(const string &)const;
	void output_Datafile(const string data_file)const;										//输入到样本数据文件
};


//---------------------------------------------------------------------------
#endif
//===========================================================================
