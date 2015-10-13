//===========================================================================
// EllipseMade.h
// ����������������ͷ�ļ����ų�����
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
	//��Ա����
	EllipseMade() {};
	void ellip_generation(ifstream &infile,string output_file,string data_file, double mod=0);	//mod=0��ʾ����tecplotͼ
private:
	//-----------------------------------------------------------------
	//���ݳ�Ա
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
	//��Ա����
	string Get_Line(ifstream &)const;																//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
	double ellipse_volume(const struct elliparam &,const double&);					//���������ڵ����е����
	void draw_tecplot(const string&,const double&)const;								//
	void output_EllipseData(const string &)const;
	void output_Datafile(const string data_file)const;										//���뵽���������ļ�
};


//---------------------------------------------------------------------------
#endif
//===========================================================================
