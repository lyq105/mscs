
#ifndef ELLIPSEGEN_H
#define ELLIPSEGEN_H

#include <strstream>
#include <fstream>
#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "time.h"
#include <vector>
#include "Hns.h"

using namespace hns;

class EllipseGen
{
	public:
		int uniell_generation(string input_file,string out_file, string data_file);
		double realrate;
	private: 
		int no_uniform(int distribution_sign2,ifstream& input_stream,string out_file, string data_file);
		int uniform(ifstream& input_stream, string out_file);
		int compact_elli(ifstream& input_stream,int distribution_sign1);
		void firstsettley(struct parameter0 **p,int j,int i);
		void firstsettlez(struct parameter0 ***p,int k,int j,int i);
		double secondsettlez(struct parameter0 ***p0,int k,int j,int i);
		double dminz(struct parameter0 *p1,struct parameter0 *p2);
		int judgez(struct parameter0 *,struct parameter0 *,double ,double );
		double disz(struct parameter0 *,struct parameter0 *,double ,double );
		double secondsettley(struct parameter0 ***p0,int k,int j,int i);
		double dminy(struct parameter0 *p1,struct parameter0 *p2);
		int judgey(struct parameter0 *,struct parameter0 *,double ,double );
		double disy(struct parameter0 *,struct parameter0 *,double ,double );
		double secondsettlex(struct parameter0 ***p0,int k,int j,int i);
		double dminx(struct parameter0 *p1,struct parameter0 *p2);
		int judgex(struct parameter0 *,struct parameter0 *,double ,double );
		double disx(struct parameter0 *,struct parameter0 *,double ,double );
		int insertz(struct parameter0 *p0,struct parameter0 **p);
		int inserty(struct parameter0 *p0,struct parameter0 *p);
		void replace(struct parameter0 *p1,struct parameter0 *p2);
		void parametergeneration1(struct parameter0 *p1,int f_e_c_division,int distribution_sign1);
		void parametergeneration2(struct parameter0 *p0);
		double recursion(double *iseed);
		double unifrnd(double c,double d,double *iseed);//均匀分布发生器
		double gasdev(double *iseed);//标准正态分布发生器
		double min(double x0,double x1);
		double max(double x0,double x1);
		char* get_line(ifstream& input_stream, char* rline);//用来读取输入数据
		double volume(double x1,double y1,double z1,struct parameter *p0);//计算椭球的体积
		double normfunction(double t0,double t1,double x);//正态分布函数
		double expfunction(double t,double x);//指数分布函数
		double volume(double x1,double y1,double z1,struct basicparameter *p0);
		double neededrate;
		double seed();  //生成随机数的初值
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
