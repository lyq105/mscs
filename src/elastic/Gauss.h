//=========================================
//高斯类
//Gauss.h
//=========================================


#ifndef GAUSS_H
#define GAUSS_H

#include<iostream>
#include<vector>
#include<cmath>

#include"Fem.h"
#include"Vector2D.h"
using namespace std;
//----------------------------------------------------------

class Gauss{
public:
	//数据变量；
	int precision;
	vector<Node> gauss; 
    vector<double> wight;
    //构造函数；
	Gauss(){precision=2;};
	//成员函数；
int Generate_gauss(int type=3);
};//-----------------------------------------------------
#endif
//==========================================================

