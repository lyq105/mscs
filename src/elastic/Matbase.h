//===========================================
//Matbase.cpp
//材料库类的定义
//A Class of Material database
//===========================================

#ifndef MATBASE_H
#define MATBASE_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include"MatPro.h"
#include "Hns.h"
#include "stdlib.h"
using namespace hns;
//-----------------------------------------------------

class Matbase{
public:
	//数据变量
	vector<MatPro> mats_vec;
	//构造函数
	Matbase(){};
	//成员函数
	int Generate_matbase(ifstream &infile);
protected:
	//-------------------------------------------------------
	//成员函数
	//读入信息一行，跳过注释行（以%开头）
	string Get_Line(ifstream &infile)const;
};//=====================================================
#endif
//-------------------------------------------------------




