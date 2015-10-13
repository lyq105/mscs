//=======================================================
//计算均匀化系数(弹性矩阵，热膨胀系数)
//HomoPara.h
//=======================================================

#ifndef HOMOPARA_H
#define HOMOPARA_H

#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

#include"MathMatrix.h"
#include"Fem.h"
#include"MatPro.h"
#include"Gauss.h"
#include "Hns.h"
using namespace hns;
//----------------------------------------------------------


class HomoPara{
public:
	double E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13;
	double CTE1,CTE2;
	//数据变量
	//均匀化弹性矩阵
	double Homo_D[6][6];
	//均匀化热膨胀矩阵
	double Homo_alpha[3][3];
	//均匀化热弹性矩阵
	double Homo_beta[6];
	//构造函数
	HomoPara();
	//成员函数
	//生成均匀化弹性模量参数
	int Generate_Homo(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const double& Unitcell_V,const int &a1,const int &m);
	//生成alpha矩阵；均匀化热膨胀系数张量；
	int Generate_Homo_alpha(const double &Unitcell_V);
	//生成横观各向同性弹性工程常数；
	int Generate_Homo_engconst(MatPro* homoMat); 
	//四面体；
    void Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const int &cn);
	//三棱柱;
    void Threeprism_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//输出函数；
	void print(int mode = 2);
	void output_Datafile(const string data_file)const;
	double ele_Homo_Dcn[6];
protected:
	//单元热弹性常数
 double ele_Homo_beta[6];
    //单元D矩阵；
 double ele_Homo_D[6][6];
	//单元D矩阵中ij元素
   void Generate_ele_Dmatrixcn(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//单元热弹性常数；
	void Generate_ele_Homo_beta(const Element &e,const vector<MatPro> &mats_vec,const int &cn);
	//映射；
	int Mapping(int i,int j);
};//----------------------------------------------------------------------------------------------------------------------------------------------------------
#endif






