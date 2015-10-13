//===========================================
//MatPro.h
//材料属性类头文件
//A Class of Material Property
//===========================================

#ifndef MATPRO_H
#define MATPRO_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include"MathMatrix.h"
using namespace std;
//--------------------------------------------
class MatPro{
public:	
	int type_val;   //材料类型，目前定义3种，
	//					0：各向同性，
	//					1：横观各向同性
	//					2：正交各向异性
	double E_val,Mu_val;		//各向同性材料的弹性模量和泊松比
	double cte1,cte2;				//热膨胀系数；Coefficient of Thermal Expansion 
	double E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13;//正交各向异性的弹性模量，泊松比，剪切模量；
	double E1,Nu1,E2,Nu2,G2;	//横观各向同性的材料的弹性模量，泊松比，剪切模量；1为轴向，2为横观
	double elas_matrix[6][6];		//Dij 弹性矩阵
	double CTE_matrix[3][3];		//热膨胀系数矩阵 
	double aijhk[3][3][3][3];			//Dij 对应的四阶张量
	//韩非20070103增改
	double strength;						//强度极限（各向同性）
	double lon_str, tra_str;				//轴向和横向强度极限（横观各向同性）
	double str_11, str_22, str_33;	//三个方向的强度极限（正交异性）


	//构造函数；
	MatPro(int itype=0);
	//---------------------------------------------
	//成员函数；
	//读入文档数据确定各向同性材料弹性模量和泊松比
	void set_ela_para(double E,double Mu);
	//读入文档数据确定cte1,cte2;
	void set_expan_para(double a1,double a2);
	//读入文档数据确定正交各向异性材料的E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13；
	void set_ela_para(double iE11,double iE22,double iE33,double iNu12,double iNu23,double iNu13,double iG12,double iG23,double iG13);
	//读入文档数据确定横观各向同性材料的
	void set_ela_para(double iE1,double iNu1,double iE2,double iNu2,double iG2);
	//生成Dij弹性矩阵，此弹性矩阵对应的应变向量为（e11,e22,e33,2*e12,2*e23,2*e31）转置
	int Generate_elas_matrix();
	//生成热膨胀系数矩阵；
	int Generate_CTE_matrix();
	//由Dij弹性矩阵，映射与之对应的aijhk;
	int Generate_aijhk();
	//输入aijhk的四个下标，得到其所对应弹性矩阵Dij中的那个值；
	double Trans_aijhk2Dij(int ,int ,int ,int );
	//输出弹性矩阵Dij，aijhk四阶张量；
	//韩非20070103增改
    void set_str_para(double s){ strength = s;}	//设置强度参数（各向同性） 
    void set_str_para(double ls, double ts)			//设置强度参数（横观各向同性）
	{
		lon_str = ls;
		tra_str = ts;
    }                                         
    void set_str_para(double s1, double s2, double s3)	//设置强度参数（正交异性）
	{
		str_11 = s1 ;
		str_22 = s2 ;
		str_33 = s3 ;           
	}
	void print()const;
protected:
	//Dij和aijhk中下标对应关系；
	//输入参数为aijhk前两个或后两个下标，输出参数是Dij的行下标或列下标；
	int Mapping(int i,int j);
};//------------------------------------------------

#endif
//==============================================

