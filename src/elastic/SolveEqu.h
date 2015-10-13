//===========================================================================
// SolveEqu.h
// 解线性方程组头文件
//===========================================================================
#ifndef SOVLEEQU_H
#define SOVLEEQU_H

#include<assert.h>
#include<iostream>
#include<cmath>
#include"Fem.h"
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
//定义解方程类
class SolveEqu
{
public:
	//紧缩存储刚度矩阵（动态申请空间+izig函数）
	int Gen_izig(const vector<Node> &, const vector<Element> &, int* &, int* &, const int &);

	// 处理位移边界条件函数
	void tacp(	int M14, int N0, int *IZ, int *IG, int *IP, double *VP, double *F, double *AK );

	// 求解线性方程组
	void sol( double EX, int M14, int N, int *IZ, int *IG, int *IP, double *VP, double *A, 
		double *B, double *P, double *R, double *S, double *V, double *X );

private:
	//紧缩存储刚度矩阵
	void izig(	const int &, const int &, int* &, int* &, int* &, int* &,const int & );

	//在函数izig中使用
	void famnss( int INOD, int *IEN, int *IEM, int *IA, int *KN, int NEL);

	//在函数tacp中被调用
	void tocp2( int N0, int *IZ, int *IG, double *F, double *AK );

	//在函数sol中被调用
	void tcnoz( int KG, int KG1, int M14, int N, int *IP, double *VP, double *W, double *G );

	//在函数sol中被调用
	void mabvm( int , int , int *, int *, double *, double *, double * );

	int Mt,LKt;
	int Jt1,Jt2,Jt3,Jt4,Jt5,Jt6,Jt7,Jt8,Jt9;

};
//---------------------------------------------------------------------------
#endif
//===========================================================================
