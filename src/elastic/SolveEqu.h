//===========================================================================
// SolveEqu.h
// �����Է�����ͷ�ļ�
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
//����ⷽ����
class SolveEqu
{
public:
	//�����洢�նȾ��󣨶�̬����ռ�+izig������
	int Gen_izig(const vector<Node> &, const vector<Element> &, int* &, int* &, const int &);

	// ����λ�Ʊ߽���������
	void tacp(	int M14, int N0, int *IZ, int *IG, int *IP, double *VP, double *F, double *AK );

	// ������Է�����
	void sol( double EX, int M14, int N, int *IZ, int *IG, int *IP, double *VP, double *A, 
		double *B, double *P, double *R, double *S, double *V, double *X );

private:
	//�����洢�նȾ���
	void izig(	const int &, const int &, int* &, int* &, int* &, int* &,const int & );

	//�ں���izig��ʹ��
	void famnss( int INOD, int *IEN, int *IEM, int *IA, int *KN, int NEL);

	//�ں���tacp�б�����
	void tocp2( int N0, int *IZ, int *IG, double *F, double *AK );

	//�ں���sol�б�����
	void tcnoz( int KG, int KG1, int M14, int N, int *IP, double *VP, double *W, double *G );

	//�ں���sol�б�����
	void mabvm( int , int , int *, int *, double *, double *, double * );

	int Mt,LKt;
	int Jt1,Jt2,Jt3,Jt4,Jt5,Jt6,Jt7,Jt8,Jt9;

};
//---------------------------------------------------------------------------
#endif
//===========================================================================
