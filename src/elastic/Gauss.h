//=========================================
//��˹��
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
	//���ݱ�����
	int precision;
	vector<Node> gauss; 
    vector<double> wight;
    //���캯����
	Gauss(){precision=2;};
	//��Ա������
int Generate_gauss(int type=3);
};//-----------------------------------------------------
#endif
//==========================================================

