//===========================================
//Matbase.cpp
//���Ͽ���Ķ���
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
	//���ݱ���
	vector<MatPro> mats_vec;
	//���캯��
	Matbase(){};
	//��Ա����
	int Generate_matbase(ifstream &infile);
protected:
	//-------------------------------------------------------
	//��Ա����
	//������Ϣһ�У�����ע���У���%��ͷ��
	string Get_Line(ifstream &infile)const;
};//=====================================================
#endif
//-------------------------------------------------------




