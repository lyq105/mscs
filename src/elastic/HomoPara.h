//=======================================================
//������Ȼ�ϵ��(���Ծ���������ϵ��)
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
	//���ݱ���
	//���Ȼ����Ծ���
	double Homo_D[6][6];
	//���Ȼ������;���
	double Homo_alpha[3][3];
	//���Ȼ��ȵ��Ծ���
	double Homo_beta[6];
	//���캯��
	HomoPara();
	//��Ա����
	//���ɾ��Ȼ�����ģ������
	int Generate_Homo(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const double& Unitcell_V,const int &a1,const int &m);
	//����alpha���󣻾��Ȼ�������ϵ��������
	int Generate_Homo_alpha(const double &Unitcell_V);
	//���ɺ�۸���ͬ�Ե��Թ��̳�����
	int Generate_Homo_engconst(MatPro* homoMat); 
	//�����壻
    void Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const int &cn);
	//������;
    void Threeprism_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//���������
	void print(int mode = 2);
	void output_Datafile(const string data_file)const;
	double ele_Homo_Dcn[6];
protected:
	//��Ԫ�ȵ��Գ���
 double ele_Homo_beta[6];
    //��ԪD����
 double ele_Homo_D[6][6];
	//��ԪD������ijԪ��
   void Generate_ele_Dmatrixcn(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//��Ԫ�ȵ��Գ�����
	void Generate_ele_Homo_beta(const Element &e,const vector<MatPro> &mats_vec,const int &cn);
	//ӳ�䣻
	int Mapping(int i,int j);
};//----------------------------------------------------------------------------------------------------------------------------------------------------------
#endif






