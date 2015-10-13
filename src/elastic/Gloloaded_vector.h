//============================================================
//�ܺ���������ͷ�ļ���
//Gloloaded_vector.h
//=============================================================

#ifndef GLOLOADED_VECTOR_H
#define GLOLOADED_VECTOR_H

#include<iostream>
#include<vector>
#include<cmath>

#include"Fem.h"
#include"MatPro.h"
#include"time.h"
#include"Gauss.h"
using namespace std;
//----------------------------------------------------------

class Gloloaded_vector{
public:
	//���캯����
	Gloloaded_vector(){};
	//��Ա������
	//����a1�ܺ�������;
	int Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &m,double* &equright);
	//����a1a2�ܺ�������;
    int Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &alfa2,const int &m,double * &equright,const vector<vector<double> > &Na1_vec,const double Homo_D[][6]);
	//����a1��Ԫ����������
	int Generate_eleloaded(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//����a1a2��Ԫ����������
    int Generate_eleloaded1(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const vector<vector<double> > &Na1_vec,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//��Ԫ����������ӵ��ܺ���������
	void Insert_eletoglo(const Element &e,double* &equright);
	//������Ϊ���Բ�ֵʱ����a1��Ԫ����������
	void Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const int &cn);
	//����������a1a2��Ԫ����������
	void Tetrahedron_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
    //����������a1��Ԫ����������
	void Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//����������a1a2��Ԫ����������
	void Threeprism_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//�Ľ�����aijhk�ĺ����������뵯�Ծ���D�кŵĶ�Ӧ��ϵ��
	int Mapping(const int &i,const int &j)const;
	vector<double> eleloaded_vec;
	double homoelas[6][6];
protected:

};//----------------------------------------------------------

#endif
//========================================================================

