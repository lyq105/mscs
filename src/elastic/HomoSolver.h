//===========================================================================
// HomoSolver.h
// �����Ȼ�ϵ����ͷ�ļ�
// A Class of the Homogeneous Solver
//===========================================================================
#ifndef HOMOSOLVER_H
#define HOMOSOLVER_H

#include<vector>
#include"Fem.h"
#include"Gloloaded_vector.h"
#include"GloStiffMatrix.h"
#include"HomoPara.h"
#include"MatPro.h"
#include"SolveEqu.h"
#include"time.h"
using namespace std;

#define NDIM 3
#define IRN 1					         // the number of element group
#define MAXIG 11 
#define NGAP 1
#define Ex 0.00000005

//---------------------------------------------------------------------------
//���Ȼ����������
class HomoSolver
{
public:
	//--------------------------------------------
	//���ݱ���
	string data_file;
	int Na1a2_key;
	double Unitcell_V;
	vector<vector<double> > Na1_vec;			//Na1
	vector<vector<double> > Na1a2_vec;		//Na1a2
	double Homo_D[6][6];
	MatPro homoMat;										//���Ȼ���Ĳ��ϵĵ�Ч����
	//--------------------------------------------
	//��Ա����

	//���캯��
	HomoSolver(){};
	HomoSolver(int elas_ana_only, string datafile){ data_file = datafile; Na1a2_key = elas_ana_only; };
	//�����Ȼ�����;
	int Solve(const vector<Node> &nodes_vec,const vector<int> &bnodes_vec, const vector<Element> &elements_vec, const vector<MatPro> &mats_vec, const double unitcellV, int mod=0);
	//Na1��Na1a2�Ķ��������ݵ�����Ͷ�ȡ
	int Na_BinaryData(int mod, const vector<Node> &nodes_vec, string data_file, int CNum);			
protected:
	//--------------------------------------------
	//����Na1ϵ��
	int Cal_Na1(double* &, int* &, int* &, const vector<Node> &, const vector<int> &, const vector<Element> &, const vector<MatPro> &);
	//����Na1a2ϵ��
	int Cal_Na1a2(double* &, int* &, int* &, const vector<Node> &, const vector<int> &, const vector<Element> &, const vector<MatPro> &);
};
//---------------------------------------------------------------------------
#endif
//====================================
