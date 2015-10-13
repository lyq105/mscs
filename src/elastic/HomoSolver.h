//===========================================================================
// HomoSolver.h
// 求解均匀化系数类头文件
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
//均匀化参数求解类
class HomoSolver
{
public:
	//--------------------------------------------
	//数据变量
	string data_file;
	int Na1a2_key;
	double Unitcell_V;
	vector<vector<double> > Na1_vec;			//Na1
	vector<vector<double> > Na1a2_vec;		//Na1a2
	double Homo_D[6][6];
	MatPro homoMat;										//均匀化后的材料的等效参数
	//--------------------------------------------
	//成员函数

	//构造函数
	HomoSolver(){};
	HomoSolver(int elas_ana_only, string datafile){ data_file = datafile; Na1a2_key = elas_ana_only; };
	//求解均匀化参数;
	int Solve(const vector<Node> &nodes_vec,const vector<int> &bnodes_vec, const vector<Element> &elements_vec, const vector<MatPro> &mats_vec, const double unitcellV, int mod=0);
	//Na1和Na1a2的二进制数据的输出和读取
	int Na_BinaryData(int mod, const vector<Node> &nodes_vec, string data_file, int CNum);			
protected:
	//--------------------------------------------
	//计算Na1系数
	int Cal_Na1(double* &, int* &, int* &, const vector<Node> &, const vector<int> &, const vector<Element> &, const vector<MatPro> &);
	//计算Na1a2系数
	int Cal_Na1a2(double* &, int* &, int* &, const vector<Node> &, const vector<int> &, const vector<Element> &, const vector<MatPro> &);
};
//---------------------------------------------------------------------------
#endif
//====================================
