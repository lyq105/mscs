//===========================================================================
// GloStiffMatrix.h
// 求解总刚阵类头文件
// A Class of the Global Stiff Matrix
//===========================================================================
#ifndef GloStiffMatrixH
#define GloStiffMatrixH

#include"Fem.h"
#include"MatPro.h"
#include"Gauss.h"
#include"time.h"
//---------------------------------------------------------------------------
class GloStiffMatrix
{
public:
	//---------------------------------------------------------------------------
	//构造函数
	GloStiffMatrix( ){};
	//生成总体刚度矩阵
	int Gen_gsmatrix(double* &, int* &, int* &, const vector<Node> &, const vector<Element> &, const vector<MatPro> &);
	//生成单刚；
	int Generate_elestiff(const Element &,const vector<Node> &,const vector<MatPro> &,const vector<Node> &,const vector<double> &);
	//四面体
	void Tetrahedron_line(const vector<Node> &,double ele_elas[][6]);
	//三棱柱
	void Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const vector<Node> &gauss,const vector<double> &wight);
	//将单刚添加到总刚
	int Add_to_gsmatrix(double* &, int* &, int* &, const Element &);
	vector<vector<double> > elestiff;
protected:
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
