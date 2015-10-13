//============================================================
//总荷载向量类头文件；
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
	//构造函数；
	Gloloaded_vector(){};
	//成员函数；
	//生成a1总荷载向量;
	int Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &m,double* &equright);
	//生成a1a2总荷载向量;
    int Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &alfa2,const int &m,double * &equright,const vector<vector<double> > &Na1_vec,const double Homo_D[][6]);
	//生成a1单元荷载向量；
	int Generate_eleloaded(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//生成a1a2单元荷载向量；
    int Generate_eleloaded1(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const vector<vector<double> > &Na1_vec,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//单元荷载向量添加到总荷载向量；
	void Insert_eletoglo(const Element &e,double* &equright);
	//四面体为线性插值时生成a1单元荷载向量；
	void Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const int &cn);
	//四面体生成a1a2单元荷载向量；
	void Tetrahedron_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
    //三棱柱生成a1单元荷载向量；
	void Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//三棱柱生成a1a2单元荷载向量；
	void Threeprism_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight);
	//四阶张量aijhk的后两个坐标与弹性矩阵D列号的对应关系；
	int Mapping(const int &i,const int &j)const;
	vector<double> eleloaded_vec;
	double homoelas[6][6];
protected:

};//----------------------------------------------------------

#endif
//========================================================================

