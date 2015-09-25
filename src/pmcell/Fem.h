//===========================================================================
// Fem.h
// 有限单元方法类头文件
// A class of Finite Element Method
//===========================================================================

#ifndef HEADER_FEM
#define HEADER_FEM
#include<iostream>
#include<cmath>
#include<vector>
#include"Vector3D.h"
#include "Geometry.h"
#include "Hns.h"
using namespace hns;
//=================================================
//节点类头文件；
//--------------------------------------------------------
class Node : public Vector3D{	

public:
	int flag ;//0:内点；1：本质边界点；2：自然边界点；
	//构造函数；
	Node(double ix=0,double iy=0,double iz=0);
	Node(const Vector3D& d);
	//成员函数；
	//距离；
	double	Distance_to(const Node& n); 

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//以下是为了应用剖分程序所添加的
//-----------------------------------------------
	vector<int> neigh_nodes_vec;
	vector<int> relative_eles_vec; //相关单元向量（单元编号）
	int real_id;

	Node( Point pt );
	double distance_to( Node *node );
	double distance_to( Point *point );
    
	void move_to( Node *node );
	void move_to( Point *point );

	Node operator+ ( Node node );
	Node operator- ( Node node );
	Node operator* ( Node node );
	Node operator/ ( Node node );
	Node operator* ( double m );
	Node operator/ ( double d );
	bool operator< ( Node node );
	bool operator<= ( Node node );
	bool operator> ( Node node );
	bool operator>= ( Node node );
	bool operator== ( Node node );
	bool operator!= ( Node node );
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

};//-----------------------------------------------

//单元类头文件；
class Element{

public:	
	int type; //表示单元的形状的类型；0：三角形；1：四边形；2：四面体
	int mat;  //表示单元的材料物性，0：纤维；1-边界层层数：单胞边界层；边界层层数+1：基体；边界层层数+2：纤维丝外部；
	vector<int> nodes_id;
	//构造函数；
	Element( int node1_num=-1, int node2_num=-1, int node3_num=-1, int node4_num=-1);
	//成员函数；
	double Quality( vector<Node> &nodes_vec)const;
	friend ostream& operator<<(ostream& o, const Element& e);
protected:

};
inline ostream& operator<<(ostream& o,const Element& e){
	return o<<"("<<e.nodes_id[0]<<", "<<e.nodes_id[1]<<", "<<e.nodes_id[2]<<")"<<endl;
}//=========================================



//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//以下是为了应用剖分程序所添加的

//---------------------------------------------------------------------------
//计算3个节点坐标组成的3×3行列式的值
double cal_det( Node *node1, Node *node2, Node *node3 );
//---------------------------------------------------------------------------
//计算4节点四面体单元的体积(节点顺序满足右手定则）
double cal_volume( Node *node1, Node *node2, Node *node3, Node *node4 );
//---------------------------------------------------------------------------
class Bar
{
	public:
        Bar(int n1=0, int n2=0)
		{
			nodesId[0] = n1;
			nodesId[1] = n2;
        }
        int nodesId[2];  //the 2 numbers of the four nodes;
};
//---------------------------------------------------------------------------
class Tri_2D_ele 
{
	public:
        int nodesId[3];  //the 3 numbers of the nodes;
        int edgesId[3];  //the 3 numbers of the edges;
        Tri_2D_ele( int node1_num=0, int node2_num=0, int node3_num=0 )
		{
			nodesId[0] = node1_num;
			nodesId[1] = node2_num;
			nodesId[2] = node3_num;
			neighbor.clear();
        }
        Tri_2D_ele( int node1_num, int node2_num, int node3_num,
                    int edge1_num, int edge2_num, int edge3_num=0 )
		{
			nodesId[0] = node1_num;
			nodesId[1] = node2_num;
			nodesId[2] = node3_num;
			edgesId[0] = edge1_num;
			edgesId[1] = edge2_num;
			edgesId[2] = edge3_num;
			neighbor.clear();
        }
        vector<int> neighbor;
};
//---------------------------------------------------------------------------
//定义四面体单元
class Tetrahedron
{
	public:
        int nodesId[4];  //the 4 numbers of the four nodes;
        int materialId ;
        Tetrahedron( );
        Tetrahedron( int node1_num, int node2_num, int node3_num, int node4_num );
        Tetrahedron( Tri_2D_ele& tri, int node1_num );
        int is_contain( int node_num );
        Coors coors;                  //单元坐标系
};
//---------------------------------------------------------------------------
//定义六面体单元
class Hexahedron
{
	public:
        int nodesId[8];	 //the 8 numbers of the four nodes;
        int materialId;
        int facesId[6];
        Hexahedron( int n1, int n2, int n3, int n4, int n5, int n6, int n7, int n8 );
        Hexahedron( int *nn=NULL );
};
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//===============================================================================
#endif
