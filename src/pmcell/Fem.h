//===========================================================================
// Fem.h
// ���޵�Ԫ������ͷ�ļ�
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
//�ڵ���ͷ�ļ���
//--------------------------------------------------------
class Node : public Vector3D{	

public:
	int flag ;//0:�ڵ㣻1�����ʱ߽�㣻2����Ȼ�߽�㣻
	//���캯����
	Node(double ix=0,double iy=0,double iz=0);
	Node(const Vector3D& d);
	//��Ա������
	//���룻
	double	Distance_to(const Node& n); 

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//������Ϊ��Ӧ���ʷֳ�������ӵ�
//-----------------------------------------------
	vector<int> neigh_nodes_vec;
	vector<int> relative_eles_vec; //��ص�Ԫ��������Ԫ��ţ�
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

//��Ԫ��ͷ�ļ���
class Element{

public:	
	int type; //��ʾ��Ԫ����״�����ͣ�0�������Σ�1���ı��Σ�2��������
	int mat;  //��ʾ��Ԫ�Ĳ������ԣ�0����ά��1-�߽������������߽�㣻�߽�����+1�����壻�߽�����+2����ά˿�ⲿ��
	vector<int> nodes_id;
	//���캯����
	Element( int node1_num=-1, int node2_num=-1, int node3_num=-1, int node4_num=-1);
	//��Ա������
	double Quality( vector<Node> &nodes_vec)const;
	friend ostream& operator<<(ostream& o, const Element& e);
protected:

};
inline ostream& operator<<(ostream& o,const Element& e){
	return o<<"("<<e.nodes_id[0]<<", "<<e.nodes_id[1]<<", "<<e.nodes_id[2]<<")"<<endl;
}//=========================================



//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//������Ϊ��Ӧ���ʷֳ�������ӵ�

//---------------------------------------------------------------------------
//����3���ڵ�������ɵ�3��3����ʽ��ֵ
double cal_det( Node *node1, Node *node2, Node *node3 );
//---------------------------------------------------------------------------
//����4�ڵ������嵥Ԫ�����(�ڵ�˳���������ֶ���
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
//���������嵥Ԫ
class Tetrahedron
{
	public:
        int nodesId[4];  //the 4 numbers of the four nodes;
        int materialId ;
        Tetrahedron( );
        Tetrahedron( int node1_num, int node2_num, int node3_num, int node4_num );
        Tetrahedron( Tri_2D_ele& tri, int node1_num );
        int is_contain( int node_num );
        Coors coors;                  //��Ԫ����ϵ
};
//---------------------------------------------------------------------------
//���������嵥Ԫ
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
