//=======================================
//三维节点类定义实现；
//Fem.cpp
//=======================================
#include"Fem.h"
//==============================================================
//Node类构造函数；
Node::Node(double ix, double iy,double iz){
	x=ix;
	y=iy;
	z=iz;
	flag=0;
}//---------------------------------------

Node::Node(const Vector3D& d){
	x = d.x;
	y = d.y;
	z = d.z;
	flag=0;
}//---------------------------------------
//距离；
double Node::Distance_to(const Node& n){
	return(sqrt( (x-n.x)*(x-n.x)+(y-n.y)*(y-n.y)+(z-n.z)*(z-n.z)) );
}//----------------------------------------
//=================================================================
//Element 构造函数；
Element::Element( int node1_num, int node2_num, int node3_num, int node4_num)
{
	if( node1_num != -1) nodes_id.push_back(node1_num);
	if( node2_num != -1) nodes_id.push_back(node2_num);
	if( node3_num != -1) nodes_id.push_back(node3_num);
	if( node4_num != -1) nodes_id.push_back(node4_num);
	mat=0;
	type=0;
}
double Element::Quality( vector<Node> &nodes_vec)const{
	double s ;
	if(nodes_id.size()==3){
		double a = nodes_vec[nodes_id[0]].Distance_to(nodes_vec[nodes_id[1]]);
		double b = nodes_vec[nodes_id[0]].Distance_to(nodes_vec[nodes_id[2]]);
		double c = nodes_vec[nodes_id[1]].Distance_to(nodes_vec[nodes_id[2]]);
		double s1 = 1.0/2.0*(a+b+c);
		double area = sqrt(s1*(s1-a)*(s1-b)*(s1-c));
		s=sqrt(3.0)/6*(a*a+b*b+c*c)/area/2.0;
	}
	else if(nodes_id.size()==4){
		double a = nodes_vec[nodes_id[0]].Distance_to(nodes_vec[nodes_id[1]]);
		double b = nodes_vec[nodes_id[1]].Distance_to(nodes_vec[nodes_id[2]]);
		double c = nodes_vec[nodes_id[2]].Distance_to(nodes_vec[nodes_id[3]]);
		double d = nodes_vec[nodes_id[3]].Distance_to(nodes_vec[nodes_id[0]]);
		double e = nodes_vec[nodes_id[0]].Distance_to(nodes_vec[nodes_id[2]]);
		double f = nodes_vec[nodes_id[1]].Distance_to(nodes_vec[nodes_id[3]]);
		double s1 = 1.0/2.0*(a+d+f);
		double s2 = 1.0/2.0*(b+c+e);
		double area1 = sqrt(s1*(s1-a)*(s1-d)*(s1-f));
		double area2 = sqrt(s2*(s2-b)*(s2-c)*(s2-e));
		s = (a*a+b*b+c*c+d*d+2*e*e+2*f*f)/(6*(area1+area2))/2.0;
	}
	else{
		hout<<"There is a error in Quality()"<<endl;
	}
	return s;
}//---------------------------------------------------------------------------



//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//以下是为了应用剖分程序所添加的
//---------------------------------------------------------------------------
//计算3个节点坐标组成的3×3行列式的值
double cal_det( Node *node1, Node *node2, Node *node3 )
{
	double det = node1->x * ( node2->y * node3->z - node2->z * node3->y ) +
				 node1->y * ( node2->z * node3->x - node2->x * node3->z ) +
				 node1->z * ( node2->x * node3->y - node2->y * node3->x );
	return det ;
}
//---------------------------------------------------------------------------
//计算4节点四面体单元的体积(节点顺序满足右手定则）
double cal_volume( Node *node1, Node *node2, Node *node3, Node *node4 )
{
	double volume;
	double aai = cal_det( node2, node3, node4 );
	double aaj = cal_det( node3, node4, node1 );
	double aak = cal_det( node4, node1, node2 );
	double aal = cal_det( node1, node2, node3 );

	volume = 1.0/6.0*(aai-aaj+aak-aal) ;
	return volume ;
}
//---------------------------------------------------------------------------
Node::Node( Point pt )
{
	x = pt.x ;
	y = pt.y ;
	z = pt.z ; 
	flag = 0;
}
//---------------------------------------------------------------------------
void Node::move_to( Node *node )
{
	x = node->x;
	y = node->y;
	z = node->z;
}
//---------------------------------------------------------------------------
void Node::move_to( Point *point )
{
	x = point->x;
	y = point->y;
	z = point->z;
}
//---------------------------------------------------------------------------
double Node::distance_to( Node *node )
{
	return sqrt( (x - node->x)*(x - node->x)+
				 (y - node->y)*(y - node->y)+
				 (z - node->z)*(z - node->z) );
}   
//---------------------------------------------------------------------------
double Node::distance_to( Point *point )
{
	return sqrt( (x - point->x)*(x - point->x)+
				 (y - point->y)*(y - point->y)+
				 (z - point->z)*(z - point->z) );
}
//---------------------------------------------------------------------------
Node Node::operator+ ( Node node )
{    
	return Node( x + node.x, y + node.y, z + node.z );
}   
//---------------------------------------------------------------------------
Node Node::operator- ( Node node )
{
	return Node( x - node.x, y - node.y, z - node.z );
}  
//---------------------------------------------------------------------------
Node Node::operator* ( Node node )
{
	return Node( x * node.x, y * node.y, z * node.z );
}  
//---------------------------------------------------------------------------
Node Node::operator/ ( Node node )
{
	return Node( x / node.x, y / node.y, z / node.z );
}
//---------------------------------------------------------------------------
Node Node::operator* ( double m )
{
	return Node( x * m, y * m, z * m );
}
//---------------------------------------------------------------------------
Node Node::operator/ ( double d )
{
	return Node( x / d, y / d, z / d );
}  
//---------------------------------------------------------------------------
bool Node::operator< ( Node node )
{
	return ( x < node.x && y < node.y && z < node.z );
} 
//---------------------------------------------------------------------------
bool Node::operator<= ( Node node )
{
	return ( x <= node.x && y <= node.y && z <= node.z );
}  
//---------------------------------------------------------------------------
bool Node::operator> ( Node node )
{
	return ( x > node.x && y > node.y && z > node.z );
}  
//---------------------------------------------------------------------------
bool Node::operator>= ( Node node )
{
	return ( x >= node.x && y >= node.y && z >= node.z );
}
//---------------------------------------------------------------------------
bool Node::operator== ( Node node )
{
	return ( x == node.x && y == node.y && z == node.z );
} 
//---------------------------------------------------------------------------
bool Node::operator!= ( Node node )
{
	return ( x != node.x || y != node.y || z != node.z );
}
//---------------------------------------------------------------------------
Tetrahedron::Tetrahedron(){}
//---------------------------------------------------------------------------
Tetrahedron::Tetrahedron( int node1_num, int node2_num,
                          int node3_num, int node4_num )
{
        nodesId[0] = node1_num;
        nodesId[1] = node2_num;
        nodesId[2] = node3_num;
        nodesId[3] = node4_num;
        Coors crs;
        coors = crs;
}
//---------------------------------------------------------------------------
Tetrahedron::Tetrahedron( Tri_2D_ele& tri, int node1_num )
{
        nodesId[0] = tri.nodesId[0];
        nodesId[1] = tri.nodesId[1];
        nodesId[2] = tri.nodesId[2];
        nodesId[3] = node1_num;
        Coors crs;
        coors = crs;
}
//---------------------------------------------------------------------------

int Tetrahedron::is_contain( int node_num )
{
        for( int i=0; i<4; i++ )
		{
			if( nodesId[i] == node_num )
			{
				return i;
			}
		}
        return -1;
}
//---------------------------------------------------------------------------
//六面体单元构造函数
Hexahedron::Hexahedron( int n1, int n2, int n3, int n4,
                        int n5, int n6, int n7, int n8 )
{
	nodesId[0] = n1; nodesId[1] = n2;
	nodesId[2] = n3; nodesId[3] = n4;
	nodesId[4] = n5; nodesId[5] = n6;
	nodesId[6] = n7; nodesId[7] = n8;
}
//---------------------------------------------------------------------------
Hexahedron::Hexahedron( int *nn )
{
	if( nn == NULL )
	{
		for( int i=0; i<8; i++ )
		{
			nodesId[i]=0;
		}
	}
	else
	{
		for( int i=0; i<8; i++ )
		{
			nodesId[i] = nn[i];
		}
	}
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

