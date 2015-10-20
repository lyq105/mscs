//===========================================================================
// Mesher.h
// 四面体网格剖分类头文件
// A class of generating tetrahedron elemments
//===========================================================================
#ifndef MESHER_H
#define MESHER_H

#include<iomanip>
#include<iostream>
#include<fstream>
#include<sstream>
#include "time.h"
#include "CMCell.h"
#include "Fem.h"
#include "Geometry.h"
#include "Hns.h"
using namespace hns;

//---------------------------------------------------------------------------
class EllipsGeometry									//用于生成界面层的椭球几何类
{
public:
	vector<int> tet;										//用于记录椭球中的四面体单元
	vector<vector<int> > tet_nod;				//用于记录四面体单元中那个点是外点
	vector<vector<int> > tri;						//用于记录椭球的表面三角形面片
	vector<int> inell_nodes;						//用于记录椭球中的内点
	vector<int> outell_nodes;						//用于记录椭球中的外点
};
//---------------------------------------------------------------------------
class Mesher 
{
	public:
        vector<Node> nodes_vec;				//存储节点
        vector <int> bnodes_vec;					//边界节点
        vector<Tetrahedron> eles_vec;		//存储四面体单元
		vector<Element> elements_vec;			//单元类向量；

		Mesher(){};										//构造函数
		Mesher(CMCell* cell); 
		int Mesh_generate(ifstream &infile);	//网格生成
		int Mesh_data(int mod);						//网格数据的输出和读取
		int Mesh_BinaryData(int mod, string data_file, int CNum);			//二进制网格数据的输出和读取

	protected:
        CMCell* thecell;			//定义单胞
        double global_length;		//the global length of the edges of the elments;

        int mesh_tet_map(double glength=0.0);       //映射法生成四面体网格

		int gen_coating_mesh(double thick_ratio, int layer_num);			//韩非20060908增改
																										//韩非20061216又改

		//输出tecplot网格数据
		int export_tecplot_data();  //韩非20060913增改
		//输出生成边界层前的tecplot网格数据
		int export_tecplot_data_before_coating();  //韩非20060914增改
		//输出三维空间可视化Ensight网格数据
		int export_Ensight_data();  //韩非20070702增改

		//在不生成边界层的情况下，将四面体单元类型倒为混合单元类型
		int change_elements_to_blend();		//韩非20060915增改
    private:
        int material_num ; //材料数
        double x_min,x_max,y_min,y_max,z_min,z_max;
        double attach_length;

        //标题
        string mtitle;
		string mdata_file;

        //控制节点合法性的参数
        //两线段的夹角小于angle_min,直接连接两端节点，
        //大于angle_min 小于angle_max，在两线段中间生成一节点，到两端节点距离＝？？（见后）
        //大于angle_max时在第一线段法线上寻找一合适点
        double angle_min ;
        double angle_max ;

        vector< vector<int> > gnodes_z0_vec; //底面（z=0）上的给定节点（基体与纤维的边界点）
        vector< vector<int> > gnodes_z1_vec;
        vector< vector<int> > gnodes_x0_vec;   
        vector< vector<int> > gnodes_x1_vec;
        vector< vector<int> > gnodes_y0_vec;
        vector< vector<int> > gnodes_y1_vec;
        vector<int> gnodes_z0;               //底面上（z＝0）上的边界点（所有边界）
        vector<int> gnodes_z1;
        vector<int> gnodes_x0;
        vector<int> gnodes_x1;
        vector<int> gnodes_y0;
        vector<int> gnodes_y1;
        //给定节点（边界点）向量向量，大小等于纤维或颗粒（增强项）数目
        vector<vector<int> > gnodes_vec;
        //存放内部节点，进行光滑处理时使用（底面生成和三维剖分时都使用这个向量存放内部节点）

        //郁闷：六面体节点的局部编号失误，没办法，先用个转换数组顶着
        //理想编号：先x 再y 最后z，目前编号，绕z旋转编号
        int ideal_nn[8];     //由现编号查理想编号
   //   int rideal_nn[8];    //由理想编号查现编号(跟上面一样的）

        void initial();		//单胞初始化

        //判断一个节点是否在边界上，0－－否，非零－－所在的增强颗粒上，
        //mod控制要不要检查纤维和基体交界面
        //给定节点号，检查节点是否在边界上，不使用where_is,直接取is_on_nodes中的值
        int is_on_edge( int node_num, int mod=1 );
        //判断两个节点是否在同一条边界上，mod控制要不要检查增强颗粒和基体交界面
        int is_on_same_edge(int node1_num, int node2_num, int mod = 1);

        int attach_to_edge(Node* node, double a_length);   //检查给定节点和边界的距离，小于a_length则吸附到边界上

        //确定给定节点的位置（基体、第n个椭球）
        int where_is( Node *node,  int &is_on, int node_count = -1 );

        void gen_1tet(int node1_num, int node2_num, int node3_num, int node4_num, int mod);
        void gen_1tet(int node_num, Tri_2D_ele* tri, int mod);

        int is_on_face(Node* node, int mod = 1);        //检测一个节点是否在单胞面上，mod控制要不要检查纤维和基体交界面（很费时间）
        int is_on_face(int node_num, int mod = 1);
        int is_on_face(Tetrahedron* tet, int mod = 1);  //检测四面体是否有面在边界面上
        int is_on_edge(Tetrahedron* tet, int mod = 1);  //检测四面体是否有边在边界线上

        //检测两个节点是否在同一个单胞边界面上，mod控制要不要检查纤维和基体交界面（很费时间）
        int is_on_same_face(int node1_num, int node2_num, int mod = 1);

        //光滑处理三维节点 （拉普拉欣），返回值为是否调整过节点（是否仍需要光顺）
        int smoothing_3d(int mod);

        //快速排序
        void quick_sort( vector<int> *int_vec, int min_n=-1, int max_n=-1 );

        //消除整型向量中连续出现的重复元素
        void unique( vector<int> *int_vec );

        //计算四面体单元的体积
        double cal_tet_volume( Tetrahedron* tet );

        //计算和节点相关的单元的体积最小值
        double cal_relative_min_volume( int node_num );

        //检查每个相关单元（四面体），并在节点Node1_num 和 Node2_num之间插入节点node3_num
        void harmonize_tets( int node1_num, int node2_num, int node3_num,
                             double alength, vector<Bar> &new_edges, int ele_kind );

        //生成背景网格（六面体）
        void gen_bg_mesh(int x_sec_num, int y_sec_num, int z_sec_num, vector<Hexahedron> &hexes_v);

        //检查给定节点的所有相关单元是否四个节点全在边界面上
        int is_re_all_on_face( int node_num, int wh=-1 );

        //检查边界面上的单元在边界面上的节点是否离边界很近
        int check_bft_node( Tetrahedron *tet, double alength );

        //生成边界线上的节点
        int gen_edge_bn( double glength );

        //判断一个单元是否有边在单胞边界线上跨边界
        int has_cross_edge_edge( Tetrahedron *act_tet );

        //判断一个单元是否有边在单胞边界面上跨边界
        int has_cross_edge_edge_on_face( Tetrahedron *act_tet );

        //确定应该调整的边
        //首先选择边界线上跨边界的边，其次选择边界面上跨边界的边，再其次选择两点都能移动的边
        //最后选择一点在边界面上，一点在单胞内部的边
        //最后两个参数标志两端节点是否能够移动
        int deter_oper_edge( Tetrahedron *act_tet, int& left_movable, int& right_movable);

        //根据给定的边的序号，确定节点序列，进行单元拆分的时候使用
        int* deter_nodes_id(int edge_num, int* nodes_id);

        //判断一个节点是否为角点（三个面的交点）
        int is_corner( int node_num );

        //把指定单元分类放入待处理集合（跨边界的单元）
        int put_into_cbev( int ele_num, double alength, int type );

        bool debuging_3d;
        bool debuging_2d;

        double min_volume ;  //容许的最小体积
                                                 
        vector<int> where_is_nodes;  //存放每个节点属于基体还是哪个椭球
        vector<int> is_on_nodes;      //存放节点是否在椭球和基体的交界面上

        vector<int> nodes_gen_later;  //调整四面体以保证其在基体或颗粒内时操作过的节点
                                      //以备以后检查有无四个节点都在椭球面上的单元时使用
        vector<int> thin_tets;        //扁元

        vector<int> deleting_ele;     //删除的单元

        //跨边界的单元
        vector<int> *eles_across_boundary;
                                                   
        
        //边界线上的单元(在边界线上跨边界的单元）
        vector<int> eles_across_boundary_be;
        //边界面上的单元(在边界面上跨边界的单元）
        vector<int> eles_across_boundary_bf;
        //边界面上的单元(在边界面上，但不跨边界的单元）
        vector<int> eles_across_boundary_nbf;
        //两个端点都可以移动的跨边界单元
        vector<int> eles_across_boundary_dmb;
        //两个端点中有一个不可以移动的跨边界单元
        vector<int> eles_across_boundary_ndmb;

        double volume_ratio();   //计算纤维体积比

        double min_volume_of_cell(int &ele_num);    //计算所有单元中的最小体积

        double max_volume_of_cell(int &ele_num);    //计算所有单元中的最小体积

		//韩非20060915增改
        void deter_boundary_nodes(int coat=0);                //确定单胞边界上的节点

        int deter_mat_ele(int ele_num);            //确定一个单元的材料编号,返回材料编号
        int deter_mat_ele(Tetrahedron* tet,int mod =0);            //确定一个单元的材料编号,返回材料编号

        //处理边界线上的单元(在边界线上跨边界的单元）
        int deal_with_eles_across_boundary_be(vector<Bar> &new_edges,int debuging=0);

        //处理边界面上的单元(在边界面上跨边界的单元）
        int deal_with_eles_across_boundary_bf(vector<Bar> &new_edges,int debuging=0);

        //处理两个端点都可以移动的跨边界单元
        int deal_with_eles_across_boundary_dmb(vector<Bar> &new_edges,int debuging=0);

        //处理两个端点中有一个不可以移动的跨边界单元
        int deal_with_eles_across_boundary_ndmb(vector<Bar> &new_edges,int debuging=0);

        //处理跨边界单元
        int deal_with_ele_across_boundary(int ele_num, int ele_kind,
                                          vector<Bar> &new_edges, int mod);

        //把一个六面体分成6个四面体（分法根据六个面的具体情况而定）
        int hex2_6tet(int hex[14], vector<Tetrahedron > &tets,double volume=0.0);
        
        //把一个六面体分成12个四面体（在中心插入一个新节点，用来处理不能成功分成6个四面体的情况）
        int hex2_12tet(int hex[14], vector<Tetrahedron > &tets, double volume=0.0);
                            
        //按给定参数，把六面体（长方体）的6个面分成三角面片，并确定其相互关系（每条边线两边的单元）
        int deter_rela_tri_face(int *hex,vector<Tri_2D_ele> &tri_faces);   

        //根据给定的两个节点号，判断是在六面体的哪个侧面上
        int deter_face_num(int n1, int n2, int *style=NULL);

        //给定两个节点，判断是否为六面体的一条边
        int is_edge_of_hex(int n1, int n2);

        //给定一个三角面片和一个节点，判断它们是否能构成一个四面体，主要考虑边界协调问题
        int can_form_1tet(Tri_2D_ele &tri, int nn, int *hex, int &diag_line, int ls[3]);

        //调整生成的六面体，使所有的六面体都能拆分成四面体，而没有四面体跨边界
        int adjust_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces);

        //根据给定的节点编号（六面体的局部编号），确定相邻的3个节点号（仍然是局部编号）
        int deter_nei(int i, int *nei_i);

        //确定需要检查的4个节点（三角片面对面的四个节点，并把最近的节点放在最后）
        //节点的顺序为：长对角线节点、第一条边对应的斜边节点、第二条边对应的斜边节点、中点对应的节点（构成直边的节点）
        int deter_check_nodes(int i,int n,int check_nn[4]);

        //根据给定的节点号（或者边号，只能指定一个），确定六面体的剖分样式（给定节点（边）在颗粒内的边)
        int deter_hex_split_style(int node_num, int edge_num, int face_style[8]);

        //拆分六面体成四面体
        int split_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces);

        //确定每个四面体的材料
        int deter_eles_mat();

        //根据已有的参数，尽量提供一个能拆分成6个四面体的方案（只修改split_info中等于-1的元素）
        int best_split_hex2_tet(int* hex_info,vector<Tetrahedron> &tets);
        
        //计算两个三角面片之间的夹角,node1_num node2_num 分别为不在公共边上的节点
        double angle_of_2tri(Tri_2D_ele* tri1, Tri_2D_ele* tri2, int node1_num=-1, int node2_num=-1);

        //计算给定节点到给定三角片面的夹角(拆分六面体成6个四面体时使用）
        double angle_point2_tri(int j, Tri_2D_ele* tri, int* hex, int *check_nn);

		//韩非20060908增改
        //根据给定的i j k以及最大i_max j_max k_max，决定给定节点的位置（角点、边界线、边界面、内部）
        //0-5: 边界面上（返回值）
        //6-17: 边界线
        //18: 角点 （所有角点都不能移动，归为一类）
        //产生背景网格时使用
        int deter_node_flag(int i, int j, int k, int i_max, int j_max, int k_max);

		//韩非20060908增改
		//根据给定的节点的编号和单胞的x,y,z方向的最大最小值，
		// 决定给定节点的位置（角点、边界线、边界面、内部）
		//1-3: 边界面上 （6个面归为3类）
		//4-6: 边界线  (12条边界线归为3类)
		//7: 角点（所有角点都不能移动，归为1类）

		int deter_node_flag(int node_num);

        //计算给定节点到椭球面(sur)的投影点（注意单胞侧面上的节点投影也要在侧面上）
        Node project2_elli_surf(Node &point,Surface* sur,int *error=NULL);

        //给定两条直边号，返回构成的三角面片号
        int deter_face_num_from_2bar(int bn1, int bn2);

        //根据给定的侧面编号（六面体的侧面局部编号），以及斜边方向，确定一个节点编号数组
        //数组内容：斜边节点1、斜边节点2、非斜边节点1、非斜边节点2
        int deter_nn_seq(int fn, int style, int nnseq[4]);

        //计算一点到另外两点连线的距离
        double dis2_line(Node *n1, Node *n2, Node *n3);

        //恢复移动过的边界节点（单胞侧面上的点）
        int recover_bnodes(vector<int> &bnodes);

        //检查所有的单元，如果有两个面的夹角非常大（接近180），则在对边插入一个节点
        int rectify_tet(int ele_num, vector<int> &new_ele_num);

        //在两个节点中间插入节点
        int insert_node(int n1, int n2, int n3, vector<int> &new_ele_num);

        //处理四个节点都在界面上的单元（在长边中间插入节点）
        int deal_grov_eles(vector<int> &grov_eles, vector<int> &new_eles_num);

        //检查有无四个节点都在界面上的单元，有则拆之
        int split_grovelling_eles();

        //计算给定两个节点相关单元的材料号（同时包含两个节点的所有单元颗粒编号,0表示基体,1表示颗粒）
        int face_particle(int n1,int n2);

        //判断一个节点的所有相连节点是否有在颗粒内（检查节点所有的相关单元的节点，有节点在颗粒内，返回0，否则（全部在基体内或在颗粒表面）返回1）
        int can_be_extracted(int n1);

        //判断一个单元是否为扁元（节点到对面三角形面片的距离小于global_length/5.0);
        int is_thin_tet(int ele_num);

        //找出一个四面体四个面的最大夹角，并返回相应的节点序列（0-1-2  3-2-1为夹角最大的两面）
        double max_angle_of_tet(int ele_num, int *nseq=NULL);

        //检查一个节点的相关单元中有几个扁元
        int num_of_thin_tet_re(int nn);

        //计算所有单元最小体积和最大体积（数组min_max_en[2]保存最小最大体积的单元号，min_max_vl[2]保存最小最大体积值）
        int cal_min_max_tet_volume(int min_max_en[2], double min_max_vl[2]);

		//计算所有单元的体积分布情况（最小体积向下取整，最大体积向上取整，记录长度）
		//将这段长度分成n份（一般情况n取10），记录每份中所占单元数与总单元数的比值
		//mod=1输出结果，mod=0不输出结果
		void tet_vol_distrib( double min_max_vl[2], string str="", int n=10, int mod=1);

		//计算所有四面体单元的质量度量分布情况（最小值0，最大值取1，记录长度）
		//将这段长度分成n份（一般情况n取10），记录每份中所占单元数与总单元数的比值；
		void tet_quality( string str="", int n=10, int mod=1 );

		//用于计算四面体四个面的面积
		void tet_tri_area( Tetrahedron* tet , double tria[] );

		//计算三棱柱体积
		double Threeprism_volume(const vector<Node> &);

		//用于计算四面体对边边长之积
		void side_length_pro( Tetrahedron* tet , double silp[] );

		//确定所有节点的相关单元信息
		void deter_nodes_relative_eles();

		//人为增加的体分比
		void Reincrese_vol_raito();

		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
