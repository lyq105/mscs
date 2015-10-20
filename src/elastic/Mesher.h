//===========================================================================
// Mesher.h
// �����������ʷ���ͷ�ļ�
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
class EllipsGeometry									//�������ɽ��������򼸺���
{
public:
	vector<int> tet;										//���ڼ�¼�����е������嵥Ԫ
	vector<vector<int> > tet_nod;				//���ڼ�¼�����嵥Ԫ���Ǹ��������
	vector<vector<int> > tri;						//���ڼ�¼����ı�����������Ƭ
	vector<int> inell_nodes;						//���ڼ�¼�����е��ڵ�
	vector<int> outell_nodes;						//���ڼ�¼�����е����
};
//---------------------------------------------------------------------------
class Mesher 
{
	public:
        vector<Node> nodes_vec;				//�洢�ڵ�
        vector <int> bnodes_vec;					//�߽�ڵ�
        vector<Tetrahedron> eles_vec;		//�洢�����嵥Ԫ
		vector<Element> elements_vec;			//��Ԫ��������

		Mesher(){};										//���캯��
		Mesher(CMCell* cell); 
		int Mesh_generate(ifstream &infile);	//��������
		int Mesh_data(int mod);						//�������ݵ�����Ͷ�ȡ
		int Mesh_BinaryData(int mod, string data_file, int CNum);			//�������������ݵ�����Ͷ�ȡ

	protected:
        CMCell* thecell;			//���嵥��
        double global_length;		//the global length of the edges of the elments;

        int mesh_tet_map(double glength=0.0);       //ӳ�䷨��������������

		int gen_coating_mesh(double thick_ratio, int layer_num);			//����20060908����
																										//����20061216�ָ�

		//���tecplot��������
		int export_tecplot_data();  //����20060913����
		//������ɱ߽��ǰ��tecplot��������
		int export_tecplot_data_before_coating();  //����20060914����
		//�����ά�ռ���ӻ�Ensight��������
		int export_Ensight_data();  //����20070702����

		//�ڲ����ɱ߽�������£��������嵥Ԫ���͵�Ϊ��ϵ�Ԫ����
		int change_elements_to_blend();		//����20060915����
    private:
        int material_num ; //������
        double x_min,x_max,y_min,y_max,z_min,z_max;
        double attach_length;

        //����
        string mtitle;
		string mdata_file;

        //���ƽڵ�Ϸ��ԵĲ���
        //���߶εļн�С��angle_min,ֱ���������˽ڵ㣬
        //����angle_min С��angle_max�������߶��м�����һ�ڵ㣬�����˽ڵ���룽����������
        //����angle_maxʱ�ڵ�һ�߶η�����Ѱ��һ���ʵ�
        double angle_min ;
        double angle_max ;

        vector< vector<int> > gnodes_z0_vec; //���棨z=0���ϵĸ����ڵ㣨��������ά�ı߽�㣩
        vector< vector<int> > gnodes_z1_vec;
        vector< vector<int> > gnodes_x0_vec;   
        vector< vector<int> > gnodes_x1_vec;
        vector< vector<int> > gnodes_y0_vec;
        vector< vector<int> > gnodes_y1_vec;
        vector<int> gnodes_z0;               //�����ϣ�z��0���ϵı߽�㣨���б߽磩
        vector<int> gnodes_z1;
        vector<int> gnodes_x0;
        vector<int> gnodes_x1;
        vector<int> gnodes_y0;
        vector<int> gnodes_y1;
        //�����ڵ㣨�߽�㣩������������С������ά���������ǿ���Ŀ
        vector<vector<int> > gnodes_vec;
        //����ڲ��ڵ㣬���й⻬����ʱʹ�ã��������ɺ���ά�ʷ�ʱ��ʹ�������������ڲ��ڵ㣩

        //���ƣ�������ڵ�ľֲ����ʧ��û�취�����ø�ת�����鶥��
        //�����ţ���x ��y ���z��Ŀǰ��ţ���z��ת���
        int ideal_nn[8];     //���ֱ�Ų�������
   //   int rideal_nn[8];    //�������Ų��ֱ��(������һ���ģ�

        void initial();		//������ʼ��

        //�ж�һ���ڵ��Ƿ��ڱ߽��ϣ�0�����񣬷��㣭�����ڵ���ǿ�����ϣ�
        //mod����Ҫ��Ҫ�����ά�ͻ��彻����
        //�����ڵ�ţ����ڵ��Ƿ��ڱ߽��ϣ���ʹ��where_is,ֱ��ȡis_on_nodes�е�ֵ
        int is_on_edge( int node_num, int mod=1 );
        //�ж������ڵ��Ƿ���ͬһ���߽��ϣ�mod����Ҫ��Ҫ�����ǿ�����ͻ��彻����
        int is_on_same_edge(int node1_num, int node2_num, int mod = 1);

        int attach_to_edge(Node* node, double a_length);   //�������ڵ�ͱ߽�ľ��룬С��a_length���������߽���

        //ȷ�������ڵ��λ�ã����塢��n������
        int where_is( Node *node,  int &is_on, int node_count = -1 );

        void gen_1tet(int node1_num, int node2_num, int node3_num, int node4_num, int mod);
        void gen_1tet(int node_num, Tri_2D_ele* tri, int mod);

        int is_on_face(Node* node, int mod = 1);        //���һ���ڵ��Ƿ��ڵ������ϣ�mod����Ҫ��Ҫ�����ά�ͻ��彻���棨�ܷ�ʱ�䣩
        int is_on_face(int node_num, int mod = 1);
        int is_on_face(Tetrahedron* tet, int mod = 1);  //����������Ƿ������ڱ߽�����
        int is_on_edge(Tetrahedron* tet, int mod = 1);  //����������Ƿ��б��ڱ߽�����

        //��������ڵ��Ƿ���ͬһ�������߽����ϣ�mod����Ҫ��Ҫ�����ά�ͻ��彻���棨�ܷ�ʱ�䣩
        int is_on_same_face(int node1_num, int node2_num, int mod = 1);

        //�⻬������ά�ڵ� ������������������ֵΪ�Ƿ�������ڵ㣨�Ƿ�����Ҫ��˳��
        int smoothing_3d(int mod);

        //��������
        void quick_sort( vector<int> *int_vec, int min_n=-1, int max_n=-1 );

        //���������������������ֵ��ظ�Ԫ��
        void unique( vector<int> *int_vec );

        //���������嵥Ԫ�����
        double cal_tet_volume( Tetrahedron* tet );

        //����ͽڵ���صĵ�Ԫ�������Сֵ
        double cal_relative_min_volume( int node_num );

        //���ÿ����ص�Ԫ�������壩�����ڽڵ�Node1_num �� Node2_num֮�����ڵ�node3_num
        void harmonize_tets( int node1_num, int node2_num, int node3_num,
                             double alength, vector<Bar> &new_edges, int ele_kind );

        //���ɱ������������壩
        void gen_bg_mesh(int x_sec_num, int y_sec_num, int z_sec_num, vector<Hexahedron> &hexes_v);

        //�������ڵ��������ص�Ԫ�Ƿ��ĸ��ڵ�ȫ�ڱ߽�����
        int is_re_all_on_face( int node_num, int wh=-1 );

        //���߽����ϵĵ�Ԫ�ڱ߽����ϵĽڵ��Ƿ���߽�ܽ�
        int check_bft_node( Tetrahedron *tet, double alength );

        //���ɱ߽����ϵĽڵ�
        int gen_edge_bn( double glength );

        //�ж�һ����Ԫ�Ƿ��б��ڵ����߽����Ͽ�߽�
        int has_cross_edge_edge( Tetrahedron *act_tet );

        //�ж�һ����Ԫ�Ƿ��б��ڵ����߽����Ͽ�߽�
        int has_cross_edge_edge_on_face( Tetrahedron *act_tet );

        //ȷ��Ӧ�õ����ı�
        //����ѡ��߽����Ͽ�߽�ıߣ����ѡ��߽����Ͽ�߽�ıߣ������ѡ�����㶼���ƶ��ı�
        //���ѡ��һ���ڱ߽����ϣ�һ���ڵ����ڲ��ı�
        //�������������־���˽ڵ��Ƿ��ܹ��ƶ�
        int deter_oper_edge( Tetrahedron *act_tet, int& left_movable, int& right_movable);

        //���ݸ����ıߵ���ţ�ȷ���ڵ����У����е�Ԫ��ֵ�ʱ��ʹ��
        int* deter_nodes_id(int edge_num, int* nodes_id);

        //�ж�һ���ڵ��Ƿ�Ϊ�ǵ㣨������Ľ��㣩
        int is_corner( int node_num );

        //��ָ����Ԫ�������������ϣ���߽�ĵ�Ԫ��
        int put_into_cbev( int ele_num, double alength, int type );

        bool debuging_3d;
        bool debuging_2d;

        double min_volume ;  //�������С���
                                                 
        vector<int> where_is_nodes;  //���ÿ���ڵ����ڻ��廹���ĸ�����
        vector<int> is_on_nodes;      //��Žڵ��Ƿ�������ͻ���Ľ�������

        vector<int> nodes_gen_later;  //�����������Ա�֤���ڻ���������ʱ�������Ľڵ�
                                      //�Ա��Ժ��������ĸ��ڵ㶼���������ϵĵ�Ԫʱʹ��
        vector<int> thin_tets;        //��Ԫ

        vector<int> deleting_ele;     //ɾ���ĵ�Ԫ

        //��߽�ĵ�Ԫ
        vector<int> *eles_across_boundary;
                                                   
        
        //�߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
        vector<int> eles_across_boundary_be;
        //�߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
        vector<int> eles_across_boundary_bf;
        //�߽����ϵĵ�Ԫ(�ڱ߽����ϣ�������߽�ĵ�Ԫ��
        vector<int> eles_across_boundary_nbf;
        //�����˵㶼�����ƶ��Ŀ�߽絥Ԫ
        vector<int> eles_across_boundary_dmb;
        //�����˵�����һ���������ƶ��Ŀ�߽絥Ԫ
        vector<int> eles_across_boundary_ndmb;

        double volume_ratio();   //������ά�����

        double min_volume_of_cell(int &ele_num);    //�������е�Ԫ�е���С���

        double max_volume_of_cell(int &ele_num);    //�������е�Ԫ�е���С���

		//����20060915����
        void deter_boundary_nodes(int coat=0);                //ȷ�������߽��ϵĽڵ�

        int deter_mat_ele(int ele_num);            //ȷ��һ����Ԫ�Ĳ��ϱ��,���ز��ϱ��
        int deter_mat_ele(Tetrahedron* tet,int mod =0);            //ȷ��һ����Ԫ�Ĳ��ϱ��,���ز��ϱ��

        //����߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
        int deal_with_eles_across_boundary_be(vector<Bar> &new_edges,int debuging=0);

        //����߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
        int deal_with_eles_across_boundary_bf(vector<Bar> &new_edges,int debuging=0);

        //���������˵㶼�����ƶ��Ŀ�߽絥Ԫ
        int deal_with_eles_across_boundary_dmb(vector<Bar> &new_edges,int debuging=0);

        //���������˵�����һ���������ƶ��Ŀ�߽絥Ԫ
        int deal_with_eles_across_boundary_ndmb(vector<Bar> &new_edges,int debuging=0);

        //�����߽絥Ԫ
        int deal_with_ele_across_boundary(int ele_num, int ele_kind,
                                          vector<Bar> &new_edges, int mod);

        //��һ��������ֳ�6�������壨�ַ�����������ľ������������
        int hex2_6tet(int hex[14], vector<Tetrahedron > &tets,double volume=0.0);
        
        //��һ��������ֳ�12�������壨�����Ĳ���һ���½ڵ㣬���������ܳɹ��ֳ�6��������������
        int hex2_12tet(int hex[14], vector<Tetrahedron > &tets, double volume=0.0);
                            
        //�������������������壨�����壩��6����ֳ�������Ƭ����ȷ�����໥��ϵ��ÿ���������ߵĵ�Ԫ��
        int deter_rela_tri_face(int *hex,vector<Tri_2D_ele> &tri_faces);   

        //���ݸ����������ڵ�ţ��ж�������������ĸ�������
        int deter_face_num(int n1, int n2, int *style=NULL);

        //���������ڵ㣬�ж��Ƿ�Ϊ�������һ����
        int is_edge_of_hex(int n1, int n2);

        //����һ��������Ƭ��һ���ڵ㣬�ж������Ƿ��ܹ���һ�������壬��Ҫ���Ǳ߽�Э������
        int can_form_1tet(Tri_2D_ele &tri, int nn, int *hex, int &diag_line, int ls[3]);

        //�������ɵ������壬ʹ���е������嶼�ܲ�ֳ������壬��û���������߽�
        int adjust_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces);

        //���ݸ����Ľڵ��ţ�������ľֲ���ţ���ȷ�����ڵ�3���ڵ�ţ���Ȼ�Ǿֲ���ţ�
        int deter_nei(int i, int *nei_i);

        //ȷ����Ҫ����4���ڵ㣨����Ƭ�������ĸ��ڵ㣬��������Ľڵ�������
        //�ڵ��˳��Ϊ�����Խ��߽ڵ㡢��һ���߶�Ӧ��б�߽ڵ㡢�ڶ����߶�Ӧ��б�߽ڵ㡢�е��Ӧ�Ľڵ㣨����ֱ�ߵĽڵ㣩
        int deter_check_nodes(int i,int n,int check_nn[4]);

        //���ݸ����Ľڵ�ţ����߱ߺţ�ֻ��ָ��һ������ȷ����������ʷ���ʽ�������ڵ㣨�ߣ��ڿ����ڵı�)
        int deter_hex_split_style(int node_num, int edge_num, int face_style[8]);

        //����������������
        int split_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces);

        //ȷ��ÿ��������Ĳ���
        int deter_eles_mat();

        //�������еĲ����������ṩһ���ܲ�ֳ�6��������ķ�����ֻ�޸�split_info�е���-1��Ԫ�أ�
        int best_split_hex2_tet(int* hex_info,vector<Tetrahedron> &tets);
        
        //��������������Ƭ֮��ļн�,node1_num node2_num �ֱ�Ϊ���ڹ������ϵĽڵ�
        double angle_of_2tri(Tri_2D_ele* tri1, Tri_2D_ele* tri2, int node1_num=-1, int node2_num=-1);

        //��������ڵ㵽��������Ƭ��ļн�(����������6��������ʱʹ�ã�
        double angle_point2_tri(int j, Tri_2D_ele* tri, int* hex, int *check_nn);

		//����20060908����
        //���ݸ�����i j k�Լ����i_max j_max k_max�����������ڵ��λ�ã��ǵ㡢�߽��ߡ��߽��桢�ڲ���
        //0-5: �߽����ϣ�����ֵ��
        //6-17: �߽���
        //18: �ǵ� �����нǵ㶼�����ƶ�����Ϊһ�ࣩ
        //������������ʱʹ��
        int deter_node_flag(int i, int j, int k, int i_max, int j_max, int k_max);

		//����20060908����
		//���ݸ����Ľڵ�ı�ź͵�����x,y,z����������Сֵ��
		// ���������ڵ��λ�ã��ǵ㡢�߽��ߡ��߽��桢�ڲ���
		//1-3: �߽����� ��6�����Ϊ3�ࣩ
		//4-6: �߽���  (12���߽��߹�Ϊ3��)
		//7: �ǵ㣨���нǵ㶼�����ƶ�����Ϊ1�ࣩ

		int deter_node_flag(int node_num);

        //��������ڵ㵽������(sur)��ͶӰ�㣨ע�ⵥ�������ϵĽڵ�ͶӰҲҪ�ڲ����ϣ�
        Node project2_elli_surf(Node &point,Surface* sur,int *error=NULL);

        //��������ֱ�ߺţ����ع��ɵ�������Ƭ��
        int deter_face_num_from_2bar(int bn1, int bn2);

        //���ݸ����Ĳ����ţ�������Ĳ���ֲ���ţ����Լ�б�߷���ȷ��һ���ڵ�������
        //�������ݣ�б�߽ڵ�1��б�߽ڵ�2����б�߽ڵ�1����б�߽ڵ�2
        int deter_nn_seq(int fn, int style, int nnseq[4]);

        //����һ�㵽�����������ߵľ���
        double dis2_line(Node *n1, Node *n2, Node *n3);

        //�ָ��ƶ����ı߽�ڵ㣨���������ϵĵ㣩
        int recover_bnodes(vector<int> &bnodes);

        //������еĵ�Ԫ�������������ļнǷǳ��󣨽ӽ�180�������ڶԱ߲���һ���ڵ�
        int rectify_tet(int ele_num, vector<int> &new_ele_num);

        //�������ڵ��м����ڵ�
        int insert_node(int n1, int n2, int n3, vector<int> &new_ele_num);

        //�����ĸ��ڵ㶼�ڽ����ϵĵ�Ԫ���ڳ����м����ڵ㣩
        int deal_grov_eles(vector<int> &grov_eles, vector<int> &new_eles_num);

        //��������ĸ��ڵ㶼�ڽ����ϵĵ�Ԫ�������֮
        int split_grovelling_eles();

        //������������ڵ���ص�Ԫ�Ĳ��Ϻţ�ͬʱ���������ڵ�����е�Ԫ�������,0��ʾ����,1��ʾ������
        int face_particle(int n1,int n2);

        //�ж�һ���ڵ�����������ڵ��Ƿ����ڿ����ڣ����ڵ����е���ص�Ԫ�Ľڵ㣬�нڵ��ڿ����ڣ�����0������ȫ���ڻ����ڻ��ڿ������棩����1��
        int can_be_extracted(int n1);

        //�ж�һ����Ԫ�Ƿ�Ϊ��Ԫ���ڵ㵽������������Ƭ�ľ���С��global_length/5.0);
        int is_thin_tet(int ele_num);

        //�ҳ�һ���������ĸ�������нǣ���������Ӧ�Ľڵ����У�0-1-2  3-2-1Ϊ�н��������棩
        double max_angle_of_tet(int ele_num, int *nseq=NULL);

        //���һ���ڵ����ص�Ԫ���м�����Ԫ
        int num_of_thin_tet_re(int nn);

        //�������е�Ԫ��С�����������������min_max_en[2]������С�������ĵ�Ԫ�ţ�min_max_vl[2]������С������ֵ��
        int cal_min_max_tet_volume(int min_max_en[2], double min_max_vl[2]);

		//�������е�Ԫ������ֲ��������С�������ȡ��������������ȡ������¼���ȣ�
		//����γ��ȷֳ�n�ݣ�һ�����nȡ10������¼ÿ������ռ��Ԫ�����ܵ�Ԫ���ı�ֵ
		//mod=1��������mod=0��������
		void tet_vol_distrib( double min_max_vl[2], string str="", int n=10, int mod=1);

		//�������������嵥Ԫ�����������ֲ��������Сֵ0�����ֵȡ1����¼���ȣ�
		//����γ��ȷֳ�n�ݣ�һ�����nȡ10������¼ÿ������ռ��Ԫ�����ܵ�Ԫ���ı�ֵ��
		void tet_quality( string str="", int n=10, int mod=1 );

		//���ڼ����������ĸ�������
		void tet_tri_area( Tetrahedron* tet , double tria[] );

		//�������������
		double Threeprism_volume(const vector<Node> &);

		//���ڼ���������Ա߱߳�֮��
		void side_length_pro( Tetrahedron* tet , double silp[] );

		//ȷ�����нڵ����ص�Ԫ��Ϣ
		void deter_nodes_relative_eles();

		//��Ϊ���ӵ���ֱ�
		void Reincrese_vol_raito();

		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
