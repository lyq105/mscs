//===========================================================================
// Mesher.cpp
// �����������ʷ����Ա����
// Member Functions in a Class of generating tetrahedron elemments
//===========================================================================
#include "Mesher.h"
#define PI 3.141592654
#define ZERO 1e-8

//---------------------------------------------------------------------------
//���캯��
Mesher::Mesher(CMCell* cmccell)
{
	thecell=cmccell;
	initial();
}
//--------------------------------------------------------------------------
//��ʼ��
void Mesher::initial(  )
{
	attach_length = global_length / 2.0 ;
	for( int i=0; i<(int)thecell->surfaces_vec.size(); i++ )
	{
		vector<int> nodes;
		gnodes_vec.push_back(nodes); //��ʼ��
		vector<int> nodes_z0;
		gnodes_z0_vec.push_back(nodes_z0);        
		vector<int> nodes_z1;
		gnodes_z1_vec.push_back(nodes_z1);
		vector<int> nodes_y0;
		gnodes_y0_vec.push_back(nodes_y0);
		vector<int> nodes_y1;
		gnodes_y1_vec.push_back(nodes_y1);
		vector<int> nodes_x0;
		gnodes_x0_vec.push_back(nodes_x0);
		vector<int> nodes_x1;
		gnodes_x1_vec.push_back(nodes_x1);    
	}

	ideal_nn[0] = 0;
	ideal_nn[1] = 1;
	ideal_nn[2] = 3;
	ideal_nn[3] = 2;
	ideal_nn[4] = 4;
	ideal_nn[5] = 5;
	ideal_nn[6] = 7;
	ideal_nn[7] = 6;    

	debuging_3d = true;
	debuging_2d = false;
}
//---------------------------------------------------------------------------
//��������
int Mesher::Mesh_generate(ifstream &infile)
{
	//��ȡ�����ʷֳߴ�
	istringstream in(Get_Line(infile));
	double mesh_size;
	in >> mesh_size;
	//�����йؽ������Ϣ
	istringstream in1(Get_Line(infile));
	int coat, layer_num;
	double thick_ratio;
	in1 >> coat >> thick_ratio >> layer_num;

	if(mesh_tet_map(mesh_size) == 0 )
	{
		hout << "~_~ �ʷ�ʧ�ܣ�����ԭ����ǿ���ص�������С�����Ԫ����������Ԫ���̫��!" << endl;
		hout << "****************************************************************************************" <<endl;	
		return 0;
	}
	
	//������ɱ߽��ǰ��tecplot��������
	//if (export_tecplot_data_before_coating() == 0)
	//{
	//	hout << "~_~ ������ɱ߽��ǰ��tecplot��������ʧ�ܣ�" << endl;
	//	hout << "*************************************************" << endl;
	//	return 0;
	//}

	if (coat==1)
	{
		//����ȷ��materialId ==0 �ǻ��嵥Ԫ��materialId ==1 ������Ԫ
		clock_t ct0,ct1;
		ct0 = clock();
		hout << "-_- ��ʼ���ɽ��������... " << endl;
		if(gen_coating_mesh(thick_ratio, layer_num) == 0)
		{
			hout << "~_~ ���ɽ��������ʧ�ܣ�" << endl;
			hout << "*******************************************" <<endl;	
			return 0;
		}

		hout << "    ���tecplot��������....." << endl;
		if (export_tecplot_data() == 0)
		{
			hout << "~_~ ���tecplot��������ʧ�ܣ�" << endl;
			hout << "*******************************************" << endl;
			return 0;
		}

		//hout<<"   �����ά�ռ���ӻ�Ensight��������"
		//if (export_Ensight_data() == 0)
		//{
		//	hout << "~_~ �����ά�ռ���ӻ�Ensight��������ʧ�ܣ�" << endl;
		//	hout << "*******************************************" << endl;
		//	return 0;
		//}

		ct1 = clock();
		hout << "    ���ɽ�����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��" << endl;
		hout << "^_^ ����������������!" << endl<<endl;
	}
	else
	{
		//�ڲ����ɱ߽�������£��������嵥Ԫ���͵�Ϊ��ϵ�Ԫ����
		clock_t ct0,ct1;
		ct0 = clock();
		hout << "-_- ������������ת��...... " << endl;
		if (change_elements_to_blend() == 0)
		{
			hout << " �������嵥Ԫ���͵�Ϊ��ϵ�Ԫ���Ͳ���ʧ�ܣ�" << endl;
			hout << "*******************************************" << endl;
			return 0;
		}
		ct1 = clock();
		hout << "    ת�������������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��" << endl;
		hout << "^_^ ת���������������!" << endl << endl;
	}
	
	//ȷ�����нڵ����ص�Ԫ��Ϣ
	deter_nodes_relative_eles();

	return 1;
}
//---------------------------------------------------------------------------
//������ȡ��������
int Mesher::Mesh_data(int mod)
{
	if(mod==0)			//�������
	{
		//---------------------------------------------------------------------
		//�����������
		ofstream oute("Data_Elements.dat");
		int nume = int(elements_vec.size());
		oute << nume << endl;
		for(int i=0; i<nume; i++)
		{
			int num_enod = int(elements_vec[i].nodes_id.size());
			oute << num_enod << "  ";
			for(int j=0; j<num_enod; j++)
			{
				oute << elements_vec[i].nodes_id[j] << "  ";
			}
			oute << elements_vec[i].type << "  " << elements_vec[i].mat << "  " << endl;
		}
		oute.close();
		//---------------------------------------------------------------------
		//����ڵ�����
		ofstream outn("Data_Nodes.dat");
		outn << int(nodes_vec.size()) << endl;
		for(int i=0; i<int(nodes_vec.size()); i++)
		{
			outn << nodes_vec[i].flag << "  " << setprecision(12)
				    << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
		}
		outn.close();
		//---------------------------------------------------------------------
		//����߽�ڵ�����
		ofstream outb("Data_Bnodes.dat");
		outb << int(bnodes_vec.size()) << endl;
		for(int i=0; i<int(bnodes_vec.size()); i++)
		{
			outb << bnodes_vec[i] << endl;
		}
		outb.close();
	}
	else if(mod==1)	//��ȡ����
	{
		//---------------------------------------------------------------------
		//��ȡ��������
		ifstream ine("Data_Elements.dat");
		int nume;
		ine >> nume;
		for(int i=0; i<nume; i++)
		{
			Element temp_ele;
			elements_vec.push_back(temp_ele);
			int num_enod;
			ine >> num_enod;
			for(int j=0; j<num_enod; j++)
			{
				int nodes_id;
				ine >> nodes_id;
				elements_vec.back().nodes_id.push_back(nodes_id);
			}
			ine >> elements_vec.back().type;
			ine >> elements_vec.back().mat;
		}
		ine.close();
		//---------------------------------------------------------------------
		//��ȡ�ڵ�����
		ifstream inn("Data_Nodes.dat");
		int numn;
		inn >> numn;
		for(int i=0; i<numn; i++)
		{
			Node temp_node;
			nodes_vec.push_back(temp_node);
			inn >> nodes_vec.back().flag >> nodes_vec.back().x >> nodes_vec.back().y >> nodes_vec.back().z;
		}
		inn.close();
		//---------------------------------------------------------------------
		//��ȡ�߽�ڵ�����
		ifstream inb("Data_Bnodes.dat");
		int numb;
		inb >> numb;
		for(int i=0; i<numb; i++)
		{
			int bnode;
			inb >> bnode;
			bnodes_vec.push_back(bnode);
		}
		inb.close();		
	}

	return 1;
}
//---------------------------------------------------------------------------
//������ȡ��������������
int Mesher::Mesh_BinaryData(int mod, string data_file, int CNum)
{
	int num1 = CNum/10;
	int num2 = CNum - num1*10;
	char ch[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	if(mod==0)			//�������
	{
		//---------------------------------------------------------------------
		//�����������
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Mesh";
		string MeshName = title+str+ch[num1]+ch[num2]+".dat";
		ofstream out(MeshName.c_str(),ios::binary);
		if(!out) { hout << "���ܴ����������ļ�" << MeshName << "!" << endl;  exit(0); }
		//---------------------------------------------------------------------
		//��Ԫ����
		int nume = int(elements_vec.size());
		out.write((char *)&nume, sizeof(int));
		for(int i=0; i<nume; i++)
		{
			int num_enod = int(elements_vec[i].nodes_id.size());
			out.write((char *)&num_enod, sizeof(int));
			for(int j=0; j<num_enod; j++)	out.write((char *)&elements_vec[i].nodes_id[j], sizeof(int));
			out.write((char *)&elements_vec[i].type, sizeof(int));
			out.write((char *)&elements_vec[i].mat, sizeof(int));
		}
		//---------------------------------------------------------------------
		//�ڵ�����
		int numn = int(nodes_vec.size());
		out.write((char *)&numn, sizeof(int));
		for(int i=0; i<numn; i++)
		{
			out.write((char *)&nodes_vec[i].flag, sizeof(int));
			out.write((char *)&nodes_vec[i].x, sizeof(double));
			out.write((char *)&nodes_vec[i].y, sizeof(double));
			out.write((char *)&nodes_vec[i].z, sizeof(double));
			int num_relative = int(nodes_vec[i].relative_eles_vec.size());
			out.write((char *)&num_relative, sizeof(int));
			for(int j=0; j<num_relative; j++)	out.write((char *)&nodes_vec[i].relative_eles_vec[j], sizeof(int));
		}
		//---------------------------------------------------------------------
		//�߽ڵ�����
		int numb = int(bnodes_vec.size());
		out.write((char *)&numb, sizeof(int));
		for(int i=0; i<numb; i++)	out.write((char *)&bnodes_vec[i],sizeof(int));
		//---------------------------------------------------------------------
		out.close();
	}
	else if(mod==1)	//��ȡ����
	{
		//---------------------------------------------------------------------
		//�����������
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Mesh";
		string MeshName = title+str+ch[num1]+ch[num2]+".dat";
		ifstream in(MeshName.c_str(),ios::binary);
		if(!in) { hout << "���ܴ����������ļ�" << MeshName << "!" << endl;  exit(0); }
		//---------------------------------------------------------------------
		//��Ԫ����
		int nume;
		in.read((char *)&nume, sizeof(int));
		Element temp_ele;
		elements_vec.assign(nume,temp_ele);
		for(int i=0; i<nume; i++)
		{
			int num_enod;
			in.read((char *)&num_enod, sizeof(int));
			elements_vec[i].nodes_id.assign(num_enod, 0);
			for(int j=0; j<num_enod; j++)	in.read((char *)&elements_vec[i].nodes_id[j], sizeof(int));
			in.read((char *)&elements_vec[i].type, sizeof(int));
			in.read((char *)&elements_vec[i].mat, sizeof(int));
		}
		//---------------------------------------------------------------------
		//�ڵ�����
		int numn;
		in.read((char *)&numn, sizeof(int));
		Node temp_node;
		nodes_vec.assign(numn,temp_node);
		for(int i=0; i<numn; i++)
		{
			in.read((char *)&nodes_vec[i].flag, sizeof(int));
			in.read((char *)&nodes_vec[i].x, sizeof(double));
			in.read((char *)&nodes_vec[i].y, sizeof(double));
			in.read((char *)&nodes_vec[i].z, sizeof(double));
			int num_relative;
			in.read((char *)&num_relative, sizeof(int));
			nodes_vec[i].relative_eles_vec.assign(num_relative, 0);
			for(int j=0; j<num_relative; j++)		in.read((char *)&nodes_vec[i].relative_eles_vec[j], sizeof(int));
		}
		//---------------------------------------------------------------------
		//�߽ڵ�����
		int numb;
		in.read((char *)&numb, sizeof(int));
		bnodes_vec.assign(numb,0);
		for(int i=0; i<numb; i++)	in.read((char *)&bnodes_vec[i],sizeof(int));
		//---------------------------------------------------------------------
		in.close();	
	}
	return 1;
}
//---------------------------------------------------------------------------
//ӳ�䷨����������������Ҫ����Σ�
int Mesher::mesh_tet_map( double glength )
{
	clock_t ct_mesh_begin = clock();

	bool debuging_map_mesh = false;

	//ȫ��ʹ���ʷֳ߶�
	global_length = glength;
	//��Ԫ��С���
	min_volume = 0.0 ;
                                
	//�����������С����
	x_min=thecell->origin_x ;
	x_max=thecell->origin_x + thecell->clength ;
	y_min=thecell->origin_y ;
	y_max=thecell->origin_y + thecell->cwidth ;
	z_min=thecell->origin_z ;
	z_max=thecell->origin_z + thecell->cheight ;

	nodes_vec.clear();

	hout<< "-_- ��ʼ�����������ʷ�" << endl<<endl;
	//---------------------------------------------------------------------------
	hout<< "-_- ��ʼ���ɱ������������壩......" << endl;
	clock_t ct_background_begin = clock();

	//���ɱ�������
	//����ÿ�����ϷֵĶ�����ceil(num)����ѧ����ĺ��������ز�С��num������(������ֵ��ʾ)
	int x_sec_num = (int)ceil((x_max-x_min)/glength) ;
	int y_sec_num = (int)ceil((y_max-y_min)/glength) ;
	int z_sec_num = (int)ceil((z_max-z_min)/glength) ;
         
	vector< Hexahedron > hexes_vec;
	//���ɱ������������壩
	gen_bg_mesh(x_sec_num,y_sec_num,z_sec_num, hexes_vec);

	where_is_nodes.assign(nodes_vec.size(), -2);	//���ÿ���ڵ����ڻ��廹���ĸ�����
	is_on_nodes.assign(nodes_vec.size(), 0);			//��Žڵ��Ƿ�������ͻ���Ľ�������
	where_is_nodes.reserve(nodes_vec.size()*2);
	is_on_nodes.reserve(nodes_vec.size()*2);

	//�����ı�����Ƭ
	int num_rec_xfaces=(x_sec_num+1)*y_sec_num*z_sec_num;
	int num_rec_yfaces=(y_sec_num+1)*x_sec_num*z_sec_num;
	int num_rec_zfaces=(z_sec_num+1)*x_sec_num*y_sec_num;
	int num_rec_faces=num_rec_xfaces+num_rec_yfaces+num_rec_zfaces;
	vector<int> sign_rec_faces(num_rec_faces,-1);

	clock_t ct_background_end = clock();
	hout << "    ���ɣ������壩���������ʱ��" << (double)( ct_background_end - ct_background_begin )/CLOCKS_PER_SEC << " �롣" << endl;
	hout << "^_^ ���ɱ������������壩��ϡ�" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- ��ʼ�������ɵ�������������Ͻ���......" << endl;
	clock_t ct_adj_begin = clock();

	//�������ɵ������壬ʹ���е������嶼�ܲ�ֳ������壬��û���������߽�
	if(adjust_hexes(hexes_vec,sign_rec_faces) == 0)
	{
		hout << " �������ɵ����������(adjust_hexes)ʧ�ܣ�" << endl;
		hout << "*********************************************" << endl;
		return 0;
	}

	clock_t ct_adj_end = clock();
	hout << "    �������������������ʱ��" << (double)( ct_adj_end - ct_adj_begin )/CLOCKS_PER_SEC << " �롣" << endl;
	hout << "^_^ �������ɵ�������������ϡ�" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- ��ʼ���������Ϊ�����壬��ϱ߽�......" << endl;  
	clock_t ct_split_begin = clock();

	//���������Ϊ������
	if(split_hexes(hexes_vec,sign_rec_faces) == 0) 
	{
		hout << "���������Ϊ���������ʧ�ܣ�" << endl;
		hout << "***********************************" << endl;
		return 0;
	}

	hout << "    ���ɽڵ�nodes: " << (int)nodes_vec.size() << "    ���ɵ�Ԫelements: " << (int)eles_vec.size() << endl;
        
	clock_t ct_split_end = clock();
	hout << "    ���������Ϊ�������ʱ��" << (double)( ct_split_end - ct_split_begin )/CLOCKS_PER_SEC << " �롣" << endl;
	hout << "^_^ ���������Ϊ��������ϡ�" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- ȷ��ÿ��������Ĳ�������......" << endl;
 	clock_t ct_mat_begin = clock();

	//ȷ��ÿ��������Ĳ���
	deter_eles_mat();

 	clock_t ct_mat_end = clock();;
	hout << "    ȷ��ÿ��������Ĳ������ʺ�ʱ��"<< (double)( ct_mat_end - ct_mat_begin )/CLOCKS_PER_SEC << " �롣" << endl;
	hout << "^_^ ȷ��ÿ��������Ĳ���������ϡ�" << endl<<endl;

	//---------------------------------------------------------------------------
	hout<< "-_- ��������浥Ԫ��ʹ���е�Ԫ����̴��߽磬��������Ԫ��" << endl ;
	clock_t ct_bound_begin = clock();
	double alength = glength*0.5 ;

	//���ÿ����Ԫ�ı��Ƿ�ᴩ�������ǿ����
	//�ҳ���߽�ĵ�Ԫ�������ȴ������߽��Ͽ�߽�ĵ�Ԫ
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		int ele_mat = eles_vec[i].materialId;
		if( ele_mat == -1 )
		{
			Tetrahedron *act_tet = &eles_vec[i] ;

			double alength1 = alength /1.0;

			int pic = put_into_cbev( i, alength, 0 );
			if( pic == -1 )
			{
				if( check_bft_node( act_tet, alength1 ) == 1 )
				{
					i-- ;
				}
			}
		}
	}

	//�������
	if( debuging_map_mesh ) hout << " eles_across_boundary_be.size(): " << (int)eles_across_boundary_be.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_bf.size(): " << (int)eles_across_boundary_bf.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_nbf.size(): " << (int)eles_across_boundary_nbf.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_dmb.size(): " << (int)eles_across_boundary_dmb.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_ndmb.size(): " << (int)eles_across_boundary_ndmb.size() << endl;
        
  
	hout << "    ��߽絥Ԫ����" <<(int)(	eles_across_boundary_be.size()+
																	eles_across_boundary_bf.size()+
																	eles_across_boundary_dmb.size()+
																	eles_across_boundary_ndmb.size() ) << endl;

	//��¼��������Ϣ��ÿ��Ԫ������������int�����ݣ�
	//ǰ����Ϊ������ڵ���߶ζ˵�ڵ�ţ�����������Ľڵ�ţ�
	vector<vector<int> > inserted_point;

	eles_across_boundary = &eles_across_boundary_be ;

	vector<Bar> new_edges;


	//��������ĸ��ڵ㶼�ڽ����ϵĵ�Ԫ�������֮
	if( split_grovelling_eles() == 1 )
	{
		//�����ֱ߽絥Ԫ����ѭ��
		int debuging = 0;
		if( deal_with_eles_across_boundary_be(new_edges,debuging) == 0 ) return 0;
		if( deal_with_eles_across_boundary_bf(new_edges,debuging) == 0 ) return 0;
		if( deal_with_eles_across_boundary_dmb(new_edges,debuging) == 0 ) return 0;
		if( deal_with_eles_across_boundary_ndmb(new_edges,debuging) == 0 ) return 0;
	}
	else
	{
		return 0;
	}

	//for(int i=0; i < (int)eles_vec.size(); i++)
	//{
	//	if (eles_vec[i].materialId == -1)
	//	{
	//		hout << "����-1���֣�" << endl;
	//	}
	//}

	//�����ɾ��һЩ������Ԫ
	quick_sort(&deleting_ele);
	for( int i=(int)deleting_ele.size()-1; i>=0; i-- )
	{
		eles_vec.erase(eles_vec.begin()+deleting_ele[i]);
		if( debuging_map_mesh ) hout << "erase element: " << deleting_ele[i] << endl;
	}

	//��������ɵı������޺�С�ıߣ�����������ϲ��ڵ㣩
//  double min_dis=global_length/4.0;
//	check_tiny_edges(new_edges,min_dis); ���˹�����δ��ɣ�

	clock_t ct_bound_end = clock();
	hout<< "    ��������浥Ԫ��ʱ��" << (double)( ct_bound_end - ct_bound_begin )/CLOCKS_PER_SEC << " �롣" << endl;
	hout<< "^_^ ��������浥Ԫ�ɹ����." << endl;

	int min_max_en[2];
	double min_max_vl[2];
	cal_min_max_tet_volume(min_max_en,min_max_vl);
	hout << "      ��С�����" << min_max_vl[0] << " ��Ԫ��" << min_max_en[0] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[0]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[0]].materialId << endl;

	hout << "      ��������" << min_max_vl[1] << " ��Ԫ��" << min_max_en[1] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[1]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[1]].materialId << endl;
	hout << endl;
	
	//�������������ֲ����
//	tet_vol_distrib(min_max_vl,"_a");
	//��������������ֲ����
//	tet_quality("_a");
	//---------------------------------------------------------------------------
	hout<< "-_- ��ʼ���й�˳����......" << endl;
	clock_t ct_smooth_begin = clock();

	//��˳����
	//Ϊÿ���ڵ�������ڽڵ㣨��ص�Ԫ��
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		nodes_vec[i].neigh_nodes_vec.clear();
		nodes_vec[i].relative_eles_vec.clear();
	}
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		for( int j=0; j<4; j++ )
		{
			for( int k=0; k<4; k++ )
			{
				if( k != j )
				{
					nodes_vec[eles_vec[i].nodesId[j]].neigh_nodes_vec.push_back( eles_vec[i].nodesId[k] );
				}
			}
			nodes_vec[eles_vec[i].nodesId[j]].relative_eles_vec.push_back(i);
		}
	}

	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		quick_sort( &nodes_vec[i].neigh_nodes_vec );
		unique( &nodes_vec[i].neigh_nodes_vec );
	}

	int circle_num = 0 ;
	//��˳�Ǹ���εĹ���
	while( smoothing_3d(1) && smoothing_3d(0) && circle_num < 15 )
	{
		circle_num ++ ;
		hout << "    ��˳���� NO: " << circle_num << "��" <<endl;
	}

	//���һ����û�е�Ԫ��������ûȷ����ֻ�������������Ƿ�������������Ϊ������˺ܶ࣬���ұ�֤��
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		if( eles_vec[i].materialId == -1 )
		{
			hout <<"Element "	<< i ;
			hout << " Nodes: "	<< eles_vec[i].nodesId[0] << " "
											<< eles_vec[i].nodesId[1] << " "
											<< eles_vec[i].nodesId[2] << " "
											<< eles_vec[i].nodesId[3] ;
			hout << " material id: -1." <<endl;

			//ǿ����Щ��Ԫ����Ϊ��������
			eles_vec[i].materialId = 1;
			hout << "ǿ��Element" << i << "materialidΪ1." << endl;
		}
	}

	clock_t ct_smooth_end = clock();
	hout<< "    ��˳�����ʱ��" << (double)( ct_smooth_end - ct_smooth_begin )/CLOCKS_PER_SEC << " ��" << endl;
	hout<< "    ��Ԫ��"<< (int)eles_vec.size() <<" �ڵ㣺" << (int)nodes_vec.size() << endl;
	hout<< "^_^ ��˳����ɹ����." << endl;


	cal_min_max_tet_volume(min_max_en,min_max_vl);
	hout << "      ��С�����" << min_max_vl[0] << " ��Ԫ��" << min_max_en[0] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[0]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[0]].materialId << endl;

	hout << "      ��������" << min_max_vl[1] << " ��Ԫ��" << min_max_en[1] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[1]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[1]].materialId << endl;

	double  min_max_vl_radio = min_max_vl[0]/min_max_vl[1];
	hout << "      ��С�����������֮��Ϊ��" << min_max_vl_radio << endl;
	if(  min_max_vl_radio < 0.01 )
	{
		hout << "      ��ע��! �ʷֵ�Ԫ����С��������������̫�󣡣�" << endl;
//		return 0;
	}

	//�������������ֲ����
//	tet_vol_distrib(min_max_vl,"_b");
	//��������������ֲ����
//	tet_quality("_b");

	//��Ϊ���ӵ���ֱ�
//	Reincrese_vol_raito();
	
	double vol_ratio = volume_ratio();
	hout << "      ���������������֮�ȣ�" << vol_ratio << endl <<endl;

	 //ȷ���߽�ڵ�͵�Ԫ
	deter_boundary_nodes();

	clock_t ct_mesh_end = clock();
	hout << "      ���ʷֽڵ㣺" << (int)nodes_vec.size() << "  ��Ԫ��" << (int)eles_vec.size() << endl;
	hout << "      �����ʷֹ��̹���ʱ��" << (double)( ct_mesh_end - ct_mesh_begin )/CLOCKS_PER_SEC << " ��" << endl;
	hout << "^_^ �����������ʷֳɹ�!" << endl<<endl; 

	return 1;
}

//---------------------------------------------------------------------------
//���ɱ�������
void Mesher::gen_bg_mesh(int x_sec_num, int y_sec_num, int z_sec_num, vector<Hexahedron> &hexes_v)
{
	double dx = (x_max-x_min)/x_sec_num;
	double dy = (y_max-y_min)/y_sec_num;
	double dz = (z_max-z_min)/z_sec_num;
	x_sec_num ++ ;  //x����ȷֵ���
	y_sec_num ++ ;  //y����ȷֵ���
	z_sec_num ++ ;  //z����ȷֵ���

	if( dx < global_length ) global_length = dx;
	if( dy < global_length ) global_length = dy;
	if( dz < global_length ) global_length = dz;

	double hex_volume = dx*dy*dz;
        
	//���ɽڵ�
	for( int i=0; i<z_sec_num; i++ )
	{
		double z = z_min + i * dz ;
		for( int j=0; j<y_sec_num; j++ )
		{
			double y = y_min + j * dy ;
			for( int k=0; k<x_sec_num; k++ )
			{
				double x = x_min + k * dx ;

				Node nd(x,y,z);
				nd.flag = deter_node_flag(i,j,k,z_sec_num-1,y_sec_num-1,x_sec_num-1);	//��ע�ڵ��λ��
				if( nd.flag >= 0 ) bnodes_vec.push_back( (int)nodes_vec.size() );				//��ע�߽�ڵ����ڽڵ������еı��
				nodes_vec.push_back(nd);
			}
		}
	}

	//��������Ŀ
	int num_eles=(x_sec_num-1)*(y_sec_num-1)*(z_sec_num-1);
	//������Ƭ����Ŀ
	int num_rec_xfaces=x_sec_num*(y_sec_num-1)*(z_sec_num-1);
	int num_rec_yfaces=y_sec_num*(x_sec_num-1)*(z_sec_num-1);
	int num_rec_zfaces=z_sec_num*(x_sec_num-1)*(y_sec_num-1);
	int num_rec_faces=num_rec_xfaces+num_rec_yfaces+num_rec_zfaces;

	//���������嵥Ԫ(8���ڵ�ţ�6��������Ƭ�ţ�
	int now_i=0;
	for( int i=0; i<z_sec_num-1; i++ )
	{
		for( int j=0; j<y_sec_num-1; j++ )
		{
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//������İ˸�����
				int node1_num = i * x_sec_num * y_sec_num + j * x_sec_num + k ;
				int node2_num = i * x_sec_num * y_sec_num + j * x_sec_num + k + 1 ;
				int node3_num = i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1 ;
				int node4_num = i * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k ;
				int node5_num = ( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k ;
				int node6_num = ( i + 1 ) * x_sec_num * y_sec_num + j * x_sec_num + k + 1 ;
				int node7_num = ( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k + 1 ;
				int node8_num = ( i + 1 ) * x_sec_num * y_sec_num + ( j + 1 ) * x_sec_num + k ;

				Hexahedron hex(	node1_num, node2_num, node3_num, node4_num,
											node5_num, node6_num, node7_num, node8_num	);
				hex.facesId[0] = i * x_sec_num *(y_sec_num-1) + j * x_sec_num + k;
				hex.facesId[1] = hex.facesId[0]+1;
				hex.facesId[2] = i * (x_sec_num-1) * y_sec_num + j * (x_sec_num-1) + k+num_rec_xfaces;
				hex.facesId[3] = i * (x_sec_num-1) * y_sec_num + (j+1) * (x_sec_num-1) + k+num_rec_xfaces;
				hex.facesId[4] = i * (x_sec_num-1) * (y_sec_num-1) + j * (x_sec_num-1) + k+num_rec_xfaces+num_rec_yfaces;
				hex.facesId[5] = (i+1) * (x_sec_num-1) * (y_sec_num-1) + j * (x_sec_num-1) + k+num_rec_xfaces+num_rec_yfaces;

				hexes_v.push_back(hex);
			}
		}
	}
}

//---------------------------------------------------------------------------
//���ݸ�����i j k�Լ����i_max j_max k_max�����������ڵ��λ�ã��ǵ㡢�߽��ߡ��߽��桢�ڲ���
//0-5: �߽����ϣ�����ֵ��
//6-17: �߽���
//18: �ǵ㣨���нǵ㶼�����ƶ�����Ϊһ�ࣩ
//������������ʱʹ��
//i--z����
//j--y����
//k--x����
int Mesher::deter_node_flag(int i, int j, int k, int i_max, int j_max, int k_max)
{
	int flag;
	vector<int> faces_num;   //�����ĸ���
	if( i == 0 )          faces_num.push_back(4);
	else if( i == i_max ) faces_num.push_back(5);

	if( j == 0 )          faces_num.push_back(2);
	else if( j == j_max ) faces_num.push_back(3);

	if( k == 0 )          faces_num.push_back(0);
	else if( k == k_max ) faces_num.push_back(1);

	int bfn = (int)faces_num.size();		//��¼�ڼ�������
	if( bfn == 0 )									//�ڲ��ڵ�
	{
		flag = -1;	
	}
	else if( bfn == 1 )							//�߽����ϵĽڵ�	
	{
		flag = faces_num[0];	
	}
	else if( bfn == 2 )							//�ڱ߽�����
	{            
		int min_bfn = min(faces_num[0],faces_num[1]);
		int max_bfn = max(faces_num[0],faces_num[1]);
		if( min_bfn == 0 )
		{
			if( max_bfn == 2 ) flag = 6+8;
			if( max_bfn == 3 ) flag = 6+11;
			if( max_bfn == 4 ) flag = 6+3;
			if( max_bfn == 5 ) flag = 6+7;
		}
		else if( min_bfn == 1 )
		{
			if( max_bfn == 2 ) flag = 6+9;
			if( max_bfn == 3 ) flag = 6+10;
			if( max_bfn == 4 ) flag = 6+1;
			if( max_bfn == 5 ) flag = 6+5;
		}
		else if( min_bfn == 2 )
		{
			if( max_bfn == 4 ) flag = 6+0;
			if( max_bfn == 5 ) flag = 6+4;
		}
		else if( min_bfn == 3 )
		{
			if( max_bfn == 4 ) flag = 6+2;
			if( max_bfn == 5 ) flag = 6+6;
		}
	}
	else if( bfn == 3 )
	{
		//�ǵ㣬�����������ĸ��ǵ��ˣ����нǵ�Ҳ�����ƶ�
		flag = 18;
	}
	return flag;
}

//---------------------------------------------------------------------------
//����20060908����
//���ݸ����Ľڵ�ı�ź͵�����x,y,z����������Сֵ��
// ���������ڵ��λ�ã��ǵ㡢�߽��ߡ��߽��桢�ڲ���
//1-3: �߽����� ��6�����Ϊ3�ࣩ
//4-6: �߽���  (12���߽��߹�Ϊ3��)
//7: �ǵ㣨���нǵ㶼�����ƶ�����Ϊ1�ࣩ

int Mesher::deter_node_flag(int node_num)
{
	int flag;
	vector<int> faces_num;   //����������
	Node *node = &nodes_vec[node_num];

	if (fabs(node->x - x_min) <= ZERO || fabs(node->x - x_max) <= ZERO)
	{
		faces_num.push_back(1);		//��x��Сֵ�����ֵ�����ڣ�x���ܶ���y,z���Զ�
	}

	if (fabs(node->y - y_min) <= ZERO|| fabs(node->y - y_max) <= ZERO)
	{
		faces_num.push_back(2);		//��y��Сֵ�����ֵ�����ڣ�y���ܶ���x,z���Զ�
	}

	if( fabs(node->z - z_min) <= ZERO|| fabs(node->z - z_max) <= ZERO)
	{
		faces_num.push_back(3);		//��z��Сֵ�����ֵ�����ڣ�z���ܶ���x,y���Զ�
	}

	int bfn = (int)faces_num.size();		//��¼�ڼ�������
	if( bfn == 0 )									//�ڲ��ڵ�
	{
		flag = 0;	
	}
	else if( bfn == 1 )							//�߽����ϵĽڵ�	
	{
		flag = faces_num[0];	
	}
	else if( bfn == 2 )							//�ڱ߽�����
	{            
		int min_bfn = min(faces_num[0],faces_num[1]);
		int max_bfn = max(faces_num[0],faces_num[1]);
		if( min_bfn == 1 )
		{
			if( max_bfn == 2 ) flag = 4;	//��x��y��Ľ����ϣ�ֻ��z�ɶ�
			if( max_bfn == 3 ) flag = 5;   //��x��z��Ľ����ϣ�ֻ��y�ɶ�
		}
		else if( min_bfn == 2 )
		{
			if( max_bfn == 3 ) flag =6;	//��y��z��Ľ����ϣ�ֻ��x�ɶ�
		}
	}
	else if( bfn == 3 )
	{
		//�ǵ㣬�����������ĸ��ǵ��ˣ����нǵ�Ҳ�����ƶ�
		flag = 7;
	}

	return flag;
}

//---------------------------------------------------------------------------
//�������ɵ������壬ʹ���е������嶼�ܲ�ֳ������壬��û���������߽�
int Mesher::adjust_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces)
{
	bool debuging = false;

//   hout << "global_length: " << global_length << endl;

	//��С���뷧ֵ��С��������뽫�ڵ��ƶ������棬��������ڵ�ľ�����ò�Ҫ�ƶ�
	//��������Ƭ���㵽�Ա߾���С����������ʾ��Ƭ���Ϸ�����������         
	double min_dis = global_length*3.0/8.0;

	//ȷ��ÿ���ڵ��λ�ã��ڻ����ڻ�����ǿ���ڣ����ǽ����ϣ�
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		int is_on;
		where_is_nodes[i] = where_is( &nodes_vec[i],  is_on, i );
		if( where_is_nodes[i] == -3 ) 
		{
			hout << " ȷ��ÿ���ڵ��λ�ò���(where_is)ʧ��! " << endl;
			hout << "*******************************************" << endl;
			return 0;
		}
		is_on_nodes[i] = is_on;
	}

	//��־�ڵ��Ƿ�����ƶ�
	//���Ѿ��ƶ����������ϵĲ��������ƶ���
	//������˽ڵ���������������У��Ѿ�������һ�����������ĸ��ڵ��ڽ����ϣ��˽ڵ�Ҳ�����ƶ�)
	//	0: ��ʾ�������ƶ�     1: ��ʾ�����ƶ�
	vector<int> movable_ns( nodes_vec.size(), 1 );
	int pv=0;
	int pl=0;
	//��ÿ����������п���
	for( int i=0; i<(int)hexes.size(); i++ )
	{
		if( debuging )	//���������Ԫ���ڵ��źͽڵ�����λ�����
		{
			hout << "-----------------------------------------------------------" << endl;
			hout << "hexahedron " << i << "(" << (int)hexes.size() << "): ";
			for( int j=0; j<8; j++ )
			{
				hout << "      " << hexes[i].nodesId[j];
			}
			hout << endl;

			int show = 0;
			for( int j=0; j<8; j++ )
			{
				int nn = hexes[i].nodesId[j];
				int win = where_is_nodes[nn];
				if( win != -1 && is_on_nodes[nn] == 0 )
				{
					show=1;
					break;
				}
			}
			if( show == 1)
			{
				hout << "hex No. " << i << endl;
				hout << "    where is the nodes: " ;
				hout.setf(ios::right,ios::adjustfield);
				for( int j=0; j<8; j++ )
				{
					hout << "    " << where_is_nodes[hexes[i].nodesId[j]];
				}
				hout << endl;
				hout << "    is on the surface : " ;
				for( int j=0; j<8; j++ )
				{
					hout << "    " << is_on_nodes[hexes[i].nodesId[j]];
				}
				hout << endl;
			}
		}

		//����8���ڵ��Ƿ��ڽ����ϵ���Ϣ������ʱ���
		vector<int> is_on_ln;
		for( int j=0; j<8; j++ )
		{
			int nn = hexes[i].nodesId[j];
			if( is_on_nodes[nn] == 1 )
			{
				is_on_ln.push_back( j );
			}
		}

		for( int j=0; j<8; j++ )
		{
			int nn = hexes[i].nodesId[j];
			int win = where_is_nodes[nn];
			if( win != -1 && is_on_nodes[nn] == 0 )			//�����ڿ����ڵĽڵ�
			{   
				if( debuging )
				{
					hout << "       node: " << j << endl;
				}

				int nei_i[3];    //�������ڵĽڵ��ţ��ֲ���ţ�
				deter_nei(j,nei_i);
				int nei_num[3],nei_win[3],nei_is_on[3];
				for( int k=0; k<3; k++ )		//ͨ��������ֲ�����ҵ���ؽڵ��ȫ�ֱ�Ų�ȷ��λ��
				{
					nei_num[k] = hexes[i].nodesId[nei_i[k]];
					nei_win[k] = where_is_nodes[nei_num[k]];
					nei_is_on[k] = is_on_nodes[nei_num[k]];
				}
				//������ڵ������ڵ㶼�ͱ��ڵ���ͬһ�������������������ƶ��������ϣ���֤�������ڽڵ���ͬһ�������ڻ���������ϣ���
				if( nei_win[0] == win && nei_win[1] == win && nei_win[2] == win ) continue;

				//���ڵ��Լ��ڽڵ㵽�������ͶӰ�㣨ע�ⵥ�������ϵĽڵ�ͶӰҲҪ�ڲ����ϣ�
				Node zerop;
				vector<Node> pnodes(4,zerop);
				int error[4]={-1,-1,-1,-1};     //��־�Ƿ�ɹ�ͶӰ   0��ʧ��    1���ɹ�
				pnodes[0] = project2_elli_surf(nodes_vec[nn],thecell->surfaces_vec[win],&error[0]);
				pv=pv+1;
                //�����뱾�ڵ㲻��ͬһ�����ڵĽڵ㵽�ÿ����Ľ����ϣ������ڽڵ㣬���ߵ������ڵ㣩
				//����ȷ���ǵ����ڽڵ㻹�ǵ������ڵ㣬��������ڽڵ㣬���뱾�ڵ㲻��ͬһ�����ڵ��ڽڵ㶼Ҫ�����������ϣ�����������ڵ㣬�������ڵ���Բ��ٿ��ǣ�����ѭ��
				//���ݵ�������ľ������ж�
				double dis0;

				if( error[0] == 0 )
				{
					//ͶӰʧ�ܣ�����dis0һ���㹻���������ֹ�ж��ڽڵ�ʱ��ѡ��˽ڵ�
					dis0 = 2.0*global_length;
					pl=pl+1;
				}
				else
				{
					dis0 = nodes_vec[nn].distance_to(&pnodes[0]);
				}

				if( debuging ) hout << "         dis0: " << dis0 << endl;

				double nei_dis[3];

				if( debuging ) hout << "         nei dis: " ;

				int move_nei=1;			//�Ƿ��ܹ��ƶ����ڽڵ�ı�־������0�������ƶ�������1�����ƶ�

				//�����ڽڵ㵽�������ͶӰ����(nei_dis)
				for( int k=0; k<3; k++ )
				{     
					//������ڽڵ��뱾�ڵ���ͬһ������������
					if( nei_win[k] == win ) continue;    
					if( nei_win[k] != -1 && nei_is_on[k] == 0 )
					{
						//�����������ڣ�����ǿ���ƶ�����������(���ĸ���������ƶ����ĸ��ϣ�ע���п��ܿ����ص���
						Point pp0 = {nodes_vec[nn].x,nodes_vec[nn].y,nodes_vec[nn].z};
						Point pp1 = {nodes_vec[nei_num[k]].x,nodes_vec[nei_num[k]].y,nodes_vec[nei_num[k]].z};
						Point ip0,ip1;
						//��pp0��pp1�������߻��ӳ������win����nei_win[k]����������Ľ��㣬ipo���ؽ��㣬
						//perr1==1��ʾ�󽻵�ɹ���perr1==0��ʾ�󽻵�ʧ��
						int perr1 = thecell->surfaces_vec[win]->intersect(pp0,pp1,ip0);
						int perr2 = thecell->surfaces_vec[nei_win[k]]->intersect(pp0,pp1,ip1);

						//�ٶ��������ҵ�����
						//���k�����ڵ㵽����ľ���
						double dis1 = nodes_vec[nei_num[k]].distance_to(&ip0);
						double dis2 = nodes_vec[nei_num[k]].distance_to(&ip1);

						if( debuging )
						{
							hout << "nn: " << nn << " nn1: " << nei_num[k] ;
							hout << " dis1: " << dis1 << " dis2: " << dis2 << endl;
						}

						if( dis1 > dis2 )			//������ڽ��ڵ�ֱ�Ӵ������ڵ�����ֱ����������ڵ���������
						{
							nodes_vec[nei_num[k]].move_to(&ip1);
							is_on_nodes[nei_num[k]] = 1;
							movable_ns[nei_num[k]] = 0;      //�������ƶ���
							nei_is_on[k] = 1;
						}
						else
						{
							hout << " dis2 > dis1���������ص��������ݲ����� " << endl;
							hout << "*******************************************" << endl;
							return 0;
						}
					}
                   //����Ҫ�������㴦��Ҫ���ձ任���������Ϸ���ͶӰ��
					pnodes[k+1] = project2_elli_surf(nodes_vec[nei_num[k]],thecell->surfaces_vec[win],&error[k+1]);
                                             
					if( error[k+1] == 0 )
					{
						//ͶӰʧ�ܣ�����dis0һ���㹻���������ֹ�ж��ڽڵ�ʱ��ѡ��˽ڵ�
						nei_dis[k] = 2.0*global_length;
					}
					else
					{
						nei_dis[k] = nodes_vec[nei_num[k]].distance_to(&pnodes[k+1]);
					}

					if( debuging ) hout << "       " << nei_num[k] << "        " << nei_dis[k] ;

				}

				if( debuging ) hout << endl;

				//����Ƿ�Ӧ���ƶ��ڽڵ㵽��������
				for( int k=0; k<3; k++ )
				{
					//������ڽڵ��뱾�ڵ���ͬһ�������Ծ�����
					if( nei_win[k] == win ) continue;

					 if( nei_dis[k] > dis0 || (movable_ns[nei_num[k]] == 0&&nei_dis[k]>min_dis) )
					{
						if(debuging) hout << "movable_ns["<<nei_num[k]<<"]: " << movable_ns[nei_num[k]]<< " nei_dis["<<k<<"]: " << nei_dis[k] << endl;
						move_nei = 0;		//ֻҪ��һ�������ƶ�����������ڵ����ƶ�Ҳ����
						break;
					}
				}

				if( debuging )
				{
					hout << "Project result: ";
					for( int k=0; k<4; k++ )
					{
						hout << "    " << error[k] ;
					}
					hout << endl;
					hout << "move_nei: " << move_nei << endl;
				}

				if( move_nei == 0 )
				{
					//��鱾�ڵ��Ƿ����ƶ����Ƿ��Ѿ��ƶ�����
					if( error[0] == 0 ) continue;				//ͶӰʧ�ܣ�����������������һЩ���߿�Խ����߽���
					if( movable_ns[nn] == 0 ) continue;	//�ƶ������������ƶ�

					if( dis0 > global_length/2.0 ) continue;	//�ƶ�����̫Զ�����ƶ�

					if( debuging )
					{
						hout << "Moved Node " << nn << "(1) to the surface." << endl;
						hout << "    Before: " << "      " << nodes_vec[nn].x ;
						hout <<  "      " << nodes_vec[nn].y ;
						hout <<  "      " << nodes_vec[nn].z ;
						hout << endl;
					}
                                        
					//�ƶ����ڵ�󣬱��ڵ�������߽����ϣ�����Ӧ�ü���is_on_ln������		
					is_on_ln.push_back(j);

					//�ƶ����ڵ�
					nodes_vec[nn].move_to(&pnodes[0]);
					is_on_nodes[nn] = 1;
					movable_ns[nn] = 0;          

					if( debuging )
					{
						hout << "    After : " << "      " << nodes_vec[nn].x ;
						hout <<  "     " << nodes_vec[nn].y ;
						hout <<  "     " << nodes_vec[nn].z ;
						hout << endl;
					}
				}
				else
				{
					for( int k=0; k<3; k++ )
					{
						if( debuging )
						{
							hout << "Moving No." << k << " neighbor (node " << nei_num[k] << ")" <<endl;
							hout << "    nei_win: " << nei_win[k] << " project error: " << error[k+1] ;
							hout << " movable: " << movable_ns[nei_num[k]] << endl;
						}
						//������ڽڵ��뱾�ڵ���ͬһ�������Ծ�����
						if( nei_win[k] == win ) continue;
						if( error[k+1] == 0 ) continue;   //ͶӰʧ�ܣ�����������������һЩ���߿�Խ����߽���
						//���ڵ��Ƿ����ƶ����Ƿ��Ѿ��ƶ�����
						if( movable_ns[nei_num[k]] == 0 )
						{
							int is_skip=1;
							//�Ѿ��ڱ�Ŀ������棬���ַǳ��������Ǻ���ǿ����Ϊ�Ѿ��ڱ߽����ˣ�������������������
							if( debuging )
							{
								hout << "nei_dis[" << k <<"]: " << nei_dis[k] ;
								hout << " min_dis: " << min_dis << endl;
							}
							if( nei_dis[k] == 0 )
							{
								//�Ѿ�����������
								//�����������������Ϊ���ʼ��ʱ��Ѹýڵ��ƶ������������ϣ�
								//�����ֱ�ǿ����Ϊ�ڱ�������ϣ���������Ҫ���ڸ������ϣ������鷳��
								continue;
							}
							else if( nei_dis[k] < min_dis )
							{
								//����������ǳ�����������ڱ�������ϣ���ǿ����ΪҲ�ڵ�ǰ�����ϣ�ճ��һ���ˣ�
								//������ڻ����ϣ������������ƶ�����ֹ�����ĸ��ڵ㶼���������ϵĵ�Ԫ�����õ��ģ��⻬����������״��)
								if( where_is_nodes[nei_num[k]] != -1 )
								{
									where_is_nodes[nei_num[k]] = win;
									is_on_nodes[nei_num[k]] = 1;
									continue;
								}
								else
								{
									is_skip = 0;
								}
							}
							if( is_skip == 1 ) continue;
						}
                                                                                             
						if( nei_dis[k] > global_length/2.0 ) continue;
                                                  
						if( debuging )
						{
							hout << "Moved Node " << nei_num[k] << "(2) to the surface." << endl;
							hout << "    Before: " << "      " << nodes_vec[nei_num[k]].x ;
							hout <<  "      " << nodes_vec[nei_num[k]].y ;
							hout <<  "      " << nodes_vec[nei_num[k]].z ;
							hout << endl;
						}

						is_on_ln.push_back(nei_i[k]);
                                                
						nodes_vec[nei_num[k]].move_to(&pnodes[k+1]);
						where_is_nodes[nei_num[k]] = win;
						is_on_nodes[nei_num[k]] = 1;

						movable_ns[nei_num[k]] = 0;
                                                
						if( debuging )
						{
							hout << "    After : " << "      " << nodes_vec[nei_num[k]].x ;
							hout << "      " << nodes_vec[nei_num[k]].y ;
							hout << "      " << nodes_vec[nei_num[k]].z ;
							hout << endl;
						}                                          
					}
				}
			}
		}

		//����������Ľڵ��Ƿ����ƶ���������������Ѿ���4���ڵ��ڽ����ϣ����нڵ㲻�����ƶ���
		if( is_on_ln.size() >= 4 )
		{
			for( int j=0; j<8; j++ )
			{
				movable_ns[hexes[i].nodesId[j]] = 0;
			}
		}

		//ȷ����ǰ������Ĳ��·��
		//����ȷ����Щ�ڵ��ڿ�����
		vector<int> wb;						//��Խ�Ŀ������
		vector<vector<int> > wnodes;	//����ÿ�������Ľڵ���
		wb.clear();
		wnodes.clear();
		for( int j=0; j<8; j++ )
		{
			int nn = hexes[i].nodesId[j];
			int win = where_is_nodes[nn];

			if( debuging )
			{
				hout << "checking node " << j << "(" << nn << ").............." << endl;
				hout << "    win: " << win << endl;
			}
                        
			if( win != -1 && is_on_nodes[nn] == 0 )
			{
				int nw=-1;
				for( int k=0; k<(int)wb.size(); k++ )
				{
					if( win == wb[k] )
					{
						nw = k;
						break;
					}
				}
				if( nw == -1 )
				{
					wb.push_back(win);
					vector<int> temp_v;
					wnodes.push_back(temp_v);
					wnodes.back().push_back(j);
				}
				else
				{
					wnodes[nw].push_back(j);
				}

				//������п��ܸ��˵㹹��������Ľڵ㣬
				//�������������ı��棬���鵽�ÿ����ľ��룬
				//�����С����ǿ����Ϊ�Ǵ˿��������ϵĵ㣬
				//��Ҫ����Ϊ�˷�ֹ���������ܽ�ʱ���ֺ�С�ĵ�Ԫ
				for( int k=0; k<8; k++ )
				{
					int nn1 = hexes[i].nodesId[k];
					int win1 = where_is_nodes[nn1];
					if( win1 != win && is_on_nodes[nn1] == 1 )
					{
						Point pp1 = {nodes_vec[nn1].x,nodes_vec[nn1].y,nodes_vec[nn1].z};
						double dis1 = thecell->surfaces_vec[win]->dis_to(&pp1);

						if(debuging)
						{
							hout << "    Node " << k << "(" << nn1 << "): ";
							hout << " win1: " << win1 << " dis1: " << dis1 << endl;
						}

						if( dis1 < ZERO ) continue; //�Ѿ��ڽ�������
						if( dis1 < min_dis )
						{
							//ǿ����Ϊ�˽ڵ���win������
							where_is_nodes[nn1]=win;
						}
					}
				}
			}
		}

		//��ǰ������Ĳ����ʽ���������б�߷���
		int face_style[6];
		if( debuging ) hout << "face_style initial: " ;
		for( int j=0; j<6; j++ )
		{
			int sn = hexes[i].facesId[j];
			face_style[j] = sign_faces[sn];
			if( debuging ) hout << "    " << sign_faces[sn];
		}
		if( debuging ) hout << endl;
                
		for( int j=0; j<(int)wb.size(); j++ )
		{
			if( debuging )
			{
				hout << "Belongs to particle: " << wb[j] << ": ";
				for( int k=0; k<(int)wnodes[j].size(); k++ )
				{
					hout << "    " << wnodes[j][k] ;
				}
				hout << endl;
			}

			if( wnodes[j].size() == 1 )
			{
				//ֻ��һ���ڵ���ͬһ������
				int ann = wnodes[j][0];
				deter_hex_split_style(ann,-1,face_style);
			}
			else if( wnodes[j].size() == 2 )
			{
				//ֻ�������ڵ���ͬһ������
				int edge_num = is_edge_of_hex(wnodes[j][0],wnodes[j][1]);
				deter_hex_split_style(-1,edge_num,face_style);
			}
			else if( wnodes[j].size() == 3 || wnodes[j].size() == 4 )
			{
				//ֻ��һ���ڵ��ڿ�����
				//�����ҳ��˽ڵ�(�п���û�д˽ڵ㣬�����ĸ��ڵ��ڽ����ϣ�
				int ann = -1;
				for( int k=0; k<8; k++ )
				{
					int nn = hexes[i].nodesId[k];
					int win = where_is_nodes[nn];
					if( win != wb[j] )
					{
						ann = k;
					}
				}
				if( ann != -1 )
				{
					deter_hex_split_style(ann,-1,face_style);
				}
			}
		}
		for( int k=0; k<6; k++ )
		{
			if( face_style[k] == -1 )
			{
				continue;
			}
			int sn = hexes[i].facesId[k];
			sign_faces[sn] = face_style[k];
		}
		if( debuging )
		{
			int show = 0;
			for( int j=0; j<8; j++ )
			{
				int nn = hexes[i].nodesId[j];
				int win = where_is_nodes[nn];
				if( win != -1 && is_on_nodes[nn] == 0 )
				{
					show=1;
					break;
				}
			}
			if( show == 1)
			{
				hout << "hex No. " << i << endl;
				hout << "    where is the nodes: " ;
				hout.setf(ios::right,ios::adjustfield);
				for( int j=0; j<8; j++ )
				{
					hout << "     " << where_is_nodes[hexes[i].nodesId[j]];
				}
				hout << endl;
				hout << "    is on the surface : " ;
				for( int j=0; j<8; j++ )
				{
					hout << "     " << is_on_nodes[hexes[i].nodesId[j]];
				}
				hout << endl;

				hout << "    surface signal    : " ;
				for( int j=0; j<6; j++ )
				{
					hout << "     " << hexes[i].facesId[j];
				}
				hout << endl;
				hout << "                        " ;
				for( int j=0; j<6; j++ )
				{
					hout << "     " << sign_faces[hexes[i].facesId[j]];
				}
				hout << endl;
			}
		}
	}

//	��֤���������ϵĵ㻹Ҫ�ڲ����ϣ������ϵĵ㻹Ҫ�ڱ����ϣ��ǵ㱣�ֲ���
	recover_bnodes(bnodes_vec);

	hout << "    ��ͶӰ��" << pv <<"�Σ�" ;
	hout << "    ʧ�ܣ�" << pl << "�Σ�" << endl;
/*
	string file="adjusted_hexes.dat";
	string oldfile=hout.output_file;
	ofstream cout_hex( file.c_str() ) ;
	if ( !cout_hex )  hout << "�޷����ļ���" << file << endl;
	cout_hex << "TITLE = draw" << endl;
	cout_hex << "VARIABLES = X, Y, Z" << endl;
	cout_hex << "ZONE N=" << (int)nodes_vec.size() <<",E=" << (int)hexes.size() << ", F=FEPOINT, ET=BRICK" << endl;

	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		cout_hex << nodes_vec[i].x <<"  " << nodes_vec[i].y <<"  " << nodes_vec[i].z << endl;
	}
	for( int i=0; i<(int)hexes.size(); i++ )
	{
		for( int j=0; j<8 ; j++ )
		{
			cout_hex << hexes[i].nodesId[j]+1 <<"  ";
		}
		cout_hex << endl;
	}	
	if( oldfile.size() > 0 ) open_deffo_stream( (char*)oldfile.c_str() , ios::app );
*/
	return 1;
}

//---------------------------------------------------------------------------
//ȷ�������ڵ��λ�ã����塢��n������
int Mesher::where_is( Node* node,  int &is_on, int node_count )
{
	int rvalue = -1;
	is_on = 0;
	Point thepoint={ node->x, node->y, node->z };
	int already_contained =-1;
	for( int i=0; i<(int)thecell->surfaces_vec.size(); i++)
	{
		int is_contain = thecell->surfaces_vec[i]->is_contain(&thepoint);

		if(is_contain < 0 )
		{
			rvalue = i;

			if( already_contained == -1 )
			{
				already_contained = i;
			}
			else
			{
				if( node_count >= 0 )
				{
					hout << node_count << "�Žڵ�ͬʱ��" << already_contained << "��" << rvalue << "�����������ڻ�����ϣ���where is����is_contain < 0�У� " << endl;
				}
				else
				{
					hout << "�����нڵ�ͬʱ��" << already_contained << "��" << rvalue << "�����������ڻ�����ϣ���where is����is_contain < 0�У� " << endl;
				}
				rvalue = -3;
			}
//			break;		//ֱ������
		}
		else if( is_contain == 0)
		{
			is_on = 1;
			rvalue = i;

			if( already_contained == -1 )
			{
				already_contained = i;
			}
			else
			{
				if( node_count >= 0 )
				{
					hout << node_count << "�Žڵ�ͬʱ��" << already_contained << "��" << rvalue << "�����������ڻ�����ϣ���where is����is_contain == 0�У� " << endl;
				}
				else
				{
					hout << "�����нڵ�ͬʱ��" << already_contained << "��" << rvalue << "�����������ڻ�����ϣ���where is����is_contain == 0�У� " << endl;
				}
				rvalue = -3;
			}
//			break;		//ֱ������
		}
	}
	return rvalue;
}

//---------------------------------------------------------------------------
//��������ڵ㵽������(sur)��ͶӰ�㣨ע�ⵥ�������ϵĽڵ�ͶӰҲҪ�ڲ����ϣ�
//��ŷ������ҵ��õ����������ĵ�������������Ľ��㣨project����������������Ϊ��ʼ�㣬
//						Ӧ���ݶȷ����������õ����������ϵķ���ͶӰ�㣬����ǲ��桢�����ϵ�
//						Ϊ�˱������ѣ����ܻ��÷���ͶӰ��任�����桢�����ϵĵ��뱾������������
//						������淨��ͶӰ��
//						���������ϵĵ㻹Ҫ�ڲ����ϣ������ϵĵ㻹Ҫ�ڱ����ϣ��ǵ㱣�ֲ���������
//						��adjust_hexes��������recover_bnodes(bnodes_vec)�Ӻ���ʵ��
//						ֻ��project_normal�������ʱ��ɢ�����*error==0
Node Mesher::project2_elli_surf(Node &node,Surface* sur, int *error)
{
        bool debuging = false;
//      if( error != NULL && *error == 2 ) debuging = true;
        if( debuging ) hout << "Node: " << node.x << " " << node.y << " " << node.z << " flag: " << node.flag << endl;


        int lerror;
        if( error != NULL ) lerror = *error;

        Point np = {node.x, node.y, node.z};
        Point rp = sur->project_normal(&np,&lerror);
        if( error != NULL ) *error = lerror;

        if( debuging )
		{
			hout << "First projection(" << lerror << "): " ;
			hout << rp.x << " " << rp.y << " " << rp.z ;
			hout << " ff: " << sur->fvalue(rp.x, rp.y, rp.z ) << endl;
		}

        if( node.flag < 0 ) return Node(rp);   //�ڲ��ڵ�

        Point rp1=rp;
        int lerror1=0;
        if( node.flag < 6 )
		{
                //�߽���ڵ�
                if( node.flag == 0 || node.flag == 1 )
				{
                        //x�����ϵĵ�
                        rp1.x = node.x;
                }
                else if( node.flag == 2 || node.flag == 3 )
				{
                        //y�����ϵĵ�
                        rp1.y = node.y;
                }
                else if( node.flag == 4 || node.flag == 5 )
				{
                        //z�����ϵĵ�
                        rp1.z = node.z;
                }
                TDVector pvec(np,rp1);
                rp1=sur->project(&np,&pvec,&lerror1);
        }
        else if( node.flag < 18 && node.flag >= 6 )
		{
                //�߽��߽ڵ�
                if( node.flag == 6 || node.flag == 8 || node.flag == 10 || node.flag == 12 )
				{
                        //x�᷽��߽����ϵĵ�
                        rp1.x = node.x + 1;
                        rp1.y = node.y;
                        rp1.z = node.z;
                }
                else if( node.flag == 7 || node.flag == 9 || node.flag == 11 || node.flag == 13 )
				{
                        //y�᷽��߽����ϵĵ�
                        rp1.x = node.x;
                        rp1.y = node.y + 1;
                        rp1.z = node.z;
                }
                else if( node.flag == 14 || node.flag == 15 || node.flag == 16 || node.flag == 17 )
				{
                        //z�᷽��߽����ϵĵ�
                        rp1.x = node.x;
                        rp1.y = node.y;
                        rp1.z = node.z + 1;
                }
                TDVector pvec(np,rp1);
                rp1=sur->project(&np,&pvec,&lerror1);
        }
             
        if( debuging )
		{
                hout << "Second projection(" << lerror1 << "): " ;
                hout << rp1.x << " " << rp1.y << " " << rp1.z << endl;
        }

        if( debuging ) hout << "lerror1: " << lerror1 << endl;

        if( lerror1 == 0 )
		{
                //�ڲ����ڣ���߽����ϣ�ͶӰʧ�ܣ����ط���ͶӰ
                return Node(rp);
        }
        else
		{
                //�ڲ����ڣ���߽����ϣ�ͶӰ�ɹ��������룬
                //�������ͶӰ�������Է���ͶӰ���Զ�����Ƿ��ط���ͶӰ
                //��Ҫ�Ƿ�ֹ�����Ͻڵ㼷��һ��
                double dis1 = rp.distance_to(np);
                double dis2 = rp1.distance_to(np);
                if( dis2 > dis1*2.0 )
				{
                        return Node(rp);
                }
                
                if( error != NULL ) *error = lerror1;
                return Node(rp1);
        }
        
        return Node(rp);
}

//---------------------------------------------------------------------------
//���ݸ����Ľڵ��ţ�������ľֲ���ţ���ȷ�����ڵ�3���ڵ�ţ���Ȼ�Ǿֲ���ţ�
int Mesher::deter_nei(int i, int *nei_i)
{
        switch(i)
		{
        case 0:
                nei_i[0]=1;
                nei_i[1]=4;
                nei_i[2]=3;
                break;
        case 1:
                nei_i[0]=0;
                nei_i[1]=2;
                nei_i[2]=5;
                break;
        case 2:
                nei_i[0]=1;
                nei_i[1]=3;
                nei_i[2]=6;
                break;
        case 3:
                nei_i[0]=0;
                nei_i[1]=7;
                nei_i[2]=2;
                break;
        case 4:
                nei_i[0]=0;
                nei_i[1]=5;
                nei_i[2]=7;
                break;
        case 5:
                nei_i[0]=1;
                nei_i[1]=6;
                nei_i[2]=4;
                break;
        case 6:
                nei_i[0]=5;
                nei_i[1]=2;
                nei_i[2]=7;
                break;
        case 7:
                nei_i[0]=4;
                nei_i[1]=6;
                nei_i[2]=3;
                break;
        default:
                break;
        }
        return 1;
}

//---------------------------------------------------------------------------
//���ݸ����Ľڵ�ţ����߱ߺţ�ֻ��ָ��һ������ȷ����������ʷ���ʽ�������ڵ㣨�ߣ��ڿ����ڵı�)
int Mesher::deter_hex_split_style(int node_num, int edge_num, int face_style[8])
{
	if( node_num != -1 )
	{
		if( node_num == 0 || node_num == 2 || node_num == 5 || node_num == 7 )
		{
			for( int k=0; k<6; k++ )
			{
				if( k%2 == 0 )
				{
					face_style[k] = 1;
				}
				else
				{
					face_style[k] = 0;
				}
			}
		}
		else
		{
			for( int k=0; k<6; k++ )
			{
				if( k%2 == 0 )
				{
					face_style[k] = 0;
				}
				else
				{
					face_style[k] = 1;
				}
			}
		}
		return 1;
	}
	if( edge_num != -1 )
	{
		if( edge_num == 0 || edge_num == 6 )
		{
			face_style[0] = 1;
			face_style[1] = 1;
		}
		else if( edge_num == 1 || edge_num == 7 )
		{
			face_style[2] = 0;
			face_style[3] = 0;
		}
		else if( edge_num == 2 || edge_num == 4 )
		{
			face_style[0] = 0;
			face_style[1] = 0;
		}
		else if( edge_num == 3 || edge_num == 5 )
		{
			face_style[2] = 1;
			face_style[3] = 1;
		}
		else if( edge_num == 8 || edge_num == 10 )
		{
			face_style[3] = 0;
			face_style[4] = 1;
			face_style[5] = 1;
		}
		else if( edge_num == 9 || edge_num == 11 )
		{
			face_style[3] = 1;
			face_style[4] = 0;
			face_style[5] = 0;
		}
		return 1;
	}
	return 0;
}

//---------------------------------------------------------------------------
//���������ڵ㣬�ж��Ƿ�Ϊ�������һ����
int Mesher::is_edge_of_hex(int n1, int n2)
{
        int min_n = min(n1,n2);
        int max_n = max(n1,n2);
        if( min_n == 0 )
		{
			if( max_n == 1 ) return 0;
			if( max_n == 3 ) return 3;
			if( max_n == 4 ) return 8;
        }
        else if( min_n == 1 )
		{
			if( max_n == 2 ) return 1;
			if( max_n == 5 ) return 9;
        }
        else if( min_n == 2 )
		{
			if( max_n == 3 ) return 2;
			if( max_n == 6 ) return 10;
        }
        else if( min_n == 3 )
		{
			if( max_n == 7 ) return 11;
        }
        else if( min_n == 4 )
		{
			if( max_n == 5 ) return 4;
			if( max_n == 7 ) return 7;
        }
        else if( min_n == 5 )
		{
			if( max_n == 6 ) return 5;
        }
        else if( min_n == 6 )
		{
			if( max_n == 7 ) return 6;
        }
        return -1;
}

//---------------------------------------------------------------------------
//�ָ��ƶ����ı߽�ڵ㣨���������ϵĵ㣩
int Mesher::recover_bnodes(vector<int> &bnodes)
{
        bool debuging = false;

        //�ָ��ƶ����ı߽�ڵ�
        for( int i=0; i<(int)bnodes.size(); i++ )
		{                                   
                Node vertexes[8];
                double diss[8];
                double dis_min = 1.0e200;
                int j,nn;

                Node *node = &nodes_vec[bnodes_vec[i]];
                if( debuging )
				{
                        hout << "------------------------------------------------------" << endl;
                        hout << "i: " << i << " node(" << bnodes[i] << "): " ;
                        hout << "     " << node->x ;
                        hout << "     " << node->y ;
                        hout << "     " << node->z ;
                        hout << "     " << node->flag ;
                        hout << endl;
                }
                switch( node->flag )
				{
                        case 0:
                                node->x = x_min;
                                break;
                        case 1:
                                node->x = x_max;
                                break;
                        case 2:
                                node->y = y_min;
                                break;
                        case 3:
                                node->y = y_max;
                                break;
                        case 4:
                                node->z = z_min;
                                break;
                        case 5:
                                node->z = z_max;
                                break;
                        case 6:
                                node->y = y_min;
                                node->z = z_min;
                                break;
                        case 7:
                                node->x = x_max;
                                node->z = z_min;
                                break;   
                        case 8:
                                node->y = y_max;
                                node->z = z_min;
                                break;    
                        case 9:
                                node->x = x_min;
                                node->z = z_min;
                                break;
                        case 10:
                                node->y = y_min;
                                node->z = z_max;
                                break;
                        case 11:
                                node->x = x_max;
                                node->z = z_max;
                                break;   
                        case 12:
                                node->y = y_max;
                                node->z = z_max;
                                break;    
                        case 13:
                                node->x = x_min;
                                node->z = z_max;
                                break;
                        case 14:
                                node->x = x_min;
                                node->y = y_min;
                                break;
                        case 15:
                                node->x = x_max;
                                node->y = y_min;
                                break;   
                        case 16:
                                node->x = x_max;
                                node->y = y_max;
                                break;    
                        case 17:
                                node->x = x_min;
                                node->y = y_max;
                                break;
                        case 18:
                                //�����˸��ǵ㣬�ƶ������������һ��
                                vertexes[0]=Node(x_min,y_min,z_min);
                                vertexes[1]=Node(x_max,y_min,z_min);
                                vertexes[2]=Node(x_max,y_max,z_min);
                                vertexes[3]=Node(x_min,y_max,z_min);
                                vertexes[4]=Node(x_min,y_min,z_max);
                                vertexes[5]=Node(x_max,y_min,z_max);
                                vertexes[6]=Node(x_max,y_max,z_max);
                                vertexes[7]=Node(x_min,y_max,z_max);
                                for( j=0; j<8; j++ )
								{
                                        diss[j] = node->distance_to(&vertexes[j]);
                                        if( diss[j] < dis_min )
										{
                                                nn = j;
                                                dis_min = diss[j];
                                        }
                                }
                                node->x = vertexes[nn].x;
                                node->y = vertexes[nn].y;
                                node->z = vertexes[nn].z;
                                break;
                        default:
                                break;
                }
                if( debuging )
				{
                        hout << "After recovery: " ;
                        hout << "      " << node->x ;
                        hout << "      " << node->y ;
                        hout << "      " << node->z ;
                        hout << "      " << node->flag ;
                        hout << endl;
                }
        }
        return 1;
}

//---------------------------------------------------------------------------
//����������������
int Mesher::split_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces)
{
	bool debuging = false;

	double min_dis=global_length/4.0;

	for( int i=0; i<(int)hexes.size(); i++ )
	{
		if( debuging )
		{
			hout << "-------------------------------------------------------" << endl;
			hout << "Hexahedron  " << "    " << i << "(" << "    " << (int)hexes.size() << "): " ;
			for( int j=0; j<8; j++ )
			{
				hout << "    " << hexes[i].nodesId[j] ;
			}
			hout << endl;
			hout << "where is nodes        : " ;
			for( int j=0; j<8; j++ )
			{
				hout << "      " << where_is_nodes[hexes[i].nodesId[j]] ;
			}
			hout << endl;
			hout << "in on face            : " ;
			for( int j=0; j<8; j++ )
			{
				hout << "      " << is_on_nodes[hexes[i].nodesId[j]] ;
			}
			hout << endl;
			hout << "surfaces signal       : " ;
			for( int j=0; j<6; j++ )
			{
				hout << "      " << hexes[i].facesId[j] ;
			}
			hout << endl;
			hout << "                        " ;
			for( int j=0; j<6; j++ )
			{
				hout << "      " << sign_faces[hexes[i].facesId[j]] ;
			}
			hout << endl;
		}

		int hex_info[14];
		for( int j=0; j<8; j++ )
		{
			hex_info[j] = hexes[i].nodesId[j];
		}

		int face_style[6];
		for( int j=0; j<6; j++ )
		{
			face_style[j] = sign_faces[hexes[i].facesId[j]];
		}

		//�޸��ı�����Ƭ��б�߷���
		//�ȼ���Ƿ����ֵ㵽б�ߵľ�����������������취����
		//������֤������ֵ㵽б�ߵľ�����������     
		for( int j=0; j<6; j++ )
		{
			int nnseq[4];
			if( face_style[j] != -1 )
			{
				//б�߷����Ѿ�ȷ������Ȼ���
				deter_nn_seq(j,face_style[j],nnseq);
				double dis1 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[2]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[0]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[1]]] );
				double dis2 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[3]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[0]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[1]]] );

				if( dis1 < min_dis || dis2 < min_dis )
				{
					face_style[j] = (face_style[j]+1)%2;
				}
			}
			else
			{
				//б�߷�����δȷ����ȷ���������������
				deter_nn_seq(j,0,nnseq);
				double dis1 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[2]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[0]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[1]]] );
				double dis2 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[3]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[0]]],
													 &nodes_vec[hexes[i].nodesId[nnseq[1]]] );

				if( dis1 < min_dis || dis2 < min_dis )
				{
					face_style[j] = 1;
					continue;
				}

				deter_nn_seq(j,1,nnseq);
				dis1 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[2]]],
										  &nodes_vec[hexes[i].nodesId[nnseq[0]]],
										  &nodes_vec[hexes[i].nodesId[nnseq[1]]] );
				dis2 = dis2_line( &nodes_vec[hexes[i].nodesId[nnseq[3]]],
										  &nodes_vec[hexes[i].nodesId[nnseq[0]]],
										  &nodes_vec[hexes[i].nodesId[nnseq[1]]] );

				if( dis1 < min_dis || dis2 < min_dis )
				{
					face_style[j] = 0;
				}                         
			}
		}

		if( debuging )
		{
			hout << "surfaces signal (new) : " ;
			for( int j=0; j<6; j++ )
			{
				hout << "      " << face_style[j] ;
			}
			hout << endl;
		}

		for( int j=0; j<6; j++ )
		{
			hex_info[j+8] = face_style[j];
		}

		vector<Tetrahedron > tets;
		if( best_split_hex2_tet(hex_info,tets) == 0 )
		{
/*
			hout << "���������(" << i << ")��������ʧ�ܣ���ֳ�12��������!!" << endl;
			hout << "  hex_info: " ;
			for( int i=0; i<14; i++ )
			{
				hout << hex_info[i] << " " ;
			}
			hout  << endl;
			hout  << "  Created elements: " << (int)eles_vec.size()
				<< " : " << (int)eles_vec.size() + 11 << endl;
*/
			if( hex2_12tet(hex_info,tets) == 0 ) return 0;
		}
		//���¾�����Ƭб�߱�־
		for( int j=0; j<6; j++ )
		{
			int fn = hexes[i].facesId[j];
			if( sign_faces[fn] == -1 )
			{
				sign_faces[fn] = hex_info[j+8];
			}
		}

		if(debuging) hout << "splited into elements: " << endl;
		for( int j=0; j<(int)tets.size(); j++ )
		{
			int ele_num = (int)eles_vec.size();
			eles_vec.push_back(tets[j]);
			nodes_vec[tets[j].nodesId[0]].relative_eles_vec.push_back(ele_num);
			nodes_vec[tets[j].nodesId[1]].relative_eles_vec.push_back(ele_num);
			nodes_vec[tets[j].nodesId[2]].relative_eles_vec.push_back(ele_num);
			nodes_vec[tets[j].nodesId[3]].relative_eles_vec.push_back(ele_num);

			if( debuging )
			{
				hout << "Element id: " << ele_num << " Nodes: " ;
				for(int k=0; k<4; k++ )
				{
					hout << "      " << tets[j].nodesId[k];
				}
				hout << endl;
			}
		}
	}
	return 1;        
}

//---------------------------------------------------------------------------
//���ݸ����Ĳ����ţ�������Ĳ���ֲ���ţ����Լ�б�߷���ȷ��һ���ڵ�������
//�������ݣ�б�߽ڵ�1��б�߽ڵ�2����б�߽ڵ�1����б�߽ڵ�2
//Ҳ���ǲ����ı�����η������Σ���֤nnseq[4]��0,1,2�Žڵ��0,1,3�Žڵ�����һ�������Σ�
int Mesher::deter_nn_seq(int fn, int style, int nnseq[4])
{
	if( fn == 0 )
	{
		nnseq[0] = 0;
		nnseq[1] = 7;
		nnseq[2] = 3;
		nnseq[3] = 4;
	}
	else if( fn == 1 )
	{
		nnseq[0] = 1;
		nnseq[1] = 6;
		nnseq[2] = 2;
		nnseq[3] = 5;
	}
	else if( fn == 2 )
	{
		nnseq[0] = 0;
		nnseq[1] = 5;
		nnseq[2] = 1;
		nnseq[3] = 4;
	}
	else if( fn == 3 )
	{
		nnseq[0] = 3;
		nnseq[1] = 6;
		nnseq[2] = 2;
		nnseq[3] = 7;
	}
	else if( fn == 4 )
	{
		nnseq[0] = 0;
		nnseq[1] = 2;
		nnseq[2] = 1;
		nnseq[3] = 3;
	}
	else if( fn == 5 )
	{
		nnseq[0] = 4;
		nnseq[1] = 6;
		nnseq[2] = 5;
		nnseq[3] = 7;
	}
	if( style == 1 )
	{
		int temp_int = nnseq[0];
		nnseq[0] = nnseq[2];
		nnseq[2] = nnseq[1];
		nnseq[1] = nnseq[3];
		nnseq[3] = temp_int;
	}
	return 1;
}

//---------------------------------------------------------------------------
//����һ�㵽�����������ߵľ���
double Mesher::dis2_line( Node *n1, Node *n2, Node *n3 )
{
	//�ɵ�p1�Լ�p2ָ��p3����������ƽ��(p2ָ��p3��������ƽ��ķ�����)
	//���ƽ����p2��p3�����ߵĽ��㣬Ȼ����p1��˽���ľ���
	TDVector vec23(n3->x - n2->x, n3->y - n2->y, n3->z - n2->z);
	TDVector vec21(n1->x - n2->x, n1->y - n2->y, n1->z - n2->z);

	double len=vec23.length();
	if( len == 0 ) return 0.0;

	double tt = vec21.dot_product(&vec23)/(len*len);
	Node n4 = *n3-*n2;
	Node n5 = n4*tt;
	n5 = n5 + *n2;

	return n5.distance_to(n1);
}

//---------------------------------------------------------------------------
//�������еĲ����������ṩһ���ܲ�ֳ�6��������ķ�����ֻ�޸�split_info�е���-1��Ԫ�أ�
int Mesher::best_split_hex2_tet(int* hex_info,vector<Tetrahedron> &tets)
{
	bool debuging = false;
	if( debuging )
	{
		hout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		hout << "split info: " ;
		for( int i=0; i<6; i++ )
		{
			hout << "    " << hex_info[i+8] ;
		}
		hout << endl;
	}
	for( int i=0; i<6; i++ )
	{
		if( hex_info[i+8] == -1 )
		{
			//����ѡ�������б���෴�ķ���
			int ci ;
			if( i%2 == 0 )
			{
				ci = i + 1;
			}
			else
			{
				ci = i - 1;
			}
			if( hex_info[ci+8] != -1 )
			{
				hex_info[i+8] = (hex_info[ci+8]+1)%2;
			}
			else
			{
				hex_info[i+8] = 0;
			}
			if( best_split_hex2_tet(hex_info,tets) == 1 ) return 1;
			hex_info[i+8] = (hex_info[i+8]+1)%2;
			return best_split_hex2_tet(hex_info,tets);
		}
	}
	return hex2_6tet(hex_info,tets);
}

//---------------------------------------------------------------------------
//��һ��������ֳ�6�������壨�ַ�����������ľ������������
int Mesher::hex2_6tet(int hex[14], vector<Tetrahedron > &tets,double volume)
{
	bool debuging = false;
        
	tets.clear();
          
	if(debuging)
	{                
		hout << "==========================================" << endl;
		hout << "hex: " ;
		for( int i=0; i<14; i++ )
		{
			hout << hex[i] << " " ;
		}
		hout << endl;
	}
	vector<Tetrahedron> local_tets;  //�Ծֲ��ڵ��ű�ʾ��������

	//�����Խ���
	//1��1��7
	//2��2��8
	//3��3��5
	//4��4��6
	int diag_line = 0;
	int vertexes[8]={0,0,0,0,0,0,0,0};

	//���12��������Ƭ
	vector<Tri_2D_ele> tri_faces;
	tri_faces.reserve(20);

	//��־12���������Ƿ��Ѿ���������Ѿ��γ��������壩
	int sign_tri_face[12]={0,0,0,0,0,0,0,0,0,0,0,0};

	//12�������ڵ�������Ƭ��ÿ��������������12��2�����飬��һ��Ԫ��ָ������Ƭ��ţ��ڶ���Ԫ��ָ����Ľڵ�ţ�
	deter_rela_tri_face(hex,tri_faces);

	vector<vector<int> > temp_v1;
	vector<vector<vector<int> > > node_for_face(12,temp_v1);  //��¼ÿ�����Ը�������Ƭ����������Ľڵ�

	int sign_diag_line = 0;       //��¼�ڴ����ĸ�������Ƭʱ����ĶԽ���

	//��¼�Ƿ��Ѿ�������������Ƭ�Ŀ��ýڵ㣨����֮����������Ľڵ㣩
	int sign_dealed_tri_faces[12];
	for( int i=0; i<12; i++ )
	{
		sign_dealed_tri_faces[i] = 0;
	}

	int succeeded = 1;
	vector<int> i_vec; //������һ����������������Ƭ���Ա�����ʱʹ��
	for( int i=0; i<12&&i>=0; i++)
	{
		if( sign_tri_face[i] > 0 ) continue;		//��ʾ�Ѿ���������
		if(debuging)
		{
			hout << "-------------------------------------" << endl;
			hout << "tri face " << i+1 << endl;
			hout << "sign_dealed_tri_faces: " << sign_dealed_tri_faces[i] << endl;
		}
		//�ռ������i��������Ƭ����������Ľڵ�
		if( sign_dealed_tri_faces[i] == 0 )
		{
			i_vec.push_back(i);  //��¼���������ýڵ��������Ƭ��ţ��Ա�����ʱʹ��
			if( debuging )
			{
				hout << "vertexes signal: " ;
				for( int j=0; j<8; j++ )
				{
					hout << vertexes[j] << " " ;
				}
				hout << endl;
				hout << "Surfaces signal: " ;
				for( int j=0; j<8; j++ )
				{
					hout << sign_dealed_tri_faces[j] << " " ;
				}
				hout << endl;
			}
			//ȷ����Ҫ����4���ڵ㣨������Ƭ������ĸ��ڵ㣬��������Ľڵ�������
			int check_nn[4];
			deter_check_nodes(i,tri_faces[i].nodesId[1],check_nn);
			for( int jj=0; jj<4; jj++ )
			{
				int j=check_nn[jj];
				int ls[3];
				if( vertexes[j] == 0 )
				{
					//���Э���ԣ������б߲��ܽ��棩
					int cf1t = can_form_1tet(tri_faces[i],j,hex,diag_line,ls);
					//�����Ҫ�õ���������Ƭ�Ƿ��Ѿ�ʹ�ù�
					if( ls[0] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[0],tri_faces[i].edgesId[0]);
						if( fn < i ) continue;  //��Ž�С��������Ƭһ���Ѿ��������
						if( sign_tri_face[fn] > 0 ) continue;
					}
					else if( ls[1] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[1],tri_faces[i].edgesId[0]);
						if( fn < i ) continue;  //��Ž�С��������Ƭһ���Ѿ��������
						if( sign_tri_face[fn] > 0 ) continue;

						fn = deter_face_num_from_2bar(ls[1],tri_faces[i].edgesId[1]);
						if( fn < i ) continue;  //��Ž�С��������Ƭһ���Ѿ��������
						if( sign_tri_face[fn] > 0 ) continue;

					}
					else if( ls[2] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[2],tri_faces[i].edgesId[1]);
						if( fn < i ) continue;  //��Ž�С��������Ƭһ���Ѿ��������
						if( sign_tri_face[fn] > 0 ) continue;
					}

					if( cf1t == 1 )
					{
						//���Ƕȣ�����������Ƭ�ļнǲ��ܹ���
						Tri_2D_ele act_tri(	hex[tri_faces[i].nodesId[0]],
														hex[tri_faces[i].nodesId[1]],
														hex[tri_faces[i].nodesId[2]]	);

						double angle = angle_point2_tri(jj,&act_tri,hex,check_nn);

						if( angle > 160.0*PI/180.0 ) continue;

						vector<int> temp_v;
						temp_v.push_back(j);
						temp_v.push_back(ls[0]);
						temp_v.push_back(ls[1]);
						temp_v.push_back(ls[2]);
						node_for_face[i].push_back(temp_v);
					}
				}
			}
			sign_dealed_tri_faces[i]=1;
		}

		if(debuging)
		{
			hout << "The nodes can be used to form a tet: " ;
			hout << (int)node_for_face[i].size() << endl;
			for( int j=0; j<(int)node_for_face[i].size(); j++ )
			{
				hout << "     " << node_for_face[i][j][0] ;
				hout << "     " << node_for_face[i][j][1] ;
				hout << "     " << node_for_face[i][j][2] ;
				hout << "     " << node_for_face[i][j][3] ;
				hout << endl;
			}
		}

		if( node_for_face[i].size() > 0 )
		{
			vector<int> act_node = node_for_face[i].back();
			node_for_face[i].pop_back();
			//����һ��������
			local_tets.push_back(Tetrahedron(tri_faces[i],act_node[0]));
			//��������±���������Ϊб�ߣ�����������ȡ�������һ���ǣ�����Ӧ��ǽǵ㣬��ֹ�´�������
			vector<int> tf;
			int tn = 0;

			//�Ƿ���ֱ��(���һ����
			int sn = -1;
			for( int j=0; j<3; j++ )
			{
				if( act_node[j+1] < 12 )
				{
					sn = act_node[j+1];
				}
				else if( act_node[j+1] < 18 )
				{
					tf.push_back(act_node[j+1]-12);   //��¼б�߱��
				}

				else if( act_node[j+1] == 18 )
				{
					diag_line = min(act_node[0],tri_faces[i].nodesId[j])+1;
					sign_diag_line = i+1;
					if(debuging)
					{
						hout << "diag_line: " << diag_line << endl;
					}
				}
			}
			if( tf.size() >1 && sign_diag_line != i+1 )
			{
				vertexes[tri_faces[i].nodesId[1]] = i+1;
				if( debuging )
				{
					hout << "Sign the node NO. " << tri_faces[i].nodesId[1] << " for used." << endl;
				}
			}
			//����ֱ�ߣ���Ҫ����Ѿ��������������Ƭ������������ظ���Ԫ
			if( sn > -1 )
			{
				vector<int> rela_tri_faces;
				for( int j=0; j<(int)tf.size(); j++ )
				{
					int fn = tf[j]*2;
					for( int k=0; k<3; k++ )
					{
						if( tri_faces[fn].edgesId[k] == sn )
						{
							sign_tri_face[fn] = i+1;
							if(debuging)
							{
								hout << "Sign the surface NO. " << fn << " for used." << endl;
							}
						}
					}
					fn += 1;
					for( int k=0; k<3; k++ )
					{
						if( tri_faces[fn].edgesId[k] == sn )
						{
							sign_tri_face[fn] = i+1;
							if(debuging)
							{
								hout << "Sign the surface NO. " << fn << " for used." << endl;
							}
						}
					}
				}
			}
		}
		else
		{
			if( debuging )
			{
				hout << "Not pass. Back." << endl;
			}
			//����ʧ�ܣ��������һ�������ı�ǲ����˵���һ��
			sign_dealed_tri_faces[i]=0;
			i_vec.pop_back();
			if(i_vec.size() == 0)
			{
				i=-1;
			}
			else 
			{
				i=i_vec[i_vec.size()-1];
			}
			for( int j=0; j<12; j++ )
			{
				if( sign_tri_face[j] == i+1 )
				{
					sign_tri_face[j]=0;
					if(debuging)
					{
						hout << "Unmark the surface NO. " << j << " for used." << endl;
					}
				}
			}
			if( i == -1 )
			{
				succeeded = 0;
				break;
			}
			if( sign_diag_line == i+1 )
			{
				diag_line = 0;
				sign_diag_line = 0;
				if(debuging)
				{
					hout << "Clear the diag_line signal." << endl;
				}
			}
			for( int j=0; j<8; j++ )
			{
				if( vertexes[j] == i+1 )
				{
					vertexes[j] = 0;
					if( debuging )
					{
						hout << "Unmark the node NO. " << j << " for used." << endl;
					}
				}
			}
			i--;
			local_tets.pop_back();
		}
	}

	if( succeeded == 0 )
	{
		return 0;
	}
	else
	{
		if(local_tets.size() == 4)
		{
			if( hex[8] == 0 )
			{
				local_tets.push_back(Tetrahedron(0,2,7,5));
			}
			else
			{
				local_tets.push_back(Tetrahedron(1,3,4,6));
			}
		}
	}

	for( int i=0; i<(int)local_tets.size(); i++ )
	{
		tets.push_back(Tetrahedron(	hex[local_tets[i].nodesId[0]],
													hex[local_tets[i].nodesId[1]],
													hex[local_tets[i].nodesId[2]],
													hex[local_tets[i].nodesId[3]])	);
	}

	return 1;
}

//---------------------------------------------------------------------------
//�������������������壨�����壩��6����ֳ�������Ƭ����ȷ�����໥��ϵ��ÿ���������ߵĵ�Ԫ��
int Mesher::deter_rela_tri_face(int *hex,vector<Tri_2D_ele> &tri_faces)
{
	Tri_2D_ele temp_tri;
	tri_faces.assign(12,temp_tri);
	//����12��������Ƭ
	//���еĽڵ����Լ��ߵı�ž�Ϊ���ֶ���Ĵָ�����������ڲ�����б�������
	if( hex[8] == 0 )
	{
		tri_faces[0] = Tri_2D_ele(7,4,0,7,8,12);
		tri_faces[1] = Tri_2D_ele(0,3,7,3,11,12);
	}
	else
	{
		tri_faces[0] = Tri_2D_ele(4,0,3,8,3,12);
		tri_faces[1] = Tri_2D_ele(3,7,4,11,7,12);
	}

	if( hex[9] == 0 )
	{
		tri_faces[2] = Tri_2D_ele(1,5,6,9,5,13);
		tri_faces[3] = Tri_2D_ele(6,2,1,10,1,13);
	}
	else
	{
		tri_faces[2] = Tri_2D_ele(2,1,5,1,9,13);
		tri_faces[3] = Tri_2D_ele(5,6,2,5,10,13);
	}      
	if( hex[10] == 0 )
	{
		tri_faces[4] = Tri_2D_ele(0,4,5,8,4,14);
		tri_faces[5] = Tri_2D_ele(5,1,0,9,0,14);
	}
	else
	{
		tri_faces[4] = Tri_2D_ele(1,0,4,0,8,14);
		tri_faces[5] = Tri_2D_ele(4,5,1,4,9,14);
	}

	if( hex[11] == 0 )
	{
		tri_faces[6] = Tri_2D_ele(6,7,3,6,11,15);
		tri_faces[7] = Tri_2D_ele(3,2,6,2,10,15);
	}
	else
	{
		tri_faces[6] = Tri_2D_ele(7,3,2,11,2,15);
		tri_faces[7] = Tri_2D_ele(2,6,7,10,6,15);
	}

	if( hex[12] == 0 )
	{
		tri_faces[8] = Tri_2D_ele(2,3,0,2,3,16);
		tri_faces[9] = Tri_2D_ele(0,1,2,0,1,16);
	}
	else
	{
		tri_faces[8] = Tri_2D_ele(3,0,1,3,0,16);
		tri_faces[9] = Tri_2D_ele(1,2,3,1,2,16);
	}

	if( hex[13] == 0 )
	{
		tri_faces[10] = Tri_2D_ele(4,7,6,7,6,17);
		tri_faces[11] = Tri_2D_ele(6,5,4,5,4,17);
	}
	else
	{
		tri_faces[10] = Tri_2D_ele(5,4,7,4,7,17);
		tri_faces[11] = Tri_2D_ele(7,6,5,6,5,17);
	}

	return 1;
}

//---------------------------------------------------------------------------
//ȷ����Ҫ����4���ڵ㣨������Ƭ������ĸ��ڵ㣬��������Ľڵ�������
        //�ڵ��˳��Ϊ�����Խ��߽ڵ㡢��һ���߶�Ӧ��б�߽ڵ㡢�ڶ����߶�Ӧ��б�߽ڵ㡢�е��Ӧ�Ľڵ㣨����ֱ�ߵĽڵ㣩
//i��������Ƭ�ı��
//n��������Ƭ�ڶ����ڵ�ı�ţ��ֲ���ţ�
int Mesher::deter_check_nodes(int i,int n,int check_nn[4])
{
	switch(i)
	{
	case 0:
		if( n==0 )
		{
			check_nn[0] = 6;
			check_nn[1] = 5;
			check_nn[2] = 2;
			check_nn[3] = 1;
		}
		else
		{              
			check_nn[0] = 2;
			check_nn[1] = 6;
			check_nn[2] = 1;
			check_nn[3] = 5;
		}
		break;
	case 1:
		if( n==3 )
		{         
			check_nn[0] = 5;
			check_nn[1] = 1;
			check_nn[2] = 6;
			check_nn[3] = 2;
		}
		else
		{
			check_nn[0] = 1;
			check_nn[1] = 2;
			check_nn[2] = 5;
			check_nn[3] = 6;
		}
		break;
	case 2:
		if( n==1 )
		{
			check_nn[0] = 7;
			check_nn[1] = 3;
			check_nn[2] = 4;
			check_nn[3] = 0;
		}
		else
		{
			check_nn[0] = 3;
			check_nn[1] = 0;
			check_nn[2] = 7;
			check_nn[3] = 4;
		}
		break;
	case 3:
		if( n==2 )
		{
			check_nn[0] = 4;
			check_nn[1] = 7;
			check_nn[2] = 0;
			check_nn[3] = 3;
		}
		else
		{
			check_nn[0] = 0;
			check_nn[1] = 4;
			check_nn[2] = 3;
			check_nn[3] = 7;
		}
		break;
	case 4:
		if( n==0 )
		{     
			check_nn[0] = 6;
			check_nn[1] = 2;
			check_nn[2] = 7;
			check_nn[3] = 3;
		}
		else
		{
			check_nn[0] = 2;
			check_nn[1] = 3;
			check_nn[2] = 6;
			check_nn[3] = 7;
		}
		break;
	case 5:
		if( n==1 )
		{
			check_nn[0] = 7;
			check_nn[1] = 6;
			check_nn[2] = 3;
			check_nn[3] = 2;
		}
		else
		{
			check_nn[0] = 3;
			check_nn[1] = 7;
			check_nn[2] = 2;
			check_nn[3] = 6;
		}
		break;
	case 6:
		if( n==3 )
		{
			check_nn[0] = 5;
			check_nn[1] = 4;
			check_nn[2] = 1;
			check_nn[3] = 0;
		}
		else
		{
			check_nn[0] = 1;
			check_nn[1] = 5;
			check_nn[2] = 0;
			check_nn[3] = 4;
		}
		break;
	case 7:
		if( n==2 )
		{
			check_nn[0] = 4;
			check_nn[1] = 0;
			check_nn[2] = 5;
			check_nn[3] = 1;
		}
		else
		{
			check_nn[0] = 0;
			check_nn[1] = 1;
			check_nn[2] = 4;
			check_nn[3] = 5;
		}
		break;
	case 8:
		if( n==0 )
		{
			check_nn[0] = 6;
			check_nn[1] = 7;
			check_nn[2] = 5;
			check_nn[3] = 4;
		}
		else
		{
			check_nn[0] = 5;
			check_nn[1] = 6;
			check_nn[2] = 4;
			check_nn[3] = 7;
		}
		break;
	case 9:
		if( n==1 )
		{       
			check_nn[0] = 7;
			check_nn[1] = 4;
			check_nn[2] = 6;
			check_nn[3] = 5;
		}
		else
		{
			check_nn[0] = 4;
			check_nn[1] = 5;
			check_nn[2] = 7;
			check_nn[3] = 6;
		}
		break;
	case 10:
		if( n==4 )
		{           
			check_nn[0] = 2;
			check_nn[1] = 1;
			check_nn[2] = 3;
			check_nn[3] = 0;
		}
		else
		{
			check_nn[0] = 1;
			check_nn[1] = 0;
			check_nn[2] = 2;
			check_nn[3] = 3;
		}
		break;
	case 11:
		if( n==5 )
		{
			check_nn[0] = 3;
			check_nn[1] = 2;
			check_nn[2] = 0;
			check_nn[3] = 1;
		}
		else
		{
			check_nn[0] = 0;
			check_nn[1] = 3;
			check_nn[2] = 1;
			check_nn[3] = 2;
		}
		break;
	default:
		break;
	}
	return 1;
}

//---------------------------------------------------------------------------
//����һ��������Ƭ��һ���ڵ㣬�ж������Ƿ��ܹ���һ�������壬��Ҫ���Ǳ߽�Э������
//ls[3]: �����±ߵı�ţ����������еı�ţ�
int Mesher::can_form_1tet(Tri_2D_ele &tri, int nn, int *hex, int &diag_line, int ls[3])
{
	bool debuging=false;
	//�ж������±��Ƿ�����б߽�Э�������нڵ��ž�Ϊ�ֲ���ţ�
	for( int i=0; i<3; i++ )
	{
		ls[i] = -1;
		int n0 = tri.nodesId[i];
		if( debuging ) hout << "n0: " << n0 << " nn: " << nn <<endl;
		//�Ƿ��ظ�
		if( n0 == nn ) return 0;
		//�Ƿ�Ϊֱ��
		int ieoh = is_edge_of_hex(n0,nn);
		if( ieoh >= 0 )
		{
			ls[i] = ieoh;
			continue;
		}
		//�Ƿ�Ϊб�ߣ�����ǣ�Ҫ����Ƿ�Э��
		int style;
		int is_dia_edge = deter_face_num(n0,nn,&style);
		if( debuging ) hout << "is dia edge: " << is_dia_edge << endl;
		if( is_dia_edge >= 0 )
		{
			if( hex[8+is_dia_edge] == style )
			{
				ls[i] = is_dia_edge+12;
				continue;
			}
			else
			{ 
				return 0;
			}
		}
		//�Ƿ�Ϊ�Խ���
		int min_n=min(n0,nn);
		if( diag_line == 0 )
		{
			ls[i] = 18;
			continue;
		}
		else
		{
			if( min_n+1 == diag_line )
			{
				ls[i] = 18;
				continue;
			}
			else
			{
				return 0;
			}
		}
	}
	if( ls[0] >= 0 && ls[1] >= 0 && ls[2] >= 0 ) return 1;
	return 0;
}

//---------------------------------------------------------------------------
//���ݸ����������ڵ�ţ��ж�������������ĸ�������
int Mesher::deter_face_num(int n1, int n2, int *style)
{
	int min_n=min(n1,n2);
	int max_n=max(n1,n2);
	if( min_n == 0 )
	{
		if( max_n == 7 )
		{
			if( style != NULL ) *style = 0;
			return 0;
		}
		else if( max_n == 5 )
		{
			if( style != NULL ) *style = 0;
			return 2;
		}
		else if( max_n == 2 )
		{
			if( style != NULL ) *style = 0;
			return 4;
		}
	}
	else if( min_n == 1 )
	{
		if( max_n == 6 )
		{
			if( style != NULL ) *style = 0;
			return 1;
		}
		else if( max_n == 4 )
		{
			if( style != NULL ) *style = 1;
			return 2;
		}
		else if( max_n == 3 )
		{
			if( style != NULL ) *style = 1;
			return 4;
		}
	}
	else if( min_n == 2 )
	{
		if( max_n == 5 )
		{
			if( style != NULL ) *style = 1;
			return 1;
		}
		else if( max_n == 7 )
		{
			if( style != NULL ) *style = 1;
			return 3;
		}
	}
	else if( min_n == 3 )
	{
		if( max_n == 4 )
		{
			if( style != NULL ) *style = 1;
			return 0;
		}
		else if( max_n == 6 )
		{
			if( style != NULL ) *style = 0;
			return 3;
		}
	}
	else if( min_n == 4 )
	{
		if( max_n == 6 )
		{
			if( style != NULL ) *style = 0;
			return 5;
		}
	}
	else if( min_n == 5 )
	{
		if( max_n == 7 )
		{
			if( style != NULL ) *style = 1;
			return 5;
		}
	}
	return -1;
}

//---------------------------------------------------------------------------
//��������ֱ�ߺţ����ع��ɵ�������Ƭ��
int Mesher::deter_face_num_from_2bar(int bn1, int bn2)
{
	int fn = -1;
	int min_bn = min(bn1,bn2);
	int max_bn = max(bn1,bn2);
	if( min_bn == 0 )
	{
		if( max_bn == 1 )
		{
			fn = 9;
		}
		else if( max_bn == 3 )
		{
			fn = 8;
		}
		else if( max_bn == 8 )
		{
			fn = 4;
		}
		else if( max_bn == 9 )
		{
			fn = 5;
		}
	}
	else if( min_bn == 1 )
	{
		if( max_bn == 2 )
		{
			fn = 9;
		}
		else if( max_bn == 9 )
		{
			fn = 2;
		}
		else if( max_bn == 10 )
		{
			fn = 3;
		}
	}
	else if( min_bn == 2 )
	{
		if( max_bn == 3 )
		{
			fn = 8;
		}
		else if( max_bn == 10 )
		{
			fn = 7;
		}
		else if( max_bn == 11 )
		{
			fn = 6;
		}
	} 
	else if( min_bn == 3 )
	{
		if( max_bn == 8 )
		{
			fn = 0;
		}
		else if( max_bn == 11 )
		{
			fn = 1;
		}
	}
	else if( min_bn == 4 )
	{
		if( max_bn == 5 )
		{
			fn = 11;
		}
		else if( max_bn == 7 )
		{
			fn = 10;
		}
		else if( max_bn == 8 )
		{
			fn = 4;
		}
		else if( max_bn == 9 )
		{
			fn = 5;
		}
	}
	else if( min_bn == 5 )
	{
		if( max_bn == 6 )
		{
			fn = 11;
		}
		else if( max_bn == 9 )
		{
			fn = 2;
		}
		else if( max_bn == 10 )
		{
			fn = 3;
		}
	}
	else if( min_bn == 6 )
	{
		if( max_bn == 7 )
		{
			fn = 10;
		}
		else if( max_bn == 10 )
		{
			fn = 7;
		}
		else if( max_bn == 11 )
		{
			fn = 6;
		}
	}  
	else if( min_bn == 7 )
	{
		if( max_bn == 8 )
		{
			fn = 0;
		}
		else if( max_bn == 11 )
		{
			fn = 1;
		}
	}
	return fn;
}

//---------------------------------------------------------------------------
//��������ڵ㵽��������Ƭ��ļн�(����������6��������ʱʹ�ã�
double Mesher::angle_point2_tri(int j, Tri_2D_ele* tri, int *hex, int *check_nn)
{
	bool debuging = false;
	if( debuging )
	{
		hout << "----------------------------------------------------" << endl;
		hout << "j: " << j << endl;
		hout << "tri      : " ;
		for( int i=0; i<3; i++ )
		{
			hout << "     " << tri->nodesId[i] ;
		}
		hout << endl;
		hout << "Hex      : " ;
		for( int i=0; i<14; i++ )
		{
			hout << "     " << hex[i] ;
		}
		hout << endl;
		hout << "check_nn : " ;
		for( int i=0; i<4; i++ )
		{
			hout << "     " << check_nn[i] ;
		}
		hout << endl;
	}
	if( j==0 )
	{
		int node1_num = hex[check_nn[j]];
		int node2_num = tri->nodesId[0];
		int node3_num = tri->nodesId[2];
		if( debuging )
		{
			hout << "Checked tri: ";
			hout << "     " << node1_num ;
			hout << "     " << node2_num ;
			hout << "     " << node3_num ;
			hout << endl;
		}
		Tri_2D_ele tri1(node1_num,node2_num,node3_num);
		return angle_of_2tri(tri,&tri1);
	}
	else if( j==1 )
	{
		int node1_num = hex[check_nn[j]];
		int node2_num = tri->nodesId[1];
		int node3_num = tri->nodesId[0];
		Tri_2D_ele tri1(node1_num,node2_num,node3_num);
		return angle_of_2tri(tri,&tri1);
	} 
	else if( j==2 )
	{
		int node1_num = hex[check_nn[j]];
		int node2_num = tri->nodesId[2];
		int node3_num = tri->nodesId[1];
		Tri_2D_ele tri1(node1_num,node2_num,node3_num);
		return angle_of_2tri(tri,&tri1);
	}
	else if( j==3 )
	{
		int node1_num = hex[check_nn[j]];
		int node2_num = tri->nodesId[1];
		int node3_num = tri->nodesId[0];
		Tri_2D_ele tri1(node1_num,node2_num,node3_num);
		double angle1 = angle_of_2tri(tri,&tri1);
		node1_num = hex[check_nn[j]];
		node2_num = tri->nodesId[2];
		node3_num = tri->nodesId[1];
		Tri_2D_ele tri2(node1_num,node2_num,node3_num);
		double angle2 = angle_of_2tri(tri,&tri2);
		return max(angle1,angle2);
	}
	return PI;
}

//---------------------------------------------------------------------------
//��������������Ƭ֮��ļн�,node1_num node2_num �ֱ�Ϊ���ڹ������ϵĽڵ�
double Mesher::angle_of_2tri(Tri_2D_ele* tri1, Tri_2D_ele* tri2, int node1_num, int node2_num)
{
	//���tri1 tri2����Ч�ԣ������нڵ��ظ���
	int nn10 = tri1->nodesId[0];
	int nn20 = tri2->nodesId[0];
	for( int i=1; i<3; i++ )
	{
		if( tri1->nodesId[i] == nn10 )
		{
			hout << "Error occurs while trying to compute the angle between two tri face." << endl;
			return 0.0;
		}
		if( tri2->nodesId[i] == nn20 )
		{
			hout << "Error occurs while trying to compute the angle between two tri face." << endl;
			return 0.0;
		}
	}

	//����tri1�ķ�����
	Node node21 = nodes_vec[tri1->nodesId[1]] - nodes_vec[tri1->nodesId[0]] ;
	TDVector v11 = TDVector( node21.x, node21.y, node21.z );
	Node node31 = nodes_vec[tri1->nodesId[2]] - nodes_vec[tri1->nodesId[0]] ;
	TDVector v12 = TDVector( node31.x, node31.y, node31.z );
	TDVector normal1 = v11.cro_product( &v12 );
	//����tri2�ķ�����
	node21 = nodes_vec[tri2->nodesId[1]] - nodes_vec[tri2->nodesId[0]] ;
	TDVector v21 = TDVector( node21.x, node21.y, node21.z );
	node31 = nodes_vec[tri2->nodesId[2]] - nodes_vec[tri2->nodesId[0]] ;
	TDVector v22 = TDVector( node31.x, node31.y, node31.z );
	TDVector normal2 = v21.cro_product( &v22 );

	//���������������ļн�
	double acos_v = normal1.dot_product( &normal2 ) / ( normal1.length() * normal2.length() );
	double theter1 = acos( acos_v );
	//����Ƭ�н�Ϊ����нǵ����
	double theter2 = PI - theter1 ;

	//�ж�����Ƭ���Ӵ��ǰ��λ�͹��
	//�ҳ����ڹ������ϵ�����
	Node node1, node2;
	int is_coplane = 1;
	for ( int i=0; i<3; i++ )
	{
		if( node1_num != -1 )
		{
			node1 = nodes_vec[node1_num] ;
		}
		else
		{
			node1 = nodes_vec[tri1->nodesId[i]] ;
		}
		if ( fabs(	(node1.x - nodes_vec[tri2->nodesId[0]].x)*normal2.x +
						(node1.y - nodes_vec[tri2->nodesId[0]].y)*normal2.y +
						(node1.z - nodes_vec[tri2->nodesId[0]].z)*normal2.z		) > 1e-8 )
		{
			is_coplane = 0;
			break ;
		}
		if( node1_num != -1 ) 
		{
			break ;
		}
	}
	if( is_coplane )
	{
		return theter2 ;
	}
	for ( int i=0; i<3; i++ )
	{
		if( node2_num != -1 )
		{
			node2 = nodes_vec[node2_num] ;
		}
		else
		{
			node2 = nodes_vec[tri2->nodesId[i]] ;
		}
		if ( fabs(	(node2.x - nodes_vec[tri1->nodesId[0]].x)*normal1.x +
						(node2.y - nodes_vec[tri1->nodesId[0]].y)*normal1.y +
						(node2.z - nodes_vec[tri1->nodesId[0]].z)*normal1.z		) > 1e-8 )
		{
			break;
		}
		if( node2_num != -1 )
		{
			break ;
		}
	}
	TDVector v3 = TDVector( node1.x - node2.x, node1.y - node2.y, node1.z - node2.z );

	acos_v = v3.dot_product( &normal2 ) / ( v3.length() * normal2.length() );
	if( fabs(acos_v) > 1.0+ZERO )
	{
		hout << "acos_v: " << acos_v << endl;
	}
	theter1 = acos( acos_v );

	if ( theter1 < PI/2.0 )
	{
		return theter2;
	}
	else
	{
		return 2.0 * PI - theter2 ;
	}
}

//---------------------------------------------------------------------------
//��һ��������ֳ�12�������壨�����Ĳ���һ���½ڵ㣬���������ܳɹ��ֳ�6��������������
int Mesher::hex2_12tet(int hex[14], vector<Tetrahedron > &tets, double volume)
{
	tets.clear();

	//     vector<Tetrahedron> local_tets;
	//�ֳ�12�������壨�����������Ĳ���һ���ڵ㣩
	Node cen_node = nodes_vec[hex[0]];
	for( int i=1; i<8; i++ )
	{
		cen_node = cen_node + nodes_vec[hex[i]];
	}
	cen_node = cen_node/8.0;
	int node_num = (int)nodes_vec.size();
	nodes_vec.push_back(cen_node); 
	int is_on;
	int win = where_is( &cen_node, is_on );
	if( win == -3 ) return 0;
	where_is_nodes.push_back(win);
	is_on_nodes.push_back(is_on);

	Tri_2D_ele tri_faces[12];
	if( hex[12] == 0 )
	{
		tri_faces[1] = Tri_2D_ele(hex[2],hex[3],hex[0]);
		tri_faces[0] = Tri_2D_ele(hex[0],hex[1],hex[2]);
	}
	else
	{
		tri_faces[1] = Tri_2D_ele(hex[3],hex[0],hex[1]);
		tri_faces[0] = Tri_2D_ele(hex[1],hex[2],hex[3]);
	}
	if( hex[13] == 0 )
	{
		tri_faces[2] = Tri_2D_ele(hex[6],hex[5],hex[4]);
		tri_faces[3] = Tri_2D_ele(hex[4],hex[7],hex[6]);
	}
	else
	{
		tri_faces[2] = Tri_2D_ele(hex[7],hex[6],hex[5]);
		tri_faces[3] = Tri_2D_ele(hex[5],hex[4],hex[7]);
	}
	if( hex[10] == 0 )
	{
		tri_faces[4] = Tri_2D_ele(hex[0],hex[4],hex[5]);
		tri_faces[5] = Tri_2D_ele(hex[5],hex[1],hex[0]);
	}
	else
	{
		tri_faces[4] = Tri_2D_ele(hex[1],hex[0],hex[4]);
		tri_faces[5] = Tri_2D_ele(hex[4],hex[5],hex[1]);
	}
	if( hex[9] == 0 )
	{
		tri_faces[6] = Tri_2D_ele(hex[1],hex[5],hex[6]);
		tri_faces[7] = Tri_2D_ele(hex[6],hex[2],hex[1]);
	}
	else
	{
		tri_faces[6] = Tri_2D_ele(hex[2],hex[1],hex[5]);
		tri_faces[7] = Tri_2D_ele(hex[5],hex[6],hex[2]);
	}
	if( hex[11] == 0 )
	{
		tri_faces[9] = Tri_2D_ele(hex[3],hex[2],hex[6]);
		tri_faces[8] = Tri_2D_ele(hex[6],hex[7],hex[3]);
	}
	else
	{
		tri_faces[9] = Tri_2D_ele(hex[2],hex[6],hex[7]);
		tri_faces[8] = Tri_2D_ele(hex[7],hex[3],hex[2]);
	}
	if( hex[8] == 0 )
	{
		tri_faces[10] = Tri_2D_ele(hex[0],hex[3],hex[7]);
		tri_faces[11] = Tri_2D_ele(hex[7],hex[4],hex[0]);
	}
	else
	{
		tri_faces[10] = Tri_2D_ele(hex[3],hex[7],hex[4]);
		tri_faces[11] = Tri_2D_ele(hex[4],hex[0],hex[3]);
	}
	for( int i=0; i<12; i++ )
	{
		tets.push_back(Tetrahedron(tri_faces[i],node_num));
	}

	return 1;
}

//---------------------------------------------------------------------------
//ȷ��ÿ��������Ĳ���
int Mesher::deter_eles_mat()
{
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		deter_mat_ele(i);
	}
	return 1;
}

//---------------------------------------------------------------------------
//ȷ��һ����Ԫ�Ĳ��ϱ��
int Mesher::deter_mat_ele(int ele_num)
{
	bool debuging = false;
	int is_deleted = 1;
	for( int i=1; i<4; i++ )
	{
		if( eles_vec[ele_num].nodesId[i] == eles_vec[ele_num].nodesId[0] )
		{
			return 1;
		}
	}

	if( debuging ) hout << "--------------------deter_mat_ele(int ele_num)---------------------" << endl;
	if( debuging ) hout << "element " << ele_num << "\n";

	int rv = deter_mat_ele(&eles_vec[ele_num],int(debuging));
	int itt = is_thin_tet(ele_num);
	if( debuging ) hout << "itt: " << itt << endl;
	if( itt == 1 )
	{
		thin_tets.push_back(ele_num);
	}
	if( debuging ) hout << "material ID. " << eles_vec[ele_num].materialId << endl;
	if( debuging ) hout << "rv: " << rv << endl;   
	if( debuging ) hout << "-----------------------------------------" << endl;
	return rv;
}

//---------------------------------------------------------------------------
//ȷ��һ����Ԫ�Ĳ��ϱ��
//����ֵ��-1: �����ϵĵ�Ԫ
//					0,1 �����壬��ǿ��
//					-3	: ����Ӧ���˳�����
int Mesher::deter_mat_ele(Tetrahedron* act_tet, int mod)
{
	bool debuging = false;
	if( mod == 1 ) debuging = true;
	int wws[4],ions[4]; 
	for( int i=0; i<4; i++ )
	{
		wws[i] = where_is_nodes[act_tet->nodesId[i]];
		ions[i] = is_on_nodes[act_tet->nodesId[i]];
	}

	if( debuging )
	{
		hout << "nodes of element " << ": ";
		hout << act_tet->nodesId[0] << " " ;
		hout << act_tet->nodesId[1] << " " ;
		hout << act_tet->nodesId[2] << " " ;
		hout << act_tet->nodesId[3] << endl;
		hout << "where is the nodes: " ;
		hout << wws[0] << " " ;
		hout << wws[1] << " " ;
		hout << wws[2] << " " ;
		hout << wws[3] << endl;
		hout << "is on the nodes: " ;
		hout << ions[0] << " " ;
		hout << ions[1] << " " ;
		hout << ions[2] << " " ;
		hout << ions[3] << endl;
	}

	//ȷ���˵�Ԫ�Ĳ���
	//�ڻ�����
	if( (wws[0] == -1 || ions[0] == 1) && (wws[1] == -1 || ions[1] == 1) &&
		(wws[2] == -1 || ions[2] == 1) && (wws[3] == -1 || ions[3] == 1) &&
		((wws[0] == -1 || wws[1] == -1 || wws[2] == -1 || wws[3] == -1)||
		(!(wws[0] == wws[1] && wws[0] == wws[2] && wws[0] == wws[3]))))
	{
		//�ĸ��ڵ㶼�ڻ����ڻ��߽����ϣ���������һ���ڵ��ڻ����ڻ���
		//�ĸ��ڵ��ڲ�ͬ�Ľ����ϣ���˵���˵�ԪΪ����

		//news��ע������������ɳ�����������ܳ��ֿ����ص�����
		//���������ڵĵ�Ԫ���ܻ�������е�is_on����1�����
		int is_cross_tet=-1;
		for( int i=0; i<4; i++ )
		{
			int lis_on;
			int lwh = where_is(&nodes_vec[act_tet->nodesId[i]],lis_on);
			if( lwh == -1 )
			{
				is_cross_tet = -1;
				break;
			}
			if( lis_on == 0 )
			{
				if( debuging )
				{
					hout << "    node " << act_tet->nodesId[i] << " is not on." << endl;
				}
				is_cross_tet = lwh;
			}
		}  
		if( is_cross_tet != -1 )
		{
			int mat_num = thecell->surfaces_vec[is_cross_tet]->material_num;
			act_tet->materialId = mat_num;
			return mat_num;
		}
		else
		{
			act_tet->materialId = 0 ;
			return 0;
		}
	}
	//��Ԫ����ǿ�����ڲ�
	if(wws[0] != -1 && wws[1] != -1 && wws[2] != -1 && wws[3] !=-1 )
	{

		if( wws[0] == wws[1] && wws[0] == wws[2] && wws[0] == wws[3] )
		{
			int mat_num = thecell->surfaces_vec[wws[0]]->material_num ;
			act_tet->materialId = mat_num ;
			return mat_num;
		}

		Point nodesp[4];
		for( int i=0; i<4; i++ )
		{
			Node *np = &nodes_vec[act_tet->nodesId[i]];
			Point pp = {np->x,np->y,np->z};
			nodesp[i] = pp;
		}
		int win = -1;
		for( int i=0; i<4; i++ )
		{
			win = wws[i];
			int is_same_par=1;
			for( int j=0; j<4; j++ )
			{
				if( wws[j] != win && thecell->surfaces_vec[win]->is_contain(&nodesp[j])>0 )
				{
					is_same_par = 0;
					break;
				}
			}
			if( is_same_par == 1 ) break;
			win = -1;
		}
		if( win != -1 )
		{
			int mat_num = thecell->surfaces_vec[win]->material_num ;		
			act_tet->materialId = mat_num ;
			return mat_num;
		}

	}
	//�ڽ����ϵĵ�Ԫ
	act_tet->materialId = -1 ;
//	act_tet->materialId = 0 ;	//�����ڻ����ϵ������
	return -1;
}

//---------------------------------------------------------------------------
//�ж�һ����Ԫ�Ƿ�Ϊ��Ԫ���ڵ㵽������������Ƭ�ľ���С��global_length/5.0);
int Mesher::is_thin_tet(int ele_num)
{
	bool debuging = false;

	if( debuging )
	{
		hout << "-------------------is_thin_tet()-----------------------------" << endl;
		hout << "ele_num: " << ele_num << " " ;
		hout << "  nodes: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "     " << eles_vec[ele_num].nodesId[k]  ;
		}
		hout << " mat: " << eles_vec[ele_num].materialId ;
		hout <<endl;
		hout << "  where: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "     " << where_is_nodes[eles_vec[ele_num].nodesId[k]]  ;
		}
		hout << endl;
		hout << "  is on: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "     " << is_on_nodes[eles_vec[ele_num].nodesId[k]]  ;
		}
		hout << endl;
	}

	//��鵥Ԫ��Ч�ԣ��Ƿ��Ѿ�ɾ�����ĸ��ڵ���ͬ����
	int nn0 = eles_vec[ele_num].nodesId[0];
	for( int i=1; i<4; i++ )
	{
		if( eles_vec[ele_num].nodesId[i] == nn0 )
		{
			return 0;
		}
	}

	int is_grov_tet=1;
	int n0 = eles_vec[ele_num].nodesId[0];
	for( int i=0; i<4; i++ )
	{
		int nn = eles_vec[ele_num].nodesId[i];
		if( where_is_nodes[nn] != where_is_nodes[n0] || is_on_nodes[nn] == 0 )
		{
			is_grov_tet = 0;
			break;
		}
	}
	if( debuging ) hout << " is_grov_tet: " << is_grov_tet << endl;
	if( is_grov_tet == 1 ) return 1;

	double angle_limit = 160.0*PI/180.0;
	double len_limit   = global_length/2.0;

	int nseq[4];
	double max_angle = max_angle_of_tet(ele_num,nseq);
	if( debuging ) hout << " max_angle: " << max_angle*180./PI << endl;
	if( max_angle >= angle_limit )
	{
		double dis1 = nodes_vec[nseq[0]].distance_to(&nodes_vec[nseq[3]]);

		if( dis1 > len_limit )
		{
			return 1;
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
//�ҳ�һ���������ĸ�������нǣ���������Ӧ�Ľڵ����У�0-1-2  3-2-1Ϊ�н��������棩
double Mesher::max_angle_of_tet(int ele_num, int *nseq)
{
	int nns[6][4];   //����ڵ���
	//�����������������������������Ľ��ǣ�
	//��д�ɣ��鷳
	nns[0][0] = 0;
	nns[0][1] = 1;
	nns[0][2] = 2;
	nns[0][3] = 3;
	nns[1][0] = 0;
	nns[1][1] = 3;
	nns[1][2] = 1;
	nns[1][3] = 2;
	nns[2][0] = 0;
	nns[2][1] = 2;
	nns[2][2] = 3;
	nns[2][3] = 1;
	nns[3][0] = 1;
	nns[3][1] = 2;
	nns[3][2] = 0;
	nns[3][3] = 3;
	nns[4][0] = 1;
	nns[4][1] = 0;
	nns[4][2] = 3;
	nns[4][3] = 2;
	nns[5][0] = 2;
	nns[5][1] = 0;
	nns[5][2] = 1;
	nns[5][3] = 3;

	int nsn[4];        //�ĸ��ڵ���
	//�ҳ��Ƕ�����������
	double max_ang=0.0;
	int max_aj = -1;
	for( int j=0; j<6; j++ )
	{
		for( int k=0; k<4; k++ )
		{
			nsn[k] = eles_vec[ele_num].nodesId[nns[j][k]];
		}

		Tri_2D_ele tri1(nsn[0],nsn[1],nsn[2]);
		Tri_2D_ele tri2(nsn[3],nsn[2],nsn[1]);
		double angle1 = angle_of_2tri( &tri1, &tri2, nsn[0], nsn[3]);

		if( angle1 > max_ang )
		{
			max_ang = angle1;
			max_aj = j;
		}
	}
	if( nseq != NULL )
	{
		for( int i=0; i<4; i++ )
		{
			nseq[i] = eles_vec[ele_num].nodesId[nns[max_aj][i]];
		}
	}
	return max_ang;
}

//---------------------------------------------------------------------------
//��ָ����Ԫ�������������ϣ���߽�ĵ�Ԫ��
int Mesher::put_into_cbev( int ele_num, double alength, int type )
{
	bool debuging = false;

	Tetrahedron *act_tet = &eles_vec[ele_num];
	//�����ҳ���һ�����ڵ����߽����ϵĵ�Ԫ
	if( type < 1 && has_cross_edge_edge( act_tet ) != 0 )
	{
		if( debuging ) hout << " Putting ele " << ele_num << " into be." <<endl;
		eles_across_boundary_be.push_back( ele_num );
		return 1;
	}

	//�ҳ���һ�����ڵ����߽����ϵĵ�Ԫ
	if( type < 2 )
	{
		int hceeof = has_cross_edge_edge_on_face( act_tet );
		if( hceeof >= 1 )
		{
			if( debuging ) hout << " Putting ele " << ele_num << " into bf." <<endl;
			eles_across_boundary_bf.push_back( ele_num );
			return 2;
		}
	}

	if( type < 3 )
	{
		int left_movable, right_movable;
		int edge_num = deter_oper_edge( act_tet, left_movable, right_movable );
		if( edge_num != 0 )
		{
			if( left_movable == 1 && right_movable == 1)
			{    
				if( debuging ) hout << " Putting ele " << ele_num << " into dmb." <<endl;

				eles_across_boundary_dmb.push_back(ele_num);
				return 3 ;
			}
			else 
			{                   
				if( debuging ) hout << " Putting ele " << ele_num << " into ndmb." <<endl;
				eles_across_boundary_ndmb.push_back(ele_num);
				return 4 ;
			}
		}
	}

	if( debuging ) hout << " Putting ele " << ele_num << " into ndmb." <<endl;
	eles_across_boundary_ndmb.push_back(ele_num);
	return 4;
}

//---------------------------------------------------------------------------
//�ж�һ����Ԫ�Ƿ��б��ڵ����߽����Ͽ�߽�
int Mesher::has_cross_edge_edge( Tetrahedron *act_tet )
{
	int w1 = where_is_nodes[act_tet->nodesId[0]];
	int w2 = where_is_nodes[act_tet->nodesId[1]];
	int w3 = where_is_nodes[act_tet->nodesId[2]];
	int w4 = where_is_nodes[act_tet->nodesId[3]];
	int is_on1 = is_on_nodes[act_tet->nodesId[0]];
	int is_on2 = is_on_nodes[act_tet->nodesId[1]];
	int is_on3 = is_on_nodes[act_tet->nodesId[2]];
	int is_on4 = is_on_nodes[act_tet->nodesId[3]];

	int act_ioe = is_on_edge( act_tet );

	if( act_ioe != 0 )
	{
		Node *n1 = &nodes_vec[act_tet->nodesId[0]];
		Node *n2 = &nodes_vec[act_tet->nodesId[1]];
		Node *n3 = &nodes_vec[act_tet->nodesId[2]];
		Node *n4 = &nodes_vec[act_tet->nodesId[3]];
		Point p1={n1->x, n1->y,n1->z};
		Point p2={n2->x, n2->y,n2->z};
		Point p3={n3->x, n3->y,n3->z};
		Point p4={n4->x, n4->y,n4->z};
		switch( act_ioe )
		{
		case 1:                                      
			if( w1 != w2 )
			{
				bool is_same_par1= w1!=-1 && is_on1==0 && thecell->surfaces_vec[w1]->is_contain(&p2)>0;
				bool is_same_par2= w2!=-1 && is_on2==0 && thecell->surfaces_vec[w2]->is_contain(&p1)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 1;
				}
			}
			break;
		case 2:                                        
			if( w1 != w3 )
			{
				bool is_same_par1= w1!=-1 && is_on1==0 && thecell->surfaces_vec[w1]->is_contain(&p3)>0;
				bool is_same_par2= w3!=-1 && is_on3==0 && thecell->surfaces_vec[w3]->is_contain(&p1)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 2;
				}
			}
			break;
		case 3:                                       
			if( w1 != w4 )
			{
				bool is_same_par1= w1!=-1 && is_on1==0 && thecell->surfaces_vec[w1]->is_contain(&p4)>0;
				bool is_same_par2= w4!=-1 && is_on4==0 && thecell->surfaces_vec[w4]->is_contain(&p1)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 3;
				}
			}
			break;
		case 4:                                      
			if( w2 != w3 )
			{
				bool is_same_par1= w2!=-1 && is_on2==0 && thecell->surfaces_vec[w2]->is_contain(&p3)>0;
				bool is_same_par2= w3!=-1 && is_on3==0 && thecell->surfaces_vec[w3]->is_contain(&p2)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 4;
				}
			}
			break;
		case 5:                                      
			if( w2 != w4 )
			{
				bool is_same_par1= w2!=-1 && is_on2==0 && thecell->surfaces_vec[w2]->is_contain(&p4)>0;
				bool is_same_par2= w4!=-1 && is_on4==0 && thecell->surfaces_vec[w4]->is_contain(&p2)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 5;
				}
			}
			break;
		case 6:                                      
			if( w3 != w4 )
			{
				bool is_same_par1= w3!=-1 && is_on3==0 && thecell->surfaces_vec[w3]->is_contain(&p4)>0;
				bool is_same_par2= w4!=-1 && is_on4==0 && thecell->surfaces_vec[w4]->is_contain(&p3)>0;
				if( is_same_par1 || is_same_par2 )
				{
					return 6;
				}
			}
			break;
		default:
			break;
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
int Mesher::is_on_edge( int node_num, int mod )
{
	Node *node = &nodes_vec[node_num];
	//���ڵ��Ƿ�ͬʱ���������ϣ�����ǣ���˵���˵������潻����
	int ofs = 0 ;
	if( fabs(node->x - x_min) <= ZERO )
	{
		ofs ++ ;
	}
	else if( fabs(node->x - x_max) <= ZERO )
	{
		ofs ++ ;
	}
	if( fabs(node->y - y_min) <= ZERO )
	{
		ofs ++ ;
	}
	else if( fabs(node->y - y_max) <= ZERO )
	{
		ofs ++ ;
	}
	if( fabs(node->z - z_min) <= ZERO )
	{
		ofs ++ ;
	}
	else if( fabs(node->z - z_max) <= ZERO )
	{
		ofs ++ ;
	}

	//������
	if( mod == 1)
	{
		if( is_on_nodes[node_num] == 1 ) ofs ++ ;
		int temp_int;
		int wh_node = where_is(node,temp_int);
		if( wh_node != where_is_nodes[node_num] ) ofs ++;
	}

	if( ofs > 1 ) return 1;
	return 0;
}

//---------------------------------------------------------------------------
//����������Ƿ��б��ڱ߽���
//����ֵΪ�ڱ߽����ϵıߺ�
//      0: û�б��ڱ߽�����
//      1: 0-1 ��
//      2: 0-2 ��
//      3: 0-3 ��
//      4: 1-2 ��
//      5: 1-3 ��
//      6: 2-3 ��
int Mesher::is_on_edge(Tetrahedron* act_tet, int mod)
{
	int ioe1 = is_on_same_edge( act_tet->nodesId[0],
												act_tet->nodesId[1], 0 ) ;
	if( ioe1 != 0 ) return 1;

	int ioe2 = is_on_same_edge( act_tet->nodesId[0],
												act_tet->nodesId[2], 0 ) ;
	if( ioe2 != 0 ) return 2;

	int ioe3 = is_on_same_edge( act_tet->nodesId[0],
												act_tet->nodesId[3], 0 ) ;
	if( ioe3 != 0 ) return 3;

	int ioe4 = is_on_same_edge( act_tet->nodesId[1],
												act_tet->nodesId[2], 0 ) ;
	if( ioe4 != 0 ) return 4;

	int ioe5 = is_on_same_edge( act_tet->nodesId[1],
												act_tet->nodesId[3], 0 ) ;
	if( ioe5 != 0 ) return 5;

	int ioe6 = is_on_same_edge( act_tet->nodesId[2],
												act_tet->nodesId[3], 0 ) ;
	if( ioe6 != 0 ) return 6;

	return 0;
}

//---------------------------------------------------------------------------
//�ж������ڵ��Ƿ���ͬһ���߽��ϣ�mod����Ҫ��Ҫ�����ǿ�����ͻ��彻����
int Mesher::is_on_same_edge(int node1_num, int node2_num, int mod)
{
        Node *node1 = &nodes_vec[node1_num];
        Node *node2 = &nodes_vec[node2_num];

        //����ȷ�������ڵ�ֱ�����Щ���ϣ���Щ��Ľ���
        vector<int> faces1, faces2 ;
        if( fabs(node1->x - x_min) <= ZERO ) faces1.push_back(-1) ;
        else if( fabs(node1->x - x_max) <= ZERO ) faces1.push_back(-2) ;
        if( fabs(node1->y - y_min) <= ZERO ) faces1.push_back(-3) ;
        else if( fabs(node1->y - y_max) <= ZERO ) faces1.push_back(-4) ;
        if( fabs(node1->z - z_min) <= ZERO ) faces1.push_back(-5) ;
        else if( fabs(node1->z - z_max) <= ZERO ) faces1.push_back(-6) ;
        if( mod == 1 )
		{
                int ww = where_is_nodes[node1_num] ;
                if( is_on_nodes[node1_num] == 1 ) faces1.push_back(ww + 1) ;
        }

        if( fabs(node2->x - x_min) <= ZERO ) faces2.push_back(-1) ;
        else if( fabs(node2->x - x_max) <= ZERO ) faces2.push_back(-2) ;
        if( fabs(node2->y - y_min) <= ZERO ) faces2.push_back(-3) ;
        else if( fabs(node2->y - y_max) <= ZERO ) faces2.push_back(-4) ;
        if( fabs(node2->z - z_min) <= ZERO ) faces2.push_back(-5) ;
        else if( fabs(node2->z - z_max) <= ZERO ) faces2.push_back(-6) ;
        if( mod == 1 )
		{       
                int ww = where_is_nodes[node2_num] ;
                if( is_on_nodes[node2_num] == 1 ) faces2.push_back(ww + 1) ;
        }
        //ֻҪ����������ͬ����˵����������ͬһ�߽�����
        int sf = 0 ;
        for( int i=0; i<(int)faces1.size(); i++ )
		{
                for( int j=0; j<(int)faces2.size(); j++ )
				{
                        if( faces1[i] == faces2[j] )
						{
                                sf ++ ;
                                if( sf > 1 ) return 1;
                        }
                }
        }
        return 0;
}

//---------------------------------------------------------------------------
//�ж�һ����Ԫ�Ƿ��б��ڵ����߽����Ͽ�߽�
int Mesher::has_cross_edge_edge_on_face( Tetrahedron *act_tet )
{
	int wh[4],is_on[4];
	for( int i=0; i<4; i++ )
	{
		wh[i] = where_is_nodes[act_tet->nodesId[i]];
		is_on[i] = is_on_nodes[act_tet->nodesId[i]];
	}
	Node *n1 = &nodes_vec[act_tet->nodesId[0]];
	Node *n2 = &nodes_vec[act_tet->nodesId[1]];
	Node *n3 = &nodes_vec[act_tet->nodesId[2]];
	Node *n4 = &nodes_vec[act_tet->nodesId[3]];
	Point pp[]={ {n1->x, n1->y,n1->z},
						{n2->x, n2->y,n2->z},
						{n3->x, n3->y,n3->z},
						{n4->x, n4->y,n4->z} };

	int edge_num=1;
	for( int i=0; i<3; i++ )
	{
		for( int j=i+1; j<4; j++,edge_num++ )
		{
			//���ÿ�����Ƿ��ڵ����߽����ϣ����ҿ����
			int iof = is_on_same_face( act_tet->nodesId[i],
													 act_tet->nodesId[j], 0 );
			if( iof == 0 ) continue;

			if( wh[i] != wh[j] )
			{
				bool is_same_par1= wh[i]!=-1 && is_on[i]==0 && thecell->surfaces_vec[wh[i]]->is_contain(&pp[j])>0;
				bool is_same_par2= wh[j]!=-1 && is_on[j]==0 && thecell->surfaces_vec[wh[j]]->is_contain(&pp[i])>0;
				if( is_same_par1 || is_same_par2 )
				{
					return edge_num;
				}
			}
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
//�����������ڵ��Ƿ���ͬһ�����ϣ������߽�����߻�������ǿ��Ľ��棩
//mod=0ֻ��鵥���߽��棬
//mod=1��鵥���߽���ͻ�������ǿ��Ľ���
int Mesher::is_on_same_face(int node1_num, int node2_num, int mod)
{
	Node *node1 = &nodes_vec[node1_num];
	Node *node2 = &nodes_vec[node2_num];
	if( fabs(node1->x - x_min) <= ZERO &&
		fabs(node2->x - x_min) <= ZERO ) return 1;
	else if(	fabs(node1->x - x_max) <= ZERO &&
				fabs(node2->x - x_max) <= ZERO ) return 1;
	else if(	fabs(node1->y - y_min) <= ZERO &&
				fabs(node2->y - y_min) <= ZERO ) return 1;
	else if(	fabs(node1->y - y_max) <= ZERO &&
				fabs(node2->y - y_max) <= ZERO ) return 1;
	else if(	fabs(node1->z - z_min) <= ZERO &&
				fabs(node2->z - z_min) <= ZERO ) return 1;
	else if(	fabs(node1->z - z_max) <= ZERO &&
				fabs(node2->z - z_max) <= ZERO ) return 1;

	if( mod == 1 )
	{
		Point pp = { node2->x, node2->y, node2->z };
		int ww1, is_on1;
		ww1 = where_is_nodes[node1_num];
		is_on1 = is_on_nodes[node1_num];

		if( ww1 != -1 && is_on1 == 1 )
		{
			int ww2 = where_is_nodes[node2_num];
			int is_on2 = is_on_nodes[node2_num];
			if( ww2 == ww1 && is_on2 == 1 ) return 1;

			int ic = thecell->surfaces_vec[ww1]->is_contain( &pp );

			if( ic == 0 )
			{
				return 1;
			}
		}
	}
	return 0;
}

//---------------------------------------------------------------------------
//ȷ��Ӧ�õ����ı�
//����ѡ��߽����Ͽ�߽�ıߣ����ѡ��߽����Ͽ�߽�ıߣ������ѡ�����㶼���ƶ��ı�
//���ѡ��һ���ڱ߽����ϣ�һ���ڵ����ڲ��ı�
int Mesher::deter_oper_edge( Tetrahedron *act_tet, int& left_movable, int& right_movable)
{
	bool debuging = false;

	int ww[4], is_on[4];
	ww[0] = where_is_nodes[act_tet->nodesId[0]];
	ww[1] = where_is_nodes[act_tet->nodesId[1]];
	ww[2] = where_is_nodes[act_tet->nodesId[2]];
	ww[3] = where_is_nodes[act_tet->nodesId[3]];
	is_on[0] = is_on_nodes[act_tet->nodesId[0]];
	is_on[1] = is_on_nodes[act_tet->nodesId[1]];
	is_on[2] = is_on_nodes[act_tet->nodesId[2]];
	is_on[3] = is_on_nodes[act_tet->nodesId[3]];

	int left_nn, right_nn;

	//����Ƿ����ڱ߽����Ͽ�߽�ı�
	int hcee = has_cross_edge_edge( act_tet );
	if( hcee > 0 )
	{
		switch( hcee )
		{
		case 1:
			left_nn = 0;
			right_nn = 1;
			break;
		case 2:
			left_nn = 0;
			right_nn = 2;
			break;
		case 3:
			left_nn = 0;
			right_nn = 3;
			break;
		case 4:
			left_nn = 1;
			right_nn = 2;
			break;
		case 5:
			left_nn = 1;
			right_nn = 3;
			break;
		case 6:
			left_nn = 2;
			right_nn = 3;
			break;
		default:
			break;
		}

		int left_nodeN = act_tet->nodesId[left_nn];
		int right_nodeN = act_tet->nodesId[right_nn];

		if( is_corner( left_nodeN ) != 0 ||is_on[left_nn] == 1 )
		{
			left_movable = 0;
		}
		else
		{
			left_movable = 1;
		}
		if( is_corner( right_nodeN ) != 0 || is_on[left_nn] == 1 )
		{
			right_movable = 0;
		}
		else if(ww[left_nn] != -1 && is_on[left_nn] == 0 && ww[right_nn] != -1)
		{
			right_movable = 0;
		}
		else
		{
			right_movable = 1;
		}
		return hcee ;
	}

	//����Ƿ����ڱ߽����Ͽ�߽�ı�
	int hceeof = has_cross_edge_edge_on_face( act_tet );

	Node *n1 = &nodes_vec[act_tet->nodesId[0]];
	Node *n2 = &nodes_vec[act_tet->nodesId[1]];
	Node *n3 = &nodes_vec[act_tet->nodesId[2]];
	Node *n4 = &nodes_vec[act_tet->nodesId[3]];
	Point pp[]={	{n1->x, n1->y,n1->z},
						{n2->x, n2->y,n2->z},
						{n3->x, n3->y,n3->z},
						{n4->x, n4->y,n4->z} };
	//����ѡȡ���˽ڵ㶼���ƶ��Ŀ�߽��
	//�����ҳ����п�߽�ı�
	vector< vector<int> > cross_edges ;
	for( int i=0; i<6; i++ )
	{
		//left_nn right_nn�ֲ����
		//left_ndoeN right_nodeN������
		int left_nn , right_nn ;
		if( i == 0 )			{ left_nn = 0; right_nn = 1; }
		else if( i == 1 )	{ left_nn = 0; right_nn = 2; }
		else if( i == 2 )	{ left_nn = 0; right_nn = 3; }
		else if( i == 3 )	{ left_nn = 1; right_nn = 2; }
		else if( i == 4 )	{ left_nn = 1; right_nn = 3; }
		else if( i == 5 )	{ left_nn = 2; right_nn = 3; }

		int left_nodeN = act_tet->nodesId[left_nn];
		int right_nodeN = act_tet->nodesId[right_nn];

		if( ww[left_nn] != ww[right_nn])
		{
			bool is_same_par1= ww[left_nn]!=-1 && is_on[left_nn]==0 && thecell->surfaces_vec[ww[left_nn]]->is_contain(&pp[right_nn])>0;
			bool is_same_par2= ww[right_nn]!=-1 && is_on[right_nn]==0 && thecell->surfaces_vec[ww[right_nn]]->is_contain(&pp[left_nn])>0;
			if( !is_same_par1 && !is_same_par2 )	continue;

			vector<int> temp_v;
			//ȷ�����˽ڵ��Ƿ����ƶ�
			if( hceeof >= 1 && is_on_same_face( left_nodeN, right_nodeN, 0 ) == 1 )
			{
				temp_v.push_back(i+1);
				temp_v.push_back(left_nn);
				temp_v.push_back(right_nn);

				if( is_on_edge(left_nodeN) != 0 || is_on[left_nn] == 1 )
				{
					temp_v.push_back(0);
				}
				else
				{
					temp_v.push_back(1);
				}

				if( is_on_edge(right_nodeN) != 0 || is_on[right_nn] == 1 )
				{
					temp_v.push_back(0);
				}
				else if(ww[left_nn] != -1 && is_on[left_nn] == 0 && ww[right_nn] != -1)
				{
					temp_v.push_back(0);
				}
				else
				{
					temp_v.push_back(1);
				}
			}
			else if( hceeof == 0 )
			{
				temp_v.push_back(i+1);
				temp_v.push_back(left_nn);
				temp_v.push_back(right_nn);

				if( is_on_face(left_nodeN) != 0 || is_on[left_nn] == 1 )
				{
					temp_v.push_back(0);
				}
				else
				{
					temp_v.push_back(1);
				}
				if( is_on_face(right_nodeN) != 0 || is_on[right_nn] == 1 )
				{
					temp_v.push_back(0);
				}
				else if(ww[left_nn] != -1 && is_on[left_nn] == 0 && ww[right_nn] != -1)
				{
					temp_v.push_back(0);
				}
				else
				{
					temp_v.push_back(1);
				}
			}
			if( temp_v.size() > 0 )
			{
				cross_edges.push_back( temp_v );
			}
		}
	}

	//�����û�п����ص������
	//�����б߿�Խ�������������������ص����֣�
	//���������㷨�����⣬����Ҫ�����������˰������������ҳ����⣩
	//�Ǻ� ���¸�ʦ��������Ȥ
	vector<vector<int> > cross_edges1(cross_edges.begin(),cross_edges.end());
	cross_edges.clear();
	for( int i=0; i<(int)cross_edges1.size(); i++ )
	{
		int ln = cross_edges1[i][1];
		int rn = cross_edges1[i][2];
		if( ww[ln]!=-1 && is_on[ln]==0 && ww[rn]!=-1 && is_on[rn]==0 )
		{
			//������������
			//���������������Ƿ��ص������������㣬���þ����жϣ�
			Point inter_point1, inter_point2;
			int err1 = thecell->surfaces_vec[ww[ln]]->intersect(pp[ln],pp[rn],inter_point1);
			if( err1 == 0 ) continue;
			int err2 = thecell->surfaces_vec[ww[rn]]->intersect(pp[ln],pp[rn],inter_point2);
			if( err2 == 0 ) continue;
			double dis1 = pp[ln].distance_to(inter_point1);
			double dis2 = pp[ln].distance_to(inter_point2);
			if( dis1 > dis2 )
			{
				//�ص���
				//ǿ����Ϊ���ڵ����Ҳ�ڵ����ڿ����ı߽��ϣ��Է�����ȷ��
				//�������취��־�˵�ԪΪ��߽絥Ԫ���ţ��Ժ���˵
				hout	<< "cross two particles: " << act_tet->nodesId[ln] << "("
						<< ww[ln] << ") "
						<< act_tet->nodesId[rn] << "("
						<< ww[rn] << ")" << endl;
				where_is_nodes[act_tet->nodesId[ln]] = ww[rn];
				is_on_nodes[act_tet->nodesId[ln]]    = 1;

				continue;
			}
		}
		cross_edges.push_back(cross_edges1[i]);
	}
	if( cross_edges.size() == 0 ) return 0;

	//������п�߽�ıߣ�����ѡ�����˵㶼���ƶ��ı�
	vector<vector<int> > cross_edges_dmb ;
	for( int i=0; i<(int)cross_edges.size(); i++ )
	{
//		hout << "i=" << i << endl;
		if( cross_edges[i][3] == 1 && cross_edges[i][4] == 1 )
		{
			cross_edges_dmb.push_back( cross_edges[i] );
			left_movable = 1;
			right_movable = 1;
			return cross_edges[i][0];
		}
	}
	//ѡ��ڵ��Žϴ�ı�
	if( cross_edges_dmb.size() > 0 )
	{
		vector<int> the_edge(cross_edges_dmb[0].begin(),cross_edges_dmb[0].end());
		for( int i=1; i<(int)cross_edges_dmb.size(); i++ )
		{
			if( max(	act_tet->nodesId[cross_edges_dmb[i][1]],
						act_tet->nodesId[cross_edges_dmb[i][2]]) >
				max(	act_tet->nodesId[the_edge[1]],
						act_tet->nodesId[the_edge[2]]))
			{
					the_edge.assign( cross_edges_dmb[i].begin(),
											  cross_edges_dmb[i].end() );
			}
		}
		left_movable = 1;
		right_movable = 1;
		return the_edge[0];
	}
	//������п�߽�ıߣ�����ѡ�����˵㶼������ǿ���ڵıߣ���ֹ��ѭ����
	for( int i=0; i<(int)cross_edges.size(); i++ )
	{
		if( is_on[cross_edges[i][1]] == 0 && is_on[cross_edges[i][2]] == 0 )
		{
			left_movable = cross_edges[i][3];
			right_movable = cross_edges[i][4];
			return cross_edges[i][0];
		}
	}
	left_movable = cross_edges[0][3];
	right_movable = cross_edges[0][4];
	return cross_edges[0][0];
}

//---------------------------------------------------------------------------
//�ڵ��Ƿ��ڵ����Ľ�����
int Mesher::is_on_face(Node* node, int mod)
{
	if( fabs(node->x - x_min) <= ZERO ) return -1 ;
	else if( fabs(node->x - x_max) <= ZERO ) return -2 ;
	if( fabs(node->y - y_min) <= ZERO ) return -3 ;
	else if( fabs(node->y - y_max) <= ZERO ) return -4 ;
	if( fabs(node->z - z_min) <= ZERO ) return -5 ;
	else if( fabs(node->z - z_max) <= ZERO ) return -6 ;

	if( mod == 1 )
	{
		int is_on = 0;
	}
	return 0;
}
//---------------------------------------------------------------------------
int Mesher::is_on_face(int node_num, int mod)
{
	Node *node = &nodes_vec[node_num];  
	if( fabs(node->x - x_min) <= ZERO ) return -1 ;
	else if( fabs(node->x - x_max) <= ZERO ) return -2 ;
	if( fabs(node->y - y_min) <= ZERO ) return -3 ;
	else if( fabs(node->y - y_max) <= ZERO ) return -4 ;
	if( fabs(node->z - z_min) <= ZERO ) return -5 ;
	else if( fabs(node->z - z_max) <= ZERO ) return -6 ;
                                
	if( mod == 1 )
	{
		int ww = where_is_nodes[node_num] ;
		if( is_on_nodes[node_num] == 1 ) return ww + 1;
	}

	return 0;
}
//---------------------------------------------------------------------------
//����������Ƿ������ڱ߽�����
//����ֵΪ�ڱ߽����ϵ����
//      0: û�����ڱ߽�������
//      1: 0-1-2 ��
//      2: 0-2-3 ��
//      3: 0-3-1 ��
//      4: 1-2-3 ��
int Mesher::is_on_face(Tetrahedron* act_tet, int mod)
{
        int iof1 = is_on_same_face( act_tet->nodesId[0], act_tet->nodesId[1], 0 ) ;
        int iof2 = is_on_same_face( act_tet->nodesId[0], act_tet->nodesId[2], 0 ) ;
        int iof3 = is_on_same_face( act_tet->nodesId[0], act_tet->nodesId[3], 0 ) ;
        int iof4 = is_on_same_face( act_tet->nodesId[1], act_tet->nodesId[2], 0 ) ;
        int iof5 = is_on_same_face( act_tet->nodesId[1], act_tet->nodesId[3], 0 ) ;
        int iof6 = is_on_same_face( act_tet->nodesId[2], act_tet->nodesId[3], 0 ) ;
        if( iof1 == 1 && iof2 == 1 && iof4 == 1 ) return 1 ;
        if( iof2 == 1 && iof3 == 1 && iof6 == 1 ) return 2 ;
        if( iof1 == 1 && iof3 == 1 && iof5 == 1 ) return 3 ;
        if( iof4 == 1 && iof5 == 1 && iof6 == 1 ) return 4 ;

        return 0;
}

//---------------------------------------------------------------------------
//�ж�һ���ڵ��Ƿ�Ϊ�ǵ㣨������Ľ��㣩
int Mesher::is_corner( int node_num )
{
	Node *node = &nodes_vec[node_num];
	//���ڵ��Ƿ�ͬʱ���������ϣ�����ǣ���˵���˵������潻����
	int ofs = 0 ;

	if( fabs(node->x - x_min) <= ZERO ) ofs ++ ;
	else if( fabs(node->x - x_max) <= ZERO ) ofs ++ ;
	if( fabs(node->y - y_min) <= ZERO ) ofs ++ ;
	else if( fabs(node->y - y_max) <= ZERO ) ofs ++ ;
	if( fabs(node->z - z_min) <= ZERO ) ofs ++ ;
	else if( fabs(node->z - z_max) <= ZERO ) ofs ++ ;
	if( is_on_nodes[node_num] == 1 ) ofs ++ ;

	if( ofs > 2 ) return 1;
	return 0;
}

//---------------------------------------------------------------------------
//���߽����ϵĵ�Ԫ�ڱ߽����ϵĽڵ��Ƿ���߽�ܽ�
int Mesher::check_bft_node( Tetrahedron *act_tet, double alength )
{
	//������ϵĵ��Ƿ���߽�ܽ�
	int node1_num=-1, node2_num=-1, node3_num=-1 ;
	int act_iof = is_on_face( act_tet , 0 );
	int w1 = where_is_nodes[act_tet->nodesId[0]];
	int w2 = where_is_nodes[act_tet->nodesId[1]];
	int w3 = where_is_nodes[act_tet->nodesId[2]];
	int w4 = where_is_nodes[act_tet->nodesId[3]];
	int ww = w1 ;
	if( ww == -1 )
	{
		if( w2 != -1 ) ww = w2;
		else if( w3 != -1 ) ww = w3;
		else if( w4 != -1 ) ww = w4;  
	}
	int is_on1 = is_on_nodes[act_tet->nodesId[0]];
	int is_on2 = is_on_nodes[act_tet->nodesId[1]];
	int is_on3 = is_on_nodes[act_tet->nodesId[2]];
	int is_on4 = is_on_nodes[act_tet->nodesId[3]];

	switch( act_iof )
	{
	case 1:
		if( is_on1 == 1 && is_on2 == 1 && w3 == -1 )
		{
			node1_num = act_tet->nodesId[2] ;
			node2_num = act_tet->nodesId[0] ;
			node3_num = act_tet->nodesId[1] ;
		}
		else if( is_on2 == 1 && is_on3 == 1 && w1 == -1 )
		{
			node1_num = act_tet->nodesId[0] ;
			node2_num = act_tet->nodesId[1] ;
			node3_num = act_tet->nodesId[2] ;
		}
		else if( is_on3 == 1 && is_on1 == 1 && w2 == -1 )
		{
			node1_num = act_tet->nodesId[1] ;
			node2_num = act_tet->nodesId[2] ;
			node3_num = act_tet->nodesId[0] ;
		}
		break;
	case 2:            
		if( is_on1 == 1 && is_on3 == 1 && w4 == -1 )
		{
			node1_num = act_tet->nodesId[3] ;
			node2_num = act_tet->nodesId[0] ;
			node3_num = act_tet->nodesId[2] ;
		}
		else if( is_on3 == 1 && is_on4 == 1 && w1 == -1 )
		{
			node1_num = act_tet->nodesId[0] ;
			node2_num = act_tet->nodesId[2] ;
			node3_num = act_tet->nodesId[3] ;
		}
		else if( is_on4 == 1 && is_on1 == 1 && w3 == -1 )
		{
			node1_num = act_tet->nodesId[2] ;
			node2_num = act_tet->nodesId[3] ;
			node3_num = act_tet->nodesId[0] ;
		}
		break;
	case 3:
		if( is_on1 == 1 && is_on2 == 1 && w4 == -1 )
		{
			node1_num = act_tet->nodesId[3] ;
			node2_num = act_tet->nodesId[0] ;
			node3_num = act_tet->nodesId[1] ;
		}
		else if( is_on2 == 1 && is_on4 == 1 && w1 == -1 )
		{
			node1_num = act_tet->nodesId[0] ;
			node2_num = act_tet->nodesId[1] ;
			node3_num = act_tet->nodesId[3] ;
		}
		else if( is_on4 == 1 && is_on1 == 1 && w2 == -1 )
		{
			node1_num = act_tet->nodesId[1] ;
			node2_num = act_tet->nodesId[3] ;
			node3_num = act_tet->nodesId[0] ;
		}
		break;
	case 4:               
		if( is_on2 == 1 && is_on3 == 1 && w4 == -1 )
		{
			node1_num = act_tet->nodesId[3] ;
			node2_num = act_tet->nodesId[1] ;
			node3_num = act_tet->nodesId[2] ;
		}
		else if( is_on3 == 1 && is_on4 == 1 && w2 == -1 )
		{
			node1_num = act_tet->nodesId[1] ;
			node2_num = act_tet->nodesId[2] ;
			node3_num = act_tet->nodesId[3] ;
		}
		else if( is_on4 == 1 && is_on2 == 1 && w3 == -1 )
		{
			node1_num = act_tet->nodesId[2] ;
			node2_num = act_tet->nodesId[3] ;
			node3_num = act_tet->nodesId[1] ;
		}
		break;
	default:
		break;
	}

	if( node1_num != -1 && node2_num != -1 && node3_num != -1 )
	{
		Point p1 = { nodes_vec[node1_num].x, nodes_vec[node1_num].y, nodes_vec[node1_num].z };
		Node mid_n = (nodes_vec[node2_num] + nodes_vec[node3_num])/2.0 ;
		alength = nodes_vec[node1_num].distance_to(&mid_n)/2.0;
		Point mid_p = { mid_n.x, mid_n.y, mid_n.z };
		Point rp1;
		if( thecell->surfaces_vec[ww]->intersect( p1, mid_p, rp1 ) )
		{
			Node rn1 = Node( rp1 );
			double dis = rn1.distance_to( &nodes_vec[node1_num] ) ;

			if( dis < alength )
			{
				Node temp_n = nodes_vec[node1_num] ;//�����������Ա���ִ���ʱ�ָ�
				nodes_vec[node1_num].move_to( &rn1 );
				double relative_min_volume = cal_relative_min_volume( node1_num );

				if( relative_min_volume < min_volume )
				{
					nodes_vec[node1_num].move_to( &temp_n );
				}
				else
				{
					where_is_nodes[node1_num] = ww;
					is_on_nodes[node1_num] = 1;
					return 1 ;
				}
			}
		}
	}
	return 0 ;
}

//---------------------------------------------------------------------------
//����ͽڵ���صĵ�Ԫ�������Сֵ
double Mesher::cal_relative_min_volume( int node_num )
{
	bool debuging = false;

	double min_re_volume = 1e300 ;

	for( int i=0; i<(int)nodes_vec[node_num].relative_eles_vec.size(); i++ )
	{
		Tetrahedron *act_tet = &eles_vec[nodes_vec[node_num].relative_eles_vec[i]] ;
		if( act_tet->is_contain( node_num ) == -1 ) continue ;
		double volume = cal_tet_volume( act_tet );

		if( debuging )
		{
			hout	<< "       element " << nodes_vec[node_num].relative_eles_vec[i]
					<< " : "
					<< act_tet->nodesId[0] << " "
					<< act_tet->nodesId[1] << " "
					<< act_tet->nodesId[2] << " "
					<< act_tet->nodesId[3] << " "
					<< act_tet->materialId << " : "   
					<< volume << endl;
		}

		if( volume < min_re_volume )
		{
			min_re_volume = volume ;
		}
	}

	return min_re_volume ;
}

//---------------------------------------------------------------------------
//���������嵥Ԫ�����(���������
double Mesher::cal_tet_volume( Tetrahedron* tet )
{
	return cal_volume(	 &nodes_vec[tet->nodesId[0]],
								 &nodes_vec[tet->nodesId[1]],
								 &nodes_vec[tet->nodesId[2]],
								 &nodes_vec[tet->nodesId[3]]	);
}

//---------------------------------------------------------------------------
//��������ĸ��ڵ㶼�ڽ����ϵĵ�Ԫ�������֮
int Mesher::split_grovelling_eles()
{
	bool debuing=false;
	double alength = global_length*0.8;

	vector<int> new_eles_num;
	deal_grov_eles(thin_tets,new_eles_num);
	for( int j=0; j<(int)new_eles_num.size(); j++ )
	{
		put_into_cbev( new_eles_num[j], alength, 0 );
	}
	thin_tets.clear();

	return 1;
}

//---------------------------------------------------------------------------
//�����ĸ��ڵ㶼�ڽ����ϵĵ�Ԫ���ڳ����м����ڵ㣩
int Mesher::deal_grov_eles(vector<int> &grov_eles, vector<int> &new_eles_num)
{
	for( int i=0; i<(int)grov_eles.size(); i++ )
	{
		rectify_tet(grov_eles[i],new_eles_num);
	}
	return 1;      
}

//---------------------------------------------------------------------------
//������еĵ�Ԫ�������������ļнǷǳ��󣨽ӽ�180�������ڶԱ߲���һ���ڵ�
int Mesher::rectify_tet(int ele_num, vector<int> &new_ele_num)
{
	bool debuging = false;

	double max_angle = 160.0*PI/180.0;
	double lmin_dis = global_length/2.0;

	if( debuging )
	{
		hout << "-----------------------rectify_tet()------------------------" << endl;
		hout << "ele_num: " << ele_num << " " << " size: " << (int)eles_vec.size() << endl ;
		hout << "  nodes: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "      " << eles_vec[ele_num].nodesId[k]  ;
		}
		hout << " mat: " << eles_vec[ele_num].materialId ;
		hout <<endl;
	}


	int is_deleted = 0;
	int n0 = eles_vec[ele_num].nodesId[0];
	for( int i=1; i<4; i++ )
	{
		int nn = eles_vec[ele_num].nodesId[i];
		if( nn == n0 ) is_deleted = 1;
	}
	if( debuging )
	{
		hout << "  is_deleted: " << is_deleted << endl;
	}
	if( is_deleted == 1 ) return 1;

	//�ж��Ƿ��ĸ��ڵ�ȫ���ڽ�����
	//��˳�����Ƿ��ĸ��ڵ��У�����ص�ԪΪ��Ԫ�ĸ������ܳ���2���������ֲ���������ͣ�Ĳ𣩣�
	int is_all_on_face = 1;
	int re_all_on_face = 0;
	for( int i=0; i<4; i++ )
	{
		int nn = eles_vec[ele_num].nodesId[i];
		if( debuging ) hout << "nn : " << nn << endl;

		int nn_re = num_of_thin_tet_re(nn);
		if( nn_re > re_all_on_face ) re_all_on_face = nn_re;
		if( nn_re > 1 )
		{
//			hout << "  Node " << nn << " has " << nn_re << " thin relative elements." << endl;
		}
		if( is_on_nodes[nn] == 0  )
		{
			is_all_on_face = 0;
		}
	}
	if( eles_vec[ele_num].materialId == 0 ) is_all_on_face = 0;

	if( debuging ) hout << "  is_all_on_face: " << is_all_on_face <<endl;

	if( is_all_on_face == 1 )
	{
		for( int i=0; i<4; i++ )
		{
			int nn = eles_vec[ele_num].nodesId[i];
			if( where_is_nodes[nn] == -1 ) continue;
			int cbe = can_be_extracted(nn);
			if( cbe == 1 )
			{
				//��˽ڵ����������е�Ԫ�����˵�ǰ����Ԫ�����ǻ��嵥Ԫ��ֱ�����ȥ
				where_is_nodes[nn] = -1;
				is_on_nodes[nn] = 0;
				//eles_vec[ele_num].materialId = 0;
				//�����д˽ڵ����ص�Ԫ�ĳɻ�����ϣ�˳���˳һ���£�hoho
				Node smoothed_node(0,0,0);
				int sm_num =0;
				for( int j=0; j<(int)nodes_vec[nn].relative_eles_vec.size(); j++ )
				{
					int en = nodes_vec[nn].relative_eles_vec[j];
					int ic = eles_vec[en].is_contain(nn);
					if( ic == -1 ) continue;
					eles_vec[en].materialId = 0;
					for( int k=0; k<4; k++ )
					{
						int rnn = eles_vec[en].nodesId[k];
						smoothed_node = smoothed_node + nodes_vec[rnn];
						sm_num ++;
					}
				}
				smoothed_node = smoothed_node/(double)sm_num;
				Node old_node = nodes_vec[nn];
				nodes_vec[nn].move_to(&smoothed_node);
				double rvolume = cal_relative_min_volume(nn);
				if( rvolume < 0 )
				{
					nodes_vec[nn].move_to(&old_node);
//					hout << "Smooth failed at node " << nn << "." << endl;
				}

//				hout << "  Moved node: " << nn << " into matrix." << endl;

				return 1;
			}
		}
	}
	else
	{
		//�����ǿ����ڲ��ı�Ԫ����������ӵ���������Ԫʱ���������ˣ����У������Σ�
		if( re_all_on_face > 1 ) return 0;
	}
	//��������ӵĻ��������ܶౡԪʱ���������ˣ������п��ܳ�����ѭ������
	if( re_all_on_face > 2 ) return 0;

	int itt = is_thin_tet(ele_num);
	if( itt == 0 ) return 1;
	//�ѵ�Ԫ�ֽ��4�������Σ�ɾ��
	//�����ҳ��Ƕ�����������
	int nsn[4];        //�ĸ��ڵ���
	double max_ang = max_angle_of_tet(ele_num,nsn);

	if( debuging )
	{
		hout << "  nodes: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "      " << nsn[k]  ;
		}
		hout << endl;
		hout << "  where: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "      " << where_is_nodes[nsn[k]]  ;
		}
		hout << endl;
		hout << "  is on: " ;
		for( int k=0; k<4; k++ )
		{
			hout << "      " << is_on_nodes[nsn[k]]  ;
		}
		hout << endl;
		hout << "  max_ang: " << max_ang*180.0/PI << endl;
	}

	if( max_ang < max_angle ) return 0;

	if( re_all_on_face < 2 )
	{
		//�����½ڵ�λ�ã��ĸ��ڵ��ƽ����
		Node cen_node(0,0,0);
		for( int i=0; i<4; i++ )
		{
			cen_node = cen_node + nodes_vec[nsn[i]];
		}
		cen_node = cen_node /4.0;
		int new_node_num = (int)nodes_vec.size();
		nodes_vec.push_back(cen_node);
		int is_on;
		int whnn = where_is(&cen_node,is_on);
		if( is_all_on_face == 1 )
		{
			//������ĸ��ڵ㶼��������ĵ�Ԫ�������Ƿ���Ҫ���½ڵ�ͶӰ������
			if( whnn != -1 )
			{
				int perr=-1;
				Node pnode = project2_elli_surf(cen_node,thecell->surfaces_vec[whnn],&perr);
				if( perr == 1 )
				{
					nodes_vec[new_node_num].move_to(&pnode);
					is_on = 1;
				}
			}
		}        
		where_is_nodes.push_back(whnn);
		is_on_nodes.push_back(is_on);

		if( debuging )
		{
			hout << " angle: " << max_ang*180/PI << endl;

			hout	<< "dis: " << nodes_vec[nsn[0]].distance_to(&nodes_vec[nsn[3]]) << " min: " << lmin_dis <<endl;
			hout	<< "Generated a new node: " << new_node_num
					<< " where: " << where_is_nodes[new_node_num]
					<< " is on: " << is_on_nodes[new_node_num] <<endl << endl;
		}

		//ɾ���ɵ�Ԫ
		//�ȱ�ǣ���ͷ��ɾ
		for( int i=1; i<4; i++ )
		{
			eles_vec[ele_num].nodesId[i] = eles_vec[ele_num].nodesId[0];
		}

		deleting_ele.push_back(ele_num);

		insert_node(nsn[0],nsn[3],new_node_num,new_ele_num);
		insert_node(nsn[1],nsn[2],new_node_num,new_ele_num);
	}
	else
	{
		//�ڽǶ����ĶԱ߲���ڵ�
		Node new_node = nodes_vec[nsn[0]]+nodes_vec[nsn[3]];
		new_node = new_node/2.0;

		int new_node_num = (int)nodes_vec.size();
		nodes_vec.push_back(new_node);
		int is_on;
		int wh = where_is(&new_node,is_on);

		if( eles_vec[ele_num].materialId == 0 ) wh = -1;
		if( wh != -1 )
		{
			int ww = where_is_nodes[eles_vec[ele_num].nodesId[0]];
			if( is_all_on_face == 1 )
			{
				int mn = face_particle(nsn[0],nsn[3]);
				if(debuging) hout << "mn: " << mn << endl;
				if( mn == 0 )
				{
					wh = -1;
				}
				else
				{
					wh = ww;
				}
			}
			else    wh = ww;
		}
		is_on=0;
		where_is_nodes.push_back(wh);
		is_on_nodes.push_back(is_on);

		insert_node(nsn[0],nsn[3],new_node_num,new_ele_num);
	}
	return 1;
}

//---------------------------------------------------------------------------
//���һ���ڵ����ص�Ԫ���м�����Ԫ
int Mesher::num_of_thin_tet_re(int nn)
{
	bool debuging = false;

	int te_num=0; //��ص�Ԫ�б�Ԫ��Ŀ
	for( int i=0; i<(int)nodes_vec[nn].relative_eles_vec.size(); i++ )
	{
		int ele_num = nodes_vec[nn].relative_eles_vec[i];
		int ist = is_thin_tet(ele_num);
		if( ist == 1 ) te_num ++;
	}
	return te_num;
}

//---------------------------------------------------------------------------
//�ж�һ���ڵ�����������ڵ��Ƿ����ڿ����ڣ����ڵ����е���ص�Ԫ�Ľڵ㣬�нڵ��ڿ����ڣ�����0������ȫ���ڻ����ڻ��ڿ������棩����1��
int Mesher::can_be_extracted(int n1)
{
	if( nodes_vec[n1].flag >= 0 ) return 0;  //���������ϵĵ�
	for( int i=0; i<(int)nodes_vec[n1].relative_eles_vec.size(); i++ )
	{
		int en = nodes_vec[n1].relative_eles_vec[i];
		int ic = eles_vec[en].is_contain(n1);
		if( ic == -1 ) continue;   //�Ѿ����ٰ����ڵ�n1�ˣ�û�и��µ�Ե�ʣ�
		for( int j=0; j<4; j++ )
		{
			if( j == ic ) continue;
			int nn = eles_vec[en].nodesId[j];
			if( where_is_nodes[nn] != -1 && is_on_nodes[nn] == 0 )	return 0;
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//������������ڵ���ص�Ԫ�Ĳ��Ϻţ�ͬʱ���������ڵ�����е�Ԫ�������,-1��ʾ���壩
int Mesher::face_particle(int n1,int n2)
{
	bool debuging = false;

	if( debuging ) hout << "n1: " << n1 << " n2: " << n2 << endl;
	int wn2=1;
	if( debuging ) hout << "    elements: " ;
	for( int i=0; i<(int)nodes_vec[n1].relative_eles_vec.size(); i++ )
	{
		int ele_num = nodes_vec[n1].relative_eles_vec[i];
		int isc1 = eles_vec[ele_num].is_contain(n1);
		if( isc1 == -1 ) continue;
		int isc2 = eles_vec[ele_num].is_contain(n2);
		if( isc2 == -1 ) continue;
		if( debuging ) hout << ele_num << " mat: " <<  eles_vec[ele_num].materialId << endl;
		if( eles_vec[ele_num].materialId == 0 )
		{
			wn2 = 0;
			break;
		}
	}
	return wn2;
}

//---------------------------------------------------------------------------
//�������ڵ��м����ڵ�
int Mesher::insert_node(int node1_num, int node2_num, int node3_num,
									vector<int> &new_ele_num)
{                
	bool debuging = false;

	//���ڵ�Node1_num��������ص�Ԫ
	for( int i=0; i<(int)nodes_vec[node1_num].relative_eles_vec.size(); i++ )
	{
		//���õ�Ԫ�Ƿ�����������ڵ�
		int ele_num = nodes_vec[node1_num].relative_eles_vec[i];
		int node1_i = eles_vec[ele_num].is_contain( node1_num );
		if( node1_i == -1 )		continue;
		int node2_i = eles_vec[ele_num].is_contain( node2_num );
		if( node2_i == -1 )		continue;
		if( debuging )
		{
			hout << "insert node into tet: " << ele_num ;
			hout << " nodes: " << eles_vec[ele_num].nodesId[0] << " " ;
			hout << eles_vec[ele_num].nodesId[1] << " " ;
			hout << eles_vec[ele_num].nodesId[2] << " " ;
			hout << eles_vec[ele_num].nodesId[3] << endl ;
		}

		//ȷ����������
		if( node1_i > node2_i )
		{
			int temp_i = node1_i ;
			node1_i = node2_i ;
			node2_i = temp_i ;
		}

		Tetrahedron *new_tet;
		if( node1_i == 0 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[1],
				eles_vec[ele_num].nodesId[2], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}
		else if( node1_i == 1 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[2],
				eles_vec[ele_num].nodesId[0], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}
		else if( node1_i == 2 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[0],
				eles_vec[ele_num].nodesId[1], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}

		int new_tet_num = (int)eles_vec.size();
		eles_vec.push_back( *new_tet );
		nodes_vec[node3_num].relative_eles_vec.push_back( ele_num );
		for( int k=0; k<4; k++ )
		{
			nodes_vec[eles_vec.back().nodesId[k]].relative_eles_vec.push_back( new_tet_num );
		}

		new_ele_num.push_back(ele_num);
		new_ele_num.push_back(new_tet_num);

		//ȷ����������Ĳ�������
		deter_mat_ele(ele_num);
		deter_mat_ele(new_tet_num); 
	}
	return 1;
}

//---------------------------------------------------------------------------
//����߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
int Mesher::deal_with_eles_across_boundary_be(vector<Bar> &new_edges,int debuging)
{
	bool debuging_map_mesh = false;
	if(debuging == 1)
	{
		debuging_map_mesh = true;
	}
	for( int i=0; i<(int)eles_across_boundary_be.size(); i++ )
	{
		if( debuging_map_mesh) hout	<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_map_mesh) hout	<< "deal_with_eles_across_boundary_be" << endl;
		int ele_num = eles_across_boundary_be[i] ;
		if( debuging_map_mesh) hout	<< " i: " << i << " max: " << (int)eles_across_boundary_be.size() << " ele_num: " << ele_num << endl;
		if( debuging_map_mesh) hout	<< " act_tet: " << eles_vec[ele_num].nodesId[0] << " "
														<< eles_vec[ele_num].nodesId[1] << " "
														<< eles_vec[ele_num].nodesId[2] << " "
														<< eles_vec[ele_num].nodesId[3] << endl;

		int circle_num = 0;
		while(true && circle_num < 20)
		{
			circle_num ++ ;
			if( has_cross_edge_edge( &eles_vec[ele_num] ) == 0 )
			{
				eles_across_boundary_bf.push_back( ele_num );
				if( debuging_map_mesh )
				{
					hout << " put into eles_across_boundary_bf!" << endl;
				}
				break ;
			}
			int deal_result = deal_with_ele_across_boundary(ele_num,0,new_edges,debuging);
			if( deal_result == 0 ) return 0;
			if( deter_mat_ele(ele_num) != -1 ) break;
		}

		if( thin_tets.size() > 0 )
		{
			//�˵�ԪΪ��Ԫ
			for( int j=0; j<(int)thin_tets.size(); j++ )
			{
				vector<int> new_tets;
				rectify_tet(thin_tets[j],new_tets);
				for( int k=0; k<(int)new_tets.size(); k++ )
				{
					put_into_cbev( new_tets[k], global_length*0.8, 0 );
				}
			}
			thin_tets.clear();
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//����߽����ϵĵ�Ԫ(�ڱ߽����Ͽ�߽�ĵ�Ԫ��
int Mesher::deal_with_eles_across_boundary_bf(vector<Bar> &new_edges,int debuging)
{
	bool debuging_map_mesh = false;       
	if(debuging == 1)
	{
		debuging_map_mesh = true;
	}
	for( int i=0; i<(int)eles_across_boundary_bf.size(); i++ )
	{
		if( debuging_map_mesh) hout	<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_map_mesh) hout	<< "deal_with_eles_across_boundary_bf" << endl;
		int ele_num = eles_across_boundary_bf[i] ;
		if( debuging_map_mesh) hout	<< " i: " << i << " ele_num: " << ele_num << endl;
		if( debuging_map_mesh) hout	<< " act_tet: " << eles_vec[ele_num].nodesId[0] << " "
														<< eles_vec[ele_num].nodesId[1] << " "
														<< eles_vec[ele_num].nodesId[2] << " "
														<< eles_vec[ele_num].nodesId[3] << endl;

		int circle_num = 0;
		while(true && circle_num < 20)
		{
			circle_num ++ ;
			//��鵱ǰ�������Ƿ��������ڱ߽����ϣ��ұ߽����ϵ����Խ�������ά�߽�
			int hceeof = has_cross_edge_edge_on_face( &eles_vec[ele_num] ) ;
			if( debuging_map_mesh)
			{
				hout << " hceeof: " << hceeof << endl;
			}
			if( hceeof <= 0 )
			{
				eles_across_boundary_dmb.push_back( ele_num );
				if( debuging_map_mesh)
				{
					hout << " put into eles_across_boundary_dmb!" << endl;
				}
				break;
			}
			int deal_result = deal_with_ele_across_boundary(ele_num,1,new_edges,debuging);
			if( deal_result == 0 ) return 0;
			if( deter_mat_ele(ele_num) != -1 ) break;
		}

		if( thin_tets.size() > 0 )
		{
			//�˵�ԪΪ��Ԫ
			for( int j=0; j<(int)thin_tets.size(); j++ )
			{
				vector<int> new_tets;
				rectify_tet(thin_tets[j],new_tets);
				for( int k=0; k<(int)new_tets.size(); k++ )
				{
					put_into_cbev( new_tets[k], global_length*0.8, 1 );
				}
			}
			thin_tets.clear();
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//���������˵㶼�����ƶ��Ŀ�߽絥Ԫ
int Mesher::deal_with_eles_across_boundary_dmb(vector<Bar> &new_edges,int debuging)
{
	bool debuging_map_mesh = false ;  
	if( debuging == 1 )
	{
		debuging_map_mesh = true;
	}
	for( int i=0; i<(int)eles_across_boundary_dmb.size(); i++ )
	{
		if( debuging_map_mesh) hout	<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_map_mesh) hout	<< "deal_with_eles_across_boundary_dmb" << endl;

		int ele_num = eles_across_boundary_dmb[i] ;

		if( debuging_map_mesh) hout	<< " i: " << i << " ele_num: " << ele_num << endl;
		if( debuging_map_mesh) hout	<< " act_tet: " << eles_vec[ele_num].nodesId[0] << " "
														<< eles_vec[ele_num].nodesId[1] << " "
														<< eles_vec[ele_num].nodesId[2] << " "
														<< eles_vec[ele_num].nodesId[3] << endl;

		int circle_num = 0;
		while(true && circle_num < 20)
		{
			circle_num ++ ;
			int left_movable, right_movable;
			int edge_num = deter_oper_edge( &eles_vec[ele_num], left_movable, right_movable );
			if( edge_num != 0 )
			{
				if( left_movable != 1 || right_movable != 1)
				{
					eles_across_boundary_ndmb.push_back(ele_num);
					if( debuging_map_mesh)
					{
						hout << " put into eles_across_boundary_ndmb!" << endl;
					}
					break ; 
				}
				int deal_result = deal_with_ele_across_boundary(ele_num,2,new_edges,debuging);
				if( deal_result == 0 ) return 0;           
			}
			if( deter_mat_ele(ele_num) != -1 ) break;
		}

		if( thin_tets.size() > 0 )
		{
			//�˵�ԪΪ��Ԫ
			for( int j=0; j<(int)thin_tets.size(); j++ )
			{
				vector<int> new_tets;
				rectify_tet(thin_tets[j],new_tets);
				for( int k=0; k<(int)new_tets.size(); k++ )
				{
					put_into_cbev( new_tets[k], global_length*0.8, 2 );
				}
			}
			thin_tets.clear();
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//���������˵�����һ���������ƶ��Ŀ�߽絥Ԫ
int Mesher::deal_with_eles_across_boundary_ndmb(vector<Bar> &new_edges,int debuging)
{
	bool debuging_map_mesh = false; 
	if( debuging == 1 )
	{
		debuging_map_mesh = true;
	}
	//����ƶ��ڵ�ʧ�ܣ����ָ������Ԫ������ʱ����ڴ������У�����ٴ���
	vector<int> undealed_ele;

	for( int i=0; i<(int)eles_across_boundary_ndmb.size(); i++ )
	{
		int ele_num = eles_across_boundary_ndmb[i] ;
		if( debuging_map_mesh) hout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_map_mesh) hout << "deal_with_eles_across_boundary_ndmb" << endl;
		if( debuging_map_mesh) hout << " i: " << i << " size: " << (int)eles_across_boundary_ndmb.size() << " ele_num: " << ele_num << " size: " << (int)eles_vec.size() << endl;

		int circle_num = 0;
		while(true && circle_num < 5)
		{
			circle_num ++ ;
			if( debuging_map_mesh) hout	<< "~~~~~~~~~~~~~~~~" << endl;
			if( debuging_map_mesh) hout	<< " act_tet: " << eles_vec[ele_num].nodesId[0] << " "
															<< eles_vec[ele_num].nodesId[1] << " "
															<< eles_vec[ele_num].nodesId[2] << " "
															<< eles_vec[ele_num].nodesId[3] << endl;

			int deal_result = deal_with_ele_across_boundary(ele_num,3,new_edges,debuging);
			if( deal_result == 0 ) return 0;
			if( deal_result == 2 )
			{
				undealed_ele.push_back(ele_num);
				break;
			}
			int dme = deter_mat_ele(ele_num);
			if( debuging_map_mesh) hout << "ddme: " << dme << endl;
			if( dme != -1 ) break;
		}

		if( thin_tets.size() > 0 )
		{
			//�˵�ԪΪ��Ԫ
			for( int j=0; j<(int)thin_tets.size(); j++ )
			{
				vector<int> new_tets;
				rectify_tet(thin_tets[j],new_tets);
				for( int k=0; k<(int)new_tets.size(); k++ )
				{
					put_into_cbev( new_tets[k], global_length*0.8, 3 );
				}
			}
			thin_tets.clear();
		}
	}

	debuging_map_mesh = false;
	eles_across_boundary_ndmb.assign( undealed_ele.begin(),undealed_ele.end() );

	for( int i=0; i<(int)eles_across_boundary_ndmb.size(); i++ )
	{
		int ele_num = eles_across_boundary_ndmb[i] ;
		if( debuging_map_mesh) hout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_map_mesh) hout << "deal_with_eles_across_boundary_ndmb_undealed_ele" << endl;
		if( debuging_map_mesh) hout << " i: " << i << " size: " << (int)eles_across_boundary_ndmb.size() << " ele_num: " << ele_num << " size: " << (int)eles_vec.size() << endl;

		int circle_num = 0;
		while(true && circle_num < 25)
		{
			circle_num ++ ;

			if( debuging_map_mesh) hout	<< "~~~~~~~~~~~~~~~~" << endl;
			if( debuging_map_mesh) hout	<< " act_tet: " << eles_vec[ele_num].nodesId[0] << " "
															<< eles_vec[ele_num].nodesId[1] << " "
															<< eles_vec[ele_num].nodesId[2] << " "
															<< eles_vec[ele_num].nodesId[3] << endl;

			deal_with_ele_across_boundary(ele_num,3,new_edges,debuging);
			if( debuging_map_mesh) hout << "dme: " << deter_mat_ele(ele_num) << endl;
			if( deter_mat_ele(ele_num) != -1 ) break;
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//�����߽絥Ԫ
//ele_kind:��Ԫ��������(�ڱ߽����Ͽ����ĵ�Ԫ���ڱ߽����Ͽ�߽�ĵ�Ԫ......)
//mod���������ƶ��ڵ����ָ������Ԫʱ�Ĵ�����  ��������û�ã�
//      0�������������򷵻�2��
//      1������ڵ�
int Mesher::deal_with_ele_across_boundary( int ele_num, int ele_kind, vector<Bar> &new_edges,  int mod=1 )
{
	bool debuging_map_mesh = false ;
	if( mod ==1 ) debuging_map_mesh = true;

	if( debuging_map_mesh ) hout << "++++++++++++++++++++++++++++" << endl;

	if( debuging_map_mesh) hout << "ele_num: " << ele_num << " ele type: " << ele_kind << endl;

	double alength = global_length*0.8;
	if( ele_kind == 1 ) alength = global_length*0.5;   //�����Ͽ�߽絥Ԫ
	double min_dis = global_length*3.0/8.0;
	int wh[4],is_on[4];
	for( int j=0; j<4; j++ )
	{
		wh[j] = where_is_nodes[eles_vec[ele_num].nodesId[j]];
		if( wh[j] == -3 )			//���ֽڵ�ͬʱ��������ǿ����
		{ 
			return 0;
		}
		is_on[j] = is_on_nodes[eles_vec[ele_num].nodesId[j]];
	}
	//ȷ�������߿�߽�
	int left_movable, right_movable;
	int edge_num = deter_oper_edge( &eles_vec[ele_num], left_movable, right_movable );

	if( edge_num == 0 ) return -1;     //û�б߿�߽�

	int nodes_id[4];
	deter_nodes_id(edge_num,nodes_id);

	int *pww;
	if( wh[nodes_id[0]] != -1 && is_on[nodes_id[0]] == 0 )
	{
		pww=&wh[nodes_id[0]];
	}
	else
	{
		pww=&wh[nodes_id[1]];
	}

	//������ǿ���Ͻ���Ľ���       
	if( debuging_map_mesh) hout	<< "edge_num: " << edge_num
													<< " nodes: " << eles_vec[ele_num].nodesId[nodes_id[0]] << " "
													<< eles_vec[ele_num].nodesId[nodes_id[1]] << endl;
	if( debuging_map_mesh)
	{
		hout << " where is nodes:  " << wh[0] << " " << wh[1] << " " << wh[2] << " " << wh[3] << endl;
		hout << " is on interface: " << is_on[0] << " " << is_on[1] << " " << is_on[2] << " " << is_on[3] << endl;
		hout << " left_movable: " << left_movable << " right_movable: " << right_movable << endl; 
	}

	Point p1 = {	nodes_vec[eles_vec[ele_num].nodesId[nodes_id[0]]].x,
						nodes_vec[eles_vec[ele_num].nodesId[nodes_id[0]]].y,
						nodes_vec[eles_vec[ele_num].nodesId[nodes_id[0]]].z };
	Point p2 = {	nodes_vec[eles_vec[ele_num].nodesId[nodes_id[1]]].x,
						nodes_vec[eles_vec[ele_num].nodesId[nodes_id[1]]].y,
						nodes_vec[eles_vec[ele_num].nodesId[nodes_id[1]]].z };
	Point rp ;
	if( debuging_map_mesh)
	{
		hout << " node1: " << p1.x << " " << p1.y << " " << p1.z << endl;
		hout << " node2: " << p2.x << " " << p2.y << " " << p2.z << endl;
	}
	if( thecell->surfaces_vec[*pww]->intersect( p1, p2, rp ) )
	{
		int node1_num = eles_vec[ele_num].nodesId[nodes_id[0]] ;
		int node2_num = eles_vec[ele_num].nodesId[nodes_id[1]] ;
		int node3_num = eles_vec[ele_num].nodesId[nodes_id[2]] ;
		int node4_num = eles_vec[ele_num].nodesId[nodes_id[3]] ;

		Node node( rp );
		if( debuging_map_mesh) hout << " rp: " << rp.x << "  " << rp.y << "  " << rp.z << endl;

		//��齫Ҫ�ƶ��������ڵ�����ڽڵ��Ƿ���3���ڵ㶼�ڽ������ϵ����
		int is_all_on1 = is_re_all_on_face( node1_num, *pww );
		int is_all_on2 = is_re_all_on_face( node2_num, *pww );

		if( debuging_map_mesh) hout << " is_all_on1: " << is_all_on1 << " is_all_on2: " << is_all_on2 << endl;

		//������������ڵ������ܽ�ʱ���������Ӧ���ƶ��������ϣ���ʹis_all_on==1;
		double dis2_sur1 = thecell->surfaces_vec[*pww]->dis_to(&p1);
		double dis2_sur2 = thecell->surfaces_vec[*pww]->dis_to(&p2);
		double the_dis = global_length*3.0/8.0;
		if( dis2_sur1 < the_dis ) is_all_on1 = 0;
		if( dis2_sur2 < the_dis ) is_all_on2 = 0;

		if( is_all_on1 >= 1 ) left_movable = 0;
		if( is_all_on2 >= 1 ) right_movable = 0;

		if( debuging_map_mesh) hout	<< " node1: " << node1_num
														<< " node2: " << node2_num << endl;

		double dis1 = node.distance_to( &nodes_vec[node1_num] );
		double dis2 = node.distance_to( &nodes_vec[node2_num] );

		//ȷ���ƶ��ĸ��㵽������
		int move_which_node = 0;
		if( left_movable == 1 )
		{
			if( right_movable == 1 )
			{
				//���˽ڵ㶼���ƶ���ȡ����С���ĸ�
				if( dis2 < dis1 ) move_which_node = 1;  //�ƶ��Ҷ˽ڵ�
			}                            
		}
		else
		{
			if( right_movable == 1 )
			{
				move_which_node = 1;
			}
			else
			{
				move_which_node = -1;  //�������ƶ� 
			}
		}

		if( debuging_map_mesh) hout << "dis1: " << dis1 << " dis2: " << dis2 << " alength: " << alength << endl;
		int moved_nn = -1 ;
		if( move_which_node == 0 )
		{
			if( dis1 < alength )
			{
				if( debuging_map_mesh) hout << " dis1 < dis2, move node " << node1_num << " to the new position.\n";
				Node temp_n = nodes_vec[node1_num] ;//�����������Ա���ִ���ʱ�ָ�
				nodes_vec[node1_num].move_to( &node );
				double relative_min_volume = cal_relative_min_volume( node1_num );
				if( debuging_map_mesh) hout << " relative_min_volume: " << relative_min_volume << endl;
				if( relative_min_volume < min_volume )
				{
					nodes_vec[node1_num].move_to( &temp_n );
					if( debuging_map_mesh) hout << " but failed. " << endl;
				}
				else
				{
					where_is_nodes[node1_num] = *pww;
					is_on_nodes[node1_num] = 1;
					moved_nn = node1_num;
				}
			}
		}
		else if( move_which_node == 1 ) 
		{
			if( dis2 < alength )
			{
				if( debuging_map_mesh) hout << " dis1 > dis2, move node " << node2_num << " to the new position.\n";
				Node temp_n = nodes_vec[node2_num] ;//�����������Ա���ִ���ʱ�ָ�
				nodes_vec[node2_num].move_to( &node );
				double relative_min_volume = cal_relative_min_volume( node2_num );
				if( debuging_map_mesh) hout << " relative_min_volume: " << relative_min_volume << endl;

				if( relative_min_volume < min_volume )
				{
					nodes_vec[node2_num].move_to( &temp_n );  
					if( debuging_map_mesh) hout << " but failed. " << endl;   
				}
				else
				{
					where_is_nodes[node2_num] = *pww;
					is_on_nodes[node2_num] = 1;
					moved_nn = node2_num;
				}
			}
		}

		//�������Ľڵ㣨�ƶ����������²���Ľڵ�ţ��Ա�����Ĵ���
		int oper_nn=moved_nn;

		//�ƶ��˽ڵ�
		if( debuging_map_mesh )	hout << "moved node num: " << moved_nn << endl;

		if( moved_nn != -1 )
		{
			//���������ص�Ԫ�Ĳ�������
			for( int j=0;j <(int)nodes_vec[moved_nn].relative_eles_vec.size(); j++ )
			{
				int rele_num = nodes_vec[moved_nn].relative_eles_vec[j];
				if( deter_mat_ele(rele_num) == -1 )
				{
					put_into_cbev( rele_num, alength, ele_kind );
				}
			}
		}
		else
		{
			double threshold_len = global_length/4.0;
			if( dis1 < threshold_len )
			{
				if( is_on_nodes[node1_num] == 1 && where_is_nodes[node2_num] != -1 )
				{
					where_is_nodes[node1_num] = *pww;
					return 1;
				}
			}
			else if( dis2 < threshold_len )
			{
				if( is_on_nodes[node2_num] == 1 && where_is_nodes[node1_num] != -1 )
				{
					where_is_nodes[node2_num] = *pww;
					return 1;
				}
			}

			if( debuging_map_mesh) hout << " node: " << node.x << "  " << node.y << "  " << node.z << endl;
			int node_num = (int)nodes_vec.size();    //�����ɽڵ��
			if( debuging_map_mesh) hout << " generate a node, node_num: " << node_num << endl;

			nodes_vec.push_back( node );
			where_is_nodes.push_back( *pww );
			is_on_nodes.push_back( 1 );    

			vector<int> temp_v;
			temp_v.push_back( min(node1_num, node2_num) );
			temp_v.push_back( max(node1_num, node2_num) );
			temp_v.push_back( node_num );
			Tetrahedron new_tet( node_num, node3_num,
				node4_num, eles_vec[ele_num].nodesId[3] );
			eles_vec[ele_num].nodesId[nodes_id[1]] = node_num;

			int new_tet_num = (int)eles_vec.size();
			eles_vec.push_back( new_tet );
			nodes_vec[node_num].relative_eles_vec.push_back(ele_num);
			for( int k=0; k<4; k++ )
			{
				nodes_vec[eles_vec.back().nodesId[k]].relative_eles_vec.push_back( new_tet_num );
			}

			put_into_cbev( new_tet_num, alength, ele_kind );


			//���ÿ�������˱ߵ������壬��������ѳ��������Ա�֤��Ԫ��Э��
			harmonize_tets( node1_num, node2_num, node_num, alength, new_edges, ele_kind );

			oper_nn = node_num;
		}

		if( oper_nn >= 0 )
		{
			nodes_gen_later.push_back(oper_nn);
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//�������ڵ��������ص�Ԫ�Ƿ��ĸ��ڵ�ȫ�ڱ߽�����
int Mesher::is_re_all_on_face( int node_num, int wh )
{
	bool debuging = false;

	int aof_num=0; //��ص�Ԫ���ĸ��ڵ�ȫ�ڽ����ϵĵ�Ԫ����ͨ����Ҫ����������������Ѳ��������
	for( int i=0; i<(int)nodes_vec[node_num].relative_eles_vec.size(); i++ )
	{
		int ele_num = nodes_vec[node_num].relative_eles_vec[i];

		//����Ƿ��Ѿ���ɾ�����ĸ��ڵ����ͬ������������֪�����ж��ٵط�û�м�鵥Ԫ�Ƿ��Ѿ�ɾ����
		int is_deleted = 0;
		for( int i=1; i<4; i++ )
		{
			if( eles_vec[ele_num].nodesId[i] == eles_vec[ele_num].nodesId[0] )
			{
				is_deleted = 1;
			}
		}
		if( is_deleted == 1 ) continue;

		int is_all_on=1;
		if( wh == -1 ) wh = where_is_nodes[eles_vec[ele_num].nodesId[0]];
		for( int i=0; i<4; i++ )
		{
			int nn = eles_vec[ele_num].nodesId[i];
			if( nn == node_num ) continue;
			if( where_is_nodes[nn] != wh || is_on_nodes[nn] == 0 )
			{
				is_all_on = 0;
				break;
			}
		}
		if( is_all_on == 1 )
		{
			aof_num++;
			if( debuging )  hout	<< " ele: " << ele_num  << " : "
											<< eles_vec[ele_num].nodesId[0] << " "
											<< eles_vec[ele_num].nodesId[1] << " "
											<< eles_vec[ele_num].nodesId[2] << " "
											<< eles_vec[ele_num].nodesId[3] << endl;
		}
	}
	return aof_num;
}

//---------------------------------------------------------------------------
//���ݸ����ıߵ���ţ�ȷ���ڵ����У����е�Ԫ��ֵ�ʱ��ʹ��
int* Mesher::deter_nodes_id(int edge_num, int* nodes_id)
{
	int node1_i,node2_i,node3_i,node4_i;
	switch( edge_num )
	{
	case 1:
		node1_i = 0 ;
		node2_i = 1 ;
		node3_i = 1 ;
		node4_i = 2 ;
		break;
	case 2:
		node1_i = 0;
		node2_i = 2;
		node3_i = 1;
		node4_i = 2;
		break;
	case 3:
		node1_i = 0;
		node2_i = 3;
		node3_i = 1;
		node4_i = 2;
		break;
	case 4:
		node1_i = 1;
		node2_i = 2;
		node3_i = 2;
		node4_i = 0;
		break;
	case 5:
		node1_i = 1;
		node2_i = 3;
		node3_i = 2;
		node4_i = 0;
		break;
	case 6:
		node1_i = 2;
		node2_i = 3;
		node3_i = 0;
		node4_i = 1;
		break;
	default:
		break;
	}
	nodes_id[0] = node1_i;
	nodes_id[1] = node2_i;          
	nodes_id[2] = node3_i;
	nodes_id[3] = node4_i;

	return nodes_id;
}

//---------------------------------------------------------------------------
//���ÿ����ص�Ԫ�������壩�����ڽڵ�Node1_num �� Node2_num֮�����ڵ�node3_num
void Mesher::harmonize_tets( int node1_num, int node2_num, int node3_num,
											 double alength, vector<Bar> &new_edges, int ele_kind )
{
	bool debuging = false;

	//���ڵ�Node1_num��������ص�Ԫ
	for( int i=0; i<(int)nodes_vec[node1_num].relative_eles_vec.size(); i++ )
	{
		//���õ�Ԫ�Ƿ�����������ڵ�
		int ele_num = nodes_vec[node1_num].relative_eles_vec[i];
		int node1_i = eles_vec[ele_num].is_contain( node1_num );
		if( node1_i == -1 )	 continue;
		int node2_i = eles_vec[ele_num].is_contain( node2_num );
		if( node2_i == -1 ) continue;
		if( debuging )
		{
			hout << "harmonize tet: " << ele_num ;
			hout << " nodes: " << eles_vec[ele_num].nodesId[0] << " " ;
			hout << eles_vec[ele_num].nodesId[1] << " " ;
			hout << eles_vec[ele_num].nodesId[2] << " " ;
			hout << eles_vec[ele_num].nodesId[3] << endl ;
		}     

		//���������ɵı�
		for( int k=0; k<4; k++ )
		{
			if( k == node1_i || k == node2_i ) continue;
			Bar newbar(node3_num,eles_vec[ele_num].nodesId[k]);
			new_edges.push_back(newbar);
		}

		//ȷ����������
		if( node1_i > node2_i )
		{
			int temp_i = node1_i ;
			node1_i = node2_i ;
			node2_i = temp_i ;
		}

		Tetrahedron *new_tet;
		if( node1_i == 0 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[1],
				eles_vec[ele_num].nodesId[2], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}
		else if( node1_i == 1 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[2],
				eles_vec[ele_num].nodesId[0], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}
		else if( node1_i == 2 )
		{
			Tetrahedron new_tet_t( node3_num, eles_vec[ele_num].nodesId[0],
				eles_vec[ele_num].nodesId[1], eles_vec[ele_num].nodesId[3] );
			new_tet = &new_tet_t ;
			eles_vec[ele_num].nodesId[node2_i] = node3_num;
		}

		int new_tet_num = (int)eles_vec.size();
		eles_vec.push_back( *new_tet );
		nodes_vec[node3_num].relative_eles_vec.push_back( ele_num );
		for( int k=0; k<4; k++ )
		{
			nodes_vec[eles_vec.back().nodesId[k]].relative_eles_vec.push_back( new_tet_num );
		}

		if( deter_mat_ele(ele_num) == -1 )
		{
			put_into_cbev( ele_num, alength, ele_kind );
		}
		if( deter_mat_ele(new_tet_num) == -1 )
		{
			put_into_cbev( new_tet_num, alength, ele_kind );
		}
	}
}

//---------------------------------------------------------------------------
//���ÿ������� (��С����
void Mesher::quick_sort( vector<int> *int_vec, int min_n, int max_n )
{
	if( min_n == -1 ) min_n = 0;
	if( max_n == -1 ) max_n = (int)int_vec->size()-1;
	if( max_n <= 0 ) return;
	int middle = (*int_vec)[(int)((min_n + max_n)/2)];
	int i = min_n;
	int j = max_n;
	do{
		while(((*int_vec)[i]<middle) && (i<max_n))//����ɨ�������ֵ����
		{       
			i++;
		}
		while(((*int_vec)[j]>middle) && (j>min_n))//����ɨ��С����ֵ����
		{        
			j--;
		}
		if(i<=j)//�ҵ���һ��ֵ������
		{                         
			int temp_v = (*int_vec)[i];
			(*int_vec)[i]=(*int_vec)[j];
			(*int_vec)[j]=temp_v;
			i++;
			j--;
		}
	}while( i <= j );//�������ɨ����±꽻����ֹͣ�����һ�Σ�

	//����߲�����ֵ(left<j)���ݹ�����
	if( min_n < j ) quick_sort( int_vec, min_n, j );
	//���ұ߲�����ֵ(right>i)���ݹ��Ұ��
	if( max_n > i ) quick_sort( int_vec, i, max_n );
}

//---------------------------------------------------------------------------
//�������е�Ԫ��С�����������������min_max_en[2]������С�������ĵ�Ԫ�ţ�min_max_vl[2]������С������ֵ��
int Mesher::cal_min_max_tet_volume(int min_max_en[2], double min_max_vl[2])
{
	//���ÿ����Ԫ�����
	double mmin_volume = 1.0e88 ;//��¼��С�����
	int min_vol_ele_num = -1;
	double mmax_volume = 0.0 ;		//��¼��С�����
	int max_vol_ele_num = -1;
	double volume_ele;
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		volume_ele = cal_tet_volume( &eles_vec[i] );
		if( volume_ele < mmin_volume )
		{
			mmin_volume = volume_ele;
			min_vol_ele_num = i ;
		}
		else if( volume_ele > mmax_volume )
		{
			mmax_volume = volume_ele;
			max_vol_ele_num = i ;
		}
	}
	min_max_en[0] = min_vol_ele_num ;
	min_max_en[1] = max_vol_ele_num ;
	min_max_vl[0] = mmin_volume ;
	min_max_vl[1] = mmax_volume ;
	return 1;
}

//---------------------------------------------------------------------------
//�������е�Ԫ������ֲ��������С�������ȡ��������������ȡ������¼���ȣ�
//����γ��ȷֳ�n�ݣ�һ�����nȡ10������¼ÿ������ռ��Ԫ�����ܵ�Ԫ���ı�ֵ��
void Mesher::tet_vol_distrib( double min_max_vl[2], string str, int n, int mod )
{
	//���ڼ�¼ÿ������������ĸ���
	vector<int> tetnum(n,0) ;
	double vol_min=floor(min_max_vl[0]);
	double vol_max=ceil(min_max_vl[1]);
	double vol_ele=(vol_max-vol_min)/n;
	double volume_ele;
	int m;

	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		volume_ele = cal_tet_volume( &eles_vec[i] );
		m=(int)((volume_ele-vol_min)/vol_ele);
		tetnum[m]=tetnum[m]+1;
	}

	if(mod==1)
	{
		string file="tet_vol_distrib"+str+".dat";
		ofstream disout( file.c_str() ) ;
		if ( !disout )  hout << "�޷����ļ���" << file << endl;
		disout << "����������ֲ�" << endl;
		disout << "��С�����" << min_max_vl[0] << endl;
		disout << "��������" << min_max_vl[1] << endl;
		m=(int)eles_vec.size();
		disout << "�����������" << m << endl;
		for( int i=0; i<n; i++ )
		{
			disout << i <<" :  ";
			disout << "�����" << vol_min+i*vol_ele << "~" << vol_min+(i+1)*vol_ele << "֮�䣺"; 
			disout << "�����������" << tetnum[i];
			disout << "�ٷֱȣ�" <<	100*(tetnum[i]*1.0/m) << "%" << endl;
		}
	}
}

//---------------------------------------------------------------------------
//�������������嵥Ԫ�����������ֲ��������Сֵ0�����ֵȡ1����¼���ȣ�
//����γ��ȷֳ�n�ݣ�һ�����nȡ10������¼ÿ������ռ��Ԫ�����ܵ�Ԫ���ı�ֵ��
void Mesher::tet_quality( string str, int n, int mod )
{
	//���ڼ�¼ÿ������������ĸ���
	vector<int> tetnum(n,0) ;
	double tria[4];	//���ڼ�¼�������ĸ�������(triangle_area)
	double silp[3];	//���ڼ�¼������Ա߱߳�֮��(side_length_product)
	double ir;			//������뾶
	double R;			//�����뾶
	double rou;       //��ֵrou
	double volume_ele; //���������
	int m;

	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		//�������������
		volume_ele = cal_tet_volume( &eles_vec[i] );
		//�����������ĸ�������
		tet_tri_area( &eles_vec[i] , tria );
		//����������뾶
		ir=(3.0*volume_ele)/(tria[0]+tria[1]+tria[2]+tria[3]);
		//����������Ա߱߳�֮��
		side_length_pro(&eles_vec[i] , silp );
		//���������뾶
		R=sqrt((silp[0]+silp[1]+silp[2])*(silp[0]+silp[1]-silp[2])*(silp[0]+silp[2]-silp[1])*(silp[1]+silp[2]-silp[0]))/(24*volume_ele);
		//�����ֵrou
		rou=3*ir/R;
		if(rou<0.0-ZERO||rou>1.0+ZERO)
		{
			hout << "�ڼ����" << i << "����Ԫʱ,��ֵrou=" << rou <<" С��0.0�����1.0,����" <<endl;
			hout << "ir=" << ir << "  " << "R=" << R << endl;
		}
		else
		{
			m=(int)(rou*n);
			tetnum[m]=tetnum[m]+1;
		}
	}

	if(mod==1)
	{
		string file="tet_quality"+str+".dat";
		ofstream disout( file.c_str() ) ;
		if ( !disout )  hout << "�޷����ļ���" << file << endl;
		disout << "����������ϵ���ֲ�" << endl;
		m=(int)eles_vec.size();
		disout << "�����������" << m << endl;
		for( int i=0; i<n; i++ )
		{
			disout << i <<" :  ";
			disout << "����ϵ����" << i*(1.0/n) << "~" << (i+1)*(1.0/n) << "֮�䣺"; 
			disout << "�����������" << tetnum[i];
			disout << "�ٷֱȣ�" <<	100*(tetnum[i]*1.0/m) << "%" << endl;
		}
	}
}
//---------------------------------------------------------------------------
//���ڼ����������ĸ�������
void Mesher::tet_tri_area( Tetrahedron *tet , double tria[4] )
{
	Node node1, node2, node3;
	double p,a,b,c;
	for(int i=0; i<4; i++)
	{
		switch(i)
		{
		case 0:	node1=nodes_vec[tet->nodesId[0]];
					node2=nodes_vec[tet->nodesId[1]];
					node3=nodes_vec[tet->nodesId[2]];
					break;
		case 1:	node1=nodes_vec[tet->nodesId[0]];
					node2=nodes_vec[tet->nodesId[1]];
					node3=nodes_vec[tet->nodesId[3]];
					break;
		case 2:	node1=nodes_vec[tet->nodesId[0]];
					node2=nodes_vec[tet->nodesId[2]];
					node3=nodes_vec[tet->nodesId[3]];
					break;
		case 3:	node1=nodes_vec[tet->nodesId[1]];
					node2=nodes_vec[tet->nodesId[2]];
					node3=nodes_vec[tet->nodesId[3]];
					break;
		default: break;
		}
		a=sqrt((node1.x-node2.x)*(node1.x-node2.x)+(node1.y-node2.y)*(node1.y-node2.y)+(node1.z-node2.z)*(node1.z-node2.z));
		b=sqrt((node1.x-node3.x)*(node1.x-node3.x)+(node1.y-node3.y)*(node1.y-node3.y)+(node1.z-node3.z)*(node1.z-node3.z));
		c=sqrt((node2.x-node3.x)*(node2.x-node3.x)+(node2.y-node3.y)*(node2.y-node3.y)+(node2.z-node3.z)*(node2.z-node3.z));
		p=(a+b+c)/2.0;

		tria[i]=sqrt(p*(p-a)*(p-b)*(p-c));
	}
}
//---------------------------------------------------------------------------
//���ڼ���������Ա߱߳�֮��
void Mesher::side_length_pro( Tetrahedron* tet , double silp[3] )
{
	Node node[4];
	double x,y,z;
	double len[2];
	int n[2][2];


	for(int i=0; i<4; i++)
	{
		node[i]=nodes_vec[tet->nodesId[i]];
	}
	int k;
	for(int i=1; i<4; i++)
	{
		k=0;
		//��������ֱ��ǣ�(0,1;2,3),(0,2;1,3),(0,3;1,2), ��n[2][2]��¼
		n[0][0]=0;
		n[0][1]=i;
		for(int j=1; j<4; j++)
		{
			if(j!=i)
			{
				n[1][k]=j;
				k++;
			}
		}
		for(int j=0; j<2; j++)
		{
			x=node[n[j][0]].x-node[n[j][1]].x;
			y=node[n[j][0]].y-node[n[j][1]].y;
			z=node[n[j][0]].z-node[n[j][1]].z;
			len[j]=sqrt(x*x+y*y+z*z);
		}
		silp[i-1]=len[0]*len[1];
	}
}
//---------------------------------------------------------------------------
//���������������������ֵ��ظ�Ԫ��
void Mesher::unique( vector<int> *int_vec )
{
	if( int_vec->size() < 2 ) return ;
	//���Ƹ���
	vector<int> int_vec_cop = *int_vec;
	//�������
	int_vec->clear();
	int_vec->push_back( int_vec_cop[0] ) ;
	for( int i=1; i<(int)int_vec_cop.size(); i++ )
	{
		if( int_vec_cop[i] != int_vec_cop[i-1] )
		{
			int_vec->push_back( int_vec_cop[i] );
		}
	}
}

//---------------------------------------------------------------------------
//mod=0:������˹�⻬����������һ�ڵ㣬ȡ������Ϊ��֮���ڵĽڵ�����ƽ��ֵ
//mod=1:ȡ���ڽڵ���Χ���������
int Mesher::smoothing_3d(int mod)
{
	bool debuging_smoothing = false;
	int min_vol_en;
	double min_min_volume_ori=min_volume_of_cell(min_vol_en);  //��˳ǰ��С���
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		if( debuging_smoothing ) hout	<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_smoothing ) hout	<< "node num: " << i << endl;
		if( debuging_smoothing ) hout	<< "old location: " << nodes_vec[i].x << " "
														<< nodes_vec[i].y << " "
														<< nodes_vec[i].z << endl;
		double min_volume_ori = cal_relative_min_volume( i ); //��˳ǰ�ýڵ���ص�Ԫ����С���
		if( debuging_smoothing ) hout << " min_volume_ori: " << min_volume_ori << endl;
		double sum_x=0, sum_y=0, sum_z=0; //���˽ڵ���ڽڵ��������
		int num_node=0;                   //�����ٽڵ��������
		//���ڱ߽��ϵĽڵ㣬ֻ����ͬ���߽��ϵĽڵ㣬���ܶ�����ƽ��
		//������������߽����ϵĽڵ㣬�߽����ϵĽڵ�
		int ioe =  is_on_edge( i ) ;
		int iof =  is_on_face( i ) ;
		vector<int> ioe_nodes; //�͵�ǰ�ڵ���ͬһ�߽����ϵĽڵ㼯
		double x_min=1.0e200,x_max=-1.0e200,y_min=1.0e200,y_max=-1.0e200,z_min=1.0e200,z_max=-1.0e200;

		if( debuging_smoothing ) hout << "ioe: " << ioe << " iof: " << iof << endl;
		if( debuging_smoothing ) hout << " neigh_nodes: " << endl;

		vector<int> nei_nodes;
		for( int j=0; j<(int)nodes_vec[i].neigh_nodes_vec.size(); j++ )
		{
			if( debuging_smoothing )
			{
				hout << "No. " << j << " neighbor node: " << nodes_vec[i].neigh_nodes_vec[j] << endl;
			}
			if( ioe != 0 )
			{
				if( is_on_same_edge( i, nodes_vec[i].neigh_nodes_vec[j] ) == 0 )	continue ;
				ioe_nodes.push_back( nodes_vec[i].neigh_nodes_vec[j] );
			}
			else if( iof != 0 )
			{
				if( is_on_same_face( i, nodes_vec[i].neigh_nodes_vec[j] ) == 0 )	continue ;                                   
			}
			double x_temp=nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].x;
			double y_temp=nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].y;
			double z_temp=nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].z;
			if(x_temp<x_min)x_min=x_temp;
			if( x_temp>x_max)x_max=x_temp;
			if(y_temp<y_min)y_min=y_temp;
			if( y_temp>y_max)y_max=y_temp;
			if(z_temp<z_min)z_min=z_temp;
			if( z_temp>z_max)z_max=z_temp;
			sum_x += nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].x ;
			sum_y += nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].y ;
			sum_z += nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].z ;
			num_node ++ ;
			nei_nodes.push_back(nodes_vec[i].neigh_nodes_vec[j]);
			if( debuging_smoothing )
			{
				hout	<< nodes_vec[i].neigh_nodes_vec[j] << " " 
						<< nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].x << " "
						<< nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].y << " "
						<< nodes_vec[nodes_vec[i].neigh_nodes_vec[j]].z << endl;
			}
		}
		if( num_node > 1 )
		{
			if( debuging_smoothing ) hout	<< "mod: " << mod << endl;
			if( debuging_smoothing ) hout	<< "num_node: " << num_node << endl;
			if( debuging_smoothing ) hout	<< "xyz_min_max: " << x_min << " " << x_max << " "
															<< y_min << " " << y_max << " "
															<< z_min << " " << z_max << endl;
			Point the_p = { sum_x / num_node, sum_y / num_node, sum_z / num_node };
			if( mod == 1 )
			{
				the_p.x = (x_min+x_max)/2.0;
				the_p.y = (y_min+y_max)/2.0;
				the_p.z = (z_min+z_max)/2.0;
			}
			if( debuging_smoothing )
			{
				hout << "Before projection the point: " << the_p.x << " " << the_p.y << " " << the_p.z << endl;
			}

			Node old_position = nodes_vec[i];	//�����������Ա���ִ���ʱ�ָ�
			if( ioe != 0 )
			{
				if( debuging_smoothing ) hout << " on the edge! " << endl;
				if( num_node > 2 ) continue;			//�ǵ�
				//�ҳ���˽ڵ����ڵ�����ͬ�߽����ϵĽڵ����ߵ��д������������Ľ���
				Node *node1 = &nodes_vec[ioe_nodes[0]];
				Node *node2 = &nodes_vec[ioe_nodes[1]];
				Node *node3 = &nodes_vec[i];
				if( debuging_smoothing )
				{
					hout << "node1: " << node1->x << " " << node1->y << " " << node1->z << endl;
					hout << "node2: " << node2->x << " " << node2->y << " " << node2->z << endl;
					hout << "node3: " << node3->x << " " << node3->y << " " << node3->z << endl;
				}
				int w1 = where_is_nodes[i];
				int is_on1 = is_on_nodes[i];
				if( is_on1 == 0 )
				{
					//�ڵ����߽�����
					Node temp_node = (*node1 + *node2)/2.0;			//����20061221��
					node3->x = temp_node.x;
					node3->y = temp_node.y;
					node3->z = temp_node.z;
					continue ;
				}
				//������ͻ���Ľ�������
				Node  node4 = ( *node1 + *node2 )/2.0 ;
				Node node31 = *node3 - *node1 ;
				Node node41 =  node4 - *node1 ;
				Node node21 = *node2 - *node1 ;
				double n3121 = ( node31.x * node21.x + node31.y * node21.y + node31.z * node21.z );
				if( n3121 == 0 ) continue;
				double tt = ( node41.x * node21.x + node41.y * node21.y + node41.z * node21.z )/n3121;
				Node node_31_mid = node31*tt + *node1 ;
				Node node32 = *node3 - *node2 ;
				Node node42 =  node4 - *node2 ;
				double n3221 = ( node32.x * node21.x + node32.y * node21.y + node32.z * node21.z );
				if( n3221 == 0 ) continue;
				tt = ( node42.x * node21.x + node42.y * node21.y + node42.z * node21.z )/n3221;
				Node node_32_mid = node32*tt + *node2 ;
				//����������Ľ���
				Point p1 = { node_31_mid.x,  node_31_mid.y,  node_31_mid.z };
				Point p2 = { node_32_mid.x,  node_32_mid.y,  node_32_mid.z };
				Point rp ;
				if( thecell->surfaces_vec[w1]->intersect( p1, p2, rp ) == 1 )
				{
					if( debuging_smoothing ) hout	<< "new location0: " 
						<< rp.x << " "<< rp.y << " "<< rp.z << endl;
					node3->move_to( &rp );
				}
			}
			else if( iof > 0 )
			{
				if( debuging_smoothing ) hout << " on the face No. " << iof << endl;
				if( num_node > 2 )
				{
					//����ڽڵ㳬��3�����������⻬�Ľڵ㣬���磬����ڵ�������
					//��ֻ����ͬһ���ϵ��ڽڵ�ſ������⻬����
					//�������ڽڵ㣬����ѡȡ3�������߷���ѭ�������һ���ڽڵ�
					//Ȼ������еķ���ȡƽ������ΪͶӰ����
					vector<TDVector> vecs;
					double dxx=nodes_vec[nei_nodes[1]].x-nodes_vec[nei_nodes[0]].x;
					double dyy=nodes_vec[nei_nodes[1]].y-nodes_vec[nei_nodes[0]].y;
					double dzz=nodes_vec[nei_nodes[1]].z-nodes_vec[nei_nodes[0]].z;
					TDVector vec1(dxx,dyy,dzz);
					dxx=nodes_vec[nei_nodes[2]].x-nodes_vec[nei_nodes[1]].x;
					dyy=nodes_vec[nei_nodes[2]].y-nodes_vec[nei_nodes[1]].y;
					dzz=nodes_vec[nei_nodes[2]].z-nodes_vec[nei_nodes[1]].z;
					TDVector vec2(dxx,dyy,dzz);
					TDVector vec3=vec1.cro_product(&vec2);
					vecs.push_back(vec3);
					if( num_node > 3)
					{
						for( int k=1; k<num_node; k++ )
						{
							int k1=(k+1)%num_node;
							double dxx=nodes_vec[nei_nodes[k1]].x-nodes_vec[nei_nodes[k]].x;
							double dyy=nodes_vec[nei_nodes[k1]].y-nodes_vec[nei_nodes[k]].y;
							double dzz=nodes_vec[nei_nodes[k1]].z-nodes_vec[nei_nodes[k]].z;
							TDVector vec11(dxx,dyy,dzz);
							int k2=(k+2)%num_node;
							dxx=nodes_vec[nei_nodes[k2]].x-nodes_vec[nei_nodes[k1]].x;
							dyy=nodes_vec[nei_nodes[k2]].y-nodes_vec[nei_nodes[k1]].y;
							dzz=nodes_vec[nei_nodes[k2]].z-nodes_vec[nei_nodes[k1]].z;
							TDVector vec22(dxx,dyy,dzz);
							TDVector vec33=vec1.cro_product(&vec2);

							if( debuging_smoothing )
							{
								hout << " angle: " << vec33.angle_between(&vecs[k-1]) << endl;
								hout << "       " << vec33.x << " " << vec33.y << " " << vec33.z << endl;
								hout << "       " << vecs[k-1].x << " " << vecs[k-1].y << " " << vecs[k-1].z << endl;
							}

							if(vec33.angle_between(&vecs[k-1])>PI/2.0)
							{
								TDVector ZERO_vec;
								vec33=ZERO_vec-vec33;
							}
							vecs.push_back(vec33);
						}
					}
					TDVector project_vec=vecs[0];
					if( vecs.size() > 0 )
					{
						for( int k=1; k<(int)vecs.size(); k++ )
						{
							project_vec=project_vec+vecs[k];
						}
						project_vec=project_vec/vecs.size();
					}
					if( debuging_smoothing )
					{
						hout << " project_vec: " ;
						hout << project_vec.x << " " ;
						hout << project_vec.y << " " ;
						hout << project_vec.z << endl;
					}
					Point pp = thecell->surfaces_vec[iof-1]->project( &the_p, &project_vec );
					nodes_vec[i].move_to( &pp );
				}
				else
				{
					Point pp = thecell->surfaces_vec[iof-1]->project_normal( &the_p );
					nodes_vec[i].move_to( &pp );
				}
			}
			else
			{
				if( debuging_smoothing ) hout << " neither on the edge nor on the face! " << endl;
				nodes_vec[i].move_to( &the_p );
			}
			if( debuging_smoothing ) hout << "new location: " 
				<< nodes_vec[i].x << " "<< nodes_vec[i].y << " "<< nodes_vec[i].z << endl;

			double min_volume_new = cal_relative_min_volume( i );
			if( min_volume_new < min_volume_ori )	//��ص�Ԫ��С�����С���ָ�
			{
				nodes_vec[i].move_to( &old_position );
				if( debuging_smoothing ) hout << " re_min_volume: oled: " << min_volume_ori
					<< " new: " << min_volume_new << endl ;
				if( debuging_smoothing ) hout << " smoothing failed at node " << i << endl;
			}
		}
		if( debuging_smoothing )
		{
			double min_volume_new = cal_relative_min_volume( i );
			hout << " min_volume_new: " << min_volume_new << endl;
		}
	}       
	double min_min_volume_new=min_volume_of_cell(min_vol_en);  //��˳����С���

	if( debuging_smoothing ) hout << " min_min_volume_ori: " << min_min_volume_ori << " min_min_volume_new: " << min_min_volume_new << endl;
	if( min_min_volume_new - min_min_volume_ori > fabs(min_min_volume_ori) * 0.01 ) return 1;

	return 0;
}

//---------------------------------------------------------------------------
//�������е�Ԫ�е���С���
double Mesher::min_volume_of_cell(int &ele_num)
{
	double mmin_volume = 1e80;
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		double volume = cal_tet_volume( &eles_vec[i] );
		if( volume < mmin_volume )
		{
			ele_num = i;
			mmin_volume = volume ;
		}
	}
	return mmin_volume ;
}

//---------------------------------------------------------------------------
//������ά�����
double Mesher::volume_ratio()
{
	double non_matrix_vol = 0.0;
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		if( eles_vec[i].materialId != 0 )
		{
			non_matrix_vol += cal_tet_volume( &eles_vec[i] );
		}
	}
	double total_vol = thecell->clength*thecell->cwidth*thecell->cheight;

	//hout << "    ��ǿ�������: " << non_matrix_vol << " �������: " << total_vol << endl;  
	return non_matrix_vol/total_vol;
}

//---------------------------------------------------------------------------
//ȷ�������߽��ϵĽڵ�
void Mesher::deter_boundary_nodes(int coat)
{
	bnodes_vec.clear();

	for ( int i=0; i<(int)nodes_vec.size(); i++)
	{
		int iof = is_on_face(i);   //���ڵ�λ��
		if( iof < 0 )
		{
			//�ڵ����߽���
			bnodes_vec.push_back(i);
			nodes_vec[i].flag = 1;
		}
		else if( iof == 0 )
		{
			//�ڵ����ڲ������ڵ����߽����ϣ�Ҳ���ڿ���������
			nodes_vec[i].flag = 0;
		}
		else
		{
			//�ڿ���������
			nodes_vec[i].flag = 0;			//������2��ʾ�ڿ���������
														//Ϊ�˲������ĳ����ͻ�������ȸ�Ϊ0
		}
	}
}
//---------------------------------------------------------------------------
//���ɱ߽������
int Mesher::gen_coating_mesh(double thick_ratio, int layer_num)
{
	vector<EllipsGeometry> ellgeo(thecell->surfaces_vec.size());

	//-------------------------------------------------------------------------------------------------------------
	//���������������ڽڵ㡢��Ԫ������
	int ele_size = int(eles_vec.size());
	vector<int>  eles_ellip(ele_size,-1);
	for(int i=0; i<ele_size; i++)
	{
		//if(i==566||i==571||i==1146||i==1148)
		//{
		//	hout << "i= " << i << "  mat= " <<  eles_vec[i].materialId<<endl;
		//	hout << "nodesId��";
		//	for(int k=0 ; k<4; k++)
		//	{
		//		hout << eles_vec[i].nodesId[k]+1 <<"  ";
		//	}
		//	hout << endl;
		//	hout << "where_is_nodes��";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		hout << where_is_nodes[eles_vec[i].nodesId[j]] <<"  ";
		//	}
		//	hout << endl;
		//	hout << "is_on_nodes��";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		hout << is_on_nodes[eles_vec[i].nodesId[j]] <<"  ";
		//	}
		//	hout << endl;
		//	hout << "where_is��";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		int is_on=0;
		//		hout << where_is(&nodes_vec[eles_vec[i].nodesId[j]], is_on) << "/" << is_on <<"  ";
		//	}
		//	hout << endl;
		//}
		if(eles_vec[i].materialId==1)		//��Ԫ���ڿ���
		{
			//---------------------------------------------------------------------
			//�жϴ˵�Ԫ�������ĸ�����
			int ele_num;
			int wws[4];
			for(int j=0; j<4; j++)
			{
				wws[j] = 	where_is_nodes[eles_vec[i].nodesId[j]];
			}
			if(wws[0]==wws[1]&&wws[0]==wws[2]&&wws[0]==wws[3])		//�ĸ�������������һ��
			{
				ele_num = wws[0];
			}
			else																								//�ĸ�������λ�ñ�ʶ��һ��
			{
				//����õ�Ԫ�����ĵ�����
				Point pp={0.0, 0.0, 0.0};
				for(int j=0; j<4; j++)
				{
					pp.x = pp.x + nodes_vec[eles_vec[i].nodesId[j]].x;
					pp.y = pp.y + nodes_vec[eles_vec[i].nodesId[j]].y;
					pp.z = pp.z + nodes_vec[eles_vec[i].nodesId[j]].z;
				}
				pp.x=pp.x/4.0;
				pp.y=pp.y/4.0;
				pp.z=pp.z/4.0;
				//��˵�Ԫ����������
				int count=0;
				vector<int> tem_num;
				for(int j=0; j<4; j++)
				{
					int wn=wws[j];
					if(wn==-1) continue;												//ȥ������-1�������-1��ʾ����

					int key=0;																//�ж���������Ƿ��Ѿ����жϹ���tem_num�м�¼�жϹ���������
					for(int k=0; k<int(tem_num.size()); k++)
					{
						if(wn==tem_num[k])
						{
							key=1;
							break;
						}
					}
					if(key==1)		continue;
					else	tem_num.push_back(wn);

					if(thecell->surfaces_vec[wn]->is_contain(&pp)<=0)		//�ж������������Ƿ��ڸ�������
					{
						ele_num = wn;
						count++;
					}
				}
				if(count==0)
				{
					//�ж���������Ȼ������������ڸ�������ⲿ��
					//������������ⶥ�����ڻ����У�������Ϊ�˵�Ԫ�����ڸ�����
					if(int(tem_num.size())==1)
					{
						ele_num = tem_num[0];
					}
					else if(int(tem_num.size())>=2)
					{
						eles_vec[i].materialId=0;    //�˵�Ԫ����������������֮��,��Ϊ�˵�Ԫ�ǻ������
						continue;
					}
					else
					{
						hout << "i= " << i << "  mat= " <<  eles_vec[i].materialId<<endl;
						hout << "nodesId��";
						for(int k=0 ; k<4; k++)
						{
							hout << eles_vec[i].nodesId[k]+1 <<"  ";
						}
						hout << endl;
						hout << "where_is_nodes��";
						for(int j=0 ; j<4; j++)
						{
							hout << where_is_nodes[eles_vec[i].nodesId[j]] <<"  ";
						}
						hout << endl;
						hout << "is_on_nodes��";
						for(int j=0 ; j<4; j++)
						{
							hout << is_on_nodes[eles_vec[i].nodesId[j]] <<"  ";
						}
						hout << endl;
						hout << "where_is��";
						for(int j=0 ; j<4; j++)
						{
							int is_on=0;
							hout << where_is(&nodes_vec[eles_vec[i].nodesId[j]], is_on) << "/" << is_on <<"  ";
						}
						hout << endl;
						hout << "�����и������岻�����κ�����" << endl;
						return 0;
					}
				}
				else if(count>1)
				{
					hout << "�����и�������ͬʱ���ڼ�������" << endl;
					return 0;
				}
			}
			//�궨��Ԫ�����ĸ�����
			eles_ellip[i] = ele_num;
			//����������嵥Ԫ���ݵ��������
			ellgeo[ele_num].tet.push_back(i);
		}
		else if(eles_vec[i].materialId==-1)		//��Ԫ�Ȳ��ǿ���Ҳ���ǻ���
		{
			eles_vec[i].materialId=0;
		}
	}

	//�������е�����Ԫ���ɱ���Ƭ������ڵ�
	int node_size = int(nodes_vec.size());
	vector<int> nod_ellip(node_size,-2);
	int ellgeo_size = int(ellgeo.size());
	for(int i=0; i<ellgeo_size; i++)
	{
		vector<vector<int > > line_vec;		//�����жϱ�����������Ƭ�Ƿ���
		int tet_size = int(ellgeo[i].tet.size());
		for(int j=0; j<tet_size; j++)
		{
			vector<vector<int> > temp_tri;	//��ʱ��¼��������Ƭ��������������������������Ƭ��
			int ell_tet = ellgeo[i].tet[j];			//��Ԫ���
			int nid[4];									//��Ԫ�ĸ��ڵ�ı��
			for(int k=0; k<4; k++)
			{
				nid[k] =  eles_vec[ell_tet].nodesId[k];
			}
			//��¼��������ĸ���Ƭ���ڸ�����������������Ƭ
			for(int k=0; k<2; k++)					//ѭ�����������嵥Ԫ���ĸ�����
			{
				for(int l=k+1; l<3; l++)
				{
					for(int m=l+1; m<4; m++)
					{
						int key;
						if(is_on_nodes[nid[k]]==0||is_on_nodes[nid[l]]==0||is_on_nodes[nid[m]]==0)   //�жϲ���¼������������Ƭ
						{
							if(deter_node_flag(nid[k])>0&&deter_node_flag(nid[l])>0&&deter_node_flag(nid[m])>0)		//������������ڵ�����
							{
								key=0;
							}
							else																																	//������������ڵ�������
							{
								key=1;
							}
						}
						else
						{
							key=1;
						}
						if(key==1)			//�жϴ����Ƿ�������������
						{
							vector<int> rel_tet;
							int nk_rel_size = int(nodes_vec[nid[k]].relative_eles_vec.size());					//���ýڵ���������������
							for(int n=0; n<nk_rel_size; n++)
							{
								int rel_ele_num = nodes_vec[nid[k]].relative_eles_vec[n];
								if(rel_ele_num!=ell_tet&&eles_ellip[rel_ele_num] == i)							//�õ�Ԫ�����㣬���ҵ��ĵ�ԪӦ�����ڴ�����
								{
									int nl_rel_size = int(nodes_vec[nid[l]].relative_eles_vec.size());
									for(int p=0; p<nl_rel_size; p++)
									{
										if(rel_ele_num==nodes_vec[nid[l]].relative_eles_vec[p])					
										{
											int nm_rel_size = int(nodes_vec[nid[m]].relative_eles_vec.size());
											for(int q=0; q<nm_rel_size; q++)
											{
												if(rel_ele_num==nodes_vec[nid[m]].relative_eles_vec[q])
												{
													rel_tet.push_back(rel_ele_num);											//���ִ��������Ǳ���Ԫ��������Ԫ�Ĺ�����Ƭ
													goto label_prism;
												}
											}
										}
									}
								}
							}
label_prism:			if(rel_tet.empty())		//����������Ƭû���ҵ��������������
							{
								vector<int> triangle(3,0);
								triangle[0] = k;
								triangle[1] = l;
								triangle[2] = m;
								temp_tri.push_back(triangle);
							}
						}
					}
				}
			}
			//������������ݵ��������
			int num_tri = int(temp_tri.size());
			if(num_tri==1)  //ֻ��һ�����������Σ����п�����һ���ڵ㣬Ҳ�п���һ���ڵ㶼û��
			{
				int tri[3];
				for(int k=0; k<3; k++)
				{
					tri[k] = temp_tri[0][k];
				}
				for(int k=0; k<4; k++)
				{
					if(nod_ellip[nid[k]]==-2)					//-2��ʾδʹ��
					{
						if(k==tri[0]||k==tri[1]||k==tri[2])		//���
						{
							ellgeo[i].outell_nodes.push_back(nid[k]);
							nod_ellip[nid[k]] = int(ellgeo[i].outell_nodes.size())-1;
						}
					}
				}
			}
			else if(num_tri==2||num_tri==3)
			{
				for(int k=0; k<4; k++)							//ȫ���
				{
					if(nod_ellip[nid[k]]==-2)					//-2��ʾδʹ��
					{
						ellgeo[i].outell_nodes.push_back(nid[k]);
						nod_ellip[nid[k]] = int(ellgeo[i].outell_nodes.size())-1;
					}
				}
			}
			else if(num_tri!=0)
			{
				hout << "������������ĸ���������Ƭȫ�����������������Ƭ��" << endl;
				return 0;
			}

			//������������Ƭ��Ϣ���������
			vector<vector<int> > temp_line;
			for(int k=0; k<num_tri; k++)
			{
//---------------------------------------------------------------------------------------------------------------
//<1>���¶�����������ж������ε������߶��Ƿ����߶��غ�,������
				for(int m=0; m<2; m++)
				{
					for(int n=m+1; n<3; n++)
					{
						int d[2];
						if(nid[temp_tri[k][m]]<nid[temp_tri[k][n]])
						{
							d[0] = nid[temp_tri[k][m]];
							d[1] = nid[temp_tri[k][n]];
						}
						else
						{
							d[0] = nid[temp_tri[k][n]];
							d[1] = nid[temp_tri[k][m]];
						}
						int line_size = int(temp_line.size());
						int key = -1;
						for(int p=0; p<line_size; p++)
						{
							if(temp_line[p][0]==d[0]&&temp_line[p][1]==d[1])
							{
								key = p;
								break;
							}
						}
						if(key==-1)
						{
							vector<int> tem(2);
							tem[0]=d[0];
							tem[1]=d[1];
							temp_line.push_back(tem);
						}
						else
						{
							temp_line.erase(temp_line.begin()+key);
						}
					}
				}
//---------------------------------------------------------------------------------------------------------------
				//����������Ƭ��������ĵ�Ԫ���ת��Ϊ����������������б��
				for(int l=0; l<3; l++)
				{
					temp_tri[k][l] = nod_ellip[nid[temp_tri[k][l]]];
				}
			}
//---------------------------------------------------------------------------------------------------------------
//<2>���¶�����������ж������ε������߶��Ƿ����߶��غ�,������
			int tem_size = int(temp_line.size());
			for(int k=0; k<tem_size; k++)
			{
				int key = -1;
				int line_size = int(line_vec.size());
				for(int m=0; m<line_size; m++)
				{
					if(temp_line[k][0]==line_vec[m][0]&&temp_line[k][1]==line_vec[m][1])
					{
						key = m;
						break;
					}
				}
				if(key==-1)
				{
					line_vec.push_back(temp_line[k]);
				}
				else
				{
					line_vec.erase(line_vec.begin()+key);
				}
			}
//---------------------------------------------------------------------------------------------------------------
			//������������Ƭ���������
			ellgeo[i].tri.insert(ellgeo[i].tri.end(),temp_tri.begin(),temp_tri.end());
		}
//---------------------------------------------------------------------------------------------------------------
//<3>���¶�����������ж������ε������߶��Ƿ����߶��غ�,������
		//�ж�����������Ƿ���
		if(!line_vec.empty())
		{
			//�ж��Ƿ��Ǻ͵���������ཻ�����Ĳ��֣��ⲿ�����·��Ȧ
			vector<int> poi_vec;
			int lvsize = int(line_vec.size());
			for(int j=0; j<lvsize; j++)
			{
				for(int k=0; k<2; k++)
				{
					int lvjk = line_vec[j][k];
					if(nodes_vec[lvjk].flag==1)
					{
						int key=-1;
						int pvsize = int(poi_vec.size());
						for(int m=0; m<pvsize; m++)
						{
							if(poi_vec[m]==lvjk)
							{
								key = m;
								break;
							}
						}
						if(key==-1)
						{
							poi_vec.push_back(lvjk);
						}
						else
						{
							poi_vec.erase(poi_vec.begin()+key);
						}
					}
				}
			}
			if(!poi_vec.empty())
			{
				hout << "������������������β���գ���Ҫ���Ǻ͵���������ཻ�������" << endl;
				return 0;
			}
		}

		//ȷ�������ڵ�
		for(int j=0; j<tet_size; j++)
		{
			vector<int > temp_tet_nod;		//��¼���������е��������Ӧ�ĵ�Ԫ�ֲ���ź���������������еı��
			for(int k=0; k<4; k++)
			{
				int nid =  eles_vec[ellgeo[i].tet[j]].nodesId[k];			//��Ԫ�ĸ��ڵ�ı��
				if(nod_ellip[nid]==-2)		//-2��ʾδʹ��,�����е���㶼ȷ����֮�����µľ����ڵ�
				{
					ellgeo[i].inell_nodes.push_back(nid);
					nod_ellip[nid]=-1;
				}
				else if(nod_ellip[nid]>=0)
				{
					temp_tet_nod.push_back(k);							//��¼�ڸ��������еľֲ����
					temp_tet_nod.push_back(nod_ellip[nid]);		//��¼����Ӧ������������еı��
				}
			}
			//���뵥Ԫ���������Ϣ
			ellgeo[i].tet_nod.push_back(temp_tet_nod);
		}
		//---------------------------------------------------------------------------------------------------------------
		//���¿�ʼnod_ellip����
		int insize = int(ellgeo[i].inell_nodes.size());
		for(int j=0; j<insize; j++)
		{
			nod_ellip[ellgeo[i].inell_nodes[j]] = -2;
		}
		int outsize = int(ellgeo[i].outell_nodes.size());
		for(int j=0; j<outsize; j++)
		{
			nod_ellip[ellgeo[i].outell_nodes[j]] = -2;
		}
	}

	//------------------------------------------------------------------------------------------------------
	//���������Ĺ���
	//------------------------------------------------------------------------------------------------------
	//��ԭ�������嵥Ԫ���������ϵ�Ԫ����
	for(int i=0; i<ele_size; i++)
	{
		Element temp_element;
		elements_vec.push_back(temp_element);
		elements_vec.back().mat=eles_vec[i].materialId;
		elements_vec.back().type=2;
		for(int j=0; j<4; j++)
		{
			elements_vec.back().nodes_id.push_back(eles_vec[i].nodesId[j]);
		}
	}
	//�������������񲢼ӽ����
	for(int i=0; i<int(ellgeo.size()); i++)
	{
		//------------------------------------------------------------------------------------
		//�����ڵ�����(�������嵥Ԫ�޹�)
		for(int j=0; j<int(ellgeo[i].inell_nodes.size()); j++)
		{
			int nodes_n = ellgeo[i].inell_nodes[j];
			double ratio = thick_ratio;
			
			Point point;
			point.x = 0.0;
			point.y = 0.0;
			point.z = 0.0;

			Point node_poi; //��¼�����ĳ�ʼ��
			node_poi.x = nodes_vec[nodes_n].x;
			node_poi.y = nodes_vec[nodes_n].y;
			node_poi.z = nodes_vec[nodes_n].z;

			if (thecell->surfaces_vec[i]->change_coor(&node_poi, &point, ratio)==0)
			{
				hout << i << "��������任ʱ����" << endl;
				hout << "*******************************************" << endl;
				return 0;
			}

			nodes_vec[nodes_n].x = point.x;
			nodes_vec[nodes_n].y = point.y;
			nodes_vec[nodes_n].z = point.z;
		}
		//------------------------------------------------------------------------------------
		//������������������ɽ�����������ı�ĳЩ������ڵ���
		vector<int> temp_vec(layer_num+1,0);
		vector<vector<int> > stem_node(ellgeo[i].outell_nodes.size(),temp_vec);

		for(int j=0; j<int(ellgeo[i].outell_nodes.size()); j++)
		{
			const int nodes_n = ellgeo[i].outell_nodes[j];
			stem_node[j][0]=nodes_n;

			Node new_node;

			Point point;
			point.x = 0.0;
			point.y = 0.0;
			point.z = 0.0;

			Point node_poi; //��¼�����ĳ�ʼ��
			node_poi.x = nodes_vec[nodes_n].x;
			node_poi.y = nodes_vec[nodes_n].y;
			node_poi.z = nodes_vec[nodes_n].z;

			for(int k=1; k<=layer_num; k++)
			{
				double ratio=(k*1.0/layer_num)*thick_ratio;	
				
				if (thecell -> surfaces_vec[i] -> change_coor(&node_poi, &point, ratio) == 0)
				{
					hout << i << "��������任ʱ����" << endl;
					hout << "*******************************************" << endl;
					return 0;
				}

				new_node.x = point.x;
				new_node.y = point.y;
				new_node.z = point.z; 
				new_node.flag = nodes_vec[nodes_n].flag;

				//�����µ�
				nodes_vec.push_back(new_node);
				where_is_nodes.push_back(i);
				if(new_node.flag==1) bnodes_vec.push_back(int(nodes_vec.size())-1);
				stem_node[j][k]=int(nodes_vec.size())-1;

				if (j == layer_num)	is_on_nodes.push_back(3);
				else	is_on_nodes.push_back(2);
			}
		}
		//�޸�ĳЩ������ڵ���
		for(int j=0; j<int(ellgeo[i].tet.size()); j++)
		{
			for(int k=0; k<int(ellgeo[i].tet_nod[j].size());k=k+2)
			{
				elements_vec[ellgeo[i].tet[j]].nodes_id[ellgeo[i].tet_nod[j][k]]=stem_node[ellgeo[i].tet_nod[j][k+1]][layer_num];
			}
		}
		//��������
		for(int j=0; j<int(ellgeo[i].tri.size()); j++)
		{
			for(int k=0; k<layer_num; k++)
			{
				Element temp_element;
				elements_vec.push_back(temp_element);
				elements_vec.back().mat=2;
				elements_vec.back().type=3;
				//�ж������������(Ϊ��֤����Ԫ���㣬�����������Ϊ��)
				vector <Node> elenodes_vec(6);
				int nodn[6];
				for(int n=0;n<2; n++)
				{
					for(int m=0; m<3; m++)
					{
						int nume = stem_node[ellgeo[i].tri[j][m]][k+n];
						nodn[n*3+m] = nume;
						elenodes_vec[n*3+m] = nodes_vec[nume];
					}
				}
				if(Threeprism_volume(elenodes_vec)>0)
				{
					for(int m=0; m<6; m++)
					{
						elements_vec.back().nodes_id.push_back(nodn[m]);
					}
				}
				else
				{
					for(int m=0; m<3; m++)
					{
						int nume = nodn[m];
						nodn[m] = nodn[m+3];
						nodn[m+3] = nume;

						Node elentem = elenodes_vec[m];
						elenodes_vec[m] = elenodes_vec[m+3];
						elenodes_vec[m+3] = elentem;
					}
					if(Threeprism_volume(elenodes_vec)>0)
					{
						for(int m=0; m<6; m++)
						{
							elements_vec.back().nodes_id.push_back(nodn[m]);
						}
					}
					else
					{
						cout << "������������J_val���ԡ�" <<endl;
						//for(int m=0; m<6; m++)
						//{
						//	hout << elenodes_vec[m].x << "   ";
						//	hout << elenodes_vec[m].y << "   ";
						//	hout << elenodes_vec[m].z << "   ";
						//}
						//hout << endl;
						//for(int n=0;n<2; n++)
						//{
						//	for(int m=0; m<3; m++)
						//	{
						//		hout << nodes_vec[stem_node[ellgeo[i].tri[j][m]][k+n]].x << "   ";
						//		hout << nodes_vec[stem_node[ellgeo[i].tri[j][m]][k+n]].y << "   ";
						//		hout << nodes_vec[stem_node[ellgeo[i].tri[j][m]][k+n]].z << "   ";
						//	}
						//	hout << endl;
						//}
						exit(0);
					}
				}
			}
		}
	}
	return 1;
}

//---------------------------------------------------------------------------
//�ڲ����ɱ߽�������£��������嵥Ԫ���͵�Ϊ��ϵ�Ԫ����
int Mesher::change_elements_to_blend()
{
	int num = (int)eles_vec.size();
	elements_vec.resize(num);
	for(int i=0; i < num; i++ )
	{
		elements_vec[i].mat = eles_vec[i].materialId;
		for (int j=0; j < 4; j++)
		{
			elements_vec[i].nodes_id.push_back(eles_vec[i].nodesId[j]);
		}
		elements_vec[i].type = 2;
	}

	return 1;
}

//---------------------------------------------------------------------------
//���tecplot��������
int Mesher::export_tecplot_data()
{
	//����������������������񣨰������������
	ofstream otec("mesh_in_tecplot.dat");
	otec << "TITLE = blend_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)elements_vec.size(); i++)
	{
		switch (elements_vec[i].mat)
		{
		case 0: count[0]++; break; //���嵥Ԫ
		case 1: count[1]++; break; //����Ԫ
		case 2: count[2]++; break; //�߽�㵥Ԫ
		default: hout << " error!! " << endl; break;
		}
	}
	otec << "ZONE N=" << (int)nodes_vec.size() << ", E=" << count[1] << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for (int i=0; i < (int)nodes_vec.size(); i++)
	{
		otec << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements_vec.size(); i++)
	{
		if (elements_vec[i].mat == 1)
		{
			otec	<< elements_vec[i].nodes_id[0]+1 << "  " << elements_vec[i].nodes_id[1]+1 << "  " 
					<< elements_vec[i].nodes_id[2]+1 << "  " << elements_vec[i].nodes_id[3]+1 << endl;
		}
	}
	//�����������������񣨰������������
	otec << "ZONE N=" << (int)nodes_vec.size() << ", E=" << count[2] << ", F=FEPOINT, ET=BRICK" << endl;
	for (int i=0; i < (int)nodes_vec.size(); i++)
	{
		otec << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)elements_vec.size(); i++)
	{
		if (elements_vec[i].mat == 2)
		{
//---------------------------------------------------------------------------------
//�ڼ��ʱʹ��
//			int k = elements_vec[i].nodes_id[0];
//			for (int j=1; j <6; j++)
//			{
//				int m = elements_vec[i].nodes_id[j];
//				if (where_is_nodes[k]!=where_is_nodes[m])
//				{
//					hout << "���ֲ��ȵ㣡����" << endl;
//					hout << "i=" << i << endl;
//					for (int n=0; n < 6; n++)
//					{
//						hout << n << "�ŵ�" << elements_vec[i].nodes_id[n] << "�ŵ�"<< "����������" << where_is_nodes[elements_vec[i].nodes_id[n]] << endl;
//					}
//					goto label;
//				}
//			}
//label: ;
//---------------------------------------------------------------------------------
			otec	<< elements_vec[i].nodes_id[0]+1 << "  " << elements_vec[i].nodes_id[1]+1 << "  " 
					<< elements_vec[i].nodes_id[2]+1 << "  " << elements_vec[i].nodes_id[0]+1 << "  "
					<< elements_vec[i].nodes_id[3]+1 << "  " << elements_vec[i].nodes_id[4]+1 << "  "
					<< elements_vec[i].nodes_id[5]+1 << "  " << elements_vec[i].nodes_id[3]+1 << endl;
		}
	}

	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//�����ά�ռ���ӻ�Ensight��������
int Mesher::export_Ensight_data()
{
	//�����ͷ
	//ofstream otec("Composites_Materials_Particles.geo");
	//otec << "geometry meshes" << endl;
	//otec << "composites materials with core/shell particles" << endl;
	//otec << "node id assign" << endl;
	//otec << "element id assign" << endl;
	//otec << "coordinates" << endl;
	////����ڵ���Ϣ
	//otec << setw(10) << (int)nodes_vec.size() << endl;
	//for (int i=0; i < (int)nodes_vec.size(); i++)
	//{
	//	otec <<right;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].x;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].y;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].z;
	//	otec <<endl;
	//}
	////��¼���൥Ԫ�ĸ���
	//int count[3]={0,0,0};
	//for (int i=0; i < (int)elements_vec.size(); i++)
	//{
	//	switch (elements_vec[i].mat)
	//	{
	//	case 0: count[0]++; break; //���嵥Ԫ
	//	case 1: count[1]++; break; //����Ԫ
	//	case 2: count[2]++; break; //�߽�㵥Ԫ
	//	default: hout << " error!! " << endl; break;
	//	}
	//}
	//cout<<"matrix: "<<count[0]<<"  core"<<count[1]<<"  shell"<<count[2];
	////��������˵�Ԫ
	//otec << "part 1" << endl;
	//otec << "cores of particles" << endl;
 //   otec << "tetra4" << endl;
	//otec <<setw(10)<<count[1]<<endl;
	//otec <<right;
	//for (int i=0; i <(int)elements_vec.size(); i++)
	//{
	//	if (elements_vec[i].mat == 1)
	//	{
	//		otec <<setw(8)<<elements_vec[i].nodes_id[0]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[1]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[2]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[3]+1;
	//	    otec <<endl;
	//	}
	//}
	////��������ǲ㵥Ԫ
	//otec << "part 2" << endl;
	//otec << "shells of particles" << endl;
 //   otec << "penta6" << endl;
	//otec <<setw(10)<<count[2]<<endl;
	//otec <<right;
	//for (int i=0; i <(int)elements_vec.size(); i++)
	//{
	//	if (elements_vec[i].mat == 2)
	//	{
	//		otec <<setw(8)<<elements_vec[i].nodes_id[0]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[1]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[2]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[3]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[4]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[5]+1;
	//		otec <<endl;
	//	}
	//}
	////������嵥Ԫ
	//otec << "part 3" << endl;
	//otec << "matrix" << endl;
	//otec << "tetra4" << endl;
	//otec <<setw(10)<<count[0]<<endl;
	//otec <<right;
	//for (int i=0; i <(int)elements_vec.size(); i++)
	//{
	//	if (elements_vec[i].mat == 0)
	//	{
	//		otec <<setw(8)<<elements_vec[i].nodes_id[0]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[1]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[2]+1;
	//		otec <<setw(8)<<elements_vec[i].nodes_id[3]+1;
	//		otec <<endl;
	//	}
	//}
	//otec.close();

	return 1;
}
//-----------------------------------------------------------------------------------------------
//���������
double Mesher::Threeprism_volume(const vector<Node> &elenodes_vec)
{
	//ѭ����˹�������֣�
	double volume=0;
	//����ʾ���
	//--------------------------------------------
	//�κ���N��gauss[count].x,gauss[count].y,gauss[count].z��ƫ������
	double diff[3][6];
	diff[0][0]=0.5*(1.0-0);  
	diff[0][1]=0;                         
	diff[0][2]=-0.5*(1.0-0);
	diff[0][3]=0.5*(1.0+0);
	diff[0][4]=0;
	diff[0][5]=-0.5*(1.0+0);

	diff[1][0]=0;
	diff[1][1]=diff[0][0];
	diff[1][2]=diff[0][2];
	diff[1][3]=0;
	diff[1][4]=diff[0][3];
	diff[1][5]=diff[0][5];

	diff[2][0]=-0.5*0.33333333333;
	diff[2][1]=-0.5*0.33333333333;
	diff[2][2]=-0.5*0.33333333333;
	diff[2][3]=-diff[2][0];
	diff[2][4]=-diff[2][1];
	diff[2][5]=-diff[2][2];
	//--------------------------------------------------
	//��Ԫ�ڵ��������
	double elenode[6][3];
	for(int i=0;i<6;i++){
		elenode[i][0]=elenodes_vec[i].x;
		elenode[i][1]=elenodes_vec[i].y;
		elenode[i][2]=elenodes_vec[i].z;
	}
	//--------------------------------------------------
	//J����
	double Jmatrix[3][3];
	//������������Ļ�
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			Jmatrix[i][j]=0;
			for(int k=0;k<6;k++)
			Jmatrix[i][j]=Jmatrix[i][j] + diff[i][k]*elenode[k][j];
		}
		//--------------------------------------------------
		//���J���������ʽ��
		double J_val=Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
			-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
			+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);	
		//--------------------------------------------------
		//�ɸ�˹���ֵ��������
		  volume=J_val*0.5*2.0;
	
	  return volume;
}
//---------------------------------------------------------------------------
//ȷ�����нڵ����ص�Ԫ��Ϣ
void Mesher::deter_nodes_relative_eles()
{
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		nodes_vec[i].relative_eles_vec.clear();
	}
	for( int i=0; i<(int)elements_vec.size(); i++ )
	{
		int node_size = int(elements_vec[i].nodes_id.size());
		for( int j=0; j<node_size; j++ )
		{
			nodes_vec[elements_vec[i].nodes_id[j]].relative_eles_vec.push_back(i);
		}
	}
}
//---------------------------------------------------------------------------
//����ȷ��materialId ==0 �ǻ��嵥Ԫ��materialId ==1 ������Ԫ
//������ɱ߽��ǰ��tecplot��������
int Mesher::export_tecplot_data_before_coating()
{
	//����������������������񣨰������������
	ofstream otec("mesh_in_tecplot_before_coating.dat");
	otec << "TITLE = tetrahedron_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)eles_vec.size(); i++)
	{
		switch (eles_vec[i].materialId)
		{
		case 0: count[0]++; break; //���嵥Ԫ
		case 1: count[1]++; break; //����Ԫ
		case -1: count[2]++; break;
		default: hout << " error!! " << endl; break;
		}
	}
	//����
	otec << "ZONE N=" << (int)nodes_vec.size() << ", E=" << count[0]+count[2]<< ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for (int i=0; i < (int)nodes_vec.size(); i++)
	{
		otec << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)eles_vec.size(); i++)
	{
		if (eles_vec[i].materialId != 1)
		{
			otec	<< eles_vec[i].nodesId[0]+1 << "  " << eles_vec[i].nodesId[1]+1 << "  " 
					<< eles_vec[i].nodesId[2]+1 << "  " << eles_vec[i].nodesId[3]+1 << endl;
		}
	}
	//����
	otec << "ZONE N=" << (int)nodes_vec.size() << ", E=" << count[1] << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for (int i=0; i < (int)nodes_vec.size(); i++)
	{
		otec << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
	}
	otec << endl;
	for (int i=0; i < (int)eles_vec.size(); i++)
	{
		if (eles_vec[i].materialId == 1)
		{
			otec	<< eles_vec[i].nodesId[0]+1 << "  " << eles_vec[i].nodesId[1]+1 << "  " 
					<< eles_vec[i].nodesId[2]+1 << "  " << eles_vec[i].nodesId[3]+1 << endl;
		}
	}
	otec.close();

	return 1;
}
//---------------------------------------------------------------------------
//��Ϊ���ӵ���ֱ�
void Mesher::Reincrese_vol_raito()
{
	//�����ӵ���ֱ��ļ�
	//������Ȼ�������������޷��ﵽҪ��ʱ�����ʷֺ�ı���嵥Ԫ������Ϊ��������Ϊ�����������
	ifstream in_refile;
	in_refile.open("Reincrease_vol_ratio.dat", ios::in);
	if(in_refile)
	{
		istringstream in(Get_Line(in_refile));
		double re_ratio;
		in >> re_ratio;
		hout << "re_ratio: " << re_ratio << endl;
		if(re_ratio<=0.0||re_ratio>0.1)
		{
			hout << "      ��Ҫ���ӵ�����ȷ�С�ڵ���0%���Ҵ���10%"<<endl;
		}
		else
		{
			//ȷ�����нڵ����ص�Ԫ��Ϣ(ȷ������ʹ����ȷ)
			for( int i=0; i<(int)nodes_vec.size(); i++ )
			{
				nodes_vec[i].relative_eles_vec.clear();
			}
			for( int i=0; i<(int)eles_vec.size(); i++ )
			{
				for( int j=0; j<4; j++ )
				{
					nodes_vec[eles_vec[i].nodesId[j]].relative_eles_vec.push_back(i);
				}
			}
			
			//��ʼ�������
			srand((unsigned int)time(0));

			//���������
			double total_vol = thecell->clength*thecell->cwidth*thecell->cheight;

			//����������ֱȱ���
			double increase_ratio = 0.0;
			while(increase_ratio<re_ratio)
			{
				//��ȡ�����н��ڿ������ϵĻ��嵥Ԫ
				vector<int> temp_ele(2,0);
				vector<vector<int> > matrix_ele(0,temp_ele);
				for( int i=0; i<(int)eles_vec.size(); i++ )
				{
					if( eles_vec[i].materialId == 0 )
					{
						for(int j=0; j<4; j++)
						{
							int enid = eles_vec[i].nodesId[j];
							int revs = (int)nodes_vec[enid].relative_eles_vec.size();
							for(int k=0; k<revs; k++)
							{
								if(eles_vec[nodes_vec[enid].relative_eles_vec[k]].materialId != 0)
								{
									temp_ele[0] = i;
									temp_ele[1] = eles_vec[nodes_vec[enid].relative_eles_vec[k]].materialId;
									matrix_ele.push_back(temp_ele);
									goto eles_vec_mat;
								}
							}
						}
eles_vec_mat: ;
					}
				}
				
				if((int)matrix_ele.size()==0)
				{
					hout << "matrix_ele.size()==0, ����" << endl;
					exit(0);
				}

				//����ı�һ�������������ʵĻ��嵥ԪΪ������Ԫ
				//�ʷֳߴ�global_length�������Ҫ�ı䵥Ԫ����Ϊcount = (int)(re_ratio*total_vol)/(pow(global_length,3)/6.0)
				//ÿ�θı��ܸ�����1/6
				int ele_count = (int)(re_ratio*total_vol/pow(global_length,3));
				for(int i=0; i<ele_count; i++)
				{
					int mat_ele_n = (int)(((double)(rand())/RAND_MAX)*(int)matrix_ele.size());
					if(eles_vec[matrix_ele[mat_ele_n][0]].materialId != 0) continue;
					eles_vec[matrix_ele[mat_ele_n][0]].materialId = matrix_ele[mat_ele_n][1];

					//����ı䵥Ԫ�������������ӵ����
					increase_ratio = increase_ratio + cal_tet_volume( &eles_vec[matrix_ele[mat_ele_n][0]] )/total_vol;
					if(increase_ratio>=re_ratio) break;
				}
			}
			hout << "      ��Ϊ�����������Ϊ��" << increase_ratio <<endl;
		}
	}
	else
	{
		hout << "      Reincrease_vol_ratio.dat�ļ������ڻ�򲻿�������Ƿ���Ҫ���������ӿ�����ֱȵĲ�����"<<endl;
	}
	in_refile.close();
}
//======================================================================
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Mesher::Get_Line(ifstream &infile)const{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
