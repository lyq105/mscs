//===========================================================================
// Mesher.cpp
// 四面体网格剖分类成员函数
// Member Functions in a Class of generating tetrahedron elemments
//===========================================================================
#include "Mesher.h"
#define PI 3.141592654
#define ZERO 1e-8

//---------------------------------------------------------------------------
//构造函数
Mesher::Mesher(CMCell* cmccell)
{
	thecell=cmccell;
	initial();
}
//--------------------------------------------------------------------------
//初始化
void Mesher::initial(  )
{
	attach_length = global_length / 2.0 ;
	for( int i=0; i<(int)thecell->surfaces_vec.size(); i++ )
	{
		vector<int> nodes;
		gnodes_vec.push_back(nodes); //初始化
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
//网格生成
int Mesher::Mesh_generate(ifstream &infile)
{
	//读取网格剖分尺寸
	istringstream in(Get_Line(infile));
	double mesh_size;
	in >> mesh_size;
	//读入有关界面层信息
	istringstream in1(Get_Line(infile));
	int coat, layer_num;
	double thick_ratio;
	in1 >> coat >> thick_ratio >> layer_num;

	if(mesh_tet_map(mesh_size) == 0 )
	{
		hout << "~_~ 剖分失败！可能原因：增强项重叠或者最小体积单元与最大体积单元相差太大!" << endl;
		hout << "****************************************************************************************" <<endl;	
		return 0;
	}
	
	//输出生成边界层前的tecplot网格数据
	//if (export_tecplot_data_before_coating() == 0)
	//{
	//	hout << "~_~ 输出生成边界层前的tecplot网格数据失败！" << endl;
	//	hout << "*************************************************" << endl;
	//	return 0;
	//}

	if (coat==1)
	{
		//经过确认materialId ==0 是基体单元，materialId ==1 是椭球单元
		clock_t ct0,ct1;
		ct0 = clock();
		hout << "-_- 开始生成界面层网格... " << endl;
		if(gen_coating_mesh(thick_ratio, layer_num) == 0)
		{
			hout << "~_~ 生成界面层网格失败！" << endl;
			hout << "*******************************************" <<endl;	
			return 0;
		}

		hout << "    输出tecplot网格数据....." << endl;
		if (export_tecplot_data() == 0)
		{
			hout << "~_~ 输出tecplot网格数据失败！" << endl;
			hout << "*******************************************" << endl;
			return 0;
		}

		//hout<<"   输出三维空间可视化Ensight网格数据"
		//if (export_Ensight_data() == 0)
		//{
		//	hout << "~_~ 输出三维空间可视化Ensight网格数据失败！" << endl;
		//	hout << "*******************************************" << endl;
		//	return 0;
		//}

		ct1 = clock();
		hout << "    生成界面层耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒" << endl;
		hout << "^_^ 界面层网格生成完毕!" << endl<<endl;
	}
	else
	{
		//在不生成边界层的情况下，将四面体单元类型倒为混合单元类型
		clock_t ct0,ct1;
		ct0 = clock();
		hout << "-_- 将四面体网格转换...... " << endl;
		if (change_elements_to_blend() == 0)
		{
			hout << " 将四面体单元类型倒为混合单元类型操作失败！" << endl;
			hout << "*******************************************" << endl;
			return 0;
		}
		ct1 = clock();
		hout << "    转换四面体网格耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒" << endl;
		hout << "^_^ 转换四面体网格完毕!" << endl << endl;
	}
	
	//确定所有节点的相关单元信息
	deter_nodes_relative_eles();

	return 1;
}
//---------------------------------------------------------------------------
//输出或读取网格数据
int Mesher::Mesh_data(int mod)
{
	if(mod==0)			//输出数据
	{
		//---------------------------------------------------------------------
		//输出网格数据
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
		//输出节点数据
		ofstream outn("Data_Nodes.dat");
		outn << int(nodes_vec.size()) << endl;
		for(int i=0; i<int(nodes_vec.size()); i++)
		{
			outn << nodes_vec[i].flag << "  " << setprecision(12)
				    << nodes_vec[i].x << "  " << nodes_vec[i].y << "  " << nodes_vec[i].z << endl;
		}
		outn.close();
		//---------------------------------------------------------------------
		//输出边界节点数据
		ofstream outb("Data_Bnodes.dat");
		outb << int(bnodes_vec.size()) << endl;
		for(int i=0; i<int(bnodes_vec.size()); i++)
		{
			outb << bnodes_vec[i] << endl;
		}
		outb.close();
	}
	else if(mod==1)	//读取数据
	{
		//---------------------------------------------------------------------
		//读取网格数据
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
		//读取节点数据
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
		//读取边界节点数据
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
//输出或读取二进制网格数据
int Mesher::Mesh_BinaryData(int mod, string data_file, int CNum)
{
	int num1 = CNum/10;
	int num2 = CNum - num1*10;
	char ch[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	if(mod==0)			//输出数据
	{
		//---------------------------------------------------------------------
		//输出网格数据
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Mesh";
		string MeshName = title+str+ch[num1]+ch[num2]+".dat";
		ofstream out(MeshName.c_str(),ios::binary);
		if(!out) { hout << "不能打开网格数据文件" << MeshName << "!" << endl;  exit(0); }
		//---------------------------------------------------------------------
		//单元数据
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
		//节点数据
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
		//边节点数据
		int numb = int(bnodes_vec.size());
		out.write((char *)&numb, sizeof(int));
		for(int i=0; i<numb; i++)	out.write((char *)&bnodes_vec[i],sizeof(int));
		//---------------------------------------------------------------------
		out.close();
	}
	else if(mod==1)	//读取数据
	{
		//---------------------------------------------------------------------
		//输出网格数据
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Mesh";
		string MeshName = title+str+ch[num1]+ch[num2]+".dat";
		ifstream in(MeshName.c_str(),ios::binary);
		if(!in) { hout << "不能打开网格数据文件" << MeshName << "!" << endl;  exit(0); }
		//---------------------------------------------------------------------
		//单元数据
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
		//节点数据
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
		//边节点数据
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
//映射法生成四面体网格（主要程序段）
int Mesher::mesh_tet_map( double glength )
{
	clock_t ct_mesh_begin = clock();

	bool debuging_map_mesh = false;

	//全局使用剖分尺度
	global_length = glength;
	//单元最小体积
	min_volume = 0.0 ;
                                
	//单胞的最大最小坐标
	x_min=thecell->origin_x ;
	x_max=thecell->origin_x + thecell->clength ;
	y_min=thecell->origin_y ;
	y_max=thecell->origin_y + thecell->cwidth ;
	z_min=thecell->origin_z ;
	z_max=thecell->origin_z + thecell->cheight ;

	nodes_vec.clear();

	hout<< "-_- 开始四面体网格剖分" << endl<<endl;
	//---------------------------------------------------------------------------
	hout<< "-_- 开始生成背景网格（六面体）......" << endl;
	clock_t ct_background_begin = clock();

	//生成背景网格
	//计算每条边上分的段数，ceil(num)是数学库里的函数，返回不小于num的整数(按浮点值表示)
	int x_sec_num = (int)ceil((x_max-x_min)/glength) ;
	int y_sec_num = (int)ceil((y_max-y_min)/glength) ;
	int z_sec_num = (int)ceil((z_max-z_min)/glength) ;
         
	vector< Hexahedron > hexes_vec;
	//生成背景网格（六面体）
	gen_bg_mesh(x_sec_num,y_sec_num,z_sec_num, hexes_vec);

	where_is_nodes.assign(nodes_vec.size(), -2);	//存放每个节点属于基体还是哪个椭球
	is_on_nodes.assign(nodes_vec.size(), 0);			//存放节点是否在椭球和基体的交界面上
	where_is_nodes.reserve(nodes_vec.size()*2);
	is_on_nodes.reserve(nodes_vec.size()*2);

	//定义四边形面片
	int num_rec_xfaces=(x_sec_num+1)*y_sec_num*z_sec_num;
	int num_rec_yfaces=(y_sec_num+1)*x_sec_num*z_sec_num;
	int num_rec_zfaces=(z_sec_num+1)*x_sec_num*y_sec_num;
	int num_rec_faces=num_rec_xfaces+num_rec_yfaces+num_rec_zfaces;
	vector<int> sign_rec_faces(num_rec_faces,-1);

	clock_t ct_background_end = clock();
	hout << "    生成（六面体）背景网格耗时：" << (double)( ct_background_end - ct_background_begin )/CLOCKS_PER_SEC << " 秒。" << endl;
	hout << "^_^ 生成背景网格（六面体）完毕。" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- 开始调整生成的六面体网格，拟合界面......" << endl;
	clock_t ct_adj_begin = clock();

	//调整生成的六面体，使所有的六面体都能拆分成四面体，而没有四面体跨边界
	if(adjust_hexes(hexes_vec,sign_rec_faces) == 0)
	{
		hout << " 调整生成的六面体操作(adjust_hexes)失败！" << endl;
		hout << "*********************************************" << endl;
		return 0;
	}

	clock_t ct_adj_end = clock();
	hout << "    调整背景六面体网格耗时：" << (double)( ct_adj_end - ct_adj_begin )/CLOCKS_PER_SEC << " 秒。" << endl;
	hout << "^_^ 调整生成的六面体网格完毕。" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- 开始拆分六面体为四面体，拟合边界......" << endl;  
	clock_t ct_split_begin = clock();

	//拆分六面体为四面体
	if(split_hexes(hexes_vec,sign_rec_faces) == 0) 
	{
		hout << "拆分六面体为四面体操作失败！" << endl;
		hout << "***********************************" << endl;
		return 0;
	}

	hout << "    生成节点nodes: " << (int)nodes_vec.size() << "    生成单元elements: " << (int)eles_vec.size() << endl;
        
	clock_t ct_split_end = clock();
	hout << "    拆分六面体为四面体耗时：" << (double)( ct_split_end - ct_split_begin )/CLOCKS_PER_SEC << " 秒。" << endl;
	hout << "^_^ 拆分六面体为四面体完毕。" << endl<<endl;

	//---------------------------------------------------------------------------
	hout << "-_- 确定每个四面体的材料性质......" << endl;
 	clock_t ct_mat_begin = clock();

	//确定每个四面体的材料
	deter_eles_mat();

 	clock_t ct_mat_end = clock();;
	hout << "    确定每个四面体的材料性质耗时："<< (double)( ct_mat_end - ct_mat_begin )/CLOCKS_PER_SEC << " 秒。" << endl;
	hout << "^_^ 确定每个四面体的材料性质完毕。" << endl<<endl;

	//---------------------------------------------------------------------------
	hout<< "-_- 调整跨界面单元，使所有单元不会刺穿边界，并消除薄元。" << endl ;
	clock_t ct_bound_begin = clock();
	double alength = glength*0.5 ;

	//检查每个单元的边是否贯穿基体和增强材料
	//找出跨边界的单元，并优先处理单胞边界上跨边界的单元
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

	//输出数据
	if( debuging_map_mesh ) hout << " eles_across_boundary_be.size(): " << (int)eles_across_boundary_be.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_bf.size(): " << (int)eles_across_boundary_bf.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_nbf.size(): " << (int)eles_across_boundary_nbf.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_dmb.size(): " << (int)eles_across_boundary_dmb.size() << endl;
	if( debuging_map_mesh ) hout << " eles_across_boundary_ndmb.size(): " << (int)eles_across_boundary_ndmb.size() << endl;
        
  
	hout << "    跨边界单元数：" <<(int)(	eles_across_boundary_be.size()+
																	eles_across_boundary_bf.size()+
																	eles_across_boundary_dmb.size()+
																	eles_across_boundary_ndmb.size() ) << endl;

	//记录插入点的信息（每个元素向量含三个int型数据，
	//前两个为被插入节点的线段端点节点号，第三个插入的节点号）
	vector<vector<int> > inserted_point;

	eles_across_boundary = &eles_across_boundary_be ;

	vector<Bar> new_edges;


	//检查有无四个节点都在界面上的单元，有则拆之
	if( split_grovelling_eles() == 1 )
	{
		//对四种边界单元进行循环
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
	//		hout << "警告-1出现！" << endl;
	//	}
	//}

	//排序后删除一些废弃单元
	quick_sort(&deleting_ele);
	for( int i=(int)deleting_ele.size()-1; i>=0; i-- )
	{
		eles_vec.erase(eles_vec.begin()+deleting_ele[i]);
		if( debuging_map_mesh ) hout << "erase element: " << deleting_ele[i] << endl;
	}

	//检查新生成的边中有无很小的边，有则清除（合并节点）
//  double min_dis=global_length/4.0;
//	check_tiny_edges(new_edges,min_dis); （此功能尚未完成）

	clock_t ct_bound_end = clock();
	hout<< "    调整跨界面单元耗时：" << (double)( ct_bound_end - ct_bound_begin )/CLOCKS_PER_SEC << " 秒。" << endl;
	hout<< "^_^ 调整跨界面单元成功完成." << endl;

	int min_max_en[2];
	double min_max_vl[2];
	cal_min_max_tet_volume(min_max_en,min_max_vl);
	hout << "      最小体积：" << min_max_vl[0] << " 单元：" << min_max_en[0] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[0]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[0]].materialId << endl;

	hout << "      最大体积：" << min_max_vl[1] << " 单元：" << min_max_en[1] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[1]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[1]].materialId << endl;
	hout << endl;
	
	//输出四面体体积分布情况
//	tet_vol_distrib(min_max_vl,"_a");
	//输出四面体质量分布情况
//	tet_quality("_a");
	//---------------------------------------------------------------------------
	hout<< "-_- 开始进行光顺处理......" << endl;
	clock_t ct_smooth_begin = clock();

	//光顺处理
	//为每个节点添加相邻节点（相关单元）
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
	//光顺是个多次的过程
	while( smoothing_3d(1) && smoothing_3d(0) && circle_num < 15 )
	{
		circle_num ++ ;
		hout << "    光顺次数 NO: " << circle_num << "次" <<endl;
	}

	//检查一下有没有单元材料属性没确定，只是用来检测程序是否运行正常，因为后面改了很多，不敢保证了
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

			//强制这些单元材料为颗粒材料
			eles_vec[i].materialId = 1;
			hout << "强制Element" << i << "materialid为1." << endl;
		}
	}

	clock_t ct_smooth_end = clock();
	hout<< "    光顺处理耗时：" << (double)( ct_smooth_end - ct_smooth_begin )/CLOCKS_PER_SEC << " 秒" << endl;
	hout<< "    单元："<< (int)eles_vec.size() <<" 节点：" << (int)nodes_vec.size() << endl;
	hout<< "^_^ 光顺处理成功完毕." << endl;


	cal_min_max_tet_volume(min_max_en,min_max_vl);
	hout << "      最小体积：" << min_max_vl[0] << " 单元：" << min_max_en[0] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[0]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[0]].materialId << endl;

	hout << "      最大体积：" << min_max_vl[1] << " 单元：" << min_max_en[1] << endl;
	//hout << "        " ;
	//for( int i=0; i<4; i++ )
	//{
	//	hout << eles_vec[min_max_en[1]].nodesId[i] << " " ;
	//}
	//hout << eles_vec[min_max_en[1]].materialId << endl;

	double  min_max_vl_radio = min_max_vl[0]/min_max_vl[1];
	hout << "      最小体积与最大体积之比为：" << min_max_vl_radio << endl;
	if(  min_max_vl_radio < 0.01 )
	{
		hout << "      请注意! 剖分单元的最小体积与最大体积相差太大！！" << endl;
//		return 0;
	}

	//输出四面体体积分布情况
//	tet_vol_distrib(min_max_vl,"_b");
	//输出四面体质量分布情况
//	tet_quality("_b");

	//人为增加的体分比
//	Reincrese_vol_raito();
	
	double vol_ratio = volume_ratio();
	hout << "      椭球体积与基体体积之比：" << vol_ratio << endl <<endl;

	 //确定边界节点和单元
	deter_boundary_nodes();

	clock_t ct_mesh_end = clock();
	hout << "      共剖分节点：" << (int)nodes_vec.size() << "  单元：" << (int)eles_vec.size() << endl;
	hout << "      整个剖分过程共耗时：" << (double)( ct_mesh_end - ct_mesh_begin )/CLOCKS_PER_SEC << " 秒" << endl;
	hout << "^_^ 四面体网格剖分成功!" << endl<<endl; 

	return 1;
}

//---------------------------------------------------------------------------
//生成背景网格
void Mesher::gen_bg_mesh(int x_sec_num, int y_sec_num, int z_sec_num, vector<Hexahedron> &hexes_v)
{
	double dx = (x_max-x_min)/x_sec_num;
	double dy = (y_max-y_min)/y_sec_num;
	double dz = (z_max-z_min)/z_sec_num;
	x_sec_num ++ ;  //x方向等分点数
	y_sec_num ++ ;  //y方向等分点数
	z_sec_num ++ ;  //z方向等分点数

	if( dx < global_length ) global_length = dx;
	if( dy < global_length ) global_length = dy;
	if( dz < global_length ) global_length = dz;

	double hex_volume = dx*dy*dz;
        
	//生成节点
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
				nd.flag = deter_node_flag(i,j,k,z_sec_num-1,y_sec_num-1,x_sec_num-1);	//标注节点的位置
				if( nd.flag >= 0 ) bnodes_vec.push_back( (int)nodes_vec.size() );				//标注边界节点所在节点向量中的编号
				nodes_vec.push_back(nd);
			}
		}
	}

	//六面体数目
	int num_eles=(x_sec_num-1)*(y_sec_num-1)*(z_sec_num-1);
	//矩形面片的数目
	int num_rec_xfaces=x_sec_num*(y_sec_num-1)*(z_sec_num-1);
	int num_rec_yfaces=y_sec_num*(x_sec_num-1)*(z_sec_num-1);
	int num_rec_zfaces=z_sec_num*(x_sec_num-1)*(y_sec_num-1);
	int num_rec_faces=num_rec_xfaces+num_rec_yfaces+num_rec_zfaces;

	//生成六面体单元(8个节点号，6个矩形面片号）
	int now_i=0;
	for( int i=0; i<z_sec_num-1; i++ )
	{
		for( int j=0; j<y_sec_num-1; j++ )
		{
			for( int k=0; k<x_sec_num-1; k++ )
			{
				//立方体的八个顶点
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
//根据给定的i j k以及最大i_max j_max k_max，决定给定节点的位置（角点、边界线、边界面、内部）
//0-5: 边界面上（返回值）
//6-17: 边界线
//18: 角点（所有角点都不能移动，归为一类）
//产生背景网格时使用
//i--z坐标
//j--y坐标
//k--x坐标
int Mesher::deter_node_flag(int i, int j, int k, int i_max, int j_max, int k_max)
{
	int flag;
	vector<int> faces_num;   //属于哪个面
	if( i == 0 )          faces_num.push_back(4);
	else if( i == i_max ) faces_num.push_back(5);

	if( j == 0 )          faces_num.push_back(2);
	else if( j == j_max ) faces_num.push_back(3);

	if( k == 0 )          faces_num.push_back(0);
	else if( k == k_max ) faces_num.push_back(1);

	int bfn = (int)faces_num.size();		//记录在几个面内
	if( bfn == 0 )									//内部节点
	{
		flag = -1;	
	}
	else if( bfn == 1 )							//边界面上的节点	
	{
		flag = faces_num[0];	
	}
	else if( bfn == 2 )							//在边界线上
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
		//角点，不再区分是哪个角点了，所有角点也不能移动
		flag = 18;
	}
	return flag;
}

//---------------------------------------------------------------------------
//韩非20060908增改
//根据给定的节点的编号和单胞的x,y,z方向的最大最小值，
// 决定给定节点的位置（角点、边界线、边界面、内部）
//1-3: 边界面上 （6个面归为3类）
//4-6: 边界线  (12条边界线归为3类)
//7: 角点（所有角点都不能移动，归为1类）

int Mesher::deter_node_flag(int node_num)
{
	int flag;
	vector<int> faces_num;   //属于哪种面
	Node *node = &nodes_vec[node_num];

	if (fabs(node->x - x_min) <= ZERO || fabs(node->x - x_max) <= ZERO)
	{
		faces_num.push_back(1);		//在x最小值或最大值的面内，x不能动，y,z可以动
	}

	if (fabs(node->y - y_min) <= ZERO|| fabs(node->y - y_max) <= ZERO)
	{
		faces_num.push_back(2);		//在y最小值或最大值的面内，y不能动，x,z可以动
	}

	if( fabs(node->z - z_min) <= ZERO|| fabs(node->z - z_max) <= ZERO)
	{
		faces_num.push_back(3);		//在z最小值或最大值的面内，z不能动，x,y可以动
	}

	int bfn = (int)faces_num.size();		//记录在几个面内
	if( bfn == 0 )									//内部节点
	{
		flag = 0;	
	}
	else if( bfn == 1 )							//边界面上的节点	
	{
		flag = faces_num[0];	
	}
	else if( bfn == 2 )							//在边界线上
	{            
		int min_bfn = min(faces_num[0],faces_num[1]);
		int max_bfn = max(faces_num[0],faces_num[1]);
		if( min_bfn == 1 )
		{
			if( max_bfn == 2 ) flag = 4;	//在x和y面的交线上，只有z可动
			if( max_bfn == 3 ) flag = 5;   //在x和z面的交线上，只有y可动
		}
		else if( min_bfn == 2 )
		{
			if( max_bfn == 3 ) flag =6;	//在y和z面的交线上，只有x可动
		}
	}
	else if( bfn == 3 )
	{
		//角点，不再区分是哪个角点了，所有角点也不能移动
		flag = 7;
	}

	return flag;
}

//---------------------------------------------------------------------------
//调整生成的六面体，使所有的六面体都能拆分成四面体，而没有四面体跨边界
int Mesher::adjust_hexes(vector<Hexahedron> &hexes, vector<int> &sign_faces)
{
	bool debuging = false;

//   hout << "global_length: " << global_length << endl;

	//最小距离阀值，小于这个距离将节点移动到界面，大于这个节点的距离最好不要移动
	//三角形面片顶点到对边距离小于这个距离表示面片不合法，尽量避免         
	double min_dis = global_length*3.0/8.0;

	//确定每个节点的位置（在基体内还是增强体内，还是界面上）
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		int is_on;
		where_is_nodes[i] = where_is( &nodes_vec[i],  is_on, i );
		if( where_is_nodes[i] == -3 ) 
		{
			hout << " 确定每个节点的位置操作(where_is)失败! " << endl;
			hout << "*******************************************" << endl;
			return 0;
		}
		is_on_nodes[i] = is_on;
	}

	//标志节点是否可以移动
	//（已经移动到椭球面上的不可以再移动）
	//（如果此节点所有相关六面体中，已经至少有一个六面体有四个节点在界面上，此节点也不能移动)
	//	0: 表示不可以移动     1: 表示可以移动
	vector<int> movable_ns( nodes_vec.size(), 1 );
	int pv=0;
	int pl=0;
	//对每个六面体进行考虑
	for( int i=0; i<(int)hexes.size(); i++ )
	{
		if( debuging )	//可以输出单元及节点编号和节点所处位置情况
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

		//保存8个节点是否在界面上的信息，并随时检查
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
			if( win != -1 && is_on_nodes[nn] == 0 )			//发现在颗粒内的节点
			{   
				if( debuging )
				{
					hout << "       node: " << j << endl;
				}

				int nei_i[3];    //三个相邻的节点编号（局部编号）
				deter_nei(j,nei_i);
				int nei_num[3],nei_win[3],nei_is_on[3];
				for( int k=0; k<3; k++ )		//通过六面体局部编号找到相关节点的全局编号并确定位置
				{
					nei_num[k] = hexes[i].nodesId[nei_i[k]];
					nei_win[k] = where_is_nodes[nei_num[k]];
					nei_is_on[k] = is_on_nodes[nei_num[k]];
				}
				//如果相邻的三个节点都和本节点在同一颗粒内则跳过，否则移动到界面上（保证三个相邻节点在同一个颗粒内或颗粒表面上））
				if( nei_win[0] == win && nei_win[1] == win && nei_win[2] == win ) continue;

				//本节点以及邻节点到椭球面的投影点（注意单胞侧面上的节点投影也要在侧面上）
				Node zerop;
				vector<Node> pnodes(4,zerop);
				int error[4]={-1,-1,-1,-1};     //标志是否成功投影   0：失败    1：成功
				pnodes[0] = project2_elli_surf(nodes_vec[nn],thecell->surfaces_vec[win],&error[0]);
				pv=pv+1;
                //调整与本节点不在同一颗粒内的节点到该颗粒的界面上（调整邻节点，或者调整本节点）
				//首先确定是调整邻节点还是调整本节点，如果调整邻节点，则与本节点不在同一颗粒内的邻节点都要调整到界面上，如果调整本节点，则其他节点可以不再考虑，继续循环
				//根据到椭球面的距离来判断
				double dis0;

				if( error[0] == 0 )
				{
					//投影失败，赋予dis0一个足够大的数，防止判断邻节点时再选择此节点
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

				int move_nei=1;			//是否能够移动相邻节点的标志，等于0不可以移动，等于1可以移动

				//计算邻节点到椭球面的投影距离(nei_dis)
				for( int k=0; k<3; k++ )
				{     
					//如果该邻节点与本节点在同一颗粒中则跳过
					if( nei_win[k] == win ) continue;    
					if( nei_win[k] != -1 && nei_is_on[k] == 0 )
					{
						//在其他颗粒内，首先强制移动到颗粒表面(离哪个表面近就移动到哪个上，注意有可能颗粒重叠）
						Point pp0 = {nodes_vec[nn].x,nodes_vec[nn].y,nodes_vec[nn].z};
						Point pp1 = {nodes_vec[nei_num[k]].x,nodes_vec[nei_num[k]].y,nodes_vec[nei_num[k]].z};
						Point ip0,ip1;
						//求pp0和pp1两点连线或延长线与第win（或nei_win[k]）个椭球面的交点，ipo返回交点，
						//perr1==1表示求交点成功，perr1==0表示求交点失败
						int perr1 = thecell->surfaces_vec[win]->intersect(pp0,pp1,ip0);
						int perr2 = thecell->surfaces_vec[nei_win[k]]->intersect(pp0,pp1,ip1);

						//假定总是能找到交点
						//求第k个相邻点到交点的距离
						double dis1 = nodes_vec[nei_num[k]].distance_to(&ip0);
						double dis2 = nodes_vec[nei_num[k]].distance_to(&ip1);

						if( debuging )
						{
							hout << "nn: " << nn << " nn1: " << nei_num[k] ;
							hout << " dis1: " << dis1 << " dis2: " << dis2 << endl;
						}

						if( dis1 > dis2 )			//把这个邻近节点直接从它所在的椭球粗暴的拉到本节点所在椭球
						{
							nodes_vec[nei_num[k]].move_to(&ip1);
							is_on_nodes[nei_num[k]] = 1;
							movable_ns[nei_num[k]] = 0;      //不能再移动了
							nei_is_on[k] = 1;
						}
						else
						{
							hout << " dis2 > dis1，可能有重叠的椭球，暂不处理 " << endl;
							hout << "*******************************************" << endl;
							return 0;
						}
					}
                   //不仅要拉到交点处还要最终变换到椭球面上法向投影点
					pnodes[k+1] = project2_elli_surf(nodes_vec[nei_num[k]],thecell->surfaces_vec[win],&error[k+1]);
                                             
					if( error[k+1] == 0 )
					{
						//投影失败，赋予dis0一个足够大的数，防止判断邻节点时再选择此节点
						nei_dis[k] = 2.0*global_length;
					}
					else
					{
						nei_dis[k] = nodes_vec[nei_num[k]].distance_to(&pnodes[k+1]);
					}

					if( debuging ) hout << "       " << nei_num[k] << "        " << nei_dis[k] ;

				}

				if( debuging ) hout << endl;

				//检查是否应该移动邻节点到椭球面上
				for( int k=0; k<3; k++ )
				{
					//如果该邻节点与本节点在同一颗粒中仍旧跳过
					if( nei_win[k] == win ) continue;

					 if( nei_dis[k] > dis0 || (movable_ns[nei_num[k]] == 0&&nei_dis[k]>min_dis) )
					{
						if(debuging) hout << "movable_ns["<<nei_num[k]<<"]: " << movable_ns[nei_num[k]]<< " nei_dis["<<k<<"]: " << nei_dis[k] << endl;
						move_nei = 0;		//只要有一个不能移动，其他的外节点能移动也不行
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
					//检查本节点是否能移动（是否已经移动过）
					if( error[0] == 0 ) continue;				//投影失败，这种情况造成了仍有一些边线跨越椭球边界面
					if( movable_ns[nn] == 0 ) continue;	//移动过，不能再移动

					if( dis0 > global_length/2.0 ) continue;	//移动距离太远，不移动

					if( debuging )
					{
						hout << "Moved Node " << nn << "(1) to the surface." << endl;
						hout << "    Before: " << "      " << nodes_vec[nn].x ;
						hout <<  "      " << nodes_vec[nn].y ;
						hout <<  "      " << nodes_vec[nn].z ;
						hout << endl;
					}
                                        
					//移动本节点后，本节点在椭球边界面上，所以应该加在is_on_ln向量中		
					is_on_ln.push_back(j);

					//移动本节点
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
						//如果该邻节点与本节点在同一颗粒中仍旧跳过
						if( nei_win[k] == win ) continue;
						if( error[k+1] == 0 ) continue;   //投影失败，这种情况造成了仍有一些边线跨越椭球边界面
						//检查节点是否能移动（是否已经移动过）
						if( movable_ns[nei_num[k]] == 0 )
						{
							int is_skip=1;
							//已经在别的颗粒表面，但又非常近，于是乎，强制认为已经在边界上了（把两个颗粒连起来）
							if( debuging )
							{
								hout << "nei_dis[" << k <<"]: " << nei_dis[k] ;
								hout << " min_dis: " << min_dis << endl;
							}
							if( nei_dis[k] == 0 )
							{
								//已经在椭球上了
								//（发生这种情况是因为：最开始的时候把该节点移动到该椭球面上，
								//后来又被强制认为在别的椭球上，现在又需要它在该椭球上，唉，麻烦）
								continue;
							}
							else if( nei_dis[k] < min_dis )
							{
								//距离椭球面非常近，如果是在别的椭球上，则强制认为也在当前椭球上（粘在一起了）
								//如果是在基体上，不管它，不移动，防止生成四个节点都在椭球面上的单元，不用担心，光滑处理会改善其状况)
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

		//检查此六面体的节点是否还能移动（如果此六面体已经有4个节点在界面上，所有节点不能再移动）
		if( is_on_ln.size() >= 4 )
		{
			for( int j=0; j<8; j++ )
			{
				movable_ns[hexes[i].nodesId[j]] = 0;
			}
		}

		//确定当前六面体的拆分路径
		//首先确定哪些节点在颗粒内
		vector<int> wb;						//跨越的颗粒编号
		vector<vector<int> > wnodes;	//属于每个颗粒的节点数
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

				//检查所有可能跟此点构成四面体的节点，
				//如果在另外颗粒的表面，则检查到该颗粒的距离，
				//如果很小，则强制认为是此颗粒界面上的点，
				//主要还是为了防止两个颗粒很近时出现很小的单元
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

						if( dis1 < ZERO ) continue; //已经在界面上了
						if( dis1 < min_dis )
						{
							//强制认为此节点在win颗粒上
							where_is_nodes[nn1]=win;
						}
					}
				}
			}
		}

		//当前六面体的拆分样式（六个面的斜线方向）
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
				//只有一个节点在同一颗粒内
				int ann = wnodes[j][0];
				deter_hex_split_style(ann,-1,face_style);
			}
			else if( wnodes[j].size() == 2 )
			{
				//只有两个节点在同一颗粒内
				int edge_num = is_edge_of_hex(wnodes[j][0],wnodes[j][1]);
				deter_hex_split_style(-1,edge_num,face_style);
			}
			else if( wnodes[j].size() == 3 || wnodes[j].size() == 4 )
			{
				//只有一个节点在颗粒外
				//首先找出此节点(有可能没有此节点，即有四个节点在界面上）
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

//	保证单胞侧面上的点还要在侧面上，边线上的点还要在边线上，角点保持不动
	recover_bnodes(bnodes_vec);

	hout << "    总投影：" << pv <<"次！" ;
	hout << "    失败：" << pl << "次！" << endl;
/*
	string file="adjusted_hexes.dat";
	string oldfile=hout.output_file;
	ofstream cout_hex( file.c_str() ) ;
	if ( !cout_hex )  hout << "无法打开文件：" << file << endl;
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
//确定给定节点的位置（基体、第n个椭球）
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
					hout << node_count << "号节点同时在" << already_contained << "和" << rvalue << "两个椭球体内或界面上（于where is函数is_contain < 0中） " << endl;
				}
				else
				{
					hout << "发现有节点同时在" << already_contained << "和" << rvalue << "两个椭球体内或界面上（于where is函数is_contain < 0中） " << endl;
				}
				rvalue = -3;
			}
//			break;		//直接跳出
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
					hout << node_count << "号节点同时在" << already_contained << "和" << rvalue << "两个椭球体内或界面上（于where is函数is_contain == 0中） " << endl;
				}
				else
				{
					hout << "发现有节点同时在" << already_contained << "和" << rvalue << "两个椭球体内或界面上（于where is函数is_contain == 0中） " << endl;
				}
				rvalue = -3;
			}
//			break;		//直接跳出
		}
	}
	return rvalue;
}

//---------------------------------------------------------------------------
//计算给定节点到椭球面(sur)的投影点（注意单胞侧面上的节点投影也要在侧面上）
//大概方法：找到该点与椭球中心点连线与椭球面的交点（project函数），将交点作为起始点，
//						应用梯度方法迭代求解该点在椭球面上的法向投影点，如果是侧面、边线上点
//						为了避免扎堆，可能会用法向投影点变换到侧面、边线上的点与本点连线与椭球
//						交点代替法向投影点
//						单胞侧面上的点还要在侧面上，边线上的点还要在边线上，角点保持不动的特性
//						在adjust_hexes程序最后的recover_bnodes(bnodes_vec)子函数实现
//						只有project_normal函数求解时发散的情况*error==0
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

        if( node.flag < 0 ) return Node(rp);   //内部节点

        Point rp1=rp;
        int lerror1=0;
        if( node.flag < 6 )
		{
                //边界面节点
                if( node.flag == 0 || node.flag == 1 )
				{
                        //x侧面上的点
                        rp1.x = node.x;
                }
                else if( node.flag == 2 || node.flag == 3 )
				{
                        //y侧面上的点
                        rp1.y = node.y;
                }
                else if( node.flag == 4 || node.flag == 5 )
				{
                        //z侧面上的点
                        rp1.z = node.z;
                }
                TDVector pvec(np,rp1);
                rp1=sur->project(&np,&pvec,&lerror1);
        }
        else if( node.flag < 18 && node.flag >= 6 )
		{
                //边界线节点
                if( node.flag == 6 || node.flag == 8 || node.flag == 10 || node.flag == 12 )
				{
                        //x轴方向边界线上的点
                        rp1.x = node.x + 1;
                        rp1.y = node.y;
                        rp1.z = node.z;
                }
                else if( node.flag == 7 || node.flag == 9 || node.flag == 11 || node.flag == 13 )
				{
                        //y轴方向边界线上的点
                        rp1.x = node.x;
                        rp1.y = node.y + 1;
                        rp1.z = node.z;
                }
                else if( node.flag == 14 || node.flag == 15 || node.flag == 16 || node.flag == 17 )
				{
                        //z轴方向边界线上的点
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
                //在侧面内（或边界线上）投影失败，返回法向投影
                return Node(rp);
        }
        else
		{
                //在侧面内（或边界线上）投影成功，检查距离，
                //如果侧面投影点距离相对法向投影点很远，则还是返回法向投影
                //主要是防止侧面上节点挤在一起
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
//根据给定的节点编号（六面体的局部编号），确定相邻的3个节点号（仍然是局部编号）
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
//根据给定的节点号（或者边号，只能指定一个），确定六面体的剖分样式（给定节点（边）在颗粒内的边)
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
//给定两个节点，判断是否为六面体的一条边
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
//恢复移动过的边界节点（单胞侧面上的点）
int Mesher::recover_bnodes(vector<int> &bnodes)
{
        bool debuging = false;

        //恢复移动过的边界节点
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
                                //搜索八个角点，移动到距离最近的一个
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
//拆分六面体成四面体
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

		//修改四边形面片的斜边方向
		//先检查是否会出现点到斜线的距离过近的情况，并想办法消除
		//尽量保证不会出现点到斜线的距离过近的情况     
		for( int j=0; j<6; j++ )
		{
			int nnseq[4];
			if( face_style[j] != -1 )
			{
				//斜边方向已经确定，仍然检查
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
				//斜边方向尚未确定，确保不会出现糟糕情况
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
			hout << "拆分六面体(" << i << ")成四面体失败，拆分成12个四面体!!" << endl;
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
		//更新矩形面片斜边标志
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
//根据给定的侧面编号（六面体的侧面局部编号），以及斜边方向，确定一个节点编号数组
//数组内容：斜边节点1、斜边节点2、非斜边节点1、非斜边节点2
//也就是侧面四边形如何分三角形（保证nnseq[4]中0,1,2号节点和0,1,3号节点各组成一个三角形）
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
//计算一点到另外两点连线的距离
double Mesher::dis2_line( Node *n1, Node *n2, Node *n3 )
{
	//由点p1以及p2指向p3的向量构成平面(p2指向p3的向量是平面的法向量)
	//求该平面与p2到p3的连线的交点，然后求p1与此交点的距离
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
//根据已有的参数，尽量提供一个能拆分成6个四面体的方案（只修改split_info中等于-1的元素）
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
			//优先选择与对面斜边相反的方向
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
//把一个六面体分成6个四面体（分法根据六个面的具体情况而定）
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
	vector<Tetrahedron> local_tets;  //以局部节点编号表示的四面体

	//四条对角线
	//1：1－7
	//2：2－8
	//3：3－5
	//4：4－6
	int diag_line = 0;
	int vertexes[8]={0,0,0,0,0,0,0,0};

	//存放12个三角面片
	vector<Tri_2D_ele> tri_faces;
	tri_faces.reserve(20);

	//标志12个三角面是否已经处理过（已经形成了四面体）
	int sign_tri_face[12]={0,0,0,0,0,0,0,0,0,0,0,0};

	//12条边所在的三角面片（每个边有两个）（12×2的数组，第一个元素指三角面片编号，第二个元素指对面的节点号）
	deter_rela_tri_face(hex,tri_faces);

	vector<vector<int> > temp_v1;
	vector<vector<vector<int> > > node_for_face(12,temp_v1);  //记录每个可以跟三角面片构成四面体的节点

	int sign_diag_line = 0;       //记录在处理哪个三角面片时构造的对角线

	//记录是否已经搜索过三角面片的可用节点（能与之构成四面体的节点）
	int sign_dealed_tri_faces[12];
	for( int i=0; i<12; i++ )
	{
		sign_dealed_tri_faces[i] = 0;
	}

	int succeeded = 1;
	vector<int> i_vec; //保存上一个操作过的三角面片，以备后退时使用
	for( int i=0; i<12&&i>=0; i++)
	{
		if( sign_tri_face[i] > 0 ) continue;		//表示已经经过处理
		if(debuging)
		{
			hout << "-------------------------------------" << endl;
			hout << "tri face " << i+1 << endl;
			hout << "sign_dealed_tri_faces: " << sign_dealed_tri_faces[i] << endl;
		}
		//收集能与第i个三角面片构成四面体的节点
		if( sign_dealed_tri_faces[i] == 0 )
		{
			i_vec.push_back(i);  //记录搜索过可用节点的三角面片编号，以备后退时使用
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
			//确定需要检查的4个节点（三角面片对面的四个节点，并把最近的节点放在最后）
			int check_nn[4];
			deter_check_nodes(i,tri_faces[i].nodesId[1],check_nn);
			for( int jj=0; jj<4; jj++ )
			{
				int j=check_nn[jj];
				int ls[3];
				if( vertexes[j] == 0 )
				{
					//检查协调性（与已有边不能交叉）
					int cf1t = can_form_1tet(tri_faces[i],j,hex,diag_line,ls);
					//检查需要用到的三角面片是否已经使用过
					if( ls[0] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[0],tri_faces[i].edgesId[0]);
						if( fn < i ) continue;  //编号较小的三角面片一定已经处理过了
						if( sign_tri_face[fn] > 0 ) continue;
					}
					else if( ls[1] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[1],tri_faces[i].edgesId[0]);
						if( fn < i ) continue;  //编号较小的三角面片一定已经处理过了
						if( sign_tri_face[fn] > 0 ) continue;

						fn = deter_face_num_from_2bar(ls[1],tri_faces[i].edgesId[1]);
						if( fn < i ) continue;  //编号较小的三角面片一定已经处理过了
						if( sign_tri_face[fn] > 0 ) continue;

					}
					else if( ls[2] < 12 )
					{
						int fn = deter_face_num_from_2bar(ls[2],tri_faces[i].edgesId[1]);
						if( fn < i ) continue;  //编号较小的三角面片一定已经处理过了
						if( sign_tri_face[fn] > 0 ) continue;
					}

					if( cf1t == 1 )
					{
						//检查角度（两个三角面片的夹角不能过大）
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
			//生成一个四面体
			local_tets.push_back(Tetrahedron(tri_faces[i],act_node[0]));
			//如果三条新边中有两个为斜边（四面体正好取六面体的一个角），则应标记角点，防止下次再搜索
			vector<int> tf;
			int tn = 0;

			//是否有直边(最多一条）
			int sn = -1;
			for( int j=0; j<3; j++ )
			{
				if( act_node[j+1] < 12 )
				{
					sn = act_node[j+1];
				}
				else if( act_node[j+1] < 18 )
				{
					tf.push_back(act_node[j+1]-12);   //记录斜边编号
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
			//存在直边，需要标记已经处理过的三角面片，否则会生成重复单元
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
			//构造失败，清除在这一步所做的标记并回退到上一步
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
//按给定参数，把六面体（长方体）的6个面分成三角面片，并确定其相互关系（每条边线两边的单元）
int Mesher::deter_rela_tri_face(int *hex,vector<Tri_2D_ele> &tri_faces)
{
	Tri_2D_ele temp_tri;
	tri_faces.assign(12,temp_tri);
	//生成12个三角面片
	//所有的节点编号以及边的编号均为右手定则，拇指方向朝六面体内部，且斜边在最后
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
//确定需要检查的4个节点（三角面片对面的四个节点，并把最近的节点放在最后）
        //节点的顺序为：长对角线节点、第一条边对应的斜边节点、第二条边对应的斜边节点、中点对应的节点（构成直边的节点）
//i：三角面片的编号
//n：三角面片第二个节点的编号（局部编号）
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
//给定一个三角面片和一个节点，判断它们是否能构成一个四面体，主要考虑边界协调问题
//ls[3]: 三条新边的编号（在六面体中的编号）
int Mesher::can_form_1tet(Tri_2D_ele &tri, int nn, int *hex, int &diag_line, int ls[3])
{
	bool debuging=false;
	//判断三条新边是否跟已有边界协调（所有节点编号均为局部编号）
	for( int i=0; i<3; i++ )
	{
		ls[i] = -1;
		int n0 = tri.nodesId[i];
		if( debuging ) hout << "n0: " << n0 << " nn: " << nn <<endl;
		//是否重复
		if( n0 == nn ) return 0;
		//是否为直边
		int ieoh = is_edge_of_hex(n0,nn);
		if( ieoh >= 0 )
		{
			ls[i] = ieoh;
			continue;
		}
		//是否为斜边，如果是，要检查是否协调
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
		//是否为对角线
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
//根据给定的两个节点号，判断是在六面体的哪个侧面上
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
//给定两条直边号，返回构成的三角面片号
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
//计算给定节点到给定三角片面的夹角(拆分六面体成6个四面体时使用）
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
//计算两个三角面片之间的夹角,node1_num node2_num 分别为不在公共边上的节点
double Mesher::angle_of_2tri(Tri_2D_ele* tri1, Tri_2D_ele* tri2, int node1_num, int node2_num)
{
	//检查tri1 tri2的有效性（不能有节点重复）
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

	//计算tri1的法向量
	Node node21 = nodes_vec[tri1->nodesId[1]] - nodes_vec[tri1->nodesId[0]] ;
	TDVector v11 = TDVector( node21.x, node21.y, node21.z );
	Node node31 = nodes_vec[tri1->nodesId[2]] - nodes_vec[tri1->nodesId[0]] ;
	TDVector v12 = TDVector( node31.x, node31.y, node31.z );
	TDVector normal1 = v11.cro_product( &v12 );
	//计算tri2的法向量
	node21 = nodes_vec[tri2->nodesId[1]] - nodes_vec[tri2->nodesId[0]] ;
	TDVector v21 = TDVector( node21.x, node21.y, node21.z );
	node31 = nodes_vec[tri2->nodesId[2]] - nodes_vec[tri2->nodesId[0]] ;
	TDVector v22 = TDVector( node31.x, node31.y, node31.z );
	TDVector normal2 = v21.cro_product( &v22 );

	//计算两个法向量的夹角
	double acos_v = normal1.dot_product( &normal2 ) / ( normal1.length() * normal2.length() );
	double theter1 = acos( acos_v );
	//两面片夹角为法向夹角的余角
	double theter2 = PI - theter1 ;

	//判断两面片连接处是凹形还凸形
	//找出不在公共边上的两点
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
//把一个六面体分成12个四面体（在中心插入一个新节点，用来处理不能成功分成6个四面体的情况）
int Mesher::hex2_12tet(int hex[14], vector<Tetrahedron > &tets, double volume)
{
	tets.clear();

	//     vector<Tetrahedron> local_tets;
	//分成12个四面体（在六面体中心插入一个节点）
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
//确定每个四面体的材料
int Mesher::deter_eles_mat()
{
	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		deter_mat_ele(i);
	}
	return 1;
}

//---------------------------------------------------------------------------
//确定一个单元的材料编号
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
//确定一个单元的材料编号
//返回值：-1: 界面上的单元
//					0,1 ：基体，增强项
//					-3	: 错误，应该退出程序
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

	//确定此单元的材料
	//在基体中
	if( (wws[0] == -1 || ions[0] == 1) && (wws[1] == -1 || ions[1] == 1) &&
		(wws[2] == -1 || ions[2] == 1) && (wws[3] == -1 || ions[3] == 1) &&
		((wws[0] == -1 || wws[1] == -1 || wws[2] == -1 || wws[3] == -1)||
		(!(wws[0] == wws[1] && wws[0] == wws[2] && wws[0] == wws[3]))))
	{
		//四个节点都在基体内或者界面上，而且至少一个节点在基体内或者
		//四个节点在不同的界面上，都说明此单元为基体

		//news：注意如果椭球生成程序有问题可能出现颗粒重叠区域，
		//这种区域内的单元可能会出现所有的is_on都是1的情况
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
	//单元在增强颗粒内部
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
	//在界面上的单元
	act_tet->materialId = -1 ;
//	act_tet->materialId = 0 ;	//都按在基体上的情况算
	return -1;
}

//---------------------------------------------------------------------------
//判断一个单元是否为薄元（节点到对面三角形面片的距离小于global_length/5.0);
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

	//检查单元有效性（是否已经删除（四个节点相同））
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
//找出一个四面体四个面的最大夹角，并返回相应的节点序列（0-1-2  3-2-1为夹角最大的两面）
double Mesher::max_angle_of_tet(int ele_num, int *nseq)
{
	int nns[6][4];   //六组节点编号
	//共有六种情况（四面体任意两个面的交角）
	//手写吧，麻烦
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

	int nsn[4];        //四个节点编号
	//找出角度最大的两个面
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
//把指定单元分类放入待处理集合（跨边界的单元）
int Mesher::put_into_cbev( int ele_num, double alength, int type )
{
	bool debuging = false;

	Tetrahedron *act_tet = &eles_vec[ele_num];
	//首先找出有一个边在单胞边界线上的单元
	if( type < 1 && has_cross_edge_edge( act_tet ) != 0 )
	{
		if( debuging ) hout << " Putting ele " << ele_num << " into be." <<endl;
		eles_across_boundary_be.push_back( ele_num );
		return 1;
	}

	//找出有一个面在单胞边界面上的单元
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
//判断一个单元是否有边在单胞边界线上跨边界
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
	//检查节点是否同时在两个面上，如果是，则说明此点在两面交线上
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

	//检查界面
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
//检测四面体是否有边在边界线
//返回值为在边界线上的边号
//      0: 没有边在边界线上
//      1: 0-1 边
//      2: 0-2 边
//      3: 0-3 边
//      4: 1-2 边
//      5: 1-3 边
//      6: 2-3 边
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
//判断两个节点是否在同一条边界上，mod控制要不要检查增强颗粒和基体交界面
int Mesher::is_on_same_edge(int node1_num, int node2_num, int mod)
{
        Node *node1 = &nodes_vec[node1_num];
        Node *node2 = &nodes_vec[node2_num];

        //首先确定两个节点分别在哪些面上（哪些面的交）
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
        //只要有两个面相同，就说明此两点在同一边界线上
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
//判断一个单元是否有边在单胞边界面上跨边界
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
			//检查每条边是否在单胞边界面上，并且跨界面
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
//检查给定两个节点是否在同一个面上（单胞边界面或者基体与增强项的界面）
//mod=0只检查单胞边界面，
//mod=1检查单胞边界面和基体与增强项的界面
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
//确定应该调整的边
//首先选择边界线上跨边界的边，其次选择边界面上跨边界的边，再其次选择两点都能移动的边
//最后选择一点在边界面上，一点在单胞内部的边
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

	//检查是否有在边界线上跨边界的边
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

	//检查是否有在边界面上跨边界的边
	int hceeof = has_cross_edge_edge_on_face( act_tet );

	Node *n1 = &nodes_vec[act_tet->nodesId[0]];
	Node *n2 = &nodes_vec[act_tet->nodesId[1]];
	Node *n3 = &nodes_vec[act_tet->nodesId[2]];
	Node *n4 = &nodes_vec[act_tet->nodesId[3]];
	Point pp[]={	{n1->x, n1->y,n1->z},
						{n2->x, n2->y,n2->z},
						{n3->x, n3->y,n3->z},
						{n4->x, n4->y,n4->z} };
	//优先选取两端节点都能移动的跨边界边
	//首先找出所有跨边界的边
	vector< vector<int> > cross_edges ;
	for( int i=0; i<6; i++ )
	{
		//left_nn right_nn局部编号
		//left_ndoeN right_nodeN整体编号
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
			//确定两端节点是否能移动
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

	//检查有没有颗粒重叠的情况
	//（即有边跨越两个颗粒，颗粒又有重叠部分，
	//颗粒生成算法的问题，但需要处理，唉，愁人啊，于艳真会给我出难题）
	//呵呵 余新刚师兄真有情趣
	vector<vector<int> > cross_edges1(cross_edges.begin(),cross_edges.end());
	cross_edges.clear();
	for( int i=0; i<(int)cross_edges1.size(); i++ )
	{
		int ln = cross_edges1[i][1];
		int rn = cross_edges1[i][2];
		if( ww[ln]!=-1 && is_on[ln]==0 && ww[rn]!=-1 && is_on[rn]==0 )
		{
			//在两个颗粒内
			//下面检查两个颗粒是否重叠（求两个交点，利用距离判断）
			Point inter_point1, inter_point2;
			int err1 = thecell->surfaces_vec[ww[ln]]->intersect(pp[ln],pp[rn],inter_point1);
			if( err1 == 0 ) continue;
			int err2 = thecell->surfaces_vec[ww[rn]]->intersect(pp[ln],pp[rn],inter_point2);
			if( err2 == 0 ) continue;
			double dis1 = pp[ln].distance_to(inter_point1);
			double dis2 = pp[ln].distance_to(inter_point2);
			if( dis1 > dis2 )
			{
				//重叠啦
				//强制认为左侧节点在右侧节点所在颗粒的边界上，以方便后边确定
				//最好是想办法标志此单元为跨边界单元，嗯，以后再说
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

	//检查所有跨边界的边，优先选择两端点都能移动的边
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
	//选择节点编号较大的边
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
	//检查所有跨边界的边，优先选择两端点都都在增强体内的边（防止死循环）
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
//节点是否在单胞的界面上
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
//检测四面体是否有面在边界面上
//返回值为在边界面上的面号
//      0: 没有面在边界面线上
//      1: 0-1-2 面
//      2: 0-2-3 面
//      3: 0-3-1 面
//      4: 1-2-3 面
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
//判断一个节点是否为角点（三个面的交点）
int Mesher::is_corner( int node_num )
{
	Node *node = &nodes_vec[node_num];
	//检查节点是否同时在两个面上，如果是，则说明此点在两面交线上
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
//检查边界面上的单元在边界面上的节点是否离边界很近
int Mesher::check_bft_node( Tetrahedron *act_tet, double alength )
{
	//检查面上的点是否离边界很近
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
				Node temp_n = nodes_vec[node1_num] ;//保留副本，以便出现错误时恢复
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
//计算和节点相关的单元的体积最小值
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
//计算四面体单元的体积(有向体积）
double Mesher::cal_tet_volume( Tetrahedron* tet )
{
	return cal_volume(	 &nodes_vec[tet->nodesId[0]],
								 &nodes_vec[tet->nodesId[1]],
								 &nodes_vec[tet->nodesId[2]],
								 &nodes_vec[tet->nodesId[3]]	);
}

//---------------------------------------------------------------------------
//检查有无四个节点都在界面上的单元，有则拆之
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
//处理四个节点都在界面上的单元（在长边中间插入节点）
int Mesher::deal_grov_eles(vector<int> &grov_eles, vector<int> &new_eles_num)
{
	for( int i=0; i<(int)grov_eles.size(); i++ )
	{
		rectify_tet(grov_eles[i],new_eles_num);
	}
	return 1;      
}

//---------------------------------------------------------------------------
//检查所有的单元，如果有两个面的夹角非常大（接近180），则在对边插入一个节点
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

	//判断是否四个节点全部在界面上
	//（顺便检查是否四个节点中，其相关单元为薄元的个数不能超过2个，否则拆分不收敛（不停的拆））
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
				//与此节点相连的所有单元（除了当前薄单元）都是基体单元，直接提出去
				where_is_nodes[nn] = -1;
				is_on_nodes[nn] = 0;
				//eles_vec[ele_num].materialId = 0;
				//把所有此节点的相关单元改成基体材料，顺便光顺一下下，hoho
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
		//可能是颗粒内部的薄元，当与其相接的有其他薄元时，不处理了（不行，处理，晕）
		if( re_all_on_face > 1 ) return 0;
	}
	//当与其相接的还有其他很多薄元时，不处理了，发现有可能出现死循环，烦
	if( re_all_on_face > 2 ) return 0;

	int itt = is_thin_tet(ele_num);
	if( itt == 0 ) return 1;
	//把单元分解成4个三角形，删除
	//首先找出角度最大的两个面
	int nsn[4];        //四个节点编号
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
		//计算新节点位置（四个节点的平均）
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
			//如果是四个节点都在椭球面的单元，考虑是否需要将新节点投影到球面
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

		//删除旧单元
		//先标记，回头再删
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
		//在角度最大的对边插入节点
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
//检查一个节点的相关单元中有几个薄元
int Mesher::num_of_thin_tet_re(int nn)
{
	bool debuging = false;

	int te_num=0; //相关单元中薄元数目
	for( int i=0; i<(int)nodes_vec[nn].relative_eles_vec.size(); i++ )
	{
		int ele_num = nodes_vec[nn].relative_eles_vec[i];
		int ist = is_thin_tet(ele_num);
		if( ist == 1 ) te_num ++;
	}
	return te_num;
}

//---------------------------------------------------------------------------
//判断一个节点的所有相连节点是否有在颗粒内（检查节点所有的相关单元的节点，有节点在颗粒内，返回0，否则（全部在基体内或在颗粒表面）返回1）
int Mesher::can_be_extracted(int n1)
{
	if( nodes_vec[n1].flag >= 0 ) return 0;  //单胞侧面上的点
	for( int i=0; i<(int)nodes_vec[n1].relative_eles_vec.size(); i++ )
	{
		int en = nodes_vec[n1].relative_eles_vec[i];
		int ic = eles_vec[en].is_contain(n1);
		if( ic == -1 ) continue;   //已经不再包含节点n1了（没有更新的缘故）
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
//计算给定两个节点相关单元的材料号（同时包含两个节点的所有单元颗粒编号,-1表示基体）
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
//在两个节点中间插入节点
int Mesher::insert_node(int node1_num, int node2_num, int node3_num,
									vector<int> &new_ele_num)
{                
	bool debuging = false;

	//检查节点Node1_num的所有相关单元
	for( int i=0; i<(int)nodes_vec[node1_num].relative_eles_vec.size(); i++ )
	{
		//检查该单元是否包含给定两节点
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

		//确定是哪条边
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

		//确定新四面体的材料属性
		deter_mat_ele(ele_num);
		deter_mat_ele(new_tet_num); 
	}
	return 1;
}

//---------------------------------------------------------------------------
//处理边界线上的单元(在边界线上跨边界的单元）
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
			//此单元为薄元
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
//处理边界面上的单元(在边界面上跨边界的单元）
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
			//检查当前四面体是否仍有面在边界面上，且边界面上的面跨越基体和纤维边界
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
			//此单元为薄元
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
//处理两个端点都可以移动的跨边界单元
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
			//此单元为薄元
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
//处理两个端点中有一个不可以移动的跨边界单元
int Mesher::deal_with_eles_across_boundary_ndmb(vector<Bar> &new_edges,int debuging)
{
	bool debuging_map_mesh = false; 
	if( debuging == 1 )
	{
		debuging_map_mesh = true;
	}
	//如果移动节点失败（出现负体积单元），暂时存放在此向量中，最后再处理
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
			//此单元为薄元
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
//处理跨边界单元
//ele_kind:单元处理类型(在边界线上跨界面的单元、在边界面上跨边界的单元......)
//mod：当发现移动节点后出现负体积单元时的处理方法  （废弃，没用）
//      0：不作处理，程序返回2，
//      1：插入节点
int Mesher::deal_with_ele_across_boundary( int ele_num, int ele_kind, vector<Bar> &new_edges,  int mod=1 )
{
	bool debuging_map_mesh = false ;
	if( mod ==1 ) debuging_map_mesh = true;

	if( debuging_map_mesh ) hout << "++++++++++++++++++++++++++++" << endl;

	if( debuging_map_mesh) hout << "ele_num: " << ele_num << " ele type: " << ele_kind << endl;

	double alength = global_length*0.8;
	if( ele_kind == 1 ) alength = global_length*0.5;   //侧面上跨边界单元
	double min_dis = global_length*3.0/8.0;
	int wh[4],is_on[4];
	for( int j=0; j<4; j++ )
	{
		wh[j] = where_is_nodes[eles_vec[ele_num].nodesId[j]];
		if( wh[j] == -3 )			//发现节点同时在两个增强项内
		{ 
			return 0;
		}
		is_on[j] = is_on_nodes[eles_vec[ele_num].nodesId[j]];
	}
	//确定哪条边跨边界
	int left_movable, right_movable;
	int edge_num = deter_oper_edge( &eles_vec[ele_num], left_movable, right_movable );

	if( edge_num == 0 ) return -1;     //没有边跨边界

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

	//求与增强材料界面的交点       
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

		//检查将要移动的两个节点的相邻节点是否有3个节点都在交界面上的情况
		int is_all_on1 = is_re_all_on_face( node1_num, *pww );
		int is_all_on2 = is_re_all_on_face( node2_num, *pww );

		if( debuging_map_mesh) hout << " is_all_on1: " << is_all_on1 << " is_all_on2: " << is_all_on2 << endl;

		//经验表明：当节点离界面很近时，无论如何应该移动到界面上，即使is_all_on==1;
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

		//确定移动哪个点到界面上
		int move_which_node = 0;
		if( left_movable == 1 )
		{
			if( right_movable == 1 )
			{
				//两端节点都能移动，取距离小的哪个
				if( dis2 < dis1 ) move_which_node = 1;  //移动右端节点
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
				move_which_node = -1;  //都不能移动 
			}
		}

		if( debuging_map_mesh) hout << "dis1: " << dis1 << " dis2: " << dis2 << " alength: " << alength << endl;
		int moved_nn = -1 ;
		if( move_which_node == 0 )
		{
			if( dis1 < alength )
			{
				if( debuging_map_mesh) hout << " dis1 < dis2, move node " << node1_num << " to the new position.\n";
				Node temp_n = nodes_vec[node1_num] ;//保留副本，以便出现错误时恢复
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
				Node temp_n = nodes_vec[node2_num] ;//保留副本，以便出现错误时恢复
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

		//操作过的节点（移动过，或者新插入的节点号，以备后面的处理）
		int oper_nn=moved_nn;

		//移动了节点
		if( debuging_map_mesh )	hout << "moved node num: " << moved_nn << endl;

		if( moved_nn != -1 )
		{
			//检查所有相关单元的材料属性
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
			int node_num = (int)nodes_vec.size();    //新生成节点号
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


			//检查每个包含此边的四面体，并将其分裂成两个，以保证单元的协调
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
//检查给定节点的所有相关单元是否四个节点全在边界面上
int Mesher::is_re_all_on_face( int node_num, int wh )
{
	bool debuging = false;

	int aof_num=0; //相关单元中四个节点全在界面上的单元数（通常不要超过两个，否则很难拆分消除）
	for( int i=0; i<(int)nodes_vec[node_num].relative_eles_vec.size(); i++ )
	{
		int ele_num = nodes_vec[node_num].relative_eles_vec[i];

		//检查是否已经被删除（四个节点号相同）（烦啊，不知道还有多少地方没有检查单元是否已经删除）
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
//根据给定的边的序号，确定节点序列，进行单元拆分的时候使用
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
//检查每个相关单元（四面体），并在节点Node1_num 和 Node2_num之间插入节点node3_num
void Mesher::harmonize_tets( int node1_num, int node2_num, int node3_num,
											 double alength, vector<Bar> &new_edges, int ele_kind )
{
	bool debuging = false;

	//检查节点Node1_num的所有相关单元
	for( int i=0; i<(int)nodes_vec[node1_num].relative_eles_vec.size(); i++ )
	{
		//检查该单元是否包含给定两节点
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

		//处理新生成的边
		for( int k=0; k<4; k++ )
		{
			if( k == node1_i || k == node2_i ) continue;
			Bar newbar(node3_num,eles_vec[ele_num].nodesId[k]);
			new_edges.push_back(newbar);
		}

		//确定是哪条边
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
//采用快速排序法 (从小到大）
void Mesher::quick_sort( vector<int> *int_vec, int min_n, int max_n )
{
	if( min_n == -1 ) min_n = 0;
	if( max_n == -1 ) max_n = (int)int_vec->size()-1;
	if( max_n <= 0 ) return;
	int middle = (*int_vec)[(int)((min_n + max_n)/2)];
	int i = min_n;
	int j = max_n;
	do{
		while(((*int_vec)[i]<middle) && (i<max_n))//从左扫描大于中值的数
		{       
			i++;
		}
		while(((*int_vec)[j]>middle) && (j>min_n))//从右扫描小于中值的数
		{        
			j--;
		}
		if(i<=j)//找到了一对值，交换
		{                         
			int temp_v = (*int_vec)[i];
			(*int_vec)[i]=(*int_vec)[j];
			(*int_vec)[j]=temp_v;
			i++;
			j--;
		}
	}while( i <= j );//如果两边扫描的下标交错，就停止（完成一次）

	//当左边部分有值(left<j)，递归左半边
	if( min_n < j ) quick_sort( int_vec, min_n, j );
	//当右边部分有值(right>i)，递归右半边
	if( max_n > i ) quick_sort( int_vec, i, max_n );
}

//---------------------------------------------------------------------------
//计算所有单元最小体积和最大体积（数组min_max_en[2]保存最小最大体积的单元号，min_max_vl[2]保存最小最大体积值）
int Mesher::cal_min_max_tet_volume(int min_max_en[2], double min_max_vl[2])
{
	//检查每个单元的体积
	double mmin_volume = 1.0e88 ;//记录最小的体积
	int min_vol_ele_num = -1;
	double mmax_volume = 0.0 ;		//记录最小的体积
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
//计算所有单元的体积分布情况（最小体积向下取整，最大体积向上取整，记录长度）
//将这段长度分成n份（一般情况n取10），记录每份中所占单元数与总单元数的比值；
void Mesher::tet_vol_distrib( double min_max_vl[2], string str, int n, int mod )
{
	//用于记录每份所含四面体的个数
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
		if ( !disout )  hout << "无法打开文件：" << file << endl;
		disout << "四面体体积分布" << endl;
		disout << "最小体积：" << min_max_vl[0] << endl;
		disout << "最大体积：" << min_max_vl[1] << endl;
		m=(int)eles_vec.size();
		disout << "总四面体个数" << m << endl;
		for( int i=0; i<n; i++ )
		{
			disout << i <<" :  ";
			disout << "体积在" << vol_min+i*vol_ele << "~" << vol_min+(i+1)*vol_ele << "之间："; 
			disout << "四面体个数：" << tetnum[i];
			disout << "百分比：" <<	100*(tetnum[i]*1.0/m) << "%" << endl;
		}
	}
}

//---------------------------------------------------------------------------
//计算所有四面体单元的质量度量分布情况（最小值0，最大值取1，记录长度）
//将这段长度分成n份（一般情况n取10），记录每份中所占单元数与总单元数的比值；
void Mesher::tet_quality( string str, int n, int mod )
{
	//用于记录每份所含四面体的个数
	vector<int> tetnum(n,0) ;
	double tria[4];	//用于记录四面体四个面的面积(triangle_area)
	double silp[3];	//用于记录四面体对边边长之积(side_length_product)
	double ir;			//内切球半径
	double R;			//外接球半径
	double rou;       //比值rou
	double volume_ele; //四面体体积
	int m;

	for( int i=0; i<(int)eles_vec.size(); i++ )
	{
		//计算四面体体积
		volume_ele = cal_tet_volume( &eles_vec[i] );
		//计算四面体四个面的面积
		tet_tri_area( &eles_vec[i] , tria );
		//计算内切球半径
		ir=(3.0*volume_ele)/(tria[0]+tria[1]+tria[2]+tria[3]);
		//计算四面体对边边长之积
		side_length_pro(&eles_vec[i] , silp );
		//计算外接球半径
		R=sqrt((silp[0]+silp[1]+silp[2])*(silp[0]+silp[1]-silp[2])*(silp[0]+silp[2]-silp[1])*(silp[1]+silp[2]-silp[0]))/(24*volume_ele);
		//计算比值rou
		rou=3*ir/R;
		if(rou<0.0-ZERO||rou>1.0+ZERO)
		{
			hout << "在计算第" << i << "个单元时,比值rou=" << rou <<" 小于0.0或大于1.0,出错！" <<endl;
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
		if ( !disout )  hout << "无法打开文件：" << file << endl;
		disout << "四面体质量系数分布" << endl;
		m=(int)eles_vec.size();
		disout << "总四面体个数" << m << endl;
		for( int i=0; i<n; i++ )
		{
			disout << i <<" :  ";
			disout << "质量系数在" << i*(1.0/n) << "~" << (i+1)*(1.0/n) << "之间："; 
			disout << "四面体个数：" << tetnum[i];
			disout << "百分比：" <<	100*(tetnum[i]*1.0/m) << "%" << endl;
		}
	}
}
//---------------------------------------------------------------------------
//用于计算四面体四个面的面积
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
//用于计算四面体对边边长之积
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
		//三种情况分别是：(0,1;2,3),(0,2;1,3),(0,3;1,2), 用n[2][2]记录
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
//消除整型向量中连续出现的重复元素
void Mesher::unique( vector<int> *int_vec )
{
	if( int_vec->size() < 2 ) return ;
	//复制副本
	vector<int> int_vec_cop = *int_vec;
	//清空向量
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
//mod=0:拉普拉斯光滑处理，即对任一节点，取其坐标为与之相邻的节点坐标平均值
//mod=1:取相邻节点所围区域的中心
int Mesher::smoothing_3d(int mod)
{
	bool debuging_smoothing = false;
	int min_vol_en;
	double min_min_volume_ori=min_volume_of_cell(min_vol_en);  //光顺前最小体积
	for( int i=0; i<(int)nodes_vec.size(); i++ )
	{
		if( debuging_smoothing ) hout	<< "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if( debuging_smoothing ) hout	<< "node num: " << i << endl;
		if( debuging_smoothing ) hout	<< "old location: " << nodes_vec[i].x << " "
														<< nodes_vec[i].y << " "
														<< nodes_vec[i].z << endl;
		double min_volume_ori = cal_relative_min_volume( i ); //光顺前该节点相关单元的最小体积
		if( debuging_smoothing ) hout << " min_volume_ori: " << min_volume_ori << endl;
		double sum_x=0, sum_y=0, sum_z=0; //将此节点的邻节点坐标求和
		int num_node=0;                   //将多少节点坐标求和
		//对于边界上的节点，只有在同样边界上的节点，才能对其求平均
		//有两种情况：边界线上的节点，边界面上的节点
		int ioe =  is_on_edge( i ) ;
		int iof =  is_on_face( i ) ;
		vector<int> ioe_nodes; //和当前节点在同一边界线上的节点集
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

			Node old_position = nodes_vec[i];	//保留副本，以便出现错误时恢复
			if( ioe != 0 )
			{
				if( debuging_smoothing ) hout << " on the edge! " << endl;
				if( num_node > 2 ) continue;			//角点
				//找出与此节点相邻的两个同边界线上的节点连线的中垂线与椭球表面的交点
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
					//在单胞边界线上
					Node temp_node = (*node1 + *node2)/2.0;			//韩非20061221改
					node3->x = temp_node.x;
					node3->y = temp_node.y;
					node3->z = temp_node.z;
					continue ;
				}
				//在椭球和基体的交界线上
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
				//求与椭球面的交点
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
					//如果邻节点超过3个（可用来光滑的节点，比如，如果节点在面上
					//则只有在同一面上的邻节点才可用来光滑处理）
					//对所有邻节点，依次选取3个，求法线方向，循环到最后一个邻节点
					//然后对所有的法线取平均，作为投影方向
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
			if( min_volume_new < min_volume_ori )	//相关单元最小体积减小，恢复
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
	double min_min_volume_new=min_volume_of_cell(min_vol_en);  //光顺后最小体积

	if( debuging_smoothing ) hout << " min_min_volume_ori: " << min_min_volume_ori << " min_min_volume_new: " << min_min_volume_new << endl;
	if( min_min_volume_new - min_min_volume_ori > fabs(min_min_volume_ori) * 0.01 ) return 1;

	return 0;
}

//---------------------------------------------------------------------------
//计算所有单元中的最小体积
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
//计算纤维体积比
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

	//hout << "    增强颗粒体积: " << non_matrix_vol << " 总体体积: " << total_vol << endl;  
	return non_matrix_vol/total_vol;
}

//---------------------------------------------------------------------------
//确定单胞边界上的节点
void Mesher::deter_boundary_nodes(int coat)
{
	bnodes_vec.clear();

	for ( int i=0; i<(int)nodes_vec.size(); i++)
	{
		int iof = is_on_face(i);   //检查节点位置
		if( iof < 0 )
		{
			//在单胞边界上
			bnodes_vec.push_back(i);
			nodes_vec[i].flag = 1;
		}
		else if( iof == 0 )
		{
			//在单胞内部，不在单胞边界面上，也不在颗粒界面上
			nodes_vec[i].flag = 0;
		}
		else
		{
			//在颗粒界面上
			nodes_vec[i].flag = 0;			//可以用2表示在颗粒界面上
														//为了不与后面的程序冲突，这里先改为0
		}
	}
}
//---------------------------------------------------------------------------
//生成边界层网格
int Mesher::gen_coating_mesh(double thick_ratio, int layer_num)
{
	vector<EllipsGeometry> ellgeo(thecell->surfaces_vec.size());

	//-------------------------------------------------------------------------------------------------------------
	//按椭球整理椭球内节点、单元等数据
	int ele_size = int(eles_vec.size());
	vector<int>  eles_ellip(ele_size,-1);
	for(int i=0; i<ele_size; i++)
	{
		//if(i==566||i==571||i==1146||i==1148)
		//{
		//	hout << "i= " << i << "  mat= " <<  eles_vec[i].materialId<<endl;
		//	hout << "nodesId：";
		//	for(int k=0 ; k<4; k++)
		//	{
		//		hout << eles_vec[i].nodesId[k]+1 <<"  ";
		//	}
		//	hout << endl;
		//	hout << "where_is_nodes：";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		hout << where_is_nodes[eles_vec[i].nodesId[j]] <<"  ";
		//	}
		//	hout << endl;
		//	hout << "is_on_nodes：";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		hout << is_on_nodes[eles_vec[i].nodesId[j]] <<"  ";
		//	}
		//	hout << endl;
		//	hout << "where_is：";
		//	for(int j=0 ; j<4; j++)
		//	{
		//		int is_on=0;
		//		hout << where_is(&nodes_vec[eles_vec[i].nodesId[j]], is_on) << "/" << is_on <<"  ";
		//	}
		//	hout << endl;
		//}
		if(eles_vec[i].materialId==1)		//单元属于颗粒
		{
			//---------------------------------------------------------------------
			//判断此单元在属于哪个椭球
			int ele_num;
			int wws[4];
			for(int j=0; j<4; j++)
			{
				wws[j] = 	where_is_nodes[eles_vec[i].nodesId[j]];
			}
			if(wws[0]==wws[1]&&wws[0]==wws[2]&&wws[0]==wws[3])		//四个点所在椭球编号一致
			{
				ele_num = wws[0];
			}
			else																								//四个点所在位置标识不一致
			{
				//求出该单元的重心点坐标
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
				//求此单元所在椭球编号
				int count=0;
				vector<int> tem_num;
				for(int j=0; j<4; j++)
				{
					int wn=wws[j];
					if(wn==-1) continue;												//去掉等于-1的情况，-1表示基体

					int key=0;																//判断这个椭球是否已经被判断过，tem_num中记录判断过的椭球编号
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

					if(thecell->surfaces_vec[wn]->is_contain(&pp)<=0)		//判断四面体重心是否在该椭球中
					{
						ele_num = wn;
						count++;
					}
				}
				if(count==0)
				{
					//判断特例：虽然四面体的重心在该椭球的外部，
					//但因四面体的外顶点是在基体中，可以认为此单元仍属于该椭球。
					if(int(tem_num.size())==1)
					{
						ele_num = tem_num[0];
					}
					else if(int(tem_num.size())>=2)
					{
						eles_vec[i].materialId=0;    //此单元夹在两个以上椭球之间,认为此单元是基体材料
						continue;
					}
					else
					{
						hout << "i= " << i << "  mat= " <<  eles_vec[i].materialId<<endl;
						hout << "nodesId：";
						for(int k=0 ; k<4; k++)
						{
							hout << eles_vec[i].nodesId[k]+1 <<"  ";
						}
						hout << endl;
						hout << "where_is_nodes：";
						for(int j=0 ; j<4; j++)
						{
							hout << where_is_nodes[eles_vec[i].nodesId[j]] <<"  ";
						}
						hout << endl;
						hout << "is_on_nodes：";
						for(int j=0 ; j<4; j++)
						{
							hout << is_on_nodes[eles_vec[i].nodesId[j]] <<"  ";
						}
						hout << endl;
						hout << "where_is：";
						for(int j=0 ; j<4; j++)
						{
							int is_on=0;
							hout << where_is(&nodes_vec[eles_vec[i].nodesId[j]], is_on) << "/" << is_on <<"  ";
						}
						hout << endl;
						hout << "错误！有个四面体不属于任何椭球。" << endl;
						return 0;
					}
				}
				else if(count>1)
				{
					hout << "错误！有个四面体同时属于几个椭球。" << endl;
					return 0;
				}
			}
			//标定单元属于哪个椭球
			eles_ellip[i] = ele_num;
			//插入此四面体单元数据到椭球对象
			ellgeo[ele_num].tet.push_back(i);
		}
		else if(eles_vec[i].materialId==-1)		//单元既不是颗粒也不是基体
		{
			eles_vec[i].materialId=0;
		}
	}

	//遍历所有的椭球单元生成表面片和内外节点
	int node_size = int(nodes_vec.size());
	vector<int> nod_ellip(node_size,-2);
	int ellgeo_size = int(ellgeo.size());
	for(int i=0; i<ellgeo_size; i++)
	{
		vector<vector<int > > line_vec;		//用于判断表面三角形面片是否封闭
		int tet_size = int(ellgeo[i].tet.size());
		for(int j=0; j<tet_size; j++)
		{
			vector<vector<int> > temp_tri;	//临时记录三角形面片（四面体在椭球表面的三角形面片）
			int ell_tet = ellgeo[i].tet[j];			//单元编号
			int nid[4];									//单元四个节点的编号
			for(int k=0; k<4; k++)
			{
				nid[k] =  eles_vec[ell_tet].nodesId[k];
			}
			//记录四面体的四个面片中在该椭球表面的三角形面片
			for(int k=0; k<2; k++)					//循环便利四面体单元的四个顶点
			{
				for(int l=k+1; l<3; l++)
				{
					for(int m=l+1; m<4; m++)
					{
						int key;
						if(is_on_nodes[nid[k]]==0||is_on_nodes[nid[l]]==0||is_on_nodes[nid[m]]==0)   //判断并记录表面三角形面片
						{
							if(deter_node_flag(nid[k])>0&&deter_node_flag(nid[l])>0&&deter_node_flag(nid[m])>0)		//这种情况此面在单胞内
							{
								key=0;
							}
							else																																	//这种情况此面在单胞表面
							{
								key=1;
							}
						}
						else
						{
							key=1;
						}
						if(key==1)			//判断此面是否产生椭球外表面
						{
							vector<int> rel_tet;
							int nk_rel_size = int(nodes_vec[nid[k]].relative_eles_vec.size());					//利用节点的相关四面体向量
							for(int n=0; n<nk_rel_size; n++)
							{
								int rel_ele_num = nodes_vec[nid[k]].relative_eles_vec[n];
								if(rel_ele_num!=ell_tet&&eles_ellip[rel_ele_num] == i)							//该单元本身不算，所找到的单元应该属于此椭球
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
													rel_tet.push_back(rel_ele_num);											//发现此三角形是本身单元与其他单元的公共面片
													goto label_prism;
												}
											}
										}
									}
								}
							}
label_prism:			if(rel_tet.empty())		//此三角形面片没有找到其他相关四面体
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
			//插入内外点数据到椭球对象
			int num_tri = int(temp_tri.size());
			if(num_tri==1)  //只有一个表面三角形，但有可能有一个内点，也有可能一个内点都没有
			{
				int tri[3];
				for(int k=0; k<3; k++)
				{
					tri[k] = temp_tri[0][k];
				}
				for(int k=0; k<4; k++)
				{
					if(nod_ellip[nid[k]]==-2)					//-2表示未使用
					{
						if(k==tri[0]||k==tri[1]||k==tri[2])		//外点
						{
							ellgeo[i].outell_nodes.push_back(nid[k]);
							nod_ellip[nid[k]] = int(ellgeo[i].outell_nodes.size())-1;
						}
					}
				}
			}
			else if(num_tri==2||num_tri==3)
			{
				for(int k=0; k<4; k++)							//全外点
				{
					if(nod_ellip[nid[k]]==-2)					//-2表示未使用
					{
						ellgeo[i].outell_nodes.push_back(nid[k]);
						nod_ellip[nid[k]] = int(ellgeo[i].outell_nodes.size())-1;
					}
				}
			}
			else if(num_tri!=0)
			{
				hout << "错误！四面体的四个三角形面片全是椭球表面三角形面片。" << endl;
				return 0;
			}

			//插入三角形面片信息到椭球对象
			vector<vector<int> > temp_line;
			for(int k=0; k<num_tri; k++)
			{
//---------------------------------------------------------------------------------------------------------------
//<1>以下段落程序用于判断三角形的三条线段是否有线段重合,共三段
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
				//将三角形面片三个顶点的单元编号转换为椭球对象的外点向量中编号
				for(int l=0; l<3; l++)
				{
					temp_tri[k][l] = nod_ellip[nid[temp_tri[k][l]]];
				}
			}
//---------------------------------------------------------------------------------------------------------------
//<2>以下段落程序用于判断三角形的三条线段是否有线段重合,共三段
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
			//插入三角形面片到椭球对象
			ellgeo[i].tri.insert(ellgeo[i].tri.end(),temp_tri.begin(),temp_tri.end());
		}
//---------------------------------------------------------------------------------------------------------------
//<3>以下段落程序用于判断三角形的三条线段是否有线段重合,共三段
		//判断椭球外表面是否封闭
		if(!line_vec.empty())
		{
			//判断是否是和单胞外表面相交产生的部分，这部分留下封闭圈
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
				hout << "错误！椭球外表面三角形不封闭！（要考虑和单胞外表面相交的情况）" << endl;
				return 0;
			}
		}

		//确定椭球内点
		for(int j=0; j<tet_size; j++)
		{
			vector<int > temp_tet_nod;		//记录该四面体中的外点所对应的单元局部编号和在椭球外点向量中的编号
			for(int k=0; k<4; k++)
			{
				int nid =  eles_vec[ellgeo[i].tet[j]].nodesId[k];			//单元四个节点的编号
				if(nod_ellip[nid]==-2)		//-2表示未使用,在所有的外点都确定完之后，留下的就是内点
				{
					ellgeo[i].inell_nodes.push_back(nid);
					nod_ellip[nid]=-1;
				}
				else if(nod_ellip[nid]>=0)
				{
					temp_tet_nod.push_back(k);							//记录在该四面体中的局部编号
					temp_tet_nod.push_back(nod_ellip[nid]);		//记录所对应的在外点向量中的编号
				}
			}
			//插入单元的外点编号信息
			ellgeo[i].tet_nod.push_back(temp_tet_nod);
		}
		//---------------------------------------------------------------------------------------------------------------
		//重新开始nod_ellip向量
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
	//生成棱柱的过程
	//------------------------------------------------------------------------------------------------------
	//把原有四面体单元向量导入混合单元向量
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
	//按椭球收缩网格并加界面层
	for(int i=0; i<int(ellgeo.size()); i++)
	{
		//------------------------------------------------------------------------------------
		//椭球内点收缩(与四面体单元无关)
		for(int j=0; j<int(ellgeo[i].inell_nodes.size()); j++)
		{
			int nodes_n = ellgeo[i].inell_nodes[j];
			double ratio = thick_ratio;
			
			Point point;
			point.x = 0.0;
			point.y = 0.0;
			point.z = 0.0;

			Point node_poi; //记录收缩的初始点
			node_poi.x = nodes_vec[nodes_n].x;
			node_poi.y = nodes_vec[nodes_n].y;
			node_poi.z = nodes_vec[nodes_n].z;

			if (thecell->surfaces_vec[i]->change_coor(&node_poi, &point, ratio)==0)
			{
				hout << i << "点在坐标变换时出错！" << endl;
				hout << "*******************************************" << endl;
				return 0;
			}

			nodes_vec[nodes_n].x = point.x;
			nodes_vec[nodes_n].y = point.y;
			nodes_vec[nodes_n].z = point.z;
		}
		//------------------------------------------------------------------------------------
		//椭球外点收收缩，生成界面层棱柱，改变某些四面体节点编号
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

			Point node_poi; //记录收缩的初始点
			node_poi.x = nodes_vec[nodes_n].x;
			node_poi.y = nodes_vec[nodes_n].y;
			node_poi.z = nodes_vec[nodes_n].z;

			for(int k=1; k<=layer_num; k++)
			{
				double ratio=(k*1.0/layer_num)*thick_ratio;	
				
				if (thecell -> surfaces_vec[i] -> change_coor(&node_poi, &point, ratio) == 0)
				{
					hout << i << "点在坐标变换时出错！" << endl;
					hout << "*******************************************" << endl;
					return 0;
				}

				new_node.x = point.x;
				new_node.y = point.y;
				new_node.z = point.z; 
				new_node.flag = nodes_vec[nodes_n].flag;

				//插入新点
				nodes_vec.push_back(new_node);
				where_is_nodes.push_back(i);
				if(new_node.flag==1) bnodes_vec.push_back(int(nodes_vec.size())-1);
				stem_node[j][k]=int(nodes_vec.size())-1;

				if (j == layer_num)	is_on_nodes.push_back(3);
				else	is_on_nodes.push_back(2);
			}
		}
		//修改某些四面体节点编号
		for(int j=0; j<int(ellgeo[i].tet.size()); j++)
		{
			for(int k=0; k<int(ellgeo[i].tet_nod[j].size());k=k+2)
			{
				elements_vec[ellgeo[i].tet[j]].nodes_id[ellgeo[i].tet_nod[j][k]]=stem_node[ellgeo[i].tet_nod[j][k+1]][layer_num];
			}
		}
		//生成棱柱
		for(int j=0; j<int(ellgeo[i].tri.size()); j++)
		{
			for(int k=0; k<layer_num; k++)
			{
				Element temp_element;
				elements_vec.push_back(temp_element);
				elements_vec.back().mat=2;
				elements_vec.back().type=3;
				//判断三棱柱的体积(为保证有限元计算，棱柱体积必须为正)
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
						cout << "错误！三棱柱的J_val不对。" <<endl;
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
//在不生成边界层的情况下，将四面体单元类型倒为混合单元类型
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
//输出tecplot网格数据
int Mesher::export_tecplot_data()
{
	//输出收缩后椭球四面体网格（按四面体输出）
	ofstream otec("mesh_in_tecplot.dat");
	otec << "TITLE = blend_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)elements_vec.size(); i++)
	{
		switch (elements_vec[i].mat)
		{
		case 0: count[0]++; break; //基体单元
		case 1: count[1]++; break; //椭球单元
		case 2: count[2]++; break; //边界层单元
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
	//输出界面层三棱柱网格（按六面体输出）
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
//在检测时使用
//			int k = elements_vec[i].nodes_id[0];
//			for (int j=1; j <6; j++)
//			{
//				int m = elements_vec[i].nodes_id[j];
//				if (where_is_nodes[k]!=where_is_nodes[m])
//				{
//					hout << "发现不等点！！！" << endl;
//					hout << "i=" << i << endl;
//					for (int n=0; n < 6; n++)
//					{
//						hout << n << "号点" << elements_vec[i].nodes_id[n] << "号点"<< "所在椭球是" << where_is_nodes[elements_vec[i].nodes_id[n]] << endl;
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
//输出三维空间可视化Ensight网格数据
int Mesher::export_Ensight_data()
{
	//输出题头
	//ofstream otec("Composites_Materials_Particles.geo");
	//otec << "geometry meshes" << endl;
	//otec << "composites materials with core/shell particles" << endl;
	//otec << "node id assign" << endl;
	//otec << "element id assign" << endl;
	//otec << "coordinates" << endl;
	////输出节点信息
	//otec << setw(10) << (int)nodes_vec.size() << endl;
	//for (int i=0; i < (int)nodes_vec.size(); i++)
	//{
	//	otec <<right;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].x;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].y;
	//	otec <<setw(12)<<scientific<<setprecision(4)<<nodes_vec[i].z;
	//	otec <<endl;
	//}
	////记录各类单元的个数
	//int count[3]={0,0,0};
	//for (int i=0; i < (int)elements_vec.size(); i++)
	//{
	//	switch (elements_vec[i].mat)
	//	{
	//	case 0: count[0]++; break; //基体单元
	//	case 1: count[1]++; break; //椭球单元
	//	case 2: count[2]++; break; //边界层单元
	//	default: hout << " error!! " << endl; break;
	//	}
	//}
	//cout<<"matrix: "<<count[0]<<"  core"<<count[1]<<"  shell"<<count[2];
	////输出颗粒核单元
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
	////输出颗粒壳层单元
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
	////输出基体单元
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
//三棱柱体积
double Mesher::Threeprism_volume(const vector<Node> &elenodes_vec)
{
	//循环高斯点计算积分；
	double volume=0;
	//计算Ｊ矩阵；
	//--------------------------------------------
	//形函数N对gauss[count].x,gauss[count].y,gauss[count].z的偏导矩阵
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
	//单元节点坐标矩阵
	double elenode[6][3];
	for(int i=0;i<6;i++){
		elenode[i][0]=elenodes_vec[i].x;
		elenode[i][1]=elenodes_vec[i].y;
		elenode[i][2]=elenodes_vec[i].z;
	}
	//--------------------------------------------------
	//J矩阵
	double Jmatrix[3][3];
	//以上两个矩阵的积
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			Jmatrix[i][j]=0;
			for(int k=0;k<6;k++)
			Jmatrix[i][j]=Jmatrix[i][j] + diff[i][k]*elenode[k][j];
		}
		//--------------------------------------------------
		//求出J矩阵的行列式；
		double J_val=Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
			-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
			+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);	
		//--------------------------------------------------
		//由高斯积分点求体积；
		  volume=J_val*0.5*2.0;
	
	  return volume;
}
//---------------------------------------------------------------------------
//确定所有节点的相关单元信息
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
//经过确认materialId ==0 是基体单元，materialId ==1 是椭球单元
//输出生成边界层前的tecplot网格数据
int Mesher::export_tecplot_data_before_coating()
{
	//输出收缩后椭球四面体网格（按四面体输出）
	ofstream otec("mesh_in_tecplot_before_coating.dat");
	otec << "TITLE = tetrahedron_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)eles_vec.size(); i++)
	{
		switch (eles_vec[i].materialId)
		{
		case 0: count[0]++; break; //基体单元
		case 1: count[1]++; break; //椭球单元
		case -1: count[2]++; break;
		default: hout << " error!! " << endl; break;
		}
	}
	//基体
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
	//椭球
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
//人为增加的体分比
void Mesher::Reincrese_vol_raito()
{
	//打开增加的体分比文件
	//用于自然生成椭球的体积无法达到要求时，在剖分后改变基体单元的性质为颗粒，人为增加椭球体积
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
			hout << "      需要增加的体积比分小于等于0%并且大于10%"<<endl;
		}
		else
		{
			//确定所有节点的相关单元信息(确保后面使用正确)
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
			
			//初始化随机数
			srand((unsigned int)time(0));

			//单胞总体积
			double total_vol = thecell->clength*thecell->cwidth*thecell->cheight;

			//定义增加体分比变量
			double increase_ratio = 0.0;
			while(increase_ratio<re_ratio)
			{
				//提取出所有紧邻颗粒材料的基体单元
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
					hout << "matrix_ele.size()==0, 出错！" << endl;
					exit(0);
				}

				//随机改变一个具有以上性质的基体单元为颗粒单元
				//剖分尺寸global_length，大概需要改变单元个数为count = (int)(re_ratio*total_vol)/(pow(global_length,3)/6.0)
				//每次改变总个数的1/6
				int ele_count = (int)(re_ratio*total_vol/pow(global_length,3));
				for(int i=0; i<ele_count; i++)
				{
					int mat_ele_n = (int)(((double)(rand())/RAND_MAX)*(int)matrix_ele.size());
					if(eles_vec[matrix_ele[mat_ele_n][0]].materialId != 0) continue;
					eles_vec[matrix_ele[mat_ele_n][0]].materialId = matrix_ele[mat_ele_n][1];

					//计算改变单元材料性质所增加的体积
					increase_ratio = increase_ratio + cal_tet_volume( &eles_vec[matrix_ele[mat_ele_n][0]] )/total_vol;
					if(increase_ratio>=re_ratio) break;
				}
			}
			hout << "      人为增加体积分数为：" << increase_ratio <<endl;
		}
	}
	else
	{
		hout << "      Reincrease_vol_ratio.dat文件不存在或打不开，检查是否需要做额外增加颗粒体分比的操作！"<<endl;
	}
	in_refile.close();
}
//======================================================================
//读入一行信息，并跳过注释行（以"%"开头）；
string Mesher::Get_Line(ifstream &infile)const{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
