//---------------------------------------------------------------------------

#include "analyzer.h"
#include <fstream>
#include <iostream>
#include <complex>
#include "GloStiffMatrix.h"
#include "Mesher.h"
#include "Hns.h"
#define PI 3.141592654
#define ZERO 1.e-8
using namespace std;


//---------------------------------------------------------------------------
//构造函数
Analyzer::Analyzer( vector<Node> *nodes, vector<Element> *eles, vector<MatPro> *mats ){
        nodes_vec = nodes ;
        eles_vec  = eles  ;
        mats_vec  = mats  ;
 
        //@LLi_epsilon[0] = 1.0;
        //@LLi_epsilon[1] = 1.0;
        //@LLi_epsilon[2] = 1.0;
        //@LLi_epsilon[3] = 1.0;

        //@ox  = 0.0;
        //@oy  = 0.0;
        //@oz  = 0.0;

        //@node_or_ele_strength = 1;

        //@strength_criterion_for_matrix=3;
        //@strength_criterion_for_reinforcement=0;
};
//---------------------------------------------------------------------------
//设定解析解
void Analyzer::set_ana_solution( AnalyticalSolution* as ){
        anaSolution = as ;
};
//---------------------------------------------------------------------------
//设置Na1指针
void Analyzer::set_Na1_vec(vector<vector<double> >* na1){
        Na1_vec = na1;
};  //---------------------------------------------------------------------------
//设置Na1a2指针
void Analyzer::set_Na1a2_vec(vector<vector<double> >* na1a2){
        Na1a2_vec = na1a2;
};
//---------------------------------------------------------------------------
//设置均匀化材料指针
void Analyzer::set_homo_mat(MatPro* mp){
        homoMat = mp;
};
//---------------------------------------------------------------------------
//设置样本数据文件文件名
void Analyzer::set_string_datafile(string datafile)
{
	data_file = datafile;
}
//---------------------------------------------------------------------------
//强度分析 (计算临界载荷，适用于颗粒材料)
int Analyzer::strength_analysis_par(int mod){
        bool debuging_sa = true;
        if( debuging_sa ) hout << "____________________strength_analysis_par()_______________________" <<endl;

   //     reset_cell();
   //     gsmatrix=GloStiffMatrix(nodes_vec,eles_vec,mats_vec);

        anaSolution->print();
                                  
        //初始载荷
        double T0 = 0.0;
        double T1 = anaSolution->load_value();
        
        //对应于载荷T0
        int    *max_stress_en_T0 = new int[mats_vec->size()];
        double *max_stress_vl_T0 = new double[mats_vec->size()];
        for( int i=0; i<(int)mats_vec->size(); i++ ){
                max_stress_en_T0[i] = 0  ;
                max_stress_vl_T0[i] = 0.0;
        };

        //对应于载荷T1
        int    *max_stress_en_T1 = new    int[mats_vec->size()];
        double *max_stress_vl_T1 = new double[mats_vec->size()];

        int cn = 0;
        while( cn < 30 ){
                if( debuging_sa ) hout << "~~~~~~~~~~~~计算临界载荷~~~~~~~~~~~~~~~~~" << "\n";
                if( debuging_sa ) hout << " cn: " << cn++ << "\n";
                if( debuging_sa ) hout << " 载荷：" << anaSolution->load_value() << "\n";

                if( debuging_sa ){
                        //计算均匀化后的理论最大应变向量
                        anaSolution->max_strain( strain_vmax );
                        hout << "均匀化后最大应变：\n    " ;
                        for( int i=0; i<6; i++ ){
                                hout << strain_vmax[i] << " " ;
                        };
                        hout << "\n";
                        //计算均匀化后的理论最大应力向量
                        anaSolution->max_stress( stress_vmax );
                        hout << "均匀化后最大应力：\n    " ;
                        for( int i=0; i<6; i++ ){
                                hout << stress_vmax[i] << " " ;
                        };
                        hout << "\n";
                };

                stress_analysis_par( max_stress_en_T1, max_stress_vl_T1);
                
                if( mod == 0 ) break;
                //检查材料是否已经屈服
                int is_yield = -1 ;  //-1: 没有屈服，0：临界，1：已经屈服
                double critical_load=1.0E200;
//                for( int i=0; i<(int)mats_vec->size(); i++ ){
                for( int i=1; i<=1; i++ ){
                  //      int mat_id = (*eles_vec)[max_Lstress_nn_T1[i]].materialId ;
                        MatPro *mat = &(*mats_vec)[i];
                        if( debuging_sa ) hout << " i: " << i << " mat->type(): " << mat->type_val << "\n";
                        if( mat->type_val == 0){
                                double Ystress = mat->strength;
                                //检查轴向应力
                                if( debuging_sa ) hout << "最大应力：" << max_stress_vl_T1[i] << "\n";

                                if( max_stress_vl_T1[i] < 0 )continue;
                                
                                //根据算出的基体最大等效应力，线性确定临界载荷
                               double Tc = T1-(T1-T0)*(fabs(max_stress_vl_T1[i])-Ystress)
                                                 /(fabs(max_stress_vl_T1[i])-fabs(max_stress_vl_T0[i]));
                                if( Tc < critical_load ) critical_load = Tc;
                                
                                double ds = (fabs(max_stress_vl_T1[i])-Ystress)/Ystress ;
                                if( fabs(ds) < 1e-8 ){
                                        //临界状态
                                        is_yield = 0;      
                                }
                                else if( ds > 0 ){
                                        //已经屈服
                                        is_yield = 1;  
                                };
                        };
                };
                if( is_yield == 0 ){
                        break;
                };
                T0 = T1;
                T1 = critical_load;
                anaSolution->set_load_value(T1);

                for( int i=0; i<(int)mats_vec->size(); i++ ){
                        max_stress_en_T0[i] = max_stress_en_T1[i];
                        max_stress_vl_T0[i] = max_stress_vl_T1[i];
                };

                if( debuging_sa ) hout << " is_yiled: " << is_yield << "\n";
                if( debuging_sa ) hout << " T0: " << T0 << " T1: " << T1 << "\n";
             //   anaSolution->set_load_value( (T0+T1)/2.0 ); 
        };

        hout << "临界载荷：" << anaSolution->load_value() << "\n";
        anaSolution->print();

        //计算均匀化后的理论最大应变向量
        anaSolution->max_strain( strain_vmax );
        double effective_strain = cal_effective_value(strength_criterion_for_matrix,strain_vmax,homoMat);
        strain_vmax.push_back(effective_strain);
        hout << "均匀化后最大应变：\n    " ;
        print_vec6( &strain_vmax, 0, 7 );
        hout << "\n";
        
        //计算均匀化后的理论最大应力向量
        anaSolution->max_stress( stress_vmax );
        double effective_stress = cal_effective_value(strength_criterion_for_matrix,stress_vmax,homoMat);
        stress_vmax.push_back(effective_stress);
        hout << "均匀化后最大应力：\n    " ;
        print_vec6( &stress_vmax, 0, 7 );
        hout << "\n";

        for( int i=0; i<(int)mats_vec->size(); i++ ){
                hout << i << " 号材料：\n";
                hout << "    单元 " << max_stress_en_T1[i] << " 等效应力最大。" << "\n";
                hout << "        其应力为：" << max_stress_vl_T1[i] << "\n";
                hout << "        应力分量：" ;
                print_vec6( &eles_stress[max_stress_en_T1[i]],0,6 );
                hout << "\n";
                hout << "        应变分量：" ;
                print_vec6( &eles_strain[max_stress_en_T1[i]],0,6 );
                hout << "\n";
        };

		//输出样本数据；
		output_Datafile(max_stress_en_T1, max_stress_vl_T1, data_file);

		//输出节点强度等值图
//		export_node_equ_val();

        delete []max_stress_en_T0;
        delete []max_stress_vl_T0;
        delete []max_stress_en_T1;
        delete []max_stress_vl_T1;

        return 1;   
};
//---------------------------------------------------------------------------
//应力分析（计算给定载荷下的应力分布，最大应力，适用于颗粒材料）
int Analyzer::stress_analysis_par(int *max_stress_en, double *max_stress_vl){
        bool debuging = false ;

        //计算每个节点的位移
    //    cal_displacement();

        vector<double> strain(7,0.0) ;	//第7个应变分量空缺
        vector<double> stress(7,0.0) ;	//第7个应力分量里放根据强度准则算出来的应力值

        eles_strain.assign(eles_vec->size(),strain);	
        eles_stress.assign(eles_vec->size(),stress);	

        //初始化
        for( int i=0; i<(int)mats_vec->size(); i++ ){
                max_stress_en[i] = 0;
                max_stress_vl[i] = -1.0e200;
        };

		 for( int i=0; i<(int)eles_vec->size(); i++){
            if( debuging ) hout << "---------------------------stress_analysis_par()-------------------------" << endl;
            if( debuging ) hout << "单元( " << i << " ): " << "\n";
            //材料号                                     
            int mat_id = (*eles_vec)[i].mat;
            MatPro *mat = &(*mats_vec)[mat_id];
            if( debuging ) hout << "        mat id: " << mat_id << "\n";
  //          if( mat->type_val == 0 ){
			if( true ){
                    //计算整体坐标系下的应变、应力
                    cal_glo_strain( &(*eles_vec)[i], strain );
                    cal_glo_stress( &(*eles_vec)[i], &strain, stress );
                    double eff_value = cal_effective_value(strength_criterion_for_matrix,stress,mat);
					stress[6] = eff_value;
                    //if( eff_value > max_stress_vl[mat_id] ){
                    //       max_stress_vl[mat_id] = eff_value ;
                    //       max_stress_en[mat_id] = i ;
                    //};
            };

            if( debuging ){
                    hout << "       应变：" ;
                    for( int j=0; j<6; j++ ){
                            hout << strain[j] << " " ;
                    };
                    hout << "\n";
                    hout << "       应力：" ;
                    for( int j=0; j<6; j++ ){
                            hout << stress[j] << " " ;
                    };
                    hout << "\n";
            };
            if( debuging ) hout << "        应力应变计算完毕。\n";
            eles_stress[i] = stress ;
            eles_strain[i] = strain ;
            if( debuging ) hout << "        应力应变保存完毕。\n";

        };
        
		//平均计算各节点的应力应变
		//存放节点的应力（界面上的节点可能包含多组应力值，格式为：材料号，应力风量（6个），材料号，... //在颗粒材料中可能没有，但纤维材料中有）
		stress.resize(7,0.0);
		strain.resize(7,0.0);

		nodes_stress.assign(nodes_vec->size(),stress);
		nodes_strain.assign(nodes_vec->size(),strain);

		//用于记录有多少个点参与比较危险点
		int node_count = 0;
		//用于记录前num个大的节点的应力值和编号
		int nod_num = 1500;
		vector<double> large_stress_vl(nod_num, -1.0e200);
		vector<int> large_stress_en(nod_num, 0);

		//确定每个节点的位置：内部，材料分界面或者是虽然在内部但是离分界面最近的一层节点
		vector<int> node_mark((int)nodes_vec->size(),0);
		for( int i=0; i<(int)nodes_vec->size(); i++ )
		{
			vector<int> mat_ids;
			int rel_size = int((*nodes_vec)[i].relative_eles_vec.size());
			for( int j=0; j<rel_size; j++ )
			{
				int rel_ele = (*nodes_vec)[i].relative_eles_vec[j];
				bool same_mat = false;
				for(int m=0; m<(int)mat_ids.size(); m++ )
				{
					if( mat_ids[m] == (*eles_vec)[rel_ele].mat )
					{
						same_mat = true;
						break;
					}
				}
				if(!same_mat)
				{
					mat_ids.push_back((*eles_vec)[rel_ele].mat);
					node_mark[i]++;
				}
			}
		}
		for( int i=0; i<(int)eles_vec->size(); i++ )
		{
			bool same_nod = false;
			for( int j=0; j<(int)(*eles_vec)[i].nodes_id.size(); j++ )
			{
				if(node_mark[(*eles_vec)[i].nodes_id[j]]==2)
				{
					same_nod = true;
					break;
				}
			}
			if(same_nod)
			{
				for( int j=0; j<(int)(*eles_vec)[i].nodes_id.size(); j++ )
				{
					if(node_mark[(*eles_vec)[i].nodes_id[j]]==1)
					{
						node_mark[(*eles_vec)[i].nodes_id[j]] = 3;
					}
				}
			}
		}

		for( int i=0; i<(int)nodes_vec->size(); i++ )
		{
			vector<int> mat_ids;
			int rel_size = int((*nodes_vec)[i].relative_eles_vec.size());
			for( int j=0; j<rel_size; j++ )
			{
				int rel_ele = (*nodes_vec)[i].relative_eles_vec[j];
				for( int k=0; k<7; k++)
				{
					stress[k] += eles_stress[rel_ele][k];
					strain[k] += eles_strain[rel_ele][k];
				}

				bool same_mat = false;
				for(int m=0; m<(int)mat_ids.size(); m++ )
				{
					if( mat_ids[m] == (*eles_vec)[rel_ele].mat )
					{
						same_mat = true;
						break;
					}
				}
				if(!same_mat) mat_ids.push_back((*eles_vec)[rel_ele].mat);
			}
			//计算节点平均应力
			for( int j=0; j<7; j++)
			{
				stress[j] = stress[j]/rel_size;
				strain[j] = strain[j]/rel_size;
			}

			if(node_mark[i] == 1)	//该节点属于内部节点，不在不同材料的界面上
			{
				int mat_id = mat_ids[0];
				MatPro *mat = &(*mats_vec)[mat_id];
				double eff_value = cal_effective_value(strength_criterion_for_matrix,stress,mat);
//				hout << "第" << i <<"个单元：" << "应力分量平均" <<eff_value << "  " << "强度平均" << stress[6]<<endl;
				//经过对比发现先应力平均再算强度和强度值直接平均，得到的结果十分接近，所以不用下面的赋值
//				stress[6] = eff_value;
				if((*nodes_vec)[i].flag==0)	//该节点不在单胞的边界上
				{
					//if( eff_value > max_stress_vl[mat_id] )
					//{
     //                   max_stress_vl[mat_id] = eff_value ;
     //                   max_stress_en[mat_id] = i ;
					//}

					//记录前n大的应力值和节点编号
					node_count++;
					for(int j=0; j<nod_num; j++)
					{
						if( eff_value > large_stress_vl[j] )
						{
							for(int k=nod_num-1; k>j; k--)
							{	
								large_stress_vl[k] = large_stress_vl[k-1];
								large_stress_en[k] = large_stress_en[k-1];
							}
							large_stress_vl[j] = eff_value;
							large_stress_en[j] = i;
							break;
						}
					}
	                max_stress_vl[mat_id] = large_stress_vl[31] ;
                    max_stress_en[mat_id] = large_stress_en[31] ;				
				}
			}
			for( int j=0; j<7; j++)
			{
				nodes_stress[i][j] = stress[j];
				nodes_strain[i][j] = strain[j];
			}
		}
		
		for(int i=0; i<nod_num; i++)
		{
			hout << i << "  " << large_stress_vl[i] << "  " << large_stress_en[i] << endl;
		}
		hout << "参与比较的节点个数node_count=" << node_count<<endl;
        return 1;
};
//---------------------------------------------------------------------------
//计算每个单元在整体坐标系下的应变，返回向量为6个应变分量
//Ox,Oy,Oz,Txy,Tyz,Tzx
void Analyzer::cal_glo_strain( Element* e,vector<double> &gstrain )
{
	gstrain.assign(7,0.0);
	//求出单元节点的个数；
	int node_size =  int (e->nodes_id.size());
//---------------------------------------------------------------
   //求出此单元的所有节点坐标
	vector<Node> elenodes_vec;
	for(int i=0;i<node_size;i++)
		elenodes_vec.push_back((*nodes_vec)[e->nodes_id[i]]);
//---------------------------------------------------------------
	//求出单元内形心的Na1值(单元拉格朗日插值);
	vector<double> Na1_vec_avg(27,0);
	for( int i=0; i<27; i++ ){
		for(int j=0;j<node_size;j++)
			Na1_vec_avg[i] = Na1_vec_avg[i]+(*Na1_vec)[e->nodes_id[j]][i];
		Na1_vec_avg[i] = Na1_vec_avg[i]/node_size;
	}
//----------------------------------------------------------------
//求出形心点;
	double ox=0;
	double oy=0;
	double oz=0;
	for(int i=0;i<node_size;i++){
		ox = ox+elenodes_vec[i].x;
		oy = oy+elenodes_vec[i].y;
		oz = oz+elenodes_vec[i].z;
	}
	ox =ox/node_size;
	oy =oy/node_size;
	oz =oz/node_size;
//求出形心点坐标完毕；
//----------------------------------------------------------------
//生成bdfi;
   vector<vector<double> > bdfi;
Generate_bdfi(bdfi,e,elenodes_vec);
//----------------------------------------------------------------
//求解形心点应变；
        double epsilon[3][3];      //张量形式的应变
		for( int h=0; h<3; h++ )
	{
		for( int k=0; k<3; k++ )
		{
			//第一项；u0产生的应变
			epsilon[h][k] =(	anaSolution->u_1(h,k,ox,oy,oz)+
									anaSolution->u_1(k,h,ox,oy,oz)	)/2.0;

			//第二项 (epsilon=?)
			double sum1 = 0.0;
			for( int a1=0; a1<3; a1++ )
			{
				for( int m=0; m<3; m++ )
				{
					int nn = 3*a1 + m;
					sum1 += Na1_vec_avg[3*nn+h]*
						anaSolution->u_2(m,a1,k,ox,oy,oz);
					sum1 += Na1_vec_avg[3*nn+k]*
						anaSolution->u_2(m,a1,h,ox,oy,oz);
				}
			}
			sum1 = sum1/2.0*LLi_epsilon[3];
      // 	yout << " sum1: " << sum1 ;
			epsilon[h][k] += sum1;

		//第三项 (除去三阶导数)
			sum1 = 0.0;
			for( int a1=0; a1<3; a1++ )
			{
				for( int m=0; m<3; m++ )
				{
					int nn = 3*a1 + m;
					double sum2 = 0.0;            
					for( int l=0; l<node_size; l++ )
					{
						sum2 += bdfi[k][l]*(*Na1_vec)[e->nodes_id[l]][nn*3+h];
						sum2 += bdfi[h][l]*(*Na1_vec)[e->nodes_id[l]][nn*3+k];
					}
					sum2 = sum2/2;
					sum1 += sum2*anaSolution->u_1(m,a1,ox,oy,oz);
				}
			}

			epsilon[h][k] += sum1;

        //第四项 (除去三阶导数，epsilon=?)
			sum1 = 0.0;
			for( int a1=0; a1<3; a1++ )
			{
				for( int a2=0; a2<3; a2++ )
				{
					for( int m=0; m<3; m++ )
					{
						int nn = 9*a1+3*a2 + m;
						double sum2 = 0.0;
						for( int l=0; l<node_size; l++ )
						{
							sum2 += bdfi[k][l]*(*Na1a2_vec)[e->nodes_id[l]][nn*3+h];
							sum2 += bdfi[h][l]*(*Na1a2_vec)[e->nodes_id[l]][nn*3+k];
						}
						sum2 = sum2/2;
						sum1 += sum2*anaSolution->u_2(m,a1,a2,ox,oy,oz);
					}
				}
			}
//			yout << " sum3: " << sum1 ;
//			yout << "\n";
			epsilon[h][k] += sum1*LLi_epsilon[3];
		}
	}
	gstrain[0] = epsilon[0][0];
	gstrain[1] = epsilon[1][1];
	gstrain[2] = epsilon[2][2];
	gstrain[3] = epsilon[0][1]+epsilon[1][0];
	gstrain[4] = epsilon[1][2]+epsilon[2][1];
	gstrain[5] = epsilon[0][2]+epsilon[2][0];

//	yout << " gstrain0: " ;
//	for( int i=0; i<6; i++ )
//	{
//		yout << gstrain[i] << " " ;
//	}
//	yout << "\n";

}
void Analyzer::Generate_bdfi(vector<vector<double> > &bdfi, Element *e,const vector<Node> &elenodes_vec){
     //-------------------------------------------------------------------
	 //四面体单元求解形心应变；
	if(e->type==2&&int(e->nodes_id.size())==4){
      //分配空间；
		vector<double> tem(4,0.0);
		vector<vector<double> > tem2(3,tem);
		//求单元B矩阵；
	//首先计算节点基函数Ni，Nj，Nm，Nl；（其实是计算向量a,b,c）;
	vector<double> a(4),b(4),c(4),d(4);
	int i=0;
	int j=1;
	int m=2;
	int l=3;
	for(int k=0;k<4;k++){
		a[i]= elenodes_vec[j].x*(elenodes_vec[m].y*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].y)+
			elenodes_vec[j].y*(elenodes_vec[m].z*elenodes_vec[l].x-elenodes_vec[m].x*elenodes_vec[l].z)+
			elenodes_vec[j].z*(elenodes_vec[m].x*elenodes_vec[l].y-elenodes_vec[m].y*elenodes_vec[l].x);
		b[i]= -(elenodes_vec[m].y*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].y)
			+(elenodes_vec[j].y*elenodes_vec[l].z-elenodes_vec[j].z*elenodes_vec[l].y)
			-(elenodes_vec[j].y*elenodes_vec[m].z-elenodes_vec[j].z*elenodes_vec[m].y);
		c[i]=  (elenodes_vec[m].x*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].x)
			-(elenodes_vec[j].x*elenodes_vec[l].z-elenodes_vec[j].z*elenodes_vec[l].x)
			+(elenodes_vec[j].x*elenodes_vec[m].z-elenodes_vec[j].z*elenodes_vec[m].x);
		d[i]= -(elenodes_vec[m].x*elenodes_vec[l].y-elenodes_vec[m].y*elenodes_vec[l].x)
			+(elenodes_vec[j].x*elenodes_vec[l].y-elenodes_vec[j].y*elenodes_vec[l].x)
			-(elenodes_vec[j].x*elenodes_vec[m].y-elenodes_vec[j].y*elenodes_vec[m].x);
		int n;
		n=i;
		i=j;
		j=m;
		m=l;
		l=n;
	}
	//--------------------------------------------------------------------
	//计算四面体体积；
	double volume = 1.0/6.0*(a[0]-a[1]+a[2]-a[3]);
	//--------------------------------------------------------------------
	for(int i=0;i<4;i++)
		if(i%2==1){
			b[i]=-b[i];
			c[i]=-c[i];
			d[i]=-d[i];
		}
	for(int i=0;i<4;i++){
	tem2[0][i]=b[i]/(6.0*volume);
	tem2[1][i]=c[i]/(6.0*volume);
	tem2[2][i]=d[i]/(6.0*volume);
	}
	//赋值给bdfi；
	bdfi=tem2;
 }//end 四面体；

else if(e->type==3&&int(e->nodes_id.size())==6){
	//分配变量；
	vector<double> tem(6,0.0);
	vector<vector<double> > tem2(3,tem);
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

		diff[2][0]=-0.5*0.333333333;
		diff[2][1]=-0.5*0.333333333;
		diff[2][2]=-0.5*0.333333333;
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
        //----------------------------------------------------
		 //求出J矩阵的逆矩阵；
			double Jinverse[3][3];
			Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
			Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
           Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;
            
			Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
			Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

			Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
			Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;
            
			Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
            Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
        //-------------------------------------------------------
        //求出N对x,y的偏导；
			for(int i=0;i<3;i++)
				for(int j=0;j<6;j++)
					for(int k=0;k<3;k++)
						tem2[i][j]=tem2[i][j]+Jinverse[i][k]*diff[k][j];
			bdfi = tem2;		
	}//end 三棱柱

}
//---------------------------------------------------------------------------
//计算每个单元在整体坐标系下的应力，返回向量为6个应变分量
void Analyzer::cal_glo_stress( Element* tet, vector<double> *gstrain, vector<double> &gstress ){
        gstress.assign(7,0.0);
  //      //构造弹性系数矩阵
  //      vector<vector<double> > elas_matrix(6, gstress);
        //计算弹性矩阵
 //       GloStiffMatrix gsmatrix( mats_vec );
  //      gsmatrix.gen_gdd_matrix( elas_matrix, tet );
//		(*mats_vec)[tet->mat].elas_matrix
        //计算应力
        for( int i=0; i<6; i++ ){
                for( int j=0; j<6; j++ ){
                        gstress[i] += (*mats_vec)[tet->mat].elas_matrix[i][j] * (*gstrain)[j] ;
                };
        };
};
//---------------------------------------------------------------------------
//根据强度准则计算等效值
//0：最大正应力
//1：最大正应变
//2：最大剪应力
//3：Von-Mises
double Analyzer::cal_effective_value(int cri_num, vector<double> &stress, MatPro *mat){
        double rvalue;
        switch(cri_num){
        case 0:
                rvalue = sc_max_normal_stress(stress);
                break;
        case 1:
                rvalue = sc_max_normal_strain(stress,mat);
                break;
        case 2:
                rvalue = sc_max_shear_stress(stress);
                break;
        case 3:
                rvalue = sc_von_mises(stress);
                break;
        default:
                break;

        };

        return rvalue;
};
//---------------------------------------------------------------------------
//最大正应变断裂准则
double Analyzer::sc_max_normal_strain( vector<double> &stress, MatPro* mat, double the_value, double cvalue ){
        double sigma[3];
        //计算主应力
        cal_3d_prin_stress(sigma,stress);

        double mu12 = mat->Nu12;
        return sigma[0]-mu12*(sigma[1] + sigma[2]);
};
//---------------------------------------------------------------------------
//最大剪应力屈服准则
double Analyzer::sc_max_shear_stress( vector<double> &stress, double the_value, double cvalue ){
        double sigma[3];
        //计算主应力
        cal_3d_prin_stress(sigma,stress);
        return sigma[0]-sigma[2];
};
//---------------------------------------------------------------------------
//最大正应力断裂准则
double Analyzer::sc_max_normal_stress( vector<double> &stress, double the_value, double cvalue ){
        double sigma[3];
        //计算主应力
        cal_3d_prin_stress(sigma,stress);
        //取绝对值最大的主应力
        double sigma1=fabs(sigma[0]);
        double sigma2=fabs(sigma[2]);
        return (sigma1>sigma2) ? sigma1:sigma2;
    //    return fabs(stress[2]);
};
//---------------------------------------------------------------------------
//Von-Mises等效应力屈服准则
double Analyzer::sc_von_mises( vector<double> &stress, double the_value, double cvalue ){
        double rv = (stress[0]-stress[1])*(stress[0]-stress[1])+
                    (stress[1]-stress[2])*(stress[1]-stress[2])+
                    (stress[2]-stress[0])*(stress[2]-stress[0])+
                    6.0*(stress[3]*stress[3]+
                         stress[4]*stress[4]+
                         stress[5]*stress[5]) ;
        return sqrt(0.5*rv);
};
//---------------------------------------------------------------------------
//计算三维应力状态下的主应力
void Analyzer::cal_3d_prin_stress( double sigma[3], vector<double> &stress){
        double stress_c[6];
        for(int i=0; i<6; i++)
                stress_c[i] = stress[i];
        cal_3d_prin_stress(sigma,stress_c);
};
//---------------------------------------------------------------------------
//计算三维应力状态下的主应力
void Analyzer::cal_3d_prin_stress( double sigma[3], double stress[6]){
        double aa = stress[0]+stress[1]+stress[2];
        double bb = stress[1]*stress[2]+stress[2]*stress[0]+stress[0]*stress[1]-
                    stress[3]*stress[3]-stress[4]*stress[4]-stress[5]*stress[5];
        double cc = stress[0]*stress[1]*stress[2]+2.0*stress[3]*stress[4]*stress[5]-
                    stress[0]*stress[4]*stress[4]-stress[1]*stress[5]*stress[5]-
                    stress[2]*stress[3]*stress[3];
        double dd = -36.*aa*bb+108.*cc+8.*aa*aa*aa;
        double ee = 12.*bb*bb*bb-3.*aa*aa*bb*bb-54.*aa*bb*cc+81.*cc*cc+12.*aa*aa*aa*cc;

        complex<double> ee_com(ee);

        complex<double> ff(dd+12.*sqrt(ee_com));
        complex<double> gg = pow(ff,1./3.)/12.;
        complex<double> hh = 1./3.*bb-1./9.*aa*aa;
        complex<double> mm = gg - 1./4.*hh/gg;
        complex<double> sigma1 = 2.*mm+1./3.*aa;
        complex<double> nn(0,sqrt(3.));
        complex<double> pp = (2.*gg-mm)*nn;
        complex<double> sigma2 = -mm+1./3.*aa-pp;
        complex<double> sigma3 = -mm+1./3.*aa+pp;
        if( sigma1.imag() > ZERO || sigma1.imag() > ZERO || sigma1.imag() > ZERO ){
                hout << "注意：计算三维应力状态下的主应力时出错:" << endl;
                hout << "        虚部过大：" << sigma1.imag() << " "
                                             << sigma2.imag() << " "
                                             << sigma3.imag() << endl;

        };
        sigma[0] = sigma1.real();
        sigma[1] = sigma2.real();
        sigma[2] = sigma3.real();

        //按从大到小排序
        if( sigma[1] > sigma[0] ){
                double tmp = sigma[0] ;
                sigma[0] = sigma[1];
                sigma[1] = tmp;
        };
        if( sigma[2] > sigma[0] ){
                double tmp = sigma[0] ;
                sigma[0] = sigma[2];
                sigma[2] = tmp;
        };
        if( sigma[2] > sigma[1] ){
                double tmp = sigma[1] ;
                sigma[1] = sigma[2];
                sigma[2] = tmp;
        };    
};
//---------------------------------------------------------------------------
//输出一个6分量向量
void const  Analyzer::print_vec6(vector<double> *vec, int begin, int length){
        for( int j=begin; j<(begin+length); j++ ){
                hout.width(12);
                hout.precision(6);
                hout << (*vec)[j] << " " ;
        };
};
//---------------------------------------------------------------------------
//输出节点强度等值图
void const  Analyzer::export_node_equ_val()
{
	//输出强度等值图
	ofstream otec("node_equ_val.dat");
	otec << "TITLE = blend_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z, U" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)eles_vec->size(); i++)
	{
		switch ((*eles_vec)[i].mat)
		{
		case 0: count[0]++; break; //基体单元
		case 1: count[1]++; break; //椭球单元
		case 2: count[2]++; break; //边界层单元
		default: hout << " error!! " << endl; break;
		}
	}
	//基体单元
	otec << "ZONE N=" << (int)nodes_vec->size() << ", E=" << count[0] << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for (int i=0; i < (int)nodes_vec->size(); i++)
	{
		otec << (*nodes_vec)[i].x << "  " << (*nodes_vec)[i].y << "  " << (*nodes_vec)[i].z << "  " << nodes_stress[i][6]<< endl;
	}
	otec << endl;
	for (int i=0; i < (int)eles_vec->size(); i++)
	{
		if ((*eles_vec)[i].mat == 0)
		{
			otec	<< (*eles_vec)[i].nodes_id[0]+1 << "  " << (*eles_vec)[i].nodes_id[1]+1 << "  " 
					<< (*eles_vec)[i].nodes_id[2]+1 << "  " << (*eles_vec)[i].nodes_id[3]+1 << endl;
		}
	}

	otec << "ZONE N=" << (int)nodes_vec->size() << ", E=" << count[1] << ", F=FEPOINT, ET=TETRAHEDRON" << endl;
	for (int i=0; i < (int)nodes_vec->size(); i++)
	{
		otec << (*nodes_vec)[i].x << "  " << (*nodes_vec)[i].y << "  " << (*nodes_vec)[i].z << "  " << nodes_stress[i][6]<< endl;
	}
	otec << endl;
	for (int i=0; i < (int)eles_vec->size(); i++)
	{
		if ((*eles_vec)[i].mat == 1)
		{
			otec	<< (*eles_vec)[i].nodes_id[0]+1 << "  " << (*eles_vec)[i].nodes_id[1]+1 << "  " 
					<< (*eles_vec)[i].nodes_id[2]+1 << "  " << (*eles_vec)[i].nodes_id[3]+1 << endl;
		}
	}
	//输出界面层三棱柱网格（按六面体输出）
	if(count[2]>0)
	{
		otec << "ZONE N=" << (int)nodes_vec->size() << ", E=" << count[2] << ", F=FEPOINT, ET=BRICK" << endl;
		for (int i=0; i < (int)nodes_vec->size(); i++)
		{
			otec << (*nodes_vec)[i].x << "  " << (*nodes_vec)[i].y << "  " << (*nodes_vec)[i].z << "  " << nodes_stress[i][6]<< endl;
		}
		otec << endl;
		for (int i=0; i < (int)eles_vec->size(); i++)
		{
			if ((*eles_vec)[i].mat == 2)
			{
				otec	<< (*eles_vec)[i].nodes_id[0]+1 << "  " << (*eles_vec)[i].nodes_id[1]+1 << "  " 
						<< (*eles_vec)[i].nodes_id[2]+1 << "  " << (*eles_vec)[i].nodes_id[0]+1 << "  "
						<< (*eles_vec)[i].nodes_id[3]+1 << "  " << (*eles_vec)[i].nodes_id[4]+1 << "  "
						<< (*eles_vec)[i].nodes_id[5]+1 << "  " << (*eles_vec)[i].nodes_id[3]+1 << endl;
			}
		}
	}
	otec.close();

}
//-----------------------------------------------------
//输出样本数据；
void Analyzer::output_Datafile(int *max_stress_en_T, double *max_stress_vl_T, const string data_file)const
{
	ofstream out(data_file.c_str(),ios::app);
	out <<"%临界载荷；均匀化后最大应力应变；每种材料的承受最大应力单元及应力值" << endl;
	out << anaSolution->load_value() <<"   ";
	out <<stress_vmax[6]<<" "<<strain_vmax[6]<<"   ";
    for( int i=0; i<(int)mats_vec->size(); i++ )
	{
		out << i << "  " <<  max_stress_en_T[i] << "  " << max_stress_vl_T[i] <<"   ";
    }
	out << endl;
	out.close();
}

//---------------------------------------------------------------------------
