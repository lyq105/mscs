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
//���캯��
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
//�趨������
void Analyzer::set_ana_solution( AnalyticalSolution* as ){
        anaSolution = as ;
};
//---------------------------------------------------------------------------
//����Na1ָ��
void Analyzer::set_Na1_vec(vector<vector<double> >* na1){
        Na1_vec = na1;
};  //---------------------------------------------------------------------------
//����Na1a2ָ��
void Analyzer::set_Na1a2_vec(vector<vector<double> >* na1a2){
        Na1a2_vec = na1a2;
};
//---------------------------------------------------------------------------
//���þ��Ȼ�����ָ��
void Analyzer::set_homo_mat(MatPro* mp){
        homoMat = mp;
};
//---------------------------------------------------------------------------
//�������������ļ��ļ���
void Analyzer::set_string_datafile(string datafile)
{
	data_file = datafile;
}
//---------------------------------------------------------------------------
//ǿ�ȷ��� (�����ٽ��غɣ������ڿ�������)
int Analyzer::strength_analysis_par(int mod){
        bool debuging_sa = true;
        if( debuging_sa ) hout << "____________________strength_analysis_par()_______________________" <<endl;

   //     reset_cell();
   //     gsmatrix=GloStiffMatrix(nodes_vec,eles_vec,mats_vec);

        anaSolution->print();
                                  
        //��ʼ�غ�
        double T0 = 0.0;
        double T1 = anaSolution->load_value();
        
        //��Ӧ���غ�T0
        int    *max_stress_en_T0 = new int[mats_vec->size()];
        double *max_stress_vl_T0 = new double[mats_vec->size()];
        for( int i=0; i<(int)mats_vec->size(); i++ ){
                max_stress_en_T0[i] = 0  ;
                max_stress_vl_T0[i] = 0.0;
        };

        //��Ӧ���غ�T1
        int    *max_stress_en_T1 = new    int[mats_vec->size()];
        double *max_stress_vl_T1 = new double[mats_vec->size()];

        int cn = 0;
        while( cn < 30 ){
                if( debuging_sa ) hout << "~~~~~~~~~~~~�����ٽ��غ�~~~~~~~~~~~~~~~~~" << "\n";
                if( debuging_sa ) hout << " cn: " << cn++ << "\n";
                if( debuging_sa ) hout << " �غɣ�" << anaSolution->load_value() << "\n";

                if( debuging_sa ){
                        //������Ȼ�����������Ӧ������
                        anaSolution->max_strain( strain_vmax );
                        hout << "���Ȼ������Ӧ�䣺\n    " ;
                        for( int i=0; i<6; i++ ){
                                hout << strain_vmax[i] << " " ;
                        };
                        hout << "\n";
                        //������Ȼ�����������Ӧ������
                        anaSolution->max_stress( stress_vmax );
                        hout << "���Ȼ������Ӧ����\n    " ;
                        for( int i=0; i<6; i++ ){
                                hout << stress_vmax[i] << " " ;
                        };
                        hout << "\n";
                };

                stress_analysis_par( max_stress_en_T1, max_stress_vl_T1);
                
                if( mod == 0 ) break;
                //�������Ƿ��Ѿ�����
                int is_yield = -1 ;  //-1: û��������0���ٽ磬1���Ѿ�����
                double critical_load=1.0E200;
//                for( int i=0; i<(int)mats_vec->size(); i++ ){
                for( int i=1; i<=1; i++ ){
                  //      int mat_id = (*eles_vec)[max_Lstress_nn_T1[i]].materialId ;
                        MatPro *mat = &(*mats_vec)[i];
                        if( debuging_sa ) hout << " i: " << i << " mat->type(): " << mat->type_val << "\n";
                        if( mat->type_val == 0){
                                double Ystress = mat->strength;
                                //�������Ӧ��
                                if( debuging_sa ) hout << "���Ӧ����" << max_stress_vl_T1[i] << "\n";

                                if( max_stress_vl_T1[i] < 0 )continue;
                                
                                //��������Ļ�������ЧӦ��������ȷ���ٽ��غ�
                               double Tc = T1-(T1-T0)*(fabs(max_stress_vl_T1[i])-Ystress)
                                                 /(fabs(max_stress_vl_T1[i])-fabs(max_stress_vl_T0[i]));
                                if( Tc < critical_load ) critical_load = Tc;
                                
                                double ds = (fabs(max_stress_vl_T1[i])-Ystress)/Ystress ;
                                if( fabs(ds) < 1e-8 ){
                                        //�ٽ�״̬
                                        is_yield = 0;      
                                }
                                else if( ds > 0 ){
                                        //�Ѿ�����
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

        hout << "�ٽ��غɣ�" << anaSolution->load_value() << "\n";
        anaSolution->print();

        //������Ȼ�����������Ӧ������
        anaSolution->max_strain( strain_vmax );
        double effective_strain = cal_effective_value(strength_criterion_for_matrix,strain_vmax,homoMat);
        strain_vmax.push_back(effective_strain);
        hout << "���Ȼ������Ӧ�䣺\n    " ;
        print_vec6( &strain_vmax, 0, 7 );
        hout << "\n";
        
        //������Ȼ�����������Ӧ������
        anaSolution->max_stress( stress_vmax );
        double effective_stress = cal_effective_value(strength_criterion_for_matrix,stress_vmax,homoMat);
        stress_vmax.push_back(effective_stress);
        hout << "���Ȼ������Ӧ����\n    " ;
        print_vec6( &stress_vmax, 0, 7 );
        hout << "\n";

        for( int i=0; i<(int)mats_vec->size(); i++ ){
                hout << i << " �Ų��ϣ�\n";
                hout << "    ��Ԫ " << max_stress_en_T1[i] << " ��ЧӦ�����" << "\n";
                hout << "        ��Ӧ��Ϊ��" << max_stress_vl_T1[i] << "\n";
                hout << "        Ӧ��������" ;
                print_vec6( &eles_stress[max_stress_en_T1[i]],0,6 );
                hout << "\n";
                hout << "        Ӧ�������" ;
                print_vec6( &eles_strain[max_stress_en_T1[i]],0,6 );
                hout << "\n";
        };

		//����������ݣ�
		output_Datafile(max_stress_en_T1, max_stress_vl_T1, data_file);

		//����ڵ�ǿ�ȵ�ֵͼ
//		export_node_equ_val();

        delete []max_stress_en_T0;
        delete []max_stress_vl_T0;
        delete []max_stress_en_T1;
        delete []max_stress_vl_T1;

        return 1;   
};
//---------------------------------------------------------------------------
//Ӧ����������������غ��µ�Ӧ���ֲ������Ӧ���������ڿ������ϣ�
int Analyzer::stress_analysis_par(int *max_stress_en, double *max_stress_vl){
        bool debuging = false ;

        //����ÿ���ڵ��λ��
    //    cal_displacement();

        vector<double> strain(7,0.0) ;	//��7��Ӧ�������ȱ
        vector<double> stress(7,0.0) ;	//��7��Ӧ��������Ÿ���ǿ��׼���������Ӧ��ֵ

        eles_strain.assign(eles_vec->size(),strain);	
        eles_stress.assign(eles_vec->size(),stress);	

        //��ʼ��
        for( int i=0; i<(int)mats_vec->size(); i++ ){
                max_stress_en[i] = 0;
                max_stress_vl[i] = -1.0e200;
        };

		 for( int i=0; i<(int)eles_vec->size(); i++){
            if( debuging ) hout << "---------------------------stress_analysis_par()-------------------------" << endl;
            if( debuging ) hout << "��Ԫ( " << i << " ): " << "\n";
            //���Ϻ�                                     
            int mat_id = (*eles_vec)[i].mat;
            MatPro *mat = &(*mats_vec)[mat_id];
            if( debuging ) hout << "        mat id: " << mat_id << "\n";
  //          if( mat->type_val == 0 ){
			if( true ){
                    //������������ϵ�µ�Ӧ�䡢Ӧ��
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
                    hout << "       Ӧ�䣺" ;
                    for( int j=0; j<6; j++ ){
                            hout << strain[j] << " " ;
                    };
                    hout << "\n";
                    hout << "       Ӧ����" ;
                    for( int j=0; j<6; j++ ){
                            hout << stress[j] << " " ;
                    };
                    hout << "\n";
            };
            if( debuging ) hout << "        Ӧ��Ӧ�������ϡ�\n";
            eles_stress[i] = stress ;
            eles_strain[i] = strain ;
            if( debuging ) hout << "        Ӧ��Ӧ�䱣����ϡ�\n";

        };
        
		//ƽ��������ڵ��Ӧ��Ӧ��
		//��Žڵ��Ӧ���������ϵĽڵ���ܰ�������Ӧ��ֵ����ʽΪ�����Ϻţ�Ӧ��������6���������Ϻţ�... //�ڿ��������п���û�У�����ά�������У�
		stress.resize(7,0.0);
		strain.resize(7,0.0);

		nodes_stress.assign(nodes_vec->size(),stress);
		nodes_strain.assign(nodes_vec->size(),strain);

		//���ڼ�¼�ж��ٸ������Ƚ�Σ�յ�
		int node_count = 0;
		//���ڼ�¼ǰnum����Ľڵ��Ӧ��ֵ�ͱ��
		int nod_num = 1500;
		vector<double> large_stress_vl(nod_num, -1.0e200);
		vector<int> large_stress_en(nod_num, 0);

		//ȷ��ÿ���ڵ��λ�ã��ڲ������Ϸֽ����������Ȼ���ڲ�������ֽ��������һ��ڵ�
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
			//����ڵ�ƽ��Ӧ��
			for( int j=0; j<7; j++)
			{
				stress[j] = stress[j]/rel_size;
				strain[j] = strain[j]/rel_size;
			}

			if(node_mark[i] == 1)	//�ýڵ������ڲ��ڵ㣬���ڲ�ͬ���ϵĽ�����
			{
				int mat_id = mat_ids[0];
				MatPro *mat = &(*mats_vec)[mat_id];
				double eff_value = cal_effective_value(strength_criterion_for_matrix,stress,mat);
//				hout << "��" << i <<"����Ԫ��" << "Ӧ������ƽ��" <<eff_value << "  " << "ǿ��ƽ��" << stress[6]<<endl;
				//�����Աȷ�����Ӧ��ƽ������ǿ�Ⱥ�ǿ��ֱֵ��ƽ�����õ��Ľ��ʮ�ֽӽ������Բ�������ĸ�ֵ
//				stress[6] = eff_value;
				if((*nodes_vec)[i].flag==0)	//�ýڵ㲻�ڵ����ı߽���
				{
					//if( eff_value > max_stress_vl[mat_id] )
					//{
     //                   max_stress_vl[mat_id] = eff_value ;
     //                   max_stress_en[mat_id] = i ;
					//}

					//��¼ǰn���Ӧ��ֵ�ͽڵ���
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
		hout << "����ȽϵĽڵ����node_count=" << node_count<<endl;
        return 1;
};
//---------------------------------------------------------------------------
//����ÿ����Ԫ����������ϵ�µ�Ӧ�䣬��������Ϊ6��Ӧ�����
//Ox,Oy,Oz,Txy,Tyz,Tzx
void Analyzer::cal_glo_strain( Element* e,vector<double> &gstrain )
{
	gstrain.assign(7,0.0);
	//�����Ԫ�ڵ�ĸ�����
	int node_size =  int (e->nodes_id.size());
//---------------------------------------------------------------
   //����˵�Ԫ�����нڵ�����
	vector<Node> elenodes_vec;
	for(int i=0;i<node_size;i++)
		elenodes_vec.push_back((*nodes_vec)[e->nodes_id[i]]);
//---------------------------------------------------------------
	//�����Ԫ�����ĵ�Na1ֵ(��Ԫ�������ղ�ֵ);
	vector<double> Na1_vec_avg(27,0);
	for( int i=0; i<27; i++ ){
		for(int j=0;j<node_size;j++)
			Na1_vec_avg[i] = Na1_vec_avg[i]+(*Na1_vec)[e->nodes_id[j]][i];
		Na1_vec_avg[i] = Na1_vec_avg[i]/node_size;
	}
//----------------------------------------------------------------
//������ĵ�;
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
//������ĵ�������ϣ�
//----------------------------------------------------------------
//����bdfi;
   vector<vector<double> > bdfi;
Generate_bdfi(bdfi,e,elenodes_vec);
//----------------------------------------------------------------
//������ĵ�Ӧ�䣻
        double epsilon[3][3];      //������ʽ��Ӧ��
		for( int h=0; h<3; h++ )
	{
		for( int k=0; k<3; k++ )
		{
			//��һ�u0������Ӧ��
			epsilon[h][k] =(	anaSolution->u_1(h,k,ox,oy,oz)+
									anaSolution->u_1(k,h,ox,oy,oz)	)/2.0;

			//�ڶ��� (epsilon=?)
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

		//������ (��ȥ���׵���)
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

        //������ (��ȥ���׵�����epsilon=?)
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
	 //�����嵥Ԫ�������Ӧ�䣻
	if(e->type==2&&int(e->nodes_id.size())==4){
      //����ռ䣻
		vector<double> tem(4,0.0);
		vector<vector<double> > tem2(3,tem);
		//��ԪB����
	//���ȼ���ڵ������Ni��Nj��Nm��Nl������ʵ�Ǽ�������a,b,c��;
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
	//���������������
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
	//��ֵ��bdfi��
	bdfi=tem2;
 }//end �����壻

else if(e->type==3&&int(e->nodes_id.size())==6){
	//���������
	vector<double> tem(6,0.0);
	vector<vector<double> > tem2(3,tem);
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

		diff[2][0]=-0.5*0.333333333;
		diff[2][1]=-0.5*0.333333333;
		diff[2][2]=-0.5*0.333333333;
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
        //----------------------------------------------------
		 //���J����������
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
        //���N��x,y��ƫ����
			for(int i=0;i<3;i++)
				for(int j=0;j<6;j++)
					for(int k=0;k<3;k++)
						tem2[i][j]=tem2[i][j]+Jinverse[i][k]*diff[k][j];
			bdfi = tem2;		
	}//end ������

}
//---------------------------------------------------------------------------
//����ÿ����Ԫ����������ϵ�µ�Ӧ������������Ϊ6��Ӧ�����
void Analyzer::cal_glo_stress( Element* tet, vector<double> *gstrain, vector<double> &gstress ){
        gstress.assign(7,0.0);
  //      //���쵯��ϵ������
  //      vector<vector<double> > elas_matrix(6, gstress);
        //���㵯�Ծ���
 //       GloStiffMatrix gsmatrix( mats_vec );
  //      gsmatrix.gen_gdd_matrix( elas_matrix, tet );
//		(*mats_vec)[tet->mat].elas_matrix
        //����Ӧ��
        for( int i=0; i<6; i++ ){
                for( int j=0; j<6; j++ ){
                        gstress[i] += (*mats_vec)[tet->mat].elas_matrix[i][j] * (*gstrain)[j] ;
                };
        };
};
//---------------------------------------------------------------------------
//����ǿ��׼������Чֵ
//0�������Ӧ��
//1�������Ӧ��
//2������Ӧ��
//3��Von-Mises
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
//�����Ӧ�����׼��
double Analyzer::sc_max_normal_strain( vector<double> &stress, MatPro* mat, double the_value, double cvalue ){
        double sigma[3];
        //������Ӧ��
        cal_3d_prin_stress(sigma,stress);

        double mu12 = mat->Nu12;
        return sigma[0]-mu12*(sigma[1] + sigma[2]);
};
//---------------------------------------------------------------------------
//����Ӧ������׼��
double Analyzer::sc_max_shear_stress( vector<double> &stress, double the_value, double cvalue ){
        double sigma[3];
        //������Ӧ��
        cal_3d_prin_stress(sigma,stress);
        return sigma[0]-sigma[2];
};
//---------------------------------------------------------------------------
//�����Ӧ������׼��
double Analyzer::sc_max_normal_stress( vector<double> &stress, double the_value, double cvalue ){
        double sigma[3];
        //������Ӧ��
        cal_3d_prin_stress(sigma,stress);
        //ȡ����ֵ������Ӧ��
        double sigma1=fabs(sigma[0]);
        double sigma2=fabs(sigma[2]);
        return (sigma1>sigma2) ? sigma1:sigma2;
    //    return fabs(stress[2]);
};
//---------------------------------------------------------------------------
//Von-Mises��ЧӦ������׼��
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
//������άӦ��״̬�µ���Ӧ��
void Analyzer::cal_3d_prin_stress( double sigma[3], vector<double> &stress){
        double stress_c[6];
        for(int i=0; i<6; i++)
                stress_c[i] = stress[i];
        cal_3d_prin_stress(sigma,stress_c);
};
//---------------------------------------------------------------------------
//������άӦ��״̬�µ���Ӧ��
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
                hout << "ע�⣺������άӦ��״̬�µ���Ӧ��ʱ����:" << endl;
                hout << "        �鲿����" << sigma1.imag() << " "
                                             << sigma2.imag() << " "
                                             << sigma3.imag() << endl;

        };
        sigma[0] = sigma1.real();
        sigma[1] = sigma2.real();
        sigma[2] = sigma3.real();

        //���Ӵ�С����
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
//���һ��6��������
void const  Analyzer::print_vec6(vector<double> *vec, int begin, int length){
        for( int j=begin; j<(begin+length); j++ ){
                hout.width(12);
                hout.precision(6);
                hout << (*vec)[j] << " " ;
        };
};
//---------------------------------------------------------------------------
//����ڵ�ǿ�ȵ�ֵͼ
void const  Analyzer::export_node_equ_val()
{
	//���ǿ�ȵ�ֵͼ
	ofstream otec("node_equ_val.dat");
	otec << "TITLE = blend_mesh_in_ellipse" << endl;
	otec << "VARIABLES = X, Y, Z, U" << endl;

	int count[3]={0,0,0};
	for (int i=0; i < (int)eles_vec->size(); i++)
	{
		switch ((*eles_vec)[i].mat)
		{
		case 0: count[0]++; break; //���嵥Ԫ
		case 1: count[1]++; break; //����Ԫ
		case 2: count[2]++; break; //�߽�㵥Ԫ
		default: hout << " error!! " << endl; break;
		}
	}
	//���嵥Ԫ
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
	//�����������������񣨰������������
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
//����������ݣ�
void Analyzer::output_Datafile(int *max_stress_en_T, double *max_stress_vl_T, const string data_file)const
{
	ofstream out(data_file.c_str(),ios::app);
	out <<"%�ٽ��غɣ����Ȼ������Ӧ��Ӧ�䣻ÿ�ֲ��ϵĳ������Ӧ����Ԫ��Ӧ��ֵ" << endl;
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
