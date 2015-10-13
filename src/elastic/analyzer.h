/****************************************************************************
������������ݼ������Na1 Na1a2����ڵ�λ�ƣ���һ������Ӧ�䡢Ӧ��������ʧЧ
�оݣ�����ȷ�����Լ���
ע��Ӧ������˳��xx yy zz xy yz zx
****************************************************************************/
//---------------------------------------------------------------------------  
#ifndef analyzerH
#define analyzerH
//---------------------------------------------------------------------------
#include <vector>
#include "Fem.h"
#include "MatPro.h"          
#include "AnalyticalSolution.h"
#include "Geometry.h"
#include "GloStiffMatrix.h"
//---------------------------------------------------------------------------

class Analyzer {
public:
        vector<Node>                 *nodes_vec   ;      //�ڵ��
        vector<Element>          *eles_vec    ;      //��Ԫ��
        //vector<int>                  *INodes_vec  ;      //�ڲ��ڵ�� 
        //vector<int>                  *IEles_vec   ;      //�ڲ���Ԫ��
        //vector<vector<double> >       nodes_dis   ;      //�ڵ�λ��
        vector<vector<double> >       eles_strain ;      //��ԪӦ��
        vector<vector<double> >       eles_stress ;      //��ԪӦ��
        vector<vector<double> >       nodes_strain;      //�ڵ�Ӧ�䣨��ص�Ԫ��ƽ����
        vector<vector<double> >       nodes_stress;      //�ڵ�Ӧ������ص�Ԫ��ƽ����
        vector<MatPro>               *mats_vec    ;      //���Ͽ�
        vector<vector<double> >      *Na1_vec     ;      //Na1
        vector<vector<double> >      *Na1a2_vec   ;      //Na1a2             
        vector<double> strain_vmax, stress_vmax ; //���Ȼ���Ӧ��Ӧ������

        Analyzer( vector<Node> *nodes, vector<Element> *eles, vector<MatPro> *mats );
        void set_cell_size( double Li_epsilon[4] ){
                LLi_epsilon[0] = Li_epsilon[0];
                LLi_epsilon[1] = Li_epsilon[1];
                LLi_epsilon[2] = Li_epsilon[2];
                LLi_epsilon[3] = Li_epsilon[3];
        };
                                                         //�趨�����ߴ�
        //void set_cell_origin( double xx, double yy, double zz ){
        //        ox = xx ;
        //        oy = yy ;
        //        oz = zz ;
        //};

        //void set_title( string title ){
        //        atitle = title ;
        //};
        //                                                 //�趨����λ�ã��������λ�ã�
        //                                                 
        //void set_interior_nodes( vector<int> *INodes );  //�趨�ڲ��ڵ㼯
        void set_interior_eles( vector<int> *IEles );    //�趨�ڲ���Ԫ��
        void set_ana_solution( AnalyticalSolution* );    //�趨������
        void set_Na1_vec(vector<vector<double> >*);      //����Na1ָ��
        void set_Na1a2_vec(vector<vector<double> >*);    //����Na1a2ָ��
        void set_homo_mat(MatPro*);                      //���þ��Ȼ�����ָ��
		void set_string_datafile(string);	//�������������ļ��ļ���
        //void import_Na1( char* ) ;                       //����Na1����
        //void import_Na1a2( char* ) ;                     //����Na1a2����
        //void import_homo_para( char* ) ;                 //������Ȼ�����ϵ������
        //void cal_displacement();                         //����ÿ���ڵ��λ��
        void cal_glo_strain( Element* e,vector<double> &gstrain );
                                                         //����ÿ����Ԫ��������ϵ�µ�Ӧ��
		void Generate_bdfi(vector<vector<double> > &bdfi,Element *e,const vector<Node> &elenodes_vec);
        //void cal_glo_strain_ddis( Tetrahedron* tet, vector<double> &strain );
        //                                                 //����ÿ����Ԫ��������ϵ�µ�Ӧ��(��λ����)
        //void cal_glo_strain( Tetrahedron* tet, vector<double> &strain, int mod );
        //                                                 //����ÿ����Ԫ��������ϵ�µ�Ӧ��
        //void cal_loc_strain( Tetrahedron* tet, vector<double> &strain );
        //                                                 //����ÿ����Ԫ��Ԫ����ϵ�µ�Ӧ��
        //void cal_loc_strain( Tetrahedron* tet, vector<double> &strain, int mod );
        //                                                 //����ÿ����Ԫ��Ԫ����ϵ�µ�Ӧ��
        void cal_glo_stress( Element* tet, vector<double> *strain, vector<double> &stress );
                                                           //����ÿ����Ԫ��������ϵ�µ�Ӧ������
        //void cal_loc_stress( Tetrahedron* tet, vector<double> *strain, vector<double> &stress );
        //                                                 //����ÿ����Ԫ��Ԫ����ϵ�µ�Ӧ������
        //double cal_von_stress( vector<double> *stress );
        //                                                 //����ÿ����Ԫ��Von��ЧӦ��
        //void cal_2d_prin_stress( double &sigma_1, double &sigma_2, double stress[3]);
        //                                                 //�����άӦ��״̬�µ���Ӧ��
        //void cal_2d_prin_strain( double &epsilon_1, double &epsilon_2, double strain[3]);
        //                                                 //�����άӦ��״̬�µ���Ӧ��
        void cal_3d_prin_stress( double sigma[3], double stress[6]);
                                                         //������άӦ��״̬�µ���Ӧ��
        void cal_3d_prin_stress( double sigma[3], vector<double> &stress);
                                                         //������άӦ��״̬�µ���Ӧ��

        //double cal_BiDiFi_volume( Tetrahedron *tet, double bi[4], double di[4], double fi[4] );
        //                                          //���������嵥Ԫ��bi, di, fi, ����ֵΪ����������

        //int strength_analysis(int mod = 1) ;
        //                                                 //ǿ�ȷ��� (�����ٽ��غ�)
        int strength_analysis_par(int mod = 1) ;
														   //ǿ�ȷ��� (�����ٽ��غɣ������ڿ�������)
        //int strength_analysis_par_node(int mod = 1) ;
        //                                                 //ǿ�ȷ��� (�����ٽ��غɣ������ڿ�������)(���սڵ�ƽ��Ӧ����
        //int stress_analysis(int *max_Lstrain_en, double *max_Lstrain_vl,
        //                    int *max_Tstrain_en, double *max_Tstrain_vl,
        //                    int *max_Lstress_en, double *max_Lstress_vl,
        //                    int *max_Tstress_en, double *max_Tstress_vl) ;
        //                                            //Ӧ����������������غ��µ�Ӧ���ֲ������Ӧ����
        int stress_analysis_par(int *max_stress_en, double *max_stress_vl) ;
													  //Ӧ����������������غ��µ�Ӧ���ֲ������Ӧ���������ڿ������ϣ�
        //int stress_analysis_par_node(int *max_stress_nn, double *max_stress_vl) ;
        //                                            //Ӧ����������������غ��µ�Ӧ���ֲ������Ӧ���������ڿ������ϣ�(����ڵ�ƽ��Ӧ����
        double sc_von_mises( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //Von-Mises��ЧӦ������׼��
                                                    //the_value:��ǰ��Ԫ��ֵ��cvalue:�ο��ٽ�ֵ
        double sc_max_normal_stress( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //�����Ӧ������׼��
        double sc_max_normal_strain( vector<double> &stress, MatPro *mat, double the_value=0.0, double cvalue=0.0 );
                                                    //�����Ӧ�����׼��
        double sc_max_shear_stress( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //����Ӧ������׼��

        //void export_def_iges(char* file="cell_def_fem.igs", int mod = 0 ) ;
        //                                            //�������κ�Ľڵ㣬��Ԫ��export a iges(*.igs) file��

        //void reset_cell();
        //                                
        //void write_result(int par_scale_num, int sample_num) const;
        //                                                                 
        //void write_result(char* file, char* title) const;
        //                                            //�Ѽ�����д���ļ���ͨ��ֻд��ǿ��ֵ������ͳ����ͼ��
        void set_criterions(int mat_cri, int rein_cri){
                strength_criterion_for_matrix = mat_cri;
                strength_criterion_for_reinforcement = rein_cri;
        };
                                                    //����ǿ��׼��mat_cri:����ǿ��׼��; rein_cri:��ǿ��ǿ��׼��
                                                    //0: �����Ӧ��
                                                    //1: �����Ӧ��
                                                    //2: ����Ӧ��
                                                    //3: Von-Mises

protected:
        //double ox,oy,oz;                                 //�������λ��
        double LLi_epsilon[4];                              //�����ĳߴ�
        //double homo_E, homo_Mu ;                         //���Ȼ�����ģ���Ͳ��ɱ�
        MatPro *homoMat;                                 //���Ȼ�����
        //double cal_tet_volume( Tetrahedron* tet );       //��������������

private:
//        char *Na1_file, *Na1a2_file, *homo_para_file ;	//Na1��Na1a2�Լ����Ȼ����Գ����ļ�
          AnalyticalSolution *anaSolution ;							//������
//        GloStiffMatrix gsmatrix;

        //����ʽ���Ӧ��Ӧ�����
//        void const print_node_ss(vector<double> *vec, int mat_id);
        //���һ��6��������
        void const print_vec6(vector<double> *vec, int begin, int length);

		//����������ݣ�
		void output_Datafile(int *max_stress_en_T, double *max_stress_vl_T, const string data_file)const;

		//����ڵ�ǿ�ȵ�ֵͼ
		void const  export_node_equ_val();

        //����ǿ��׼������Чֵ
        double cal_effective_value(int cri_num, vector<double> &stress, MatPro *mat);

         //int node_or_ele_strength;

//        string atitle;

        //���õ�ǿ��׼��
        //0�������Ӧ��
        //1�������Ӧ��
        //2������Ӧ��
        //3��Von-Mises
        int strength_criterion_for_matrix;
        int strength_criterion_for_reinforcement;
		string data_file;

        //����һ�����飬��־ÿ����Ԫ�����ʣ��߽絥Ԫ���Ǳ߽絥Ԫ
 //       vector<int> tets_sign;
};
//---------------------------------------------------------------------------
#endif



 
