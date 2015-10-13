/****************************************************************************
结果分析，根据计算出的Na1 Na1a2计算节点位移，进一步计算应变、应力，引入失效
判据，即可确定弹性极限
注：应力分量顺序：xx yy zz xy yz zx
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
        vector<Node>                 *nodes_vec   ;      //节点库
        vector<Element>          *eles_vec    ;      //单元库
        //vector<int>                  *INodes_vec  ;      //内部节点库 
        //vector<int>                  *IEles_vec   ;      //内部单元库
        //vector<vector<double> >       nodes_dis   ;      //节点位移
        vector<vector<double> >       eles_strain ;      //单元应变
        vector<vector<double> >       eles_stress ;      //单元应力
        vector<vector<double> >       nodes_strain;      //节点应变（相关单元的平均）
        vector<vector<double> >       nodes_stress;      //节点应力（相关单元的平均）
        vector<MatPro>               *mats_vec    ;      //材料库
        vector<vector<double> >      *Na1_vec     ;      //Na1
        vector<vector<double> >      *Na1a2_vec   ;      //Na1a2             
        vector<double> strain_vmax, stress_vmax ; //均匀化后应变应力分量

        Analyzer( vector<Node> *nodes, vector<Element> *eles, vector<MatPro> *mats );
        void set_cell_size( double Li_epsilon[4] ){
                LLi_epsilon[0] = Li_epsilon[0];
                LLi_epsilon[1] = Li_epsilon[1];
                LLi_epsilon[2] = Li_epsilon[2];
                LLi_epsilon[3] = Li_epsilon[3];
        };
                                                         //设定单胞尺寸
        //void set_cell_origin( double xx, double yy, double zz ){
        //        ox = xx ;
        //        oy = yy ;
        //        oz = zz ;
        //};

        //void set_title( string title ){
        //        atitle = title ;
        //};
        //                                                 //设定单胞位置（宏观坐标位置）
        //                                                 
        //void set_interior_nodes( vector<int> *INodes );  //设定内部节点集
        void set_interior_eles( vector<int> *IEles );    //设定内部单元集
        void set_ana_solution( AnalyticalSolution* );    //设定解析解
        void set_Na1_vec(vector<vector<double> >*);      //设置Na1指针
        void set_Na1a2_vec(vector<vector<double> >*);    //设置Na1a2指针
        void set_homo_mat(MatPro*);                      //设置均匀化材料指针
		void set_string_datafile(string);	//设置样本数据文件文件名
        //void import_Na1( char* ) ;                       //导入Na1数据
        //void import_Na1a2( char* ) ;                     //导入Na1a2数据
        //void import_homo_para( char* ) ;                 //导入均匀化弹性系数数据
        //void cal_displacement();                         //计算每个节点的位移
        void cal_glo_strain( Element* e,vector<double> &gstrain );
                                                         //计算每个单元整体坐标系下的应变
		void Generate_bdfi(vector<vector<double> > &bdfi,Element *e,const vector<Node> &elenodes_vec);
        //void cal_glo_strain_ddis( Tetrahedron* tet, vector<double> &strain );
        //                                                 //计算每个单元整体坐标系下的应变(对位移求导)
        //void cal_glo_strain( Tetrahedron* tet, vector<double> &strain, int mod );
        //                                                 //计算每个单元整体坐标系下的应变
        //void cal_loc_strain( Tetrahedron* tet, vector<double> &strain );
        //                                                 //计算每个单元单元坐标系下的应变
        //void cal_loc_strain( Tetrahedron* tet, vector<double> &strain, int mod );
        //                                                 //计算每个单元单元坐标系下的应变
        void cal_glo_stress( Element* tet, vector<double> *strain, vector<double> &stress );
                                                           //计算每个单元整体坐标系下的应力向量
        //void cal_loc_stress( Tetrahedron* tet, vector<double> *strain, vector<double> &stress );
        //                                                 //计算每个单元单元坐标系下的应力向量
        //double cal_von_stress( vector<double> *stress );
        //                                                 //计算每个单元的Von等效应力
        //void cal_2d_prin_stress( double &sigma_1, double &sigma_2, double stress[3]);
        //                                                 //计算二维应力状态下的主应力
        //void cal_2d_prin_strain( double &epsilon_1, double &epsilon_2, double strain[3]);
        //                                                 //计算二维应力状态下的主应变
        void cal_3d_prin_stress( double sigma[3], double stress[6]);
                                                         //计算三维应力状态下的主应力
        void cal_3d_prin_stress( double sigma[3], vector<double> &stress);
                                                         //计算三维应力状态下的主应力

        //double cal_BiDiFi_volume( Tetrahedron *tet, double bi[4], double di[4], double fi[4] );
        //                                          //计算四面体单元的bi, di, fi, 返回值为四面体的体积

        //int strength_analysis(int mod = 1) ;
        //                                                 //强度分析 (计算临界载荷)
        int strength_analysis_par(int mod = 1) ;
														   //强度分析 (计算临界载荷，适用于颗粒材料)
        //int strength_analysis_par_node(int mod = 1) ;
        //                                                 //强度分析 (计算临界载荷，适用于颗粒材料)(按照节电平均应力）
        //int stress_analysis(int *max_Lstrain_en, double *max_Lstrain_vl,
        //                    int *max_Tstrain_en, double *max_Tstrain_vl,
        //                    int *max_Lstress_en, double *max_Lstress_vl,
        //                    int *max_Tstress_en, double *max_Tstress_vl) ;
        //                                            //应力分析（计算给定载荷下的应力分布，最大应力）
        int stress_analysis_par(int *max_stress_en, double *max_stress_vl) ;
													  //应力分析（计算给定载荷下的应力分布，最大应力，适用于颗粒材料）
        //int stress_analysis_par_node(int *max_stress_nn, double *max_stress_vl) ;
        //                                            //应力分析（计算给定载荷下的应力分布，最大应力，适用于颗粒材料）(计算节点平均应力）
        double sc_von_mises( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //Von-Mises等效应力屈服准则
                                                    //the_value:当前单元的值，cvalue:参考临界值
        double sc_max_normal_stress( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //最大正应力断裂准则
        double sc_max_normal_strain( vector<double> &stress, MatPro *mat, double the_value=0.0, double cvalue=0.0 );
                                                    //最大正应变断裂准则
        double sc_max_shear_stress( vector<double> &stress, double the_value=0.0, double cvalue=0.0 );
                                                    //最大剪应力屈服准则

        //void export_def_iges(char* file="cell_def_fem.igs", int mod = 0 ) ;
        //                                            //导出变形后的节点，单元（export a iges(*.igs) file）

        //void reset_cell();
        //                                
        //void write_result(int par_scale_num, int sample_num) const;
        //                                                                 
        //void write_result(char* file, char* title) const;
        //                                            //把计算结果写入文件（通常只写入强度值，便于统计作图）
        void set_criterions(int mat_cri, int rein_cri){
                strength_criterion_for_matrix = mat_cri;
                strength_criterion_for_reinforcement = rein_cri;
        };
                                                    //设置强度准则，mat_cri:基体强度准则; rein_cri:增强项强度准则
                                                    //0: 最大主应力
                                                    //1: 最大主应变
                                                    //2: 最大剪应力
                                                    //3: Von-Mises

protected:
        //double ox,oy,oz;                                 //宏观坐标位置
        double LLi_epsilon[4];                              //单胞的尺寸
        //double homo_E, homo_Mu ;                         //均匀化弹性模量和泊松比
        MatPro *homoMat;                                 //均匀化材料
        //double cal_tet_volume( Tetrahedron* tet );       //计算四面体的体积

private:
//        char *Na1_file, *Na1a2_file, *homo_para_file ;	//Na1、Na1a2以及均匀化弹性常数文件
          AnalyticalSolution *anaSolution ;							//解析解
//        GloStiffMatrix gsmatrix;

        //按格式输出应力应变分量
//        void const print_node_ss(vector<double> *vec, int mat_id);
        //输出一个6分量向量
        void const print_vec6(vector<double> *vec, int begin, int length);

		//输出样本数据；
		void output_Datafile(int *max_stress_en_T, double *max_stress_vl_T, const string data_file)const;

		//输出节点强度等值图
		void const  export_node_equ_val();

        //根据强度准则计算等效值
        double cal_effective_value(int cri_num, vector<double> &stress, MatPro *mat);

         //int node_or_ele_strength;

//        string atitle;

        //采用的强度准则
        //0：最大正应力
        //1：最大正应变
        //2：最大剪应力
        //3：Von-Mises
        int strength_criterion_for_matrix;
        int strength_criterion_for_reinforcement;
		string data_file;

        //定义一个数组，标志每个单元的性质：边界单元、非边界单元
 //       vector<int> tets_sign;
};
//---------------------------------------------------------------------------
#endif



 
