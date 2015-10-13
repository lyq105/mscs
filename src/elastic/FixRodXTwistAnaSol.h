//---------------------------------------------------------------------------
//圆截面柱形杆的扭转位移场
#ifndef FixRodXTwistAnaSolH
#define FixRodXTwistAnaSolH
//---------------------------------------------------------------------------
#include <vector>
#include "AnalyticalSolution.h"
#include "BeamSecShape.h"
#include "MatPro.h"
//---------------------------------------------------------------------------
class FixRodXTwistAnaSol : public AnalyticalSolution{
public:
        FixRodXTwistAnaSol(){};                    //默认构造函数
        FixRodXTwistAnaSol( double r, double T ) ;
                                                    //构造函数
        void set_para( double r, double T ) ;
                                                    //设置参数（梁）
                                                    
        double ux( double, double, double );        //位移函数
        double uy( double, double, double );
        double uz( double, double, double );

        double ux_dx( double, double, double );     //位移函数的一阶导数
        double ux_dy( double, double, double );
        double ux_dz( double, double, double );
        double uy_dx( double, double, double );
        double uy_dy( double, double, double );
        double uy_dz( double, double, double );
        double uz_dx( double, double, double );
        double uz_dy( double, double, double );
        double uz_dz( double, double, double );

        double ux_dx_dx( double, double, double );   //位移函数的二阶导数
        double ux_dx_dy( double, double, double );
        double ux_dx_dz( double, double, double );
        double ux_dy_dy( double, double, double );
        double ux_dy_dz( double, double, double );
        double ux_dz_dz( double, double, double );
        double uy_dx_dx( double, double, double );
        double uy_dx_dy( double, double, double );
        double uy_dx_dz( double, double, double );
        double uy_dy_dy( double, double, double );
        double uy_dy_dz( double, double, double );
        double uy_dz_dz( double, double, double );
        double uz_dx_dx( double, double, double );
        double uz_dx_dy( double, double, double );
        double uz_dx_dz( double, double, double );
        double uz_dy_dy( double, double, double );
        double uz_dy_dz( double, double, double );
        double uz_dz_dz( double, double, double );

   //     void import_homo_para( char* ) ;                 //导入均匀化弹性系数数据
   //     void set_homo_mat( MatPro* ) ;                   //设置均匀化材料属性

        void max_strain( vector<double> &strain );       //计算均匀化后理论最大应变
        void max_stress( vector<double> &stress );       //计算均匀化后理论最大应力
        
        double load_value(){ return T_load; } ;
                                                         //返回载荷大小
        void   set_load_value(double T){
                T_load = T;
        } ;
                                                         //返回载荷大小
        
        void print() const;
protected:
    //    void gen_dd_matrix();                            //生成弹性矩阵D

   //     double dd[6][6];                                 //弹性矩阵  

    //    int model_id;     //模型号

   //     double E_modo ;                                 //弹性模量
   //     double Mu_pora;                                 //泊松比
        double T_load ;                                 //集中弯矩载荷大小
        double radius;                                  //杆的截面半径

        vector<double> max_strain_vec;                  //最大应变向量
        vector<double> max_stress_vec;                  //最大应力向量

};
//---------------------------------------------------------------------------
#endif
