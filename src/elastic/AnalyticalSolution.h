/****************************************************************************
给定模型的解析解
model_id:
        0:      悬臂梁端部受集中弯矩（纯弯梁）
****************************************************************************/
//---------------------------------------------------------------------------

#ifndef AnalyticalSolutionH
#define AnalyticalSolutionH
#include <vector>         
#include "MatPro.h"
//---------------------------------------------------------------------------
class AnalyticalSolution{
public:
        AnalyticalSolution(){};                      //构造函数
    //    virtual void set_para( double M, double L, double Iz );
                                                    //设置参数（梁）
                                                    
        virtual double ux( double, double, double )=0;        //位移函数
        virtual double uy( double, double, double )=0;
        virtual double uz( double, double, double )=0;
        
        virtual double u_1( int, int, double, double, double );
                                                              //位移函数的一阶导数，前两个参数确定具体导数
        virtual double ux_dx( double, double, double )=0;     //位移函数的一阶导数
        virtual double ux_dy( double, double, double )=0;
        virtual double ux_dz( double, double, double )=0;
   //     virtual double uy_dx( double, double, double )=0;
        virtual double uy_dy( double, double, double )=0;
        virtual double uy_dz( double, double, double )=0;
   //     virtual double uz_dx( double, double, double )=0;
   //     virtual double uz_dy( double, double, double )=0;
        virtual double uz_dz( double, double, double )=0;

        virtual double uy_dx( double xx, double yy, double zz ){
                return ux_dy( xx, yy, zz ) ;
        };                                                 
        virtual double uz_dx( double xx, double yy, double zz ){
                return ux_dz( xx, yy, zz );
        };                                                            
        virtual double uz_dy( double xx, double yy, double zz ){
                return uy_dz( xx, yy, zz );
        };

        virtual double u_2( int, int, int, double, double, double );
                                                              //位移函数的一阶导数，前3个参数确定具体导数
        virtual double ux_dx_dx( double, double, double )=0;   //位移函数的二阶导数
        virtual double ux_dx_dy( double, double, double )=0;
        virtual double ux_dx_dz( double, double, double )=0;
        virtual double ux_dy_dy( double, double, double )=0;
        virtual double ux_dy_dz( double, double, double )=0;
        virtual double ux_dz_dz( double, double, double )=0;
        virtual double uy_dx_dx( double, double, double )=0;
        virtual double uy_dx_dy( double, double, double )=0;
        virtual double uy_dx_dz( double, double, double )=0;
        virtual double uy_dy_dy( double, double, double )=0;
        virtual double uy_dy_dz( double, double, double )=0;
        virtual double uy_dz_dz( double, double, double )=0;
        virtual double uz_dx_dx( double, double, double )=0;
        virtual double uz_dx_dy( double, double, double )=0;
        virtual double uz_dx_dz( double, double, double )=0;
        virtual double uz_dy_dy( double, double, double )=0;
        virtual double uz_dy_dz( double, double, double )=0;
        virtual double uz_dz_dz( double, double, double )=0;

        double ux_dy_dx( double xx, double yy, double zz ){
                return ux_dx_dy( xx, yy, zz );
        };
        double ux_dz_dx( double xx, double yy, double zz ){
                return ux_dx_dz( xx, yy, zz );
        };
        double ux_dz_dy( double xx, double yy, double zz ){
                return ux_dy_dz( xx, yy, zz );
        };
        double uy_dy_dx( double xx, double yy, double zz ){
                return uy_dx_dy( xx, yy, zz );
        };
        double uy_dz_dx( double xx, double yy, double zz ){
                return uy_dx_dz( xx, yy, zz );
        };
        double uy_dz_dy( double xx, double yy, double zz ){
                return uy_dy_dz( xx, yy, zz );
        };
        double uz_dy_dx( double xx, double yy, double zz ){
                return uz_dx_dy( xx, yy, zz );
        };
        double uz_dz_dx( double xx, double yy, double zz ){
                return uz_dx_dz( xx, yy, zz );
        };
        double uz_dz_dy( double xx, double yy, double zz ){
                return uz_dy_dz( xx, yy, zz );
        };

    //    virtual void import_homo_para( char* )=0 ;                 //导入均匀化弹性系数数据
        void set_homo_mat( MatPro& ) ;                     //设置均匀化材料属性

        virtual void max_strain( vector<double> &strain )=0;       //计算均匀化后理论最大应变
        virtual void max_stress( vector<double> &stress )=0;       //计算均匀化后理论最大应力

        virtual double load_value() = 0;
                                                         //返回载荷大小
        virtual void set_load_value(double) = 0;
                                                         //返回载荷大小
        virtual void print() const = 0;
                                                         //打印
//protected:                                        
        void gen_dd_matrix();                            //生成弹性矩阵D

        double dd[6][6];                                 //弹性矩阵
        MatPro homoMat;                                //均匀化后的材料属性
        
    //    int model_id;     //模型号
   /*
        double len    ;    //梁的长度
        double E_modo ;    //弹性模量
        double Mu_pora;    //泊松比
        double M_load ;    //集中弯矩载荷大小                  
        double Iz     ;    //梁的z轴弯曲惯性矩
                   */
};
//---------------------------------------------------------------------------
#endif
