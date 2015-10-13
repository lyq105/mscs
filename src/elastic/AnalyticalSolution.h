/****************************************************************************
����ģ�͵Ľ�����
model_id:
        0:      �������˲��ܼ�����أ���������
****************************************************************************/
//---------------------------------------------------------------------------

#ifndef AnalyticalSolutionH
#define AnalyticalSolutionH
#include <vector>         
#include "MatPro.h"
//---------------------------------------------------------------------------
class AnalyticalSolution{
public:
        AnalyticalSolution(){};                      //���캯��
    //    virtual void set_para( double M, double L, double Iz );
                                                    //���ò���������
                                                    
        virtual double ux( double, double, double )=0;        //λ�ƺ���
        virtual double uy( double, double, double )=0;
        virtual double uz( double, double, double )=0;
        
        virtual double u_1( int, int, double, double, double );
                                                              //λ�ƺ�����һ�׵�����ǰ��������ȷ�����嵼��
        virtual double ux_dx( double, double, double )=0;     //λ�ƺ�����һ�׵���
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
                                                              //λ�ƺ�����һ�׵�����ǰ3������ȷ�����嵼��
        virtual double ux_dx_dx( double, double, double )=0;   //λ�ƺ����Ķ��׵���
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

    //    virtual void import_homo_para( char* )=0 ;                 //������Ȼ�����ϵ������
        void set_homo_mat( MatPro& ) ;                     //���þ��Ȼ���������

        virtual void max_strain( vector<double> &strain )=0;       //������Ȼ����������Ӧ��
        virtual void max_stress( vector<double> &stress )=0;       //������Ȼ����������Ӧ��

        virtual double load_value() = 0;
                                                         //�����غɴ�С
        virtual void set_load_value(double) = 0;
                                                         //�����غɴ�С
        virtual void print() const = 0;
                                                         //��ӡ
//protected:                                        
        void gen_dd_matrix();                            //���ɵ��Ծ���D

        double dd[6][6];                                 //���Ծ���
        MatPro homoMat;                                //���Ȼ���Ĳ�������
        
    //    int model_id;     //ģ�ͺ�
   /*
        double len    ;    //���ĳ���
        double E_modo ;    //����ģ��
        double Mu_pora;    //���ɱ�
        double M_load ;    //��������غɴ�С                  
        double Iz     ;    //����z���������Ծ�
                   */
};
//---------------------------------------------------------------------------
#endif
