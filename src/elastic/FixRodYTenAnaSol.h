//---------------------------------------------------------------------------
//�����ܼ�������������λ�Ƴ�
#ifndef FixRodYTenAnaSolH
#define FixRodYTenAnaSolH
//---------------------------------------------------------------------------
#include <vector>
#include "AnalyticalSolution.h"
#include "BeamSecShape.h"
#include "MatPro.h"
//---------------------------------------------------------------------------
class FixRodYTenAnaSol : public AnalyticalSolution{
public:
        FixRodYTenAnaSol(){};                        //Ĭ�Ϲ��캯��
        FixRodYTenAnaSol( double area, double T ) ;
                                                    //���캯��
        void set_para( double area, double T ) ;
                                                    //���ò�������������غɴ�С��
                                                    
        double ux( double, double, double );        //λ�ƺ���
        double uy( double, double, double );
        double uz( double, double, double );

        double ux_dx( double, double, double );     //λ�ƺ�����һ�׵���
        double ux_dy( double, double, double );
        double ux_dz( double, double, double );
        double uy_dx( double, double, double );
        double uy_dy( double, double, double );
        double uy_dz( double, double, double );
        double uz_dx( double, double, double );
        double uz_dy( double, double, double );
        double uz_dz( double, double, double );

        double ux_dx_dx( double, double, double );   //λ�ƺ����Ķ��׵���
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

   //     void import_homo_para( char* ) ;                 //������Ȼ�����ϵ������
   //     void set_homo_mat( MatPro& ) ;                   //���þ��Ȼ���������

        void max_strain( vector<double> &strain );       //������Ȼ����������Ӧ��
        void max_stress( vector<double> &stress );       //������Ȼ����������Ӧ��

        double load_value(){ return Tv; } ;
                                                         //�����غɴ�С
        void   set_load_value(double T){
                Tv = T;
                pv = T / Av;
        } ;
                                                         //�����غɴ�С
        void print() const ;
protected:
   //     void gen_dd_matrix();                            //���ɵ��Ծ���D

   //     double dd[6][6];                                 //���Ծ���


    //    int model_id;     //ģ�ͺ�

   //     MatPro *homoMat;                                //���Ȼ���Ĳ�������
    /*    double len    ;                                 //���ĳ���
        double M_load ;                                 //��������غɴ�С
        BeamSecShape* bSecShape;                        //���Ľ�����״
        double Iz     ;                                 //����z���������Ծ�
                                                                          */
        double Av;                                      //��������
        double Tv;                                      //��������С
        double pv;                                      //Tv/Av
        double E_modo ;                                 //����ģ��
        double Mu_pora;                                 //���ɱ�

        vector<double> max_strain_vec;                  //���Ӧ������
        vector<double> max_stress_vec;                  //���Ӧ������

}; 
//---------------------------------------------------------------------------
#endif