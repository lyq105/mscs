//--------------------------------------------------------------------------- 
#include "FixRodXTwistAnaSol.h"
#include <fstream>
#include <iostream>    
#include "Hns.h"
using namespace hns;
#define PI 3.141592654
//---------------------------------------------------------------------------
//���캯��
FixRodXTwistAnaSol::FixRodXTwistAnaSol( double r, double T ){
        radius = r ;
        T_load    = T ;
};
//---------------------------------------------------------------------------
void FixRodXTwistAnaSol::print() const{
        hout << "Բ�������θˣ�z��Ťת��" << endl;
        hout << "       ����뾶��" ;
        hout << radius << endl;
        hout << "       �غɴ�С��" << T_load << endl;
};
//---------------------------------------------------------------------------
//���ò�������������أ����ȣ����Ծ�
void FixRodXTwistAnaSol::set_para( double r, double T ){
        radius = r ;
        T_load    = T ;
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodXTwistAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };                       */
        strain.assign(6,0.0);
        strain[0] = ux_dx(0,radius,0);
        strain[1] = uy_dy(0,radius,0);
        strain[2] = uz_dz(0,radius,0);
        strain[3] = ux_dy(0,radius,0)+uy_dx(0,radius,0);
        strain[4] = uy_dz(0,radius,0)+uz_dy(0,radius,0);
        strain[5] = uz_dx(0,radius,0)+ux_dz(0,radius,0);
        max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodXTwistAnaSol::max_stress( vector<double> &stress ){
     /*   if( max_stress_vec.size() == 6 ){
                stress.assign(max_stress_vec.begin(),max_stress_vec.end());
                return ;
        };
        if( max_strain_vec.size() < 6 )     */
        max_strain( max_strain_vec );
        stress.assign(6,0.0);
        gen_dd_matrix();
        for( int i=0; i<6; i++ ){
                for( int j=0; j<6; j++ ){
                        stress[i] += dd[i][j]*max_strain_vec[j];
                };
        };
        max_stress_vec.assign(stress.begin(),stress.end());
};
//---------------------------------------------------------------------------
//λ�ƺ���
double FixRodXTwistAnaSol::ux( double x, double y, double z ){
        return -T_load*y*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12-1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy( double x, double y, double z ){ 
        return -T_load*x*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz( double x, double y, double z ){
        return T_load*x*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵���
double FixRodXTwistAnaSol::ux_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dy( double x, double y, double z ){
        return -T_load*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12-1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dz( double x, double y, double z ){
        return -T_load*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12-1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dx( double x, double y, double z ){
        return -T_load*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dz( double x, double y, double z ){
        return -T_load*x/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dx( double x, double y, double z ){
        return T_load*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dy( double x, double y, double z ){
        return T_load*x/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
//λ�ƺ����Ķ��׵���
double FixRodXTwistAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dy_dz( double x, double y, double z ){
        return -T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12-1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dx_dz( double x, double y, double z ){
        return -T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dx_dy( double x, double y, double z ){
        return T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G12+1.0/homoMat.G13) ;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTwistAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
