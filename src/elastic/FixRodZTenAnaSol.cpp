//---------------------------------------------------------------------------
#include "FixRodZTenAnaSol.h"
#include "GloStiffMatrix.h"
//---------------------------------------------------------------------------
#include <fstream>
#include <iostream>
//---------------------------------------------------------------------------
//���캯��
FixRodZTenAnaSol::FixRodZTenAnaSol( double area, double T ){
        Av = area;
        Tv = T;
        pv = T/area;
};
//---------------------------------------------------------------------------
//�����Ϣ
void FixRodZTenAnaSol::print() const {
        hout << "����ľ������죺" << endl;
        hout << "       �غɣ�" << Tv << endl;
        hout << "       ����������" << Av << endl;
};
//---------------------------------------------------------------------------
//���ò�������������أ����ȣ����Ծ�
void FixRodZTenAnaSol::set_para( double area, double T ){
        Av = area;
        Tv = T;         
        pv = T/area;
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodZTenAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };                    */
        strain.assign(6,0.0);
        strain[0] = ux_dx(0,0,0);
        strain[1] = uy_dy(0,0,0);
        strain[2] = uz_dz(0,0,0);
        strain[3] = ux_dy(0,0,0)+uy_dx(0,0,0);
        strain[4] = uy_dz(0,0,0)+uz_dy(0,0,0);
        strain[5] = uz_dx(0,0,0)+ux_dz(0,0,0);
 //@       max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodZTenAnaSol::max_stress( vector<double> &stress ){
    /*    if( max_stress_vec.size() == 6 ){
                stress.assign(max_stress_vec.begin(),max_stress_vec.end());
                return ;
        };
        if( max_strain_vec.size() < 6 )   */
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
double FixRodZTenAnaSol::ux( double x, double y, double z ){
        return -homoMat.Nu13*pv*x/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy( double x, double y, double z ){
        return -homoMat.Nu23*pv*y/homoMat.E22 ;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz( double x, double y, double z ){
        return pv*z/homoMat.E33 ;
};
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵���
double FixRodZTenAnaSol::ux_dx( double x, double y, double z ){
        return -homoMat.Nu13*pv/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dy( double x, double y, double z ){
        return -homoMat.Nu23*pv/homoMat.E22 ;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dz( double x, double y, double z ){
        return pv/homoMat.E33 ;
};
//---------------------------------------------------------------------------
//λ�ƺ����Ķ��׵���
double FixRodZTenAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTenAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
