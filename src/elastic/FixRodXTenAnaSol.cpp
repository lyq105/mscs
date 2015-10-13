//---------------------------------------------------------------------------
#include "FixRodXTenAnaSol.h"
#include "GloStiffMatrix.h"
//---------------------------------------------------------------------------
#include <fstream>
#include <iostream>
//---------------------------------------------------------------------------
//���캯��
FixRodXTenAnaSol::FixRodXTenAnaSol( double area, double T ){
        Av = area;
        Tv = T;
        pv = T/area;
};
//---------------------------------------------------------------------------
//������Ȼ�����ϵ������
/*
void FixRodXTenAnaSol::import_homo_para( char* filename ){
        ifstream import_s;
        import_s.open( filename );
        if( !import_s ) hout << "       Can not open this file: " << filename << endl;
        import_s >> E_modo ;
        import_s >> Mu_pora ;
     //   hout << " homo_E: " << homo_E << " homo_Mu: " << homo_Mu << endl;
};      */
//---------------------------------------------------------------------------
//�����Ϣ
void FixRodXTenAnaSol::print() const {
        hout << "����ľ������죺" << endl;
        hout << "       �غɣ�" << Tv << endl;
        hout << "       ����������" << Av << endl;
};
//---------------------------------------------------------------------------
//���ò�������������أ����ȣ����Ծ�
void FixRodXTenAnaSol::set_para( double area, double T ){
        Av = area;
        Tv = T;         
        pv = T/area;
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodXTenAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };     */
        strain.assign(6,0.0);
        strain[0] = ux_dx(0,0,0);
        strain[1] = uy_dy(0,0,0);
        strain[2] = uz_dz(0,0,0);
        strain[3] = ux_dy(0,0,0)+uy_dx(0,0,0);
        strain[4] = uy_dz(0,0,0)+uz_dy(0,0,0);
        strain[5] = uz_dx(0,0,0)+ux_dz(0,0,0);
        max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void FixRodXTenAnaSol::max_stress( vector<double> &stress ){
     /*   if( max_stress_vec.size() == 6 ){
                stress.assign(max_stress_vec.begin(),max_stress_vec.end());
                return ;
        };
        if( max_strain_vec.size() < 6 )       */
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
double FixRodXTenAnaSol::ux( double x, double y, double z ){
        return pv*x/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy( double x, double y, double z ){
        return -homoMat.Nu12*pv*y/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz( double x, double y, double z ){
        return -homoMat.Nu13*pv*z/homoMat.E11 ;
};
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵���
double FixRodXTenAnaSol::ux_dx( double x, double y, double z ){
        return pv/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dy( double x, double y, double z ){
        return -homoMat.Nu12*pv/homoMat.E11 ;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dz( double x, double y, double z ){
        return -homoMat.Nu13*pv/homoMat.E11 ;
};
//---------------------------------------------------------------------------
//λ�ƺ����Ķ��׵���
double FixRodXTenAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodXTenAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
