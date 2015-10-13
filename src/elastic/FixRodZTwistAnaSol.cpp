//--------------------------------------------------------------------------- 
#include "FixRodZTwistAnaSol.h"
#include <fstream>
#include <iostream>     
#include "Hns.h"
using namespace hns;
#define PI 3.141592654
//---------------------------------------------------------------------------
//构造函数
FixRodZTwistAnaSol::FixRodZTwistAnaSol( double r, double T ){
        radius = r ;
        T_load    = T ;
};
//---------------------------------------------------------------------------
void FixRodZTwistAnaSol::print() const{
        hout << "圆截面柱形杆（z向）扭转：" << endl;
        hout << "       截面半径：" ;
        hout << radius << endl;
        hout << "       载荷大小：" << T_load << endl;
};
//---------------------------------------------------------------------------
//设置参数（梁）：弯矩，长度，惯性矩
void FixRodZTwistAnaSol::set_para( double r, double T ){
        radius = r ;
        T_load    = T ;
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应变
void FixRodZTwistAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };                       */
        strain.assign(6,0.0);
        strain[0] = ux_dx(radius,0,0);
        strain[1] = uy_dy(radius,0,0);
        strain[2] = uz_dz(radius,0,0);
        strain[3] = ux_dy(radius,0,0)+uy_dx(radius,0,0);
        strain[4] = uy_dz(radius,0,0)+uz_dy(radius,0,0);
        strain[5] = uz_dx(radius,0,0)+ux_dz(radius,0,0);
        max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应力
void FixRodZTwistAnaSol::max_stress( vector<double> &stress ){
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
//位移函数
double FixRodZTwistAnaSol::ux( double x, double y, double z ){
        return -T_load*y*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy( double x, double y, double z ){
        return T_load*x*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz( double x, double y, double z ){
        return -T_load*x*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13-1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
//位移函数的一阶导数
double FixRodZTwistAnaSol::ux_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dy( double x, double y, double z ){
        return -T_load*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dz( double x, double y, double z ){
        return -T_load*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dx( double x, double y, double z ){
        return T_load*z/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dz( double x, double y, double z ){
        return T_load*x/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dx( double x, double y, double z ){
        return -T_load*y/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13-1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dy( double x, double y, double z ){
        return -T_load*x/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13-1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
//位移函数的二阶导数
double FixRodZTwistAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dy_dz( double x, double y, double z ){
        return -T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dx_dz( double x, double y, double z ){
        return T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13+1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dx_dy( double x, double y, double z ){
        return -T_load/(PI*radius*radius*radius*radius)*(1.0/homoMat.G13-1.0/homoMat.G23) ;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double FixRodZTwistAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
