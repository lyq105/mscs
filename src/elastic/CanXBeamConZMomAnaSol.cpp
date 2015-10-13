//---------------------------------------------------------------------------
#include "CanXBeamConZMomAnaSol.h"
#include <fstream>
#include <iostream>
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
//构造函数
CanXBeamConZMomAnaSol::CanXBeamConZMomAnaSol( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
        Iz        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
void CanXBeamConZMomAnaSol::print() const{
        hout << "悬臂梁（x向）纯弯曲：" << endl;
        hout << "       截面形状：" ;
        if( bSecShape->type() == 0 )
                hout << "圆形，半径：" << bSecShape->rr << endl;
        else if( bSecShape->type() == 1 )
                hout << "矩形，长、宽：" << bSecShape->hh << " "
                                         << bSecShape->bb << endl;
        hout << "       载荷大小：" << M_load << endl;
};
//---------------------------------------------------------------------------
//设置参数（梁）：弯矩，长度，惯性矩
void CanXBeamConZMomAnaSol::set_para( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
        Iz        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应变
void CanXBeamConZMomAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };   */
        strain.assign(6,0.0);
        strain[0] = ux_dx(0,bSecShape->hh/2.0,0);
        strain[1] = uy_dy(0,bSecShape->hh/2.0,0);
        strain[2] = uz_dz(0,bSecShape->hh/2.0,0);
        strain[3] = ux_dy(0,bSecShape->hh/2.0,0)+uy_dx(0,bSecShape->hh/2.0,0);
        strain[4] = uy_dz(0,bSecShape->hh/2.0,0)+uz_dy(0,bSecShape->hh/2.0,0);
        strain[5] = uz_dx(0,bSecShape->hh/2.0,0)+ux_dz(0,bSecShape->hh/2.0,0);
        max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应力
void CanXBeamConZMomAnaSol::max_stress( vector<double> &stress ){
     /*   if( max_stress_vec.size() == 6 ){
                stress.assign(max_stress_vec.begin(),max_stress_vec.end());
                return ;
        };
        if( max_strain_vec.size() < 6 )    */
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
double CanXBeamConZMomAnaSol::ux( double x, double y, double z ){
        return M_load*x*y/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy( double x, double y, double z ){  
        return -M_load*(x*x+homoMat.Nu12*y*y-homoMat.Nu13*z*z)/(2.0*homoMat.E11*Iz) ;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz( double x, double y, double z ){
        return -homoMat.Nu13*M_load*y*z/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
//位移函数的一阶导数
double CanXBeamConZMomAnaSol::ux_dx( double x, double y, double z ){
    //    hout << " ux_dx called. y: " << y << endl;
        return M_load*y/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dy( double x, double y, double z ){
        return M_load*x/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dx( double x, double y, double z ){
        return -M_load*x/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dy( double x, double y, double z ){
        return -homoMat.Nu12*M_load*y/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dz( double x, double y, double z ){
        return homoMat.Nu13*M_load*z/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dy( double x, double y, double z ){
        return -homoMat.Nu13*M_load*z/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dz( double x, double y, double z ){
        return -homoMat.Nu13*M_load*y/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
//位移函数的二阶导数
double CanXBeamConZMomAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dx_dy( double x, double y, double z ){
        return M_load/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dx_dx( double x, double y, double z ){
        return -M_load/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dy_dy( double x, double y, double z ){
        return -homoMat.Nu12*M_load/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uy_dz_dz( double x, double y, double z ){
        return homoMat.Nu13*M_load/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dy_dz( double x, double y, double z ){
        return -homoMat.Nu13*M_load/(homoMat.E11*Iz);
};
//---------------------------------------------------------------------------
double CanXBeamConZMomAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
