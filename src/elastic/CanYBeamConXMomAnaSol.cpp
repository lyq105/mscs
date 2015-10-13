//---------------------------------------------------------------------------
#include "CanYBeamConXMomAnaSol.h"
#include <fstream>
#include <iostream>
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
//构造函数
CanYBeamConXMomAnaSol::CanYBeamConXMomAnaSol( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
        Ix        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
void CanYBeamConXMomAnaSol::print() const{
        hout << "悬臂梁（y向）纯弯曲：" << endl;
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
void CanYBeamConXMomAnaSol::set_para( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
        Ix        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应变
void CanYBeamConXMomAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };    */
        strain.assign(6,0.0);
        strain[0] = ux_dx(0,0,bSecShape->hh/2.0);
        strain[1] = uy_dy(0,0,bSecShape->hh/2.0);
        strain[2] = uz_dz(0,0,bSecShape->hh/2.0);
        strain[3] = ux_dy(0,0,bSecShape->hh/2.0)+uy_dx(0,0,bSecShape->hh/2.0);
        strain[4] = uy_dz(0,0,bSecShape->hh/2.0)+uz_dy(0,0,bSecShape->hh/2.0);
        strain[5] = uz_dx(0,0,bSecShape->hh/2.0)+ux_dz(0,0,bSecShape->hh/2.0);
        max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应力
void CanYBeamConXMomAnaSol::max_stress( vector<double> &stress ){
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
//位移函数
double CanYBeamConXMomAnaSol::ux( double x, double y, double z ){
        return -homoMat.Nu12*M_load*x*z/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy( double x, double y, double z ){
        return M_load*y*z/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz( double x, double y, double z ){
        return -M_load*(y*y/homoMat.E22+homoMat.Nu23*z*z/homoMat.E22-homoMat.Nu12*x*x/homoMat.E11)/(2.0*Ix) ;
};
//---------------------------------------------------------------------------
//位移函数的一阶导数
double CanYBeamConXMomAnaSol::ux_dx( double x, double y, double z ){
        return -homoMat.Nu12*M_load*z/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dz( double x, double y, double z ){
        return -homoMat.Nu12*M_load*x/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dy( double x, double y, double z ){
        return M_load*z/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dz( double x, double y, double z ){
        return M_load*y/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dx( double x, double y, double z ){
        return homoMat.Nu12*M_load*x/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dy( double x, double y, double z ){
        return -M_load*y/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dz( double x, double y, double z ){
        return -homoMat.Nu23*M_load*z/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
//位移函数的二阶导数
double CanYBeamConXMomAnaSol::ux_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dx_dz( double x, double y, double z ){
        return -homoMat.Nu12*M_load/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::ux_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dy_dz( double x, double y, double z ){
        return M_load/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dx_dx( double x, double y, double z ){
        return homoMat.Nu12*M_load/(homoMat.E11*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dy_dy( double x, double y, double z ){
        return -M_load/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanYBeamConXMomAnaSol::uz_dz_dz( double x, double y, double z ){
        return -homoMat.Nu23*M_load/(homoMat.E22*Ix);
};
//---------------------------------------------------------------------------
