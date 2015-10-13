//---------------------------------------------------------------------------
#include "CanZBeamConYMomAnaSol.h"
#include <fstream>
#include <iostream>  
#include "Hns.h"
using namespace hns;
//---------------------------------------------------------------------------
//构造函数
CanZBeamConYMomAnaSol::CanZBeamConYMomAnaSol( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
		cout<<"M="<<M<<endl;
        Iy        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
void CanZBeamConYMomAnaSol::print() const{
        hout << "悬臂梁（z向）纯弯曲：" << endl;
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
void CanZBeamConYMomAnaSol::set_para( BeamSecShape* bss, double M ){
        bSecShape = bss ;
        M_load    = M ;
        Iy        = bSecShape->Ix();
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应变
void CanZBeamConYMomAnaSol::max_strain( vector<double> &strain ){
    /*    if( max_strain_vec.size() == 6 ){
                strain.assign(max_strain_vec.begin(),max_strain_vec.end());
                return ;
        };           */
        strain.assign(6,0.0);
        strain[0] = ux_dx(bSecShape->hh/2.0,0,0);
        strain[1] = uy_dy(bSecShape->hh/2.0,0,0);
        strain[2] = uz_dz(bSecShape->hh/2.0,0,0);
        strain[3] = ux_dy(bSecShape->hh/2.0,0,0)+uy_dx(bSecShape->hh/2.0,0,0);
        strain[4] = uy_dz(bSecShape->hh/2.0,0,0)+uz_dy(bSecShape->hh/2.0,0,0);
        strain[5] = uz_dx(bSecShape->hh/2.0,0,0)+ux_dz(bSecShape->hh/2.0,0,0);
 //@       max_strain_vec.assign(strain.begin(),strain.end());
};
//---------------------------------------------------------------------------
//计算均匀化后理论最大应力
void CanZBeamConYMomAnaSol::max_stress( vector<double> &stress ){
    /*    if( max_stress_vec.size() == 6 ){
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
double CanZBeamConYMomAnaSol::ux( double x, double y, double z ){
        return -M_load*(z*z/homoMat.E33+homoMat.Nu13*x*x/homoMat.E11-homoMat.Nu23*y*y/homoMat.E22)/(2.0*Iy) ;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy( double x, double y, double z ){
        return -homoMat.Nu23*M_load*x*y/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz( double x, double y, double z ){
        return M_load*x*z/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
//位移函数的一阶导数
double CanZBeamConYMomAnaSol::ux_dx( double x, double y, double z ){
        return -homoMat.Nu13*M_load*x/(homoMat.E11*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dy( double x, double y, double z ){
        return homoMat.Nu23*M_load*y/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dz( double x, double y, double z ){
        return -M_load*z/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dx( double x, double y, double z ){
        return -homoMat.Nu23*M_load*y/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dy( double x, double y, double z ){
        return -homoMat.Nu23*M_load*x/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dx( double x, double y, double z ){
        return M_load*z/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dz( double x, double y, double z ){
        return M_load*x/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
//位移函数的二阶导数
double CanZBeamConYMomAnaSol::ux_dx_dx( double x, double y, double z ){
        return -homoMat.Nu13*M_load/(homoMat.E11*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dy_dy( double x, double y, double z ){
        return homoMat.Nu23*M_load/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::ux_dz_dz( double x, double y, double z ){
        return -M_load/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dx_dy( double x, double y, double z ){
        return -homoMat.Nu23*M_load/(homoMat.E22*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dx_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uy_dz_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dx_dx( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dx_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dx_dz( double x, double y, double z ){
        return M_load/(homoMat.E33*Iy);
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dy_dy( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dy_dz( double x, double y, double z ){
        return 0;
};
//---------------------------------------------------------------------------
double CanZBeamConYMomAnaSol::uz_dz_dz( double x, double y, double z ){
        return 0;
};                      
//---------------------------------------------------------------------------
