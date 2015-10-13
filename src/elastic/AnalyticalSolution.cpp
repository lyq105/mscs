//---------------------------------------------------------------------------
#include "AnalyticalSolution.h"
#include "GloStiffMatrix.h"
#include <fstream>
#include <iostream>
//---------------------------------------------------------------------------
//���þ��Ȼ���������
void AnalyticalSolution::set_homo_mat( MatPro& hm ){
  /*      homoMat.set_type(2);
        if( hm.type() == 0 ){
                double E = hm.E();
                double Mu = hm.Mu();
                double G = E/(2.0*(1.0+Mu));
                homoMat.set_ela_para( E, E, E, Mu, Mu, Mu, G, G, G );
        }
        else if( hm.type() == 1 ){
                double E11 = hm.E1();
                double E22 = hm.E1();
                double E33 = hm.E2();
                double Mu12 = hm.Mu1();
                double Mu23 = E11/E33*hm.Mu2();
                double Mu13 = Mu23;
                double G12  = E11/(2.0*(1.0+Mu12));
                double G23  = hm.G2();
                double G13  = hm.G2();
                homoMat.set_ela_para( E11, E22, E33, Mu12, Mu23, Mu13, G12, G23, G13 );
        }
        else if( hm.type() == 2 ){
                homoMat = hm;
        };
        */
        homoMat = hm;
};
//---------------------------------------------------------------------------
//����20070104����
//���ɵ��Ծ���D
void AnalyticalSolution::gen_dd_matrix()
{
        homoMat.Generate_elas_matrix();
        for( int i=0; i<6; i++ ){
                for( int j=0; j<6; j++ ){
                        dd[i][j] = homoMat.elas_matrix[i][j];
                }
        }
}
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵�����ǰ��������ȷ�����嵼��
double AnalyticalSolution::u_1( int u, int x, double xx, double yy, double zz ){
        if( u == 0 ){
                if( x == 0 ){
                        return ux_dx( xx, yy, zz );
                }
                else if( x == 1 ){
                        return ux_dy( xx, yy, zz );
                }
                else if( x == 2 ){
                        return ux_dz( xx, yy, zz );
                };
        }
        else if( u == 1 ){
                if( x == 0 ){
                        return uy_dx( xx, yy, zz );
                }
                else if( x == 1 ){
                        return uy_dy( xx, yy, zz );
                }
                else if( x == 2 ){
                        return uy_dz( xx, yy, zz );
                };
        }
        else if( u == 2 ){
                if( x == 0 ){
                        return uz_dx( xx, yy, zz );
                }
                else if( x == 1 ){
                        return uz_dy( xx, yy, zz );
                }
                else if( x == 2 ){
                        return uz_dz( xx, yy, zz );
                };
        };
        return 0.0;
};
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵�����ǰ��������ȷ�����嵼��
double AnalyticalSolution::u_2( int u, int x1, int x2, double xx, double yy, double zz ){
        if( u == 0 ){
                if( x1 == 0 ){
                        if( x2 == 0 ){
                                return ux_dx_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return ux_dx_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return ux_dx_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 1 ){
                        if( x2 == 0 ){
                                return ux_dy_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return ux_dy_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return ux_dy_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 2 ){
                        if( x2 == 0 ){
                                return ux_dz_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return ux_dz_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return ux_dz_dz( xx, yy, zz );
                        };
                };
        }
        else if( u == 1 ){
                if( x1 == 0 ){
                        if( x2 == 0 ){
                                return uy_dx_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uy_dx_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uy_dx_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 1 ){
                        if( x2 == 0 ){
                                return uy_dy_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uy_dy_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uy_dy_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 2 ){
                        if( x2 == 0 ){
                                return uy_dz_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uy_dz_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uy_dz_dz( xx, yy, zz );
                        };
                };
        }
        else if( u == 2 ){
                if( x1 == 0 ){
                        if( x2 == 0 ){
                                return uz_dx_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uz_dx_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uz_dx_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 1 ){
                        if( x2 == 0 ){
                                return uz_dy_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uz_dy_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uz_dy_dz( xx, yy, zz );
                        };
                }
                else if( x1 == 2 ){
                        if( x2 == 0 ){
                                return uz_dz_dx( xx, yy, zz );
                        }
                        else if( x2 == 1 ){
                                return uz_dz_dy( xx, yy, zz );
                        }
                        else if( x2 == 2 ){
                                return uz_dz_dz( xx, yy, zz );
                        };
                };
        };
        return 0.0;
};       
//---------------------------------------------------------------------------
//���캯��
/*
AnalyticalSolution::AnalyticalSolution( int mid ){
        model_id = mid ;
};
//---------------------------------------------------------------------------
//������Ȼ�����ϵ������
void AnalyticalSolution::import_homo_para( char* filename ){
        ifstream import_s;
        import_s.open( filename );
        if( !import_s ) hout << "       Can not open this file: " << filename << endl;
        import_s >> E_modo ;
        import_s >> Mu_pora ;
     //   hout << " homo_E: " << homo_E << " homo_Mu: " << homo_Mu << endl;
};
//---------------------------------------------------------------------------
//���ò�������������أ����ȣ����Ծ�
void AnalyticalSolution::set_para( double M, double L, double I ){
        M_load  = M ;
        len     = L ;
        Iz      = I ;
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void AnalyticalSolution::max_strain( vector<double> &strain ){
        strain.assign(6,0.0);
     //   strain[0] = 
};
//---------------------------------------------------------------------------
//������Ȼ����������Ӧ��
void AnalyticalSolution::max_stress( vector<double> &stress ){

};
//---------------------------------------------------------------------------
//λ�ƺ���
double AnalyticalSolution::ux( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return M_load*x*y/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -M_load*(x*x+Mu_pora*y*y-Mu_pora*z*z)/(2.0*E_modo*Iz) ;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load*y*z/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
//λ�ƺ�����һ�׵���
double AnalyticalSolution::ux_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return M_load*y/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return M_load*x/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -M_load*x/(E_modo*Iz) ;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load*y/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return Mu_pora*M_load*z/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load*z/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load*y/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
//λ�ƺ����Ķ��׵���
double AnalyticalSolution::ux_dx_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dx_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return M_load/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dx_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dy_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dy_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::ux_dz_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dx_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -M_load/(E_modo*Iz) ;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dx_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dx_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dy_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dy_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uy_dz_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return Mu_pora*M_load/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dx_dx( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dx_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dx_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dy_dy( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dy_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return -Mu_pora*M_load/(E_modo*Iz);
                break;
        default:
                break;
        };
};
//---------------------------------------------------------------------------
double AnalyticalSolution::uz_dz_dz( double x, double y, double z ){
        switch( model_id ){
        case 1:
                return 0;
                break;
        default:
                break;
        };
};       */
//---------------------------------------------------------------------------
