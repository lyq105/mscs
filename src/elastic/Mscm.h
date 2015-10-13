//===========================================
//��Ҫ�������
//===========================================

#include "AnalyticalSolution.h"
#include "analyzer.h"
#include "BeamSecShape.h"
#include "CanXBeamConZMomAnaSol.h"
#include "CanYBeamConXMomAnaSol.h"  
#include "CanZBeamConYMomAnaSol.h"
#include "FixRodXTenAnaSol.h"
#include "FixRodYTenAnaSol.h"
#include "FixRodZTenAnaSol.h"     
#include "FixRodXTwistAnaSol.h"
#include "FixRodYTwistAnaSol.h"
#include "FixRodZTwistAnaSol.h"
#include "PCMCell.h"
#include "Matbase.h"
#include "Mesher.h"
#include "HomoSolver.h"
//#include "Nhomopara.h"
#include "Hns.h"

#include "time.h"
using namespace hns;

//---------------------------------------------------------------------------
class Mscm
{
public: 
	int Begin(string in_file, ifstream &infile, string data_file, int CNum);
protected:
	AnalyticalSolution *ana_solu;	//������
private:
	string Get_Line(ifstream &infile)const;				//������Ϣһ�У�����ע���У���%��ͷ��	
    int resolve_load_condition(ifstream &infile);		//�����غ���Ϣ
};
