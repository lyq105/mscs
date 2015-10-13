//===========================================
//主要计算过程
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
	AnalyticalSolution *ana_solu;	//解析解
private:
	string Get_Line(ifstream &infile)const;				//读入信息一行，跳过注释行（以%开头）	
    int resolve_load_condition(ifstream &infile);		//解析载荷信息
};
