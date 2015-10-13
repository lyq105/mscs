//===========================================
//主要计算过程
//===========================================
#include "Mscm.h"

int Mscm::Begin(string in_file, ifstream &infile, string data_file, int CNum)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//将执行过程状况输出在文件
	clock_t ct_begin,ct_end;
	ct_begin = clock();
	hout << endl;
	hout << "*******************************************************************" << endl;
	hout << "-_- 开始计算......"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;
	clock_t ct0,ct1;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//单胞建模；
lable_pcm: ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- 开始单胞建模......"<<endl;
	PCMCell *PCM = new PCMCell;
	if(PCM->Cell_generate(infile,data_file)==0)        //参数0：用于艳生成的椭球，1：用韩非生成的椭球，2表示读取数据
	{
		infile.close();
		infile.open(in_file.c_str());
		//跳过注释行（这样就可以跳过样本数信息）
		Get_Line(infile);
		delete PCM;
		goto lable_pcm;
	}
	ct1 = clock();
	hout << "    建模总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 单胞建模完毕！"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//生成材料库信息	
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- 开始读取材料信息......"<<endl;
	Matbase *Matrial = new Matbase;
	Matrial->Generate_matbase(infile);
//	Matrial->mats_vec[0].print();
	ct1 = clock();
	hout << "    读取材料信息总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 读取材料信息完毕！"<<endl<<endl;
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//建模并生成网格信息
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- 开始生成网格......"<<endl<<endl;
	Mesher *Mesh = new Mesher(PCM);
	if( Mesh->Mesh_generate(infile) == 0)
	{
		infile.close();
		infile.open(in_file.c_str());
		//跳过注释行（这样就可以跳过样本数信息）
		Get_Line(infile);
		delete PCM;
		delete Matrial;
		goto lable_pcm;
	}
//	Mesh->Mesh_BinaryData(0,data_file, CNum);
//	Mesh->Mesh_data(0);
	ct1 = clock();
	hout << "    生成网格总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 网格生成完毕！"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//求解Na1和Na1a2
	ct0=clock();
	hout << "======================================================" << endl;
	hout << "-_- 开始均匀化求解......"<<endl<<endl;
	//求解均匀化系数
	//读取数据，判断是否进行强度计算
	istringstream instr(Get_Line(infile));
	int elas_ana_only;
	instr >> elas_ana_only;

	HomoSolver *Hsolver = new HomoSolver(elas_ana_only,data_file);
	Hsolver->Solve(Mesh->nodes_vec, Mesh->bnodes_vec, Mesh->elements_vec, Matrial->mats_vec,PCM->Unitcell_V);
//	Hsolver->Na_BinaryData(0,Mesh->nodes_vec,data_file, CNum);
	ct1=clock();
	hout << "    均匀化求解耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 均匀化求解过程完毕！"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//求解应力应变及弹性强度极限
	if( elas_ana_only == 0 ) return 1;		//只进行均匀化弹性常数计算
    else if( elas_ana_only > 0 ) resolve_load_condition(infile);	//读取载荷信息
	else {hout <<"错误！输入的判断是否进行强度计算的标志小于0。"<<endl; return 0;}

	ct0=clock();
	hout << "======================================================" << endl;
	hout << "-_- 开始计算应力应变及弹性强度极限......"<<endl<<endl;

    double Li_epsilon[4];   //L1, L2, L3, Epsilon，取espilon为最大尺寸，Li<=1
	Li_epsilon[3] = 0.0001;
	Li_epsilon[0] = PCM->length() / Li_epsilon[3];
	Li_epsilon[1] = PCM->width() / Li_epsilon[3];
	Li_epsilon[2] = PCM->height() / Li_epsilon[3];

	ana_solu->set_homo_mat(Hsolver->homoMat);

	//结果处理，强度分析
	Analyzer str_analyzer = Analyzer(&Mesh->nodes_vec, &Mesh->elements_vec, &Matrial->mats_vec);
	//str_analyzer.set_title(title);
	//str_analyzer.set_interior_eles(&Mesh->ieles);
	//str_analyzer.set_interior_nodes(&inodes);
	str_analyzer.set_cell_size(Li_epsilon);
	//str_analyzer.set_cell_origin(ox,oy,oz);
	str_analyzer.set_ana_solution(ana_solu);
	str_analyzer.set_Na1_vec(&Hsolver->Na1_vec);
	str_analyzer.set_Na1a2_vec(&Hsolver->Na1a2_vec);
	str_analyzer.set_homo_mat(&Hsolver->homoMat);
	str_analyzer.set_criterions(3,3);
	str_analyzer.set_string_datafile(data_file);
	str_analyzer.strength_analysis_par(elas_ana_only-1);

	//写入文件强度指标
//	if( output_homo == 1 )	str_analyzer.write_result(i+1,j+1);

	ct1=clock();
	hout << "    应力应变及弹性强度极限求解耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
	hout << "^_^ 应力应变及弹性强度极限求解过程完毕！"<<endl<<endl;
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//非均匀化方法生成的热膨胀系数;
	//Nhomopara *Npara = new Nhomopara;
	//Npara->Generate_Nhomepara(*Cell,*Matrial);
	//hout<<Npara->cte1_SH<<"  "<<Npara->cte2_SH;

	delete Hsolver;
	delete PCM;
	delete Matrial;
	delete Mesh;
//	delete Npara;

    //-----------------------------------------------------------------------------------------------------------------------------------------
	ct_end = clock();
	hout << "*******************************************************************" << endl;
	hout << "    计算总耗时"<<(double)(ct_end-ct_begin)/CLOCKS_PER_SEC<<"秒。" <<endl;
	hout << "^_^ 计算完毕！"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;
	//求解完毕
	return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//解析载荷信息
int Mscm::resolve_load_condition(ifstream &infile)
{
    //处理输入的载荷信息
	int cal_Na1a2 = 1;

	//处理结构形态、载荷类型以及载荷大小数据
	istringstream  istr5(Get_Line(infile));
	int stru_load_type;
	istr5 >> stru_load_type;
	if( stru_load_type <= -1 )
	{
		hout << "错误! 读入的结构形态及载荷类型参数错误。" << endl;
		hout << "该参数的值是 " << stru_load_type << endl;
		return 0;
	}
	if( stru_load_type < 3 )
	{
		cal_Na1a2=0;           //不需要计算Na1a2
		//柱体均匀拉伸，读入柱体截面积、载荷大小（集中力）
		double area=0.0, Tv=0.0;
		istr5 >> area >> Tv;
		if( area <= 0.0 || Tv <= 0.0 )
		{
			hout << "错误! 读入的柱体截面积或载荷大小错误。" << endl;
			hout << "结构形态及载荷类型参数是 " << stru_load_type << endl;
			hout << "柱体截面积是 " << area << endl;
			hout << "载荷大小是 " << Tv << endl;
			return 0;
		}
        if( stru_load_type == 0 )
		{
			static FixRodXTenAnaSol ana_solu0( area, Tv );
			ana_solu = &ana_solu0;
        }
        else if( stru_load_type == 1 )
		{
			static FixRodYTenAnaSol ana_solu0( area, Tv );
			ana_solu = &ana_solu0;
        }
        else if( stru_load_type == 2 )
		{
			static FixRodZTenAnaSol ana_solu0( area, Tv );
			ana_solu = &ana_solu0;
        }
	}
	else if( stru_load_type < 6 )
	{
		//悬臂梁纯弯曲
		//    static CanZBeamConYMomAnaSol ana_solu0;
		int sec_type;        //截面类型：0：圆形，1：矩形
		istr5 >> sec_type ;
		static BeamSecShape bsec;
		double Mv=0.0;
		if( sec_type == 0 )
		{
			double rr=0.0;
			istr5 >> rr >> Mv;
			bsec = BeamSecShape(rr);
			//ana_solu0.set_para( &circle_sec, Mv );
		}
		else if( sec_type == 1 )
		{
			double len=0.0,wid=0.0;
			istr5 >> len >> wid >> Mv;
			bsec = BeamSecShape(len,wid);
			//ana_solu0.set_para( &rect_sec, Mv );
		}
		if( stru_load_type == 3 )
		{
			static CanXBeamConZMomAnaSol ana_solu0;
			ana_solu0.set_para(&bsec, Mv);
			ana_solu = &ana_solu0;
		}
		else if( stru_load_type == 4 )
		{
			static CanYBeamConXMomAnaSol ana_solu0;
			ana_solu0.set_para(&bsec, Mv);
			ana_solu = &ana_solu0;
		}
		else if( stru_load_type == 5 )
		{
			static CanZBeamConYMomAnaSol ana_solu0;
			ana_solu0.set_para(&bsec, Mv);
			ana_solu = &ana_solu0;
		}
	}
	else if( stru_load_type < 9 )
	{
		//圆截面柱形杆的扭转
		double r=0.0;              //截面半径
		double T=0.0;              //扭转载荷
		istr5 >> r >> T ;
		if( stru_load_type == 6 )
		{
			static FixRodXTwistAnaSol ana_solu0(r,T);
			ana_solu = &ana_solu0;
		}
		else if( stru_load_type == 7 )
		{
			static FixRodYTwistAnaSol ana_solu0(r,T);
			ana_solu = &ana_solu0;
		}
		else if( stru_load_type == 8 )
		{
			static FixRodZTwistAnaSol ana_solu0(r,T);
			ana_solu = &ana_solu0;
		}
	}
        
	return 1;
}
//---------------------------------------------------------------------------
//读入一行信息，并跳过注释行（以"%"开头）；
string Mscm::Get_Line(ifstream &infile)const
{
	string s;
	//读入信息一行
	getline(infile,s);
	//跳过注释行     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
