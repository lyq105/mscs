//===========================================
//��Ҫ�������
//===========================================
#include "Mscm.h"

int Mscm::Begin(string in_file, ifstream &infile, string data_file, int CNum)
{
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��ִ�й���״��������ļ�
	clock_t ct_begin,ct_end;
	ct_begin = clock();
	hout << endl;
	hout << "*******************************************************************" << endl;
	hout << "-_- ��ʼ����......"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;
	clock_t ct0,ct1;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//������ģ��
lable_pcm: ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ������ģ......"<<endl;
	PCMCell *PCM = new PCMCell;
	if(PCM->Cell_generate(infile,data_file)==0)        //����0�����������ɵ�����1���ú������ɵ�����2��ʾ��ȡ����
	{
		infile.close();
		infile.open(in_file.c_str());
		//����ע���У������Ϳ���������������Ϣ��
		Get_Line(infile);
		delete PCM;
		goto lable_pcm;
	}
	ct1 = clock();
	hout << "    ��ģ�ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ������ģ��ϣ�"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���ɲ��Ͽ���Ϣ	
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ��ȡ������Ϣ......"<<endl;
	Matbase *Matrial = new Matbase;
	Matrial->Generate_matbase(infile);
//	Matrial->mats_vec[0].print();
	ct1 = clock();
	hout << "    ��ȡ������Ϣ�ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ��ȡ������Ϣ��ϣ�"<<endl<<endl;
	
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//��ģ������������Ϣ
	ct0 = clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ��������......"<<endl<<endl;
	Mesher *Mesh = new Mesher(PCM);
	if( Mesh->Mesh_generate(infile) == 0)
	{
		infile.close();
		infile.open(in_file.c_str());
		//����ע���У������Ϳ���������������Ϣ��
		Get_Line(infile);
		delete PCM;
		delete Matrial;
		goto lable_pcm;
	}
//	Mesh->Mesh_BinaryData(0,data_file, CNum);
//	Mesh->Mesh_data(0);
	ct1 = clock();
	hout << "    ���������ܺ�ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ����������ϣ�"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���Na1��Na1a2
	ct0=clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ���Ȼ����......"<<endl<<endl;
	//�����Ȼ�ϵ��
	//��ȡ���ݣ��ж��Ƿ����ǿ�ȼ���
	istringstream instr(Get_Line(infile));
	int elas_ana_only;
	instr >> elas_ana_only;

	HomoSolver *Hsolver = new HomoSolver(elas_ana_only,data_file);
	Hsolver->Solve(Mesh->nodes_vec, Mesh->bnodes_vec, Mesh->elements_vec, Matrial->mats_vec,PCM->Unitcell_V);
//	Hsolver->Na_BinaryData(0,Mesh->nodes_vec,data_file, CNum);
	ct1=clock();
	hout << "    ���Ȼ�����ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ ���Ȼ���������ϣ�"<<endl<<endl;

	//-----------------------------------------------------------------------------------------------------------------------------------------
	//���Ӧ��Ӧ�估����ǿ�ȼ���
	if( elas_ana_only == 0 ) return 1;		//ֻ���о��Ȼ����Գ�������
    else if( elas_ana_only > 0 ) resolve_load_condition(infile);	//��ȡ�غ���Ϣ
	else {hout <<"����������ж��Ƿ����ǿ�ȼ���ı�־С��0��"<<endl; return 0;}

	ct0=clock();
	hout << "======================================================" << endl;
	hout << "-_- ��ʼ����Ӧ��Ӧ�估����ǿ�ȼ���......"<<endl<<endl;

    double Li_epsilon[4];   //L1, L2, L3, Epsilon��ȡespilonΪ���ߴ磬Li<=1
	Li_epsilon[3] = 0.0001;
	Li_epsilon[0] = PCM->length() / Li_epsilon[3];
	Li_epsilon[1] = PCM->width() / Li_epsilon[3];
	Li_epsilon[2] = PCM->height() / Li_epsilon[3];

	ana_solu->set_homo_mat(Hsolver->homoMat);

	//�������ǿ�ȷ���
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

	//д���ļ�ǿ��ָ��
//	if( output_homo == 1 )	str_analyzer.write_result(i+1,j+1);

	ct1=clock();
	hout << "    Ӧ��Ӧ�估����ǿ�ȼ�������ʱ"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"�롣"<<endl;
	hout << "^_^ Ӧ��Ӧ�估����ǿ�ȼ�����������ϣ�"<<endl<<endl;
	//-----------------------------------------------------------------------------------------------------------------------------------------
	//�Ǿ��Ȼ��������ɵ�������ϵ��;
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
	hout << "    �����ܺ�ʱ"<<(double)(ct_end-ct_begin)/CLOCKS_PER_SEC<<"�롣" <<endl;
	hout << "^_^ ������ϣ�"<<endl;
	hout << "*******************************************************************" << endl;
	hout << endl;
	//������
	return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//�����غ���Ϣ
int Mscm::resolve_load_condition(ifstream &infile)
{
    //����������غ���Ϣ
	int cal_Na1a2 = 1;

	//����ṹ��̬���غ������Լ��غɴ�С����
	istringstream  istr5(Get_Line(infile));
	int stru_load_type;
	istr5 >> stru_load_type;
	if( stru_load_type <= -1 )
	{
		hout << "����! ����Ľṹ��̬���غ����Ͳ�������" << endl;
		hout << "�ò�����ֵ�� " << stru_load_type << endl;
		return 0;
	}
	if( stru_load_type < 3 )
	{
		cal_Na1a2=0;           //����Ҫ����Na1a2
		//����������죬���������������غɴ�С����������
		double area=0.0, Tv=0.0;
		istr5 >> area >> Tv;
		if( area <= 0.0 || Tv <= 0.0 )
		{
			hout << "����! ����������������غɴ�С����" << endl;
			hout << "�ṹ��̬���غ����Ͳ����� " << stru_load_type << endl;
			hout << "���������� " << area << endl;
			hout << "�غɴ�С�� " << Tv << endl;
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
		//������������
		//    static CanZBeamConYMomAnaSol ana_solu0;
		int sec_type;        //�������ͣ�0��Բ�Σ�1������
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
		//Բ�������θ˵�Ťת
		double r=0.0;              //����뾶
		double T=0.0;              //Ťת�غ�
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
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Mscm::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
