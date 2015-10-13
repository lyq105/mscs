//===========================================================================
// HomoSolver.cpp
// �����Ȼ�ϵ�����Ա����
// Member Functions in a Class of the Homogeneous Solver
//===========================================================================

#include "HomoSolver.h"
//---------------------------------------------------------------------------
//��⣬mod��0 ����Na1 , 1����Na1a2
int HomoSolver::Solve(const vector<Node> &nodes_vec,const vector<int> &bnodes_vec, 
									const vector<Element> &elements_vec, const vector<MatPro> &mats_vec, const double unitcellV, int mod)
{
	clock_t ct0,ct1;
	//���������
	Unitcell_V=unitcellV;
	//-----------------------------------------------------------------------------------
	//�����洢�նȾ���    
	ct0 = clock();
	hout<<"-_- ��ʼ�����洢�նȾ���"<<endl;
	int N =	(int)nodes_vec.size();							//�ڵ����
	int LR = (int)MAXIG*N;									//Ig�������
	int *Iz = new int[N];											//��̬����Iz��Ig�洢�ռ�, Iz,Ig�Ǳ����洢�������Ϣ
	int *Ig = new int[MAXIG*N];
	for (int i=0; i<N; i++)	Iz[i] = 0;						//��ʼ����
	for (int i=0; i<MAXIG*N; i++)	Ig[i] = 0;			//��ʼ����

	SolveEqu	*Solv = new SolveEqu;						//����ⷽ������
	Solv->Gen_izig(nodes_vec, elements_vec, Iz, Ig, MAXIG*N);
	delete Solv;
	ct1 = clock();
	hout << "    �����洢�նȾ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ �����洢�նȾ������!" << endl << endl;

	//---------------------------------------------------------------------------
	//�����ܸ���
	ct0 = clock();
	hout<<"-_- ��ʼ�����ܸ���"<<endl;
	int L25 = 10*MAXIG*N;
	double *total_matrix = new double[L25];			//��̬����total_matrix�洢�ռ�
	for (int i=0; i<L25; i++)	total_matrix[i] = 0.0;		//��ʼ����

	GloStiffMatrix *gsmatrix = new GloStiffMatrix;	//��������նȾ�����
	gsmatrix->Gen_gsmatrix(total_matrix, Iz, Ig, nodes_vec, elements_vec, mats_vec);
	delete gsmatrix;
	ct1 = clock();
	hout << "    �����ܸ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout<<"^_^ �����ܸ������!"<<endl<<endl;

	//---------------------------------------------------------------------------
	//����Na1ϵ��
	ct0 = clock();
	hout<<"-_- ��ʼ����Na1ϵ��"<<endl<<endl;
	Cal_Na1(total_matrix, Iz, Ig, nodes_vec, bnodes_vec, elements_vec, mats_vec);
    ct1 = clock();
	hout << "    ����Na1ϵ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ ����Na1ϵ�����!" << endl << endl;

	 //---------------------------------------------------------------------
    //����Na1a2��
	if(Na1a2_key==1||Na1a2_key==2) Cal_Na1a2(total_matrix, Iz, Ig, nodes_vec, bnodes_vec, elements_vec, mats_vec);
	//����ǵ��͹���--�������壬λ��չ��ʽ�ж���Na1a2�ĵ���Ϊ�㣬����Ҫ��Na1a2��ֵΪ0��Ϊ�˺�������һ���ԡ�
	if(Na1a2_key==3) 
	{
		//��ʼ��Na1a2_vec����
		vector<double> temp_vec(81,0.0);
		Na1a2_vec.assign((int)nodes_vec.size(), temp_vec);
	}

	//---------------------------------------------------------------------
	//ofstream out("Na1a2.dat");
	//out << int(Na1a2_vec.size()) << endl;
	//for(int i=0; i<int(Na1a2_vec.size()); i++)
	//{
	//	int num_Na1a2 = int(Na1a2_vec[i].size());
	//	out << i << "  " << num_Na1a2 << "  ";
	//	for(int j=0; j<num_Na1a2; j++)
	//	{
	//		out << Na1a2_vec[i][j] << "  ";
	//	}
	//	out << endl;
	//}
	//out.close();

	delete [] Iz;
	delete [] Ig;
	delete [] total_matrix;

	return 1;
}
//---------------------------------------------------------------------------
//����Na1ϵ��
int HomoSolver::Cal_Na1(double* &total_matrix, int* &Iz, int* &Ig, const vector<Node> &nodes_vec, 
										const vector<int> &bnodes_vec, const vector<Element> &elements_vec, 
										const vector<MatPro> &mats_vec)
{
	clock_t ct0,ct1;
	ct0 = clock();
	hout<<"-_- ��ʼ��̬�����ڴ�ռ䲢��ʼ������"<<endl;
	//---------------------------------------------------------------------------
	//��ʼ��Na1_vec����
	vector<double> temp_vec(27,0.0);
	Na1_vec.assign((int)nodes_vec.size(), temp_vec);

	//---------------------------------------------------------------------------
    //��̬�����ڴ�ռ䲢��ʼ��
	int N	= (int)nodes_vec.size();
	int	 N_B = (int)bnodes_vec.size();

	int *ip = new int[2*N_B];
	double *s = new double[NDIM*N];	//NDIM�����ɶȵ�ά��
	double *v = new double[NDIM*N];
	double *cg = new double[NDIM*N];
	double *vp = new double[NDIM*N_B];
	double *uvw = new double[NDIM*N];
	double *remain = new double[NDIM*N];
	double *equright =new double[NDIM*N];
	//��ʼ���߽����Ϣ
	for (int i=0; i<N_B; i++)
	{
		ip[i*2]		= bnodes_vec[i]+1;
		ip[i*2+1]	= 7;			//���ڵ��ú�����SolveEqu.tacp���е�goto���
		vp[i*3]		= 0.0;
		vp[i*3+1]	= 0.0;
		vp[i*3+2]	= 0.0;
	}
	ct1 = clock();
	hout << "    ��̬�����ڴ�ռ䲢��ʼ��������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ ��̬�����ڴ�ռ䲢��ʼ���������!" <<endl<<endl;

	int count = 0;
	HomoPara Homo;
	//---------------------------------------------------------------------------
	//���ݶԳ��ԣ�������Na1m�����е�������
	for( int a1=0; a1<3; a1++ )
	{
		for( int m=a1; m<3; m++ )
		{
			clock_t t0,t1;
			t0 = clock();
			count++;
			hout << "-_- ��ʼ��" << count << "��ѭ�����" << endl;
			//---------------------------------------------------------------------------
			//��������
			for( int j=0; j<NDIM*N; j++ )
			{
				s[j]			= 0;
				v[j]			= 0;
				cg[j]			= 0;
				uvw[j]		= 0;
				remain[j]	= 0;
				equright[j]	= 0;
			}

			//---------------------------------------------------------------------------
			//����Ҷ���
			ct0 = clock();
			hout <<"      ��ʼ����Ҷ���"<<endl;
			Gloloaded_vector	*GloLoadVec = new Gloloaded_vector;		//��������Ҷ�����
			GloLoadVec->Generate_gloloaded(elements_vec, nodes_vec, mats_vec, a1, m, equright);
			delete GloLoadVec;											//ɾ������Ҷ�����
			ct1 = clock();
			hout <<"      ����Ҷ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
			hout <<"      ����Ҷ������!"<<endl;

			//---------------------------------------------------------------------------
			//����λ�Ʊ߽�����
			ct0 = clock();
			hout <<"      ��ʼ����λ�Ʊ߽�����"<<endl;
			SolveEqu	*Solv = new SolveEqu;						//����ⷽ������
			Solv->tacp(N_B, N, Iz, Ig, ip, vp, equright, total_matrix);
			ct1 = clock();
			hout <<"      ����Ҷ����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
			hout <<"      ����λ�Ʊ߽��������!"<<endl;

			//---------------------------------------------------------------------------
			//������Է�����
			ct0 = clock();
			hout <<"      ��ʼ������Է�����"<<endl;
			Solv->sol(Ex, N_B, N, Iz, Ig, ip, vp, total_matrix, equright, cg, remain, s, v, uvw);
			delete Solv;														//ɾ���ⷽ������
			ct1 = clock();
			hout <<"      ������Է������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
			hout <<"      ������Է��������!"<<endl;
			//---------------------------------------------------------------------------

			//�õ�Na1m
			int i=3*a1+m;
			for( int j=0; j<(int)nodes_vec.size(); j++ )
			{
				Na1_vec[j][3*i]   = uvw[3*j];
				Na1_vec[j][3*i+1] = uvw[3*j+1];
				Na1_vec[j][3*i+2] = uvw[3*j+2];
			}

            //�����Ȼ�ϵ��
			cout<<"��ʼ�����Ȼ�ϵ��"<<endl;
			Homo.Generate_Homo(elements_vec, nodes_vec, mats_vec, uvw, Unitcell_V, a1, m);
			cout<<"�����Ȼ�ϵ�����"<<endl;
			t1 = clock();
			hout <<"    ��" << count << "��ѭ������ʱ��" << (double)(t1 - t0)/CLOCKS_PER_SEC << "��"  << endl;
			hout << "^_^ ��" << count << "��ѭ�������ϣ�" << endl<<endl;
		}
	}
	hout<<"��ʼ�����Ȼ�������ϵ��"<<endl;
	Homo.Generate_Homo_alpha(Unitcell_V);
	Homo.Generate_Homo_engconst(&homoMat);
//	Homo.print();
	Homo.output_Datafile(data_file);
	//Ϊ���ݱ���Homo_D��ֵ��
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			Homo_D[i][j]=Homo.Homo_D[i][j];
	//---------------------------------------------------------------------------
	//���ݶԳ��ԣ�����Na1m�����е�������
	ct0 = clock();
	hout<<"-_- ��ʼ����Na1m�����е�������"<<endl;
	for( int a1=0; a1<3; a1++ )
	{
		for( int m=0; m<a1; m++ )
		{
			int i = 3*a1+m;
			int k = 3*m+a1;  
			for( int j=0; j<(int)nodes_vec.size(); j++ )
			{
				Na1_vec[j][3*i]   = Na1_vec[j][3*k];
				Na1_vec[j][3*i+1] = Na1_vec[j][3*k+1];
				Na1_vec[j][3*i+2] = Na1_vec[j][3*k+2];
			}
		}
	}
	ct1 = clock();
	hout << "    ����Na1m�����е������Ǻ�ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��"  << endl;
	hout << "^_^ ����Na1m�����е����������!" <<endl<<endl;

	//---------------------------------------------------------------------------
	//�ͷ��ڴ�
	delete [] ip;
	delete [] s;
	delete [] v;	
	delete [] cg;
	delete [] vp;
	delete [] uvw;	
	delete [] remain;
	delete [] equright;

	return 1;
}
//============================================================================
//---------------------------------------------------------------------------
//����Na1a2ϵ��
int HomoSolver::Cal_Na1a2(	double* &total_matrix, int* &Iz, int* &Ig, const vector<Node> &nodes_vec, 
											const vector<int> &bnodes_vec, const vector<Element> &elements_vec, 
											const vector<MatPro> &mats_vec)
{
	//---------------------------------------------------------------------------
	//��ʼ��Na1a2_vec����
	vector<double> temp_vec(81,0.0);
	Na1a2_vec.assign((int)nodes_vec.size(), temp_vec);

	//---------------------------------------------------------------------------
    //��̬�����ڴ�ռ䲢��ʼ��
	int N	= (int)nodes_vec.size();
	int	 N_B = (int)bnodes_vec.size();

	int *ip = new int[2*N_B];
	double *s = new double[NDIM*N];	//NDIM�����ɶȵ�ά��
	double *v = new double[NDIM*N];
	double *cg = new double[NDIM*N];
	double *vp = new double[NDIM*N_B];
	double *uvw = new double[NDIM*N];
	double *remain = new double[NDIM*N];
	double *equright =new double[NDIM*N];

	//��ʼ���߽����Ϣ
	for (int i=0; i<N_B; i++)
	{
		ip[i*2]		= bnodes_vec[i]+1;
		ip[i*2+1]	= 7;			//���ڵ��ú�����SolveEqu.tacp���е�goto���
		vp[i*3]		= 0.0;
		vp[i*3+1]	= 0.0;
		vp[i*3+2]	= 0.0;
	}

	//---------------------------------------------------------------------------
	//���ݶԳ��ԣ�������Na1m�����е�������
	for( int a1=0; a1<3; a1++ )
	{
		for( int a2=0; a2<3; a2++ )
		{
			for( int m=a1; m<3; m++ )
			{
				//---------------------------------------------------------------------------
				//��������
				for( int j=0; j<NDIM*N; j++ )
				{
					s[j]			= 0;
					v[j]			= 0;
					cg[j]			= 0;
					uvw[j]		= 0;
					remain[j]	= 0;
					equright[j]	= 0;
				}

				//---------------------------------------------------------------------------
				hout<<"��ʼ����Ҷ���"<<endl;
				Gloloaded_vector	*GloLoadVec = new Gloloaded_vector;		//��������Ҷ�����
				GloLoadVec->Generate_gloloaded(elements_vec, nodes_vec, mats_vec, a1, a2, m, equright,Na1_vec,Homo_D);
				delete GloLoadVec;											//ɾ������Ҷ�����
				hout<<"����Ҷ������"<<endl;

				//---------------------------------------------------------------------------
				//����λ�Ʊ߽�����
				SolveEqu	*Solv = new SolveEqu;						//����ⷽ������
				Solv->tacp(N_B, N, Iz, Ig, ip, vp, equright, total_matrix);

				//---------------------------------------------------------------------------
				hout<<"����λ�Ʊ߽��������"<<endl;
				//������Է�����
				Solv->sol(Ex, N_B, N, Iz, Ig, ip, vp, total_matrix, equright, cg, remain, s, v, uvw);
				delete Solv;														//ɾ���ⷽ������

				//---------------------------------------------------------------------------
				hout<<"������Է��������"<<endl;
				//�õ�Na1m
				int i=9*a1+3*a2+m;
				for( int j=0; j<(int)nodes_vec.size(); j++ )
				{
					Na1a2_vec[j][3*i]   = uvw[3*j];
					Na1a2_vec[j][3*i+1] = uvw[3*j+1];
					Na1a2_vec[j][3*i+2] = uvw[3*j+2];
				}
			}
		}
	}
	//---------------------------------------------------------------------------
	//���ݶԳ��ԣ�����Na1m�����е�������
	for( int a1=0; a1<3; a1++ )
	{
		for( int a2=0; a2<3; a2++ )
		{
			for( int m=0; m<a1; m++ )
			{
				int i = 9*a1+3*a2+m;
				int k = 9*m+3*a2+a1;  
				for( int j=0; j<(int)nodes_vec.size(); j++ )
				{
					Na1a2_vec[j][3*i]   = Na1a2_vec[j][3*k];
					Na1a2_vec[j][3*i+1] = Na1a2_vec[j][3*k+1];
					Na1a2_vec[j][3*i+2] = Na1a2_vec[j][3*k+2];
				}
			}
		}
	}
	//---------------------------------------------------------------------------
	//�ͷ��ڴ�
	delete [] ip;
	delete [] s;
	delete [] v;	
	delete [] cg;
	delete [] vp;
	delete [] uvw;	
	delete [] remain;
	delete [] equright;

	return 1;
}
//---------------------------------------------------------------------------
//Na1��Na1a2�Ķ��������ݵ�����Ͷ�ȡ
int HomoSolver::Na_BinaryData(int mod, const vector<Node> &nodes_vec, string data_file, int CNum)
{
	int num1 = CNum/10;
	int num2 = CNum - num1*10;
	char ch[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	if(mod==0)			//�������
	{
		//---------------------------------------------------------------------
		//���Na1��Na1a2����
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Na";
		string NaName = title+str+ch[num1]+ch[num2]+".dat";
		ofstream out(NaName.c_str(),ios::binary);
		if(!out) { hout << "���ܴ����������ļ�" << NaName << "!" << endl; exit(0); }
		//---------------------------------------------------------------------
		//�ڵ����
		int NodSize = int(nodes_vec.size());
		out.write((char *)&NodSize, sizeof(int));
		//---------------------------------------------------------------------
		//Na1����
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<27; j++)
			{
				out.write((char *)&Na1_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		//Na1a2����
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<81; j++)
			{
				out.write((char *)&Na1a2_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		out.close();
	}
	else if(mod==1)	//��ȡ����
	{
		//---------------------------------------------------------------------
		//��ȡNa1��Na1a2����
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Na";
		string NaName = title+str+ch[num1]+ch[num2]+".dat";
		ifstream in(NaName.c_str(),ios::binary);
		if(!in) { hout << "���ܴ����������ļ�" << NaName << "!" << endl; exit(0); }
		//---------------------------------------------------------------------
		//�ڵ����
		int NodSize;
		in.read((char *)&NodSize, sizeof(int));
		//---------------------------------------------------------------------
		//Na1����
		vector<double> temp_vec(27,0.0);
		Na1_vec.assign((int)nodes_vec.size(), temp_vec);
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<27; j++)
			{
				in.read((char *)&Na1_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		//Na1a2����
		temp_vec.assign(81,0.0);
		Na1a2_vec.assign((int)nodes_vec.size(), temp_vec);
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<81; j++)
			{
				in.read((char *)&Na1a2_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		in.close();	
	}

	return 1;
}
//============================================================================
