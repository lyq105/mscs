//===========================================================================
// HomoSolver.cpp
// 求解均匀化系数类成员函数
// Member Functions in a Class of the Homogeneous Solver
//===========================================================================

#include "HomoSolver.h"
//---------------------------------------------------------------------------
//求解，mod：0 计算Na1 , 1计算Na1a2
int HomoSolver::Solve(const vector<Node> &nodes_vec,const vector<int> &bnodes_vec, 
									const vector<Element> &elements_vec, const vector<MatPro> &mats_vec, const double unitcellV, int mod)
{
	clock_t ct0,ct1;
	//单胞体积；
	Unitcell_V=unitcellV;
	//-----------------------------------------------------------------------------------
	//紧缩存储刚度矩阵    
	ct0 = clock();
	hout<<"-_- 开始紧缩存储刚度矩阵"<<endl;
	int N =	(int)nodes_vec.size();							//节点个数
	int LR = (int)MAXIG*N;									//Ig的最大常数
	int *Iz = new int[N];											//动态申请Iz，Ig存储空间, Iz,Ig是变带宽存储的相关信息
	int *Ig = new int[MAXIG*N];
	for (int i=0; i<N; i++)	Iz[i] = 0;						//初始化零
	for (int i=0; i<MAXIG*N; i++)	Ig[i] = 0;			//初始化零

	SolveEqu	*Solv = new SolveEqu;						//定义解方程组类
	Solv->Gen_izig(nodes_vec, elements_vec, Iz, Ig, MAXIG*N);
	delete Solv;
	ct1 = clock();
	hout << "    紧缩存储刚度矩阵耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 紧缩存储刚度矩阵完毕!" << endl << endl;

	//---------------------------------------------------------------------------
	//生成总刚阵
	ct0 = clock();
	hout<<"-_- 开始生成总纲阵"<<endl;
	int L25 = 10*MAXIG*N;
	double *total_matrix = new double[L25];			//动态申请total_matrix存储空间
	for (int i=0; i<L25; i++)	total_matrix[i] = 0.0;		//初始化零

	GloStiffMatrix *gsmatrix = new GloStiffMatrix;	//定义总体刚度矩阵类
	gsmatrix->Gen_gsmatrix(total_matrix, Iz, Ig, nodes_vec, elements_vec, mats_vec);
	delete gsmatrix;
	ct1 = clock();
	hout << "    生成总纲阵耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout<<"^_^ 生成总纲阵完毕!"<<endl<<endl;

	//---------------------------------------------------------------------------
	//计算Na1系数
	ct0 = clock();
	hout<<"-_- 开始计算Na1系数"<<endl<<endl;
	Cal_Na1(total_matrix, Iz, Ig, nodes_vec, bnodes_vec, elements_vec, mats_vec);
    ct1 = clock();
	hout << "    计算Na1系数耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 计算Na1系数完毕!" << endl << endl;

	 //---------------------------------------------------------------------
    //计算Na1a2；
	if(Na1a2_key==1||Na1a2_key==2) Cal_Na1a2(total_matrix, Iz, Ig, nodes_vec, bnodes_vec, elements_vec, mats_vec);
	//如果是典型构件--拉伸柱体，位移展开式中对于Na1a2的导数为零，但是要把Na1a2赋值为0，为了后面计算的一致性。
	if(Na1a2_key==3) 
	{
		//初始化Na1a2_vec向量
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
//计算Na1系数
int HomoSolver::Cal_Na1(double* &total_matrix, int* &Iz, int* &Ig, const vector<Node> &nodes_vec, 
										const vector<int> &bnodes_vec, const vector<Element> &elements_vec, 
										const vector<MatPro> &mats_vec)
{
	clock_t ct0,ct1;
	ct0 = clock();
	hout<<"-_- 开始动态申请内存空间并初始化参数"<<endl;
	//---------------------------------------------------------------------------
	//初始化Na1_vec向量
	vector<double> temp_vec(27,0.0);
	Na1_vec.assign((int)nodes_vec.size(), temp_vec);

	//---------------------------------------------------------------------------
    //动态申请内存空间并初始化
	int N	= (int)nodes_vec.size();
	int	 N_B = (int)bnodes_vec.size();

	int *ip = new int[2*N_B];
	double *s = new double[NDIM*N];	//NDIM是自由度的维数
	double *v = new double[NDIM*N];
	double *cg = new double[NDIM*N];
	double *vp = new double[NDIM*N_B];
	double *uvw = new double[NDIM*N];
	double *remain = new double[NDIM*N];
	double *equright =new double[NDIM*N];
	//初始化边界点信息
	for (int i=0; i<N_B; i++)
	{
		ip[i*2]		= bnodes_vec[i]+1;
		ip[i*2+1]	= 7;			//用于调用函数（SolveEqu.tacp）中的goto语句
		vp[i*3]		= 0.0;
		vp[i*3+1]	= 0.0;
		vp[i*3+2]	= 0.0;
	}
	ct1 = clock();
	hout << "    动态申请内存空间并初始化参数耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 动态申请内存空间并初始化参数完毕!" <<endl<<endl;

	int count = 0;
	HomoPara Homo;
	//---------------------------------------------------------------------------
	//根据对称性，仅计算Na1m矩阵中的上三角
	for( int a1=0; a1<3; a1++ )
	{
		for( int m=a1; m<3; m++ )
		{
			clock_t t0,t1;
			t0 = clock();
			count++;
			hout << "-_- 开始第" << count << "次循环求解" << endl;
			//---------------------------------------------------------------------------
			//变量清零
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
			//求解右端项
			ct0 = clock();
			hout <<"      开始求解右端项"<<endl;
			Gloloaded_vector	*GloLoadVec = new Gloloaded_vector;		//定义求解右端项类
			GloLoadVec->Generate_gloloaded(elements_vec, nodes_vec, mats_vec, a1, m, equright);
			delete GloLoadVec;											//删除求解右端项类
			ct1 = clock();
			hout <<"      求解右端项耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
			hout <<"      求解右端项完毕!"<<endl;

			//---------------------------------------------------------------------------
			//处理位移边界条件
			ct0 = clock();
			hout <<"      开始处理位移边界条件"<<endl;
			SolveEqu	*Solv = new SolveEqu;						//定义解方程组类
			Solv->tacp(N_B, N, Iz, Ig, ip, vp, equright, total_matrix);
			ct1 = clock();
			hout <<"      求解右端项耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
			hout <<"      处理位移边界条件完毕!"<<endl;

			//---------------------------------------------------------------------------
			//求解线性方程组
			ct0 = clock();
			hout <<"      开始求解线性方程组"<<endl;
			Solv->sol(Ex, N_B, N, Iz, Ig, ip, vp, total_matrix, equright, cg, remain, s, v, uvw);
			delete Solv;														//删除解方程组类
			ct1 = clock();
			hout <<"      求解线性方程组耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
			hout <<"      求解线性方程组完毕!"<<endl;
			//---------------------------------------------------------------------------

			//得到Na1m
			int i=3*a1+m;
			for( int j=0; j<(int)nodes_vec.size(); j++ )
			{
				Na1_vec[j][3*i]   = uvw[3*j];
				Na1_vec[j][3*i+1] = uvw[3*j+1];
				Na1_vec[j][3*i+2] = uvw[3*j+2];
			}

            //求解均匀化系数
			cout<<"开始求解均匀化系数"<<endl;
			Homo.Generate_Homo(elements_vec, nodes_vec, mats_vec, uvw, Unitcell_V, a1, m);
			cout<<"求解均匀化系数完毕"<<endl;
			t1 = clock();
			hout <<"    第" << count << "次循环求解耗时：" << (double)(t1 - t0)/CLOCKS_PER_SEC << "秒"  << endl;
			hout << "^_^ 第" << count << "次循环求解完毕！" << endl<<endl;
		}
	}
	hout<<"开始求解均匀化热膨胀系数"<<endl;
	Homo.Generate_Homo_alpha(Unitcell_V);
	Homo.Generate_Homo_engconst(&homoMat);
//	Homo.print();
	Homo.output_Datafile(data_file);
	//为数据变量Homo_D赋值；
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			Homo_D[i][j]=Homo.Homo_D[i][j];
	//---------------------------------------------------------------------------
	//根据对称性，补上Na1m矩阵中的上三角
	ct0 = clock();
	hout<<"-_- 开始补上Na1m矩阵中的上三角"<<endl;
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
	hout << "    补上Na1m矩阵中的上三角耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒"  << endl;
	hout << "^_^ 补上Na1m矩阵中的上三角完毕!" <<endl<<endl;

	//---------------------------------------------------------------------------
	//释放内存
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
//计算Na1a2系数
int HomoSolver::Cal_Na1a2(	double* &total_matrix, int* &Iz, int* &Ig, const vector<Node> &nodes_vec, 
											const vector<int> &bnodes_vec, const vector<Element> &elements_vec, 
											const vector<MatPro> &mats_vec)
{
	//---------------------------------------------------------------------------
	//初始化Na1a2_vec向量
	vector<double> temp_vec(81,0.0);
	Na1a2_vec.assign((int)nodes_vec.size(), temp_vec);

	//---------------------------------------------------------------------------
    //动态申请内存空间并初始化
	int N	= (int)nodes_vec.size();
	int	 N_B = (int)bnodes_vec.size();

	int *ip = new int[2*N_B];
	double *s = new double[NDIM*N];	//NDIM是自由度的维数
	double *v = new double[NDIM*N];
	double *cg = new double[NDIM*N];
	double *vp = new double[NDIM*N_B];
	double *uvw = new double[NDIM*N];
	double *remain = new double[NDIM*N];
	double *equright =new double[NDIM*N];

	//初始化边界点信息
	for (int i=0; i<N_B; i++)
	{
		ip[i*2]		= bnodes_vec[i]+1;
		ip[i*2+1]	= 7;			//用于调用函数（SolveEqu.tacp）中的goto语句
		vp[i*3]		= 0.0;
		vp[i*3+1]	= 0.0;
		vp[i*3+2]	= 0.0;
	}

	//---------------------------------------------------------------------------
	//根据对称性，仅计算Na1m矩阵中的上三角
	for( int a1=0; a1<3; a1++ )
	{
		for( int a2=0; a2<3; a2++ )
		{
			for( int m=a1; m<3; m++ )
			{
				//---------------------------------------------------------------------------
				//变量清零
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
				hout<<"开始求解右端项"<<endl;
				Gloloaded_vector	*GloLoadVec = new Gloloaded_vector;		//定义求解右端项类
				GloLoadVec->Generate_gloloaded(elements_vec, nodes_vec, mats_vec, a1, a2, m, equright,Na1_vec,Homo_D);
				delete GloLoadVec;											//删除求解右端项类
				hout<<"求解右端项完毕"<<endl;

				//---------------------------------------------------------------------------
				//处理位移边界条件
				SolveEqu	*Solv = new SolveEqu;						//定义解方程组类
				Solv->tacp(N_B, N, Iz, Ig, ip, vp, equright, total_matrix);

				//---------------------------------------------------------------------------
				hout<<"处理位移边界条件完毕"<<endl;
				//求解线性方程组
				Solv->sol(Ex, N_B, N, Iz, Ig, ip, vp, total_matrix, equright, cg, remain, s, v, uvw);
				delete Solv;														//删除解方程组类

				//---------------------------------------------------------------------------
				hout<<"求解线性方程组完毕"<<endl;
				//得到Na1m
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
	//根据对称性，补上Na1m矩阵中的下三角
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
	//释放内存
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
//Na1和Na1a2的二进制数据的输出和读取
int HomoSolver::Na_BinaryData(int mod, const vector<Node> &nodes_vec, string data_file, int CNum)
{
	int num1 = CNum/10;
	int num2 = CNum - num1*10;
	char ch[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
	if(mod==0)			//输出数据
	{
		//---------------------------------------------------------------------
		//输出Na1和Na1a2数据
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Na";
		string NaName = title+str+ch[num1]+ch[num2]+".dat";
		ofstream out(NaName.c_str(),ios::binary);
		if(!out) { hout << "不能打开网格数据文件" << NaName << "!" << endl; exit(0); }
		//---------------------------------------------------------------------
		//节点个数
		int NodSize = int(nodes_vec.size());
		out.write((char *)&NodSize, sizeof(int));
		//---------------------------------------------------------------------
		//Na1数据
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<27; j++)
			{
				out.write((char *)&Na1_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		//Na1a2数据
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<81; j++)
			{
				out.write((char *)&Na1a2_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		out.close();
	}
	else if(mod==1)	//读取数据
	{
		//---------------------------------------------------------------------
		//读取Na1和Na1a2数据
        size_t en = data_file.rfind('.');
        string title = data_file.substr(0,en);
		string str = "Na";
		string NaName = title+str+ch[num1]+ch[num2]+".dat";
		ifstream in(NaName.c_str(),ios::binary);
		if(!in) { hout << "不能打开网格数据文件" << NaName << "!" << endl; exit(0); }
		//---------------------------------------------------------------------
		//节点个数
		int NodSize;
		in.read((char *)&NodSize, sizeof(int));
		//---------------------------------------------------------------------
		//Na1数据
		vector<double> temp_vec(27,0.0);
		Na1_vec.assign((int)nodes_vec.size(), temp_vec);
		for(int i=0; i<NodSize; i++)
			for(int j=0; j<27; j++)
			{
				in.read((char *)&Na1_vec[i][j], sizeof(double));
			}
		//---------------------------------------------------------------------
		//Na1a2数据
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
