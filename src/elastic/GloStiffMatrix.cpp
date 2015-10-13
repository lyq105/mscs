//===========================================================================
// GloStiffMatrix.cpp
// 求解总刚阵类成员函数
// Member Functions in a Class of the Global Stiff Matrix
//===========================================================================

#include "GloStiffMatrix.h"
//---------------------------------------------------------------------------
//生成总体刚度矩阵
int GloStiffMatrix::Gen_gsmatrix(double* &total_matrix, int* &Iz, int* &Ig, const vector<Node> &nodes_vec,
								 const vector<Element> &elements_vec, const vector<MatPro> &mats_vec)
{
	//---------------------------------------------------------------------------
	//取高斯点
	Gauss gau;
	gau.Generate_gauss();
	//生成总刚
	for( int i=0; i<(int)elements_vec.size(); i++ )
	{    
		//生成单刚
		Generate_elestiff(elements_vec[i],nodes_vec,mats_vec,gau.gauss,gau.wight);
		//单刚添加到总刚
		Add_to_gsmatrix(total_matrix, Iz, Ig,elements_vec[i]);
	}
	return 1;
}
//---------------------------------------------------------------------------------
//生成单刚
int GloStiffMatrix::Generate_elestiff(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const vector<Node> &gauss,const vector<double> &wight){

	//----------------------------------------------------------
	//求出此单元的所有节点坐标
	vector<Node> elenodes_vec;
	for(int i=0;i<int(e.nodes_id.size());i++)
		elenodes_vec.push_back(nodes_vec[e.nodes_id[i]]);
	//-----------------------------------------------------------------
	//读取此单元所对应的材料弹性矩阵；
	double ele_elas[6][6];
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			ele_elas[i][j] = mats_vec[e.mat].elas_matrix[i][j];
	//-------------------------------------------------------------------
	if(e.type==2&&int(e.nodes_id.size())==4){
		Tetrahedron_line(elenodes_vec,ele_elas);
	}
	else if(e.type==3&&int(e.nodes_id.size())==6){
		Threeprism_line(elenodes_vec,ele_elas,gauss,wight);
	}

	return 1;
}//-----------------------------------------------------------------------------
//四面体线性插值
void GloStiffMatrix::Tetrahedron_line(const vector<Node> &elenodes_vec,double ele_elas[][6]){
	vector<double> tem(12,0);
	vector<vector<double> > tem2(12,tem);
	//===============================================================================
	//首先计算节点基函数Ni，Nj，Nm，Nl；（其实是计算向量a,b,c）;
	vector<double> a(4),b(4),c(4),d(4);
	int i=0;
	int j=1;
	int m=2;
	int l=3;
	for(int k=0;k<4;k++){
		a[i]= elenodes_vec[j].x*(elenodes_vec[m].y*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].y)+
			elenodes_vec[j].y*(elenodes_vec[m].z*elenodes_vec[l].x-elenodes_vec[m].x*elenodes_vec[l].z)+
			elenodes_vec[j].z*(elenodes_vec[m].x*elenodes_vec[l].y-elenodes_vec[m].y*elenodes_vec[l].x);
		b[i]= -(elenodes_vec[m].y*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].y)
			+(elenodes_vec[j].y*elenodes_vec[l].z-elenodes_vec[j].z*elenodes_vec[l].y)
			-(elenodes_vec[j].y*elenodes_vec[m].z-elenodes_vec[j].z*elenodes_vec[m].y);
		c[i]= (elenodes_vec[m].x*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].x)
			-(elenodes_vec[j].x*elenodes_vec[l].z-elenodes_vec[j].z*elenodes_vec[l].x)
			+(elenodes_vec[j].x*elenodes_vec[m].z-elenodes_vec[j].z*elenodes_vec[m].x);
		d[i]= -(elenodes_vec[m].x*elenodes_vec[l].y-elenodes_vec[m].y*elenodes_vec[l].x)
			+(elenodes_vec[j].x*elenodes_vec[l].y-elenodes_vec[j].y*elenodes_vec[l].x)
			-(elenodes_vec[j].x*elenodes_vec[m].y-elenodes_vec[j].y*elenodes_vec[m].x);
		int n;
		n=i;
		i=j;
		j=m;
		m=l;
		l=n;
	}
	//--------------------------------------------------------------------
	//计算四面体体积；
	double volume = 1.0/6.0*(a[0]-a[1]+a[2]-a[3]);
	//--------------------------------------------------------------------
	//计算B矩阵
	double B[6][12];
	for(int i=0;i<6;i++)
		for(int j=0;j<12;j++)
			B[i][j]=0;

	for(int i=0;i<4;i++){
		if(i%2==1){
			b[i]=-b[i];
			c[i]=-c[i];
			d[i]=-d[i];
		}
		B[0][i*3+0]=b[i]/(6.0*volume);
		B[1][i*3+1]=c[i]/(6.0*volume);
		B[2][i*3+2]=d[i]/(6.0*volume);
		B[3][i*3+0]=c[i]/(6.0*volume);
		B[3][i*3+1]=b[i]/(6.0*volume);
		B[4][i*3+1]=d[i]/(6.0*volume);
		B[4][i*3+2]=c[i]/(6.0*volume);
		B[5][i*3+0]=d[i]/(6.0*volume);
		B[5][i*3+2]=b[i]/(6.0*volume);
	}

	//-------------------------------------------------------
	//计算B_trans;
	double B_trans[12][6];
	for(int i=0;i<12;i++)
		for(int j=0;j<6;j++)
			B_trans[i][j] = B[j][i];
	//--------------------------------------------------------------------------------------------
	//求出B_trans矩阵与ele_elas矩阵的乘积array1；
	double array1[12][6];
	for(int i=0;i<12;i++)
		for(int j=0;j<6;j++){
			array1[i][j]=0; 
			for(int k=0;k<6;k++)
				array1[i][j] = array1[i][j]+B_trans[i][k]*ele_elas[k][j];
		}
		//求出array1矩阵与B矩阵的乘积array2;
		double array2[12][12];
		for(int i=0;i<12;i++)
			for(int j=0;j<12;j++){
				array2[i][j]=0;
				for(int k=0;k<6;k++)
					array2[i][j] = array2[i][j]+array1[i][k]*B[k][j];
				tem2[i][j]=tem2[i][j]+array2[i][j]*volume;
			}
			elestiff = tem2;
}
//-----------------------------------------------------------------------------------------------
//三棱柱线性插值
void GloStiffMatrix::Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const vector<Node> &gauss,const vector<double> &wight){
	vector<double> tem(18,0);
	vector<vector<double> > tem2(18,tem);
	//--------------------------------------------------
	/*double tem7=0;
	for(int i=0;i<wight.size();i++){
	cout<<gauss[i];
	tem7=tem7+wight[i];
	cout<<tem7<<endl;
	}*/
	//循环高斯点计算积分；
	for(int count=0;count<int(wight.size());count++){
		//计算Ｊ矩阵；
		//--------------------------------------------
		//形函数N对gauss[count].x,gauss[count].y,gauss[count].z的偏导矩阵
		double diff[3][6];
		diff[0][0]=0.5*(1.0-gauss[count].z);  
		diff[0][1]=0;                         
		diff[0][2]=-0.5*(1.0-gauss[count].z);
		diff[0][3]=0.5*(1.0+gauss[count].z);
		diff[0][4]=0;
		diff[0][5]=-0.5*(1.0+gauss[count].z);

		diff[1][0]=0;
		diff[1][1]=diff[0][0];
		diff[1][2]=diff[0][2];
		diff[1][3]=0;
		diff[1][4]=diff[0][3];
		diff[1][5]=diff[0][5];

		diff[2][0]=-0.5*gauss[count].x;
		diff[2][1]=-0.5*gauss[count].y;
		diff[2][2]=-0.5*(1.0-gauss[count].x-gauss[count].y);
		diff[2][3]=-diff[2][0];
		diff[2][4]=-diff[2][1];
		diff[2][5]=-diff[2][2];
		//--------------------------------------------------
		//单元节点坐标矩阵
		double elenode[6][3];
		for(int i=0;i<6;i++){
			elenode[i][0]=elenodes_vec[i].x;
			elenode[i][1]=elenodes_vec[i].y;
			elenode[i][2]=elenodes_vec[i].z;
		}
		//--------------------------------------------------
		//J矩阵
		double Jmatrix[3][3];
		//以上两个矩阵的积
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++){
				Jmatrix[i][j]=0;
				for(int k=0;k<6;k++)
				Jmatrix[i][j]=Jmatrix[i][j] + diff[i][k]*elenode[k][j];
			}
			//--------------------------------------------------
			//求出J矩阵的行列式；
			double J_val=Jmatrix[0][0]*(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])
				-Jmatrix[0][1]*(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])
				+Jmatrix[0][2]*(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0]);
			//----------------------------------------------------
			//求出J矩阵的逆矩阵；
			double Jinverse[3][3];
			
			Jinverse[0][0]=(Jmatrix[1][1]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][1])/J_val;
			Jinverse[1][1]=(Jmatrix[0][0]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][0])/J_val;
			Jinverse[2][2]=(Jmatrix[0][0]*Jmatrix[1][1]-Jmatrix[0][1]*Jmatrix[1][0])/J_val;

			Jinverse[0][1]=-(Jmatrix[0][1]*Jmatrix[2][2]-Jmatrix[0][2]*Jmatrix[2][1])/J_val;
			Jinverse[0][2]=(Jmatrix[0][1]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][1])/J_val;

			Jinverse[1][0]=-(Jmatrix[1][0]*Jmatrix[2][2]-Jmatrix[1][2]*Jmatrix[2][0])/J_val;
			Jinverse[1][2]=-(Jmatrix[0][0]*Jmatrix[1][2]-Jmatrix[0][2]*Jmatrix[1][0])/J_val;

			Jinverse[2][0]=(Jmatrix[1][0]*Jmatrix[2][1]-Jmatrix[1][1]*Jmatrix[2][0])/J_val;
			Jinverse[2][1]=-(Jmatrix[0][0]*Jmatrix[2][1]-Jmatrix[0][1]*Jmatrix[2][0])/J_val;
			//-------------------------------------------------------
			//求出N对x,y的偏导；
			double diffxy[3][6];
			for(int i=0;i<3;i++)
				for(int j=0;j<6;j++){
					diffxy[i][j]=0;
					for(int k=0;k<3;k++)
						diffxy[i][j]=diffxy[i][j]+Jinverse[i][k]*diff[k][j];
				}
				//--------------------------------------------------------
				//求出B矩阵
				double B[6][18];
				for(int i=0;i<6;i++)
					for(int j=0;j<18;j++)
						B[i][j]=0;
				for(int i=0;i<6;i++){
					B[0][i*3+0]=diffxy[0][i];
					B[1][i*3+1]=diffxy[1][i];
					B[2][i*3+2]=diffxy[2][i];
					B[3][i*3+0]=diffxy[1][i];
					B[3][i*3+1]=diffxy[0][i];
					B[4][i*3+1]=diffxy[2][i];
					B[4][i*3+2]=diffxy[1][i];
					B[5][i*3+0]=diffxy[2][i];
					B[5][i*3+2]=diffxy[0][i];
				}
				//--------------------------------------------------------
				//计算B_trans;
				double B_trans[18][6];
				for(int i=0;i<18;i++)
					for(int j=0;j<6;j++)
						B_trans[i][j] = B[j][i];
				//--------------------------------------------------------------------------------------------
				//求出B_trans矩阵与ele_elas矩阵的乘积array1；
				double array1[18][6];
				for(int i=0;i<18;i++)
					for(int j=0;j<6;j++){
						array1[i][j]=0; 
						for(int k=0;k<6;k++)
							array1[i][j] = array1[i][j]+B_trans[i][k]*ele_elas[k][j];
					}
				//求出array1矩阵与B矩阵的乘积array2;
				double array2[18][18];
				for(int i=0;i<18;i++)
					for(int j=0;j<18;j++){
						array2[i][j]=0;
						for(int k=0;k<6;k++)
							array2[i][j] = array2[i][j]+array1[i][k]*B[k][j];
							tem2[i][j]=tem2[i][j]+array2[i][j]*J_val*0.5*wight[count];
					}
	}
	elestiff=tem2;
}
//-----------------------------------------------------------------------------------------------
//将单刚添加到总刚
int GloStiffMatrix::Add_to_gsmatrix(double* &total_matrix, int* &Iz, int* &Ig, const Element &element)
{
	int node_num = int(element.nodes_id.size());
	vector<int> enode;
	for(int i=0; i<node_num; i++)
	{
		enode.push_back(element.nodes_id[i]);
	}

	for( int i=0; i<node_num; i++ )
	{
		int ile = 3*i;
		int ni = enode[i];
		int nin = Iz[ni];
		int nint = 0;
		int nin1= 0;
		if( ni > 0 )
		{
			nin1 = Iz[ni-1];
			nint = nin - nin1;
		}

		for( int j=0; j<node_num; j++ )
		{
			int nj = enode[j];
			if( ni < nj ) continue;
			int jle = 3*j;
			int ilt = 0;
			if( ni > nj )
			{
				for( int k=0; k<nint; k++ )
				{
					if( Ig[nin1+k] != nj+1 ) continue;
					ilt = 6*ni + 9*(nin1 + k);
					break;
				}
			}
			else
			{
				ilt = 6*ni + 9*nin;
			}
			for( int k=0; k<3; k++ )
			{
				int kk = 0;
				int ilt1 = 0;

				if( ni == nj )
				{
					kk = k+1;
					ilt1 = ilt+(k+1)*k/2;
				}
				else 
				{
					kk = 3;
					ilt1 = ilt + 3*k;
				}

				for( int l=0; l<kk; l++ )
				{
					total_matrix[ilt1+l]=total_matrix[ilt1+l]+elestiff[ile+k][jle+l];    
				}      
			}
		}                     
	}

	return 1;
}


