//==================================================
//总荷载向量类的实现；
//==================================================

#include "Gloloaded_vector.h"
//-------------------------------------------------------------------------------------------------------------------------------------------------------
//生成a1总荷载向量；
int Gloloaded_vector::Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &m,double * &equright){
	//-----------------------------------------------------------
	//为之预留向量空间；
	//-----------------------------------------------------------
	//由alfa1和m，求出对应于弹性矩阵的列号；
	int cn = Mapping(alfa1,m);
	//取高斯点
	Gauss gau;
	gau.Generate_gauss();
	//-----------------------------------------------------------
	//循环开始；
	for(int ie=0;ie<int(elements_vec.size());ie++){
		Generate_eleloaded(elements_vec[ie],nodes_vec,mats_vec,cn,gau.gauss,gau.wight);
		Insert_eletoglo(elements_vec[ie],equright);
	}
	return 1;
}
//========================================================================================
//生成a1a2总荷载向量；
int Gloloaded_vector::Generate_gloloaded(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &alfa1,const int &alfa2,const int &m,double * &equright,const vector<vector<double> > &Na1_vec,const double Homo_D[][6]){
   //-----------------------------------------------------------
	//为homoelas赋值；
    for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			homoelas[i][j] = Homo_D[i][j];
	//-----------------------------------------------------------
	//由alfa1和m，求出对应于弹性矩阵的列号；
	int cn = Mapping(alfa1,m);
	//由alfa2,求出对应于弹性矩阵的行列好；
	 vector<int> ln(3);
	 ln[0] = Mapping(0,alfa2);
	 ln[1] = Mapping(1,alfa2);
	 ln[2] = Mapping(2,alfa2);
	//-----------------------------------------------------------
	//循环开始；
	for(int ie=0;ie<int(elements_vec.size());ie++){
		//取高斯点
	    Gauss gau;
		gau.Generate_gauss(elements_vec[ie].type);
		Generate_eleloaded1(elements_vec[ie],nodes_vec,mats_vec,Na1_vec,ln,cn,gau.gauss,gau.wight);
		Insert_eletoglo(elements_vec[ie],equright);
	}
	return 1;
    }

//========================================================================================
//将单元荷载向量插入总荷载向量；
void Gloloaded_vector::Insert_eletoglo(const Element &e,double* &equright){
	for(int i=0;i<int(e.nodes_id.size());i++)
		for(int j=0;j<3;j++)
			equright[e.nodes_id[i]*3+j] = equright[e.nodes_id[i]*3+j] +eleloaded_vec[i*3+j];		
}
//--------------------------------------------------------------------------------
//生成a1单元荷载向量；
int Gloloaded_vector::Generate_eleloaded(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const int &cn,const vector<Node> &gauss,const vector<double> &wight){
	//求出此单元的所有节点坐标
	vector<Node> elenodes_vec;
	for(int i=0;i<int(e.nodes_id.size());i++)
		elenodes_vec.push_back(nodes_vec[e.nodes_id[i]]);
	//-----------------------------------------------------------------------------
	//读取此单元所对应的材料弹性矩阵
	double ele_elas[6][6];
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			ele_elas[i][j] = mats_vec[e.mat].elas_matrix[i][j];
	//------------------------------------------------------------------------------
	if(e.type==2&&int(e.nodes_id.size())==4)
		Tetrahedron_line(elenodes_vec,ele_elas,cn);
	else if(e.type==3&&int(e.nodes_id.size())==6)
		Threeprism_line(elenodes_vec,ele_elas,cn,gauss,wight);
	//-------------------------------------------------------------------------------
	return 1;
}
//----------------------------------------------------------------------------------
//生成a1a2单元荷载向量；
int Gloloaded_vector::Generate_eleloaded1(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,const vector<vector<double> > &Na1_vec,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight){
   //求出此单元的所有节点坐标
	vector<Node> elenodes_vec;
	for(int i=0;i<int(e.nodes_id.size());i++)
		elenodes_vec.push_back(nodes_vec[e.nodes_id[i]]);
	//-----------------------------------------------------------------------------
	//读取此单元所对应的材料弹性矩阵
	double ele_elas[6][6];
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			ele_elas[i][j] = mats_vec[e.mat].elas_matrix[i][j];
	//-----------------------------------------------------------------------------
	//读取Na1单元节点的位移；
	int h,k;
	if(cn==0) {h=0;k=0;}
	else if(cn==1) {h=1;k=1;}
	else if(cn==2) {h=2;k=2;}
	else if(cn==3) {h=0;k=1;}
	else if(cn==4) {h=1;k=2;}
	else if(cn==5) {h=0;k=2;}
	//读出此单元所有节点的位移;
	vector<double> ele_a(int(e.nodes_id.size()*3));
	for(int i=0;i<int(e.nodes_id.size());i++)
		for(int j=0;j<3;j++)
			ele_a[3*i+j]=Na1_vec[e.nodes_id[i]][3*(3*h+k)+j];
	//------------------------------------------------------------------------------
	if(e.type==2&&int(e.nodes_id.size())==4)
		Tetrahedron_line1(elenodes_vec,ele_elas,ele_a,ln,cn,gauss,wight);
	else if(e.type==3&&int(e.nodes_id.size())==6)
		Threeprism_line1(elenodes_vec,ele_elas,ele_a,ln,cn,gauss,wight);
	//-------------------------------------------------------------------------------
	return 1;

}
//----------------------------------------------------------------------------------
//四面体为线性插值时生成a1单元荷载向量；
void Gloloaded_vector::Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const int &cn){

	vector<double> tem(12,0);
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
		c[i]=  (elenodes_vec[m].x*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].x)
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
	//求出B_trans矩阵与弹性矩阵D的cn列的乘积；		
	double array1[12];
	for(int i=0;i<12;i++){
		array1[i]=0;
		for(int j=0;j<6;j++)				
			array1[i]=array1[i]+B_trans[i][j]*ele_elas[j][cn];
		tem[i]=tem[i]-array1[i]*volume;
	}
	//生成单元荷载向量；
	eleloaded_vec=tem;		
}
//----------------------------------------------------------------------------------
//三角形为线性插值时生成a1a2单元荷载向量；
void Gloloaded_vector::Tetrahedron_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight){
    vector<double> tem(12,0);
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
		c[i]=  (elenodes_vec[m].x*elenodes_vec[l].z-elenodes_vec[m].z*elenodes_vec[l].x)
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
	double J_val = 6.0*volume;
	//--------------------------------------------------------------------
	//计算B矩阵
	double B[6][12];
	for(int i=0;i<6;i++)
		for(int j=0;j<12;j++)
			B[i][j]=0;

	for(int i=0;i<4;i++)
	{ 
		double coe;
		if(i%2==1) coe=-1.0/(6.0*volume);
		else coe=1.0/(6.0*volume);

		B[0][i*3+0]=coe*b[i];
		B[1][i*3+1]=coe*c[i];
		B[2][i*3+2]=coe*d[i];
		B[3][i*3+0]=coe*c[i];
		B[3][i*3+1]=coe*b[i];
		B[4][i*3+1]=coe*d[i];
		B[4][i*3+2]=coe*c[i];
		B[5][i*3+0]=coe*d[i];
		B[5][i*3+2]=coe*b[i];
	}
	//-------------------------------------------------------
	//计算B_trans;
	double B_trans[12][6];
	for(int i=0;i<12;i++)
		for(int j=0;j<6;j++)
			B_trans[i][j] = B[j][i];
	//--------------------------------------------------------------------------------------------
//循环高斯点计算积分；
for(int count=0;count<int(wight.size());count++){
	//计算形函数矩阵；
    double N[3][12];
	for(int i=0;i<3;i++)
		for(int j=0;j<12;j++)
			N[i][j]=0;
	for(int i=0;i<3;i++){
		N[i][i+0] = gauss[count].x;
		N[i][i+3] = gauss[count].y;
        N[i][i+6] = gauss[count].z;
		N[i][i+9] = 1.0 - gauss[count].x - gauss[count].y - gauss[count].z;
	}
    //=====================================================================================
	//N矩阵逆矩阵；
    double N_trans[12][3];
	for(int i=0;i<12;i++)
		for(int j=0;j<3;j++)
			N_trans[i][j]=N[j][i];
    //-------------------------------------------------------------------------------
    //求右端项；
    //第一块；
     vector<double> temarray1(12,0);
	 for(int i=0;i<12;i++)
	 for(int j=0;j<3;j++)
		 temarray1[i]=temarray1[i]+N_trans[i][j]*(homoelas[ln[j]][cn]-ele_elas[ln[j]][cn]);
   //第二块；
	double temarray3[12][6];
	for(int i=0;i<12;i++)
		for(int j=0;j<6;j++){
			temarray3[i][j]=0;
			for(int k=0;k<3;k++)
				temarray3[i][j]=temarray3[i][j]+N_trans[i][k]*ele_elas[ln[k]][j];
		}
    double temarray4[12][12];
		for(int i=0;i<12;i++)
		for(int j=0;j<12;j++){
			temarray4[i][j]=0;
			for(int k=0;k<6;k++)
				temarray4[i][j]=temarray4[i][j]+temarray3[i][k]*B[k][j];
		}
   vector<double> temarray2(12,0);
   for(int i=0;i<12;i++)
	   for(int j=0;j<12;j++)
	   temarray2[i]=temarray2[i]+temarray4[i][j]*ele_a[j];  
   //第三块
  double temarray5[12][3];
  for(int i=0;i<12;i++)
	 for(int j=0;j<3;j++){
		temarray5[i][j]=0;
		for(int k=0;k<6;k++)
		    temarray5[i][j]=temarray5[i][j]+B_trans[i][k]*ele_elas[k][ln[j]];
		}
  double temarray6[12][12];
  for(int i=0;i<12;i++)
	 for(int j=0;j<12;j++){
		temarray6[i][j]=0;
		for(int k=0;k<3;k++)
		    temarray6[i][j]=temarray6[i][j]+temarray5[i][k]*N[k][j];
		}
  vector<double> temarray7(12,0);
  for(int i=0;i<12;i++)
	 for(int j=0;j<12;j++)
	   temarray7[i]=temarray7[i]+temarray6[i][j]*ele_a[j];
//--------------------------------------------------
   for(int i=0;i<12;i++)
	   tem[i]=tem[i]-(temarray1[i]-temarray2[i]+temarray7[i])*wight[count]*volume;
   }
//--------------------------------------------------
//生成单元荷载向量；
   eleloaded_vec=tem;
   //for(int i=0;i<12;i++)
	//   cout<<eleloaded_vec[i]<<" ";
}
//--------------------------------------------------------------------------------------------------------------
//三棱柱生成a1单元荷载向量；
void Gloloaded_vector::Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const int &cn,const vector<Node> &gauss,const vector<double> &wight){
     vector<double> tem(18,0);
	//--------------------------------------------------	
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
		//--------------------------------------------------------
        //计算B_trans;
				double B_trans[18][6];
				for(int i=0;i<18;i++)
					for(int j=0;j<6;j++)
						B_trans[i][j] = B[j][i];
	  //--------------------------------------------------------------------------------------------
      //求出B_trans矩阵与弹性矩阵D的cn列的乘积；		
				double array1[18];
				for(int i=0;i<18;i++){
					array1[i]=0;
					for(int j=0;j<6;j++)				
						array1[i]=array1[i]+B_trans[i][j]*ele_elas[j][cn];
					tem[i]=tem[i]-array1[i]*J_val*0.5*wight[count];
				}
	}
	//生成单元荷载向量；
	eleloaded_vec=tem;    
}
//-----------------------------------------------------------------------
void Gloloaded_vector::Threeprism_line1(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const vector<int> &ln,const int &cn,const vector<Node> &gauss,const vector<double> &wight){
      vector<double> tem(18,0);
	//--------------------------------------------------	
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
      //计算N矩阵；
      double N[3][18];
      for(int i=0;i<3;i++)
		for(int j=0;j<18;j++)
				N[i][j]=0;
		for(int i=0;i<3;i++){
			N[i][i+0]=0.5*gauss[count].x*(1.0-gauss[count].z);
			N[i][i+3]=0.5*gauss[count].y*(1.0-gauss[count].z);
			N[i][i+6]=0.5*(1.0-gauss[count].x-gauss[count].y)*(1.0-gauss[count].z);
			N[i][i+9]=0.5*gauss[count].x*(1.0+gauss[count].z);
			N[i][i+12]=0.5*gauss[count].y*(1.0+gauss[count].z);
			N[i][i+6]=0.5*(1.0-gauss[count].x-gauss[count].y)*(1.0+gauss[count].z);
		}
	 //--------------------------------------------------------------------------------
	 //求N矩阵逆矩阵；
		double N_trans[18][3];
		for(int i=0;i<18;i++)
			for(int j=0;j<3;j++)
				N_trans[i][j]=N[j][i];
     //---------------------------------------------------------------------------------
     //求右端项；
	 //第一块；
     vector<double> temarray1(18,0);
	 for(int i=0;i<18;i++)
	 for(int j=0;j<3;j++)
		 temarray1[i]=temarray1[i]+N_trans[i][j]*(homoelas[ln[j]][cn]-ele_elas[ln[j]][cn]);
    //第二块；
	double temarray3[18][6];
	for(int i=0;i<18;i++)
		for(int j=0;j<6;j++){
			temarray3[i][j]=0;
			for(int k=0;k<3;k++)
				temarray3[i][j]=temarray3[i][j]+N_trans[i][k]*ele_elas[ln[k]][j];
		}
    double temarray4[18][18];
		for(int i=0;i<18;i++)
		for(int j=0;j<18;j++){
			temarray4[i][j]=0;
			for(int k=0;k<6;k++)
				temarray4[i][j]=temarray4[i][j]+temarray3[i][k]*B[k][j];
		}
   vector<double> temarray2(18,0);
   for(int i=0;i<18;i++)
	   for(int j=0;j<18;j++)
	   temarray2[i]=temarray2[i]+temarray4[i][j]*ele_a[j];
  //第三块
  double temarray5[18][3];
  for(int i=0;i<18;i++)
	 for(int j=0;j<3;j++){
		temarray5[i][j]=0;
		for(int k=0;k<6;k++)
		    temarray5[i][j]=temarray5[i][j]+B_trans[i][k]*ele_elas[k][ln[j]];
		}
  double temarray6[18][18];
  for(int i=0;i<18;i++)
	 for(int j=0;j<18;j++){
		temarray6[i][j]=0;
		for(int k=0;k<3;k++)
		    temarray6[i][j]=temarray6[i][j]+temarray5[i][k]*N[k][j];
		}
  vector<double> temarray7(18,0);
  for(int i=0;i<18;i++)
	 for(int j=0;j<18;j++)
	   temarray7[i]=temarray7[i]+temarray6[i][j]*ele_a[j];
//--------------------------------------------------------------------------------
   for(int i=0;i<18;i++)
	   tem[i]=tem[i]-(temarray1[i]-temarray2[i]+temarray7[i])*J_val*0.5*wight[count];
   }
//--------------------------------------------------
  //生成单元荷载向量
   eleloaded_vec=tem;
}
//------------------------------------------------------------------------
//四阶丈量aijhk的后两个坐标与弹性矩阵D列号的对应关系；
int Gloloaded_vector::Mapping(const int &i,const int &j)const{
	int lc;//列下标；
	if(i==0&&j==0) lc=0;
	else if(i==1&&j==1) lc=1;
	else if(i==2&&j==2) lc=2;
	else if((i==0&&j==1)||(i==1&&j==0)) lc=3;
	else if((i==1&&j==2)||(i==2&&j==1)) lc=4;
	else if((i==0&&j==2)||(i==2&&j==0)) lc=5;
	return lc;

}

//-------------------------------------------------------------------
