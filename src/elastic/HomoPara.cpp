//=============================================================================
//计算均匀化参数类的实现；
//HomoPara.cpp
//=============================================================================

#include"HomoPara.h"
//计算均匀化弹性矩阵;
HomoPara::HomoPara(){
	for(int i=0;i<6;i++){
	Homo_beta[i]=0;
	for(int j=0;j<6;j++)
    Homo_D[i][j]=0;
	}
}
//---------------------------------------------------------------------------------------------------------------------
int HomoPara::Generate_Homo(const vector<Element> &elements_vec,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const double& Unitcell_V,const int &a1,const int &m){
    //---------------------------------------------------
    int cn=Mapping(a1,m);
	//----------------------------------------------------
	//取高斯点；
	Gauss gau;
	gau.Generate_gauss();
	for(int num=0;num<int(elements_vec.size());num++){
		//计算单元均匀化D矩阵;
	    //计算Homo_D;
		Generate_ele_Dmatrixcn(elements_vec[num],nodes_vec,mats_vec,uvw,cn,gau.gauss,gau.wight);
		//---------------------------------------------------
        for(int i=0;i<6;i++)
	    Homo_D[i][cn]=Homo_D[i][cn]+ele_Homo_Dcn[i];
		//----------------------------------------------------
		//计算Homo_beta;
		Generate_ele_Homo_beta(elements_vec[num],mats_vec,cn);
		//----------------------------------------------------
		for(int i=0;i<6;i++)
		Homo_beta[i]=Homo_beta[i]+ele_Homo_beta[i];  
	}//end遍历单元
	    //---------------------------------------------------
	    //计算Homo_D;
	    for(int i=0;i<6;i++)
		Homo_D[i][cn] = Homo_D[i][cn]/Unitcell_V;	
	    //--------------------------------------------------  
return 1;
}
//-----------------------------------------------------------------------------------------------
//生成单元弹性矩阵
void HomoPara::Generate_ele_Dmatrixcn(const Element &e,const vector<Node> &nodes_vec,const vector<MatPro> &mats_vec,double * &uvw,const int &cn,const vector<Node> &gauss,const vector<double> &wight){        
	//-------------------------------------------------------------------
	//初始化
	for(int i=0;i<6;i++)
	   ele_Homo_Dcn[i]=0;
	//读出此单元的所有节点坐标;
	vector<Node> elenodes_vec;
	for(int i=0;i<int(e.nodes_id.size());i++)
		elenodes_vec.push_back(nodes_vec[e.nodes_id[i]]);
   //--------------------------------------------------------------------
	//读出此单元的弹性矩阵;
	double ele_elas[6][6];
	for(int i=0;i<6;i++)
		for(int j=0;j<6;j++)
			ele_elas[i][j] = mats_vec[e.mat].elas_matrix[i][j];
	//-------------------------------------------------------------------	
	//读出此单元所有节点的位移;
	vector<double> ele_a(int(e.nodes_id.size()*3));
	for(int i=0;i<int(e.nodes_id.size());i++)
		for(int j=0;j<3;j++)
			ele_a[3*i+j]=uvw[e.nodes_id[i]*3+j];
	//根据单元类型计算；
	//----------------------------------------------------------------------
	if(e.type==2&&int(e.nodes_id.size())==4)
		Tetrahedron_line(elenodes_vec,ele_elas,ele_a,cn);
	else if(e.type==3&&int(e.nodes_id.size())==6)
		Threeprism_line(elenodes_vec,ele_elas,ele_a,cn,gauss,wight);
	//-------------------------------------------------------------------------------
}
//-----------------------------------------------------------------------------------------------
//四面体为线性插值时生成单元荷载向量；
void HomoPara::Tetrahedron_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a, const int &cn){
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
	//--------------------------------------------------------
	//首先计算D与e_B的乘积;
	double array1[6][12];
	for(int i=0;i<6;i++)
		for(int j=0;j<12;j++){
			array1[i][j]=0;
		for(int k=0;k<6;k++)
			array1[i][j] = array1[i][j] +ele_elas[i][k]*B[k][j];
		}
	//--------------------------------------------------------
	//再计算其与ele_a的乘积;
	double array2[6];
	for(int i=0;i<6;i++){
		array2[i]=0;
	for(int j=0;j<12;j++)
		array2[i] = array2[i] + array1[i][j]*ele_a[j];
	}
	//--------------------------------------------------------
	//与D一列的和；
	for(int i=0;i<6;i++)
	ele_Homo_Dcn[i] = (ele_elas[i][cn]+array2[i])*volume;
}		
//------------------------------------------------------------------------------------------
void HomoPara::Threeprism_line(const vector<Node> &elenodes_vec,const double ele_elas[][6],const vector<double> &ele_a,const int &cn,const vector<Node> &gauss,const vector<double> &wight){  
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
				/*for(int i=0;i<6;i++){
					for(int j=0;j<18;j++)
						cout<<B[i][j]<<" ";
				    cout<<endl;
				}
				cout<<endl;
               */
   //首先计算D与e_B的乘积;
	double array1[6][18];
	for(int i=0;i<6;i++)
		for(int j=0;j<18;j++){
			array1[i][j]=0;
		for(int k=0;k<6;k++)
			array1[i][j] = array1[i][j] +ele_elas[i][k]*B[k][j];
		}
		/* for(int i=0;i<6;i++){
					for(int j=0;j<18;j++)
						cout<<array1[i][j]<<" ";
				    cout<<endl;
				}
		 cout<<endl;

		 for(int i=0;i<18;i++)
			 cout<<ele_a[i]<<" ";
		     cout<<endl<<endl;
			 */
	//--------------------------------------------------------
	//再计算其与ele_a的乘积;
	double array2[6];
	for(int i=0;i<6;i++){
		array2[i]=0;
	for(int j=0;j<18;j++)
		array2[i] = array2[i] + array1[i][j]*ele_a[j];	
	//cout<<array2[i]<<" ";
	}
	//cout<<endl<<endl;
	//--------------------------------------------------------
	//与D一列的和；
	for(int i=0;i<6;i++)
	ele_Homo_Dcn[i] =ele_Homo_Dcn[i]+(ele_elas[i][cn]+array2[i])*J_val*0.5*wight[count];
	} 
}
//---------------------------------------------------------------------------------------------------
//计算单元热弹性常数张量；
void HomoPara::Generate_ele_Homo_beta(const Element &e,const vector<MatPro> &mats_vec,const int &cn){
	//--------------------------------------------------------------------------------
	//读出此单元的热膨胀系数矩阵；
	double ele_CTE_matrix[3][3];
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			ele_CTE_matrix[i][j] = mats_vec[e.mat].CTE_matrix[i][j];     
	//转化此为一维向量；（a[1][1],a[2][2],a[3][3],a[1][2],a[2][3],a[3][1]）
	vector<double> alpha(6);
	alpha[0]=ele_CTE_matrix[0][0];
	alpha[1]=ele_CTE_matrix[1][1];
	alpha[2]=ele_CTE_matrix[2][2];
	alpha[3]=2.0*ele_CTE_matrix[0][1];
	alpha[4]=2.0*ele_CTE_matrix[1][2];
	alpha[5]=2.0*ele_CTE_matrix[2][0];	     
	//-----------------------------------------------------------------------------------------------
	//计算ele_Homo_beta;   
	for(int i=0;i<6;i++){
		ele_Homo_beta[i]=0;
		ele_Homo_beta[i] = ele_Homo_Dcn[i]*alpha[cn];	    
	}
}
//---------------------------------------------------------------------------------------------------------	
//生成热膨胀均匀化矩阵；
int HomoPara::Generate_Homo_alpha(const double &Unitcell_V){

	for(int i=0;i<6;i++)
	Homo_beta[i]=Homo_beta[i]/Unitcell_V;

	double vec[6];
	for(int i=0;i<3;i++)
		vec[i] = Homo_beta[i];
	for(int i=3;i<6;i++)
		vec[i] = 2*Homo_beta[i];
	MathMatrix matrix1(&Homo_D[0][0],6,6);
	MathMatrix	matrix2(vec,6,1);
	MathMatrix matrix3;
	matrix3=(matrix1.Inverse())*matrix2;
	Homo_alpha[0][0] = matrix3.element[0][0];
	Homo_alpha[1][1] = matrix3.element[1][0];
	Homo_alpha[2][2] = matrix3.element[2][0];
	Homo_alpha[0][1] = matrix3.element[3][0];
	Homo_alpha[1][0] = matrix3.element[3][0];
	Homo_alpha[1][2] = matrix3.element[4][0];
	Homo_alpha[2][1] = matrix3.element[4][0];
	Homo_alpha[0][2] = matrix3.element[5][0];
	Homo_alpha[2][0] = matrix3.element[5][0];	 
	return 1;
}
//-----------------------------------------------------
//生成工程常数；
int HomoPara::Generate_Homo_engconst(MatPro* homoMat){
	MathMatrix D(6,6);
	for(int i=0;i<6;i++)
		for(int j=i;j<6;j++){
			if(i==j)
		    D.element[i][j]=Homo_D[i][j];
			else{
			D.element[i][j]=(Homo_D[i][j]+Homo_D[j][i])/2;
			D.element[j][i]=D.element[i][j];
			}
	    }	
	MathMatrix S(6,6);
	S=D.Inverse();
	E11=1.0/S.element[0][0];
	E22=1.0/S.element[1][1];
	E33=1.0/S.element[2][2];
	Nu12=-S.element[0][1]*E11;
    Nu13=-S.element[0][2]*E11;
	Nu23=-S.element[1][2]*E22;
	G12=1.0/S.element[3][3];
	G23=1.0/S.element[4][4];
	G13=1.0/S.element[5][5];
	//设置均匀化后材料参数
	homoMat->set_ela_para(E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13);
	return 1;
}
//-----------------------------------------------------
//输出函数；
void HomoPara::print(int mode){
	if(mode==0){
		MathMatrix matrix1(&Homo_D[0][0],6,6);
		cout << "Homo_D:" << "均匀化弹性矩阵如下"<<endl;
		cout << matrix1 << endl;
		MathMatrix H_Alpha(&Homo_alpha[0][0],3,3);
		cout << "Homo_Alpha:" <<"均匀化热膨胀系数张量如下"<< endl;
		cout << H_Alpha << endl;
	}
	else if(mode==1){
		ofstream outfile1;
		outfile1.open("HomoPara.dat");
		outfile1 << "Homo_D:" << "均匀化弹性矩阵如下"<< endl;
		MathMatrix matrix1(&Homo_D[0][0],6,6);
		outfile1 <<  matrix1 << endl;
		MathMatrix H_Alpha(&Homo_alpha[0][0],3,3);
		outfile1 << "Homo_Alpha:"<<"均匀化热膨胀系数张量如下" << endl;
		outfile1 << H_Alpha << endl;
		outfile1.close();
	}
	else {
		//输出到文件；
		ofstream outfile1;
		outfile1.open("HomoPara.dat");
		outfile1 << "Homo_D:" << "均匀化弹性矩阵如下"<< endl;
		MathMatrix matrix1(&Homo_D[0][0],6,6);
		outfile1 << matrix1 << endl;
		MathMatrix H_Alpha(&Homo_alpha[0][0],3,3);
		outfile1 << "Homo_Alpha:"<<"均匀化热膨胀系数张量如下" << endl;
		outfile1 << H_Alpha << endl;
		//输出到屏幕；
		cout << "Homo_D:" << "均匀化弹性矩阵如下"<<endl;
		cout << matrix1 << endl;
		cout << "Homo_Alpha:" <<"均匀化热膨胀系数张量如下"<< endl;
		cout << H_Alpha << endl;
		if(mode==2){
			cout<<"杨氏模量，泊松比，剪切模量如下"<<endl;
			cout<<"E11="<<E11<<" "<<"E22="<<E22<<" "<<"E33="<<E33<<endl<<"Nu12="<<Nu12<<" "<<"Nu13="<<Nu13<<" "<<"Nu23="<<Nu23<<endl<<"G12="<<G12<<" "<<"G13="<<G13<<" "<<"G23="<<G23<<endl;
			outfile1<<"E11="<<E11<<" "<<"E22="<<E22<<" "<<"E33="<<E33<<endl<<"Nu12="<<Nu12<<" "<<"Nu13="<<Nu13<<" "<<"Nu23="<<Nu23<<endl<<"G12="<<G12<<" "<<"G13="<<G13<<" "<<"G23="<<G23<<endl;
			cout<<"alpha1="<<Homo_alpha[2][2]<<" "<<"alpah2="<<(Homo_alpha[0][0]+Homo_alpha[1][1])/2.0<<endl;
			outfile1<<"alpha1="<<Homo_alpha[2][2]<<" "<<"alpah2="<<(Homo_alpha[0][0]+Homo_alpha[1][1])/2.0<<endl;
		}
		outfile1.close();
	}

}

//-----------------------------------------------------
//输出样本数据；
void HomoPara::output_Datafile(const string data_file)const
{
	ofstream outfile1(data_file.c_str(),ios::app);
	outfile1<<"%杨氏模量：E11，E22，E33 ; 泊松比：Nu12，Nu13，Nu23 ; 剪切模量：G12，G13，G23 ; 热膨胀系数：alpha1，alpah2"<<endl;
	outfile1<<E11<<" "<<E22<<" "<<E33<<"   ";
	outfile1<<Nu12<<" "<<Nu13<<" "<<Nu23<<"   ";
	outfile1<<G12<<" "<<G13<<" "<<G23<<"   ";
	outfile1<<Homo_alpha[2][2]<<" "<<(Homo_alpha[0][0]+Homo_alpha[1][1])/2.0<<endl;
	outfile1.close();
}

//-------------------------------------------------------------------------------
//Dij和aijhk中下标对应关系；
//输入参数为aijhk前两个或后两个下标，输出参数是Dij的行下标或列下标；
int HomoPara::Mapping(int i,int j){
	int lc;//行下标或列下标；
	if(i==0&&j==0) lc=0;
	else if(i==1&&j==1) lc=1;
	else if(i==2&&j==2) lc=2;
	else if((i==0&&j==1)||(i==1&&j==0)) lc=3;
	else if((i==1&&j==2)||(i==2&&j==1)) lc=4;
	else if((i==0&&j==2)||(i==2&&j==0)) lc=5;
	return lc;
}//--------------------------------------------------------------
