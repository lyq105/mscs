//===========================================
//MatPro.cpp
//�����������Ա����
//A Class of Material Property
//===========================================
#include"MatPro.h"

//---------------------------------------------
//���캯��
MatPro::MatPro(int itype){
	type_val = itype;
}
//-----------------------------------------------
//���ø���ͬ�Բ��ϲ���������ģ���Ͳ��ɱȣ�
void MatPro::set_ela_para(double E,double Mu){
	E_val=E;
	Mu_val=Mu;
	type_val=0;
}
//--------------------------------------------------
//���ø���ͬ�Բ��ϲ�����a1��a2������ϵ����
void MatPro::set_expan_para(double a1,double a2){
	cte1=a1;
	cte2=a2;
}
//---------------------------------------------------------------------
//���ú�۸���ͬ�Բ��ϲ�����
void MatPro::set_ela_para(double iE1,double iNu1,double iE2,double iNu2,double iG2){
	E1 = iE1; Nu1 = iNu1; E2 = iE2; Nu2 = iNu2;G2 = iG2; 
	type_val=1;
}
//---------------------------------------------------------------------
//���������������Բ��ϲ�����
void MatPro::set_ela_para(double iE11,double iE22,double iE33,double iNu12,double iNu23,double iNu13,double iG12,double iG23,double iG13){
	E11 = iE11; E22 = iE22; E33 = iE33;
	Nu12 = iNu12; Nu23 = iNu23;Nu13 = iNu13;
	G12 = iG12; G23 = iG23; G13 = iG13;
	type_val=2;
}
//----------------------------------------------------------------------
//����Dij���Ծ��󣻴˵��Ծ����Ӧ��Ӧ��e����Ϊ��e11,e22,e33,2*e12,2*e23,2*e31��ת��;
int MatPro::Generate_elas_matrix(){
	//���ݲ������ͣ�������Ӧ�Ĳ��ϵ��Ծ���
	if(type_val==0){//����ͬ�Բ��ϵĵ��Ծ�������
    E11 = E_val; E22 = E_val; E33 = E_val;
	Nu12 = Mu_val; Nu23 = Mu_val;Nu13 = Mu_val;
	G12 = 0.5*E_val/(1.0+Mu_val); G23 = G12; G13 = G12;
	}//-----------------------------------------------------------------------
	//����ͬ�Բ��ϵĵ��Ծ����������;
	else if(type_val==1){
	E11 = E1; E22 = E1; E33 = E2;
	Nu12 = Nu1; Nu23 = Nu2;Nu13 = Nu2;
	G12 = E1/(2.0+2.0*Nu1); G23 = G2; G13 = G2;
	}//��۸���ͬ�Բ��ϵ��Ծ���������ϣ�
	if(type_val<0||type_val>2)
		cout<<"�������͸�����(0-2)!"<<endl;
   //����������Ⱦ���S;
     MathMatrix S(6,6);
	 S.element[0][0]=1.0/E11;
     S.element[1][1]=1.0/E22;
	 S.element[2][2]=1.0/E33;
	 S.element[0][1]=-Nu12/E11;
	 S.element[0][2]=-Nu13/E11;
	 S.element[1][2]=-Nu23/E22;
	 for(int i=0;i<3;i++)
		 for(int j=0;j<i;j++)
			 S.element[i][j]=S.element[j][i];
	 S.element[3][3]=1.0/G12;
	 S.element[4][4]=1.0/G23;
	 S.element[5][5]=1.0/G13;

	 //�����棻
     MathMatrix D(6,6);
	 D=S.Inverse();
	 for(int i=0;i<6;i++)
		 for(int j=0;j<6;j++)
			 elas_matrix[i][j]=D.element[i][j];
	
	//------------------------------------------------------------------------
	return 1;
}//=======================================================================

//---------------------------------------------------------
//����������ϵ������
int MatPro::Generate_CTE_matrix(){
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++){
			if(i==j){
				if(i==2)
					CTE_matrix[i][j]=cte1;
				else
					CTE_matrix[i][j]=cte2;
			}
			else 
				CTE_matrix[i][j]=0;
		}

		return 1;
}

//-------------------------------------------------------------------------------
//��Dij���Ծ���ӳ��������֮��Ӧ��aijhk;
int MatPro::Generate_aijhk(){
	for( int i=0; i<3; i++ )
		for( int j=0; j<3; j++ )
			for( int h=0; h<3; h++ )
				for( int k=0; k<3; k++ ){
					int ln ; //D�����кţ�
					int cn ;
					ln=Mapping(i,j);
					cn=Mapping(h,k);
					aijhk[i][j][h][k]=elas_matrix[ln][cn];
				}
				return 1;
}//------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------------
//����aijhk���ĸ��±꣬�õ�������Ӧ���Ծ���Dij�е��Ǹ�ֵ��
double MatPro::Trans_aijhk2Dij(int i,int j,int h,int k){
	int ln=Mapping(i,j);
	int cn=Mapping(h,k);
	return elas_matrix[ln][cn];
}

//-----------------------------------------------------------------------------------------
void MatPro::print()const{
	//���Ծ����������
	cout<<"���Ծ���Dij"<<endl;
	for(int i=0;i<6;i++){
		for(int j=0;j<6;j++)
			cout<<elas_matrix[i][j]<<"  ";
		cout<<endl;}
	cout<<endl;
	//--------------------------------------------------------------
	//������ϵ������������ԣ�
	cout<<"������ϵ������"<<endl;
	for( int i=0;i<3;i++){
		for(int j=0;j<3;j++)
			cout<<CTE_matrix[i][j]<<"  ";
		cout<<endl;
	}
	cout<<endl;
}//----------------------------------------------------------------------------

//-------------------------------------------------------------------------------
//Dij��aijhk���±��Ӧ��ϵ��
//�������Ϊaijhkǰ������������±꣬���������Dij�����±�����±ꣻ
int MatPro::Mapping(int i,int j){
	int lc;//���±�����±ꣻ
	if(i==0&&j==0) lc=0;
	else if(i==1&&j==1) lc=1;
	else if(i==2&&j==2) lc=2;
	else if((i==0&&j==1)||(i==1&&j==0)) lc=3;
	else if((i==1&&j==2)||(i==2&&j==1)) lc=4;
	else if((i==0&&j==2)||(i==2&&j==0)) lc=5;
	return lc;
}//--------------------------------------------------------------
