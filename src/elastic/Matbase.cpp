//===========================================
//Matbase.cpp
//���Ͽ����ʵ��
//A Class of Material database
//===========================================
#include "Matbase.h"

//======================================================================
int Matbase::Generate_matbase(ifstream &infile)
{
	//���봦���������ͣ�
	int  mats_num = 0;//��������
	istringstream in1(Get_Line(infile));
	in1>>mats_num;
	//===================================================
	if(mats_num<1)
	{
		cout<<"�ĵ��ϲ����������ˣ�"<<endl;
		exit(0);
	}
	else
	{
		mats_vec.reserve(mats_num);
		for(int i=1;i<=mats_num;i++){
			int mat_type = -1;
			istringstream in2(Get_Line(infile));
			in2>>mat_type;
			MatPro mat;
			double a1=0.0;
			double a2=0.0;
			if(mat_type==0){//����ͬ��
				double E=0.0;
				double Mu=0.0;
                double s1=0.0;
				in2>>E>>Mu>>a1>>a2>>s1;
				mat.set_ela_para(E,Mu);
                mat.set_str_para(s1);
			}
			else if(mat_type==1){//��۸���ͬ��
				double iE1,iNu1,iE2,iNu2,iG2,s1,s2;
				in2>>iE1>>iNu1>>iE2>>iNu2>>iG2>>a1>>a2>>s1>>s2;
				mat.set_ela_para(iE1,iNu1,iE2,iNu2,iG2);
				mat.set_str_para(s1, s2);
			}
			else if(mat_type==2){//������������
				double iE11,iE22,iE33,iNu12,iNu23,iNu13,iG12,iG23,iG13,s1,s2,s3;
				in2>>iE11>>iE22>>iE33>>iNu12>>iNu23>>iNu13>>iG12>>iG23>>iG13>>a1>>a2>>s1>>s2>>s3;
				mat.set_ela_para(iE11,iE22,iE33,iNu12,iNu23,iNu13,iG12,iG23,iG13);
				mat.set_str_para(s1, s2, s3);
			}
			else{
				cout<<"�������͸����ˣ�0-2������ �������ȸ��Ĳ��ϵĸ�����"<<endl;
				exit(0);
			}

			mat.set_expan_para(a1,a2);
			mat.Generate_elas_matrix();
			mat.Generate_CTE_matrix();
			mat.Generate_aijhk();
			mats_vec.push_back(mat);
		}
	}
	return 1;
}


//======================================================================
//����һ����Ϣ��������ע���У���"%"��ͷ����
string Matbase::Get_Line(ifstream &infile)const{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}//===============================================================
