//===========================================
//MatPro.h
//����������ͷ�ļ�
//A Class of Material Property
//===========================================

#ifndef MATPRO_H
#define MATPRO_H

#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<cmath>
#include"MathMatrix.h"
using namespace std;
//--------------------------------------------
class MatPro{
public:	
	int type_val;   //�������ͣ�Ŀǰ����3�֣�
	//					0������ͬ�ԣ�
	//					1����۸���ͬ��
	//					2��������������
	double E_val,Mu_val;		//����ͬ�Բ��ϵĵ���ģ���Ͳ��ɱ�
	double cte1,cte2;				//������ϵ����Coefficient of Thermal Expansion 
	double E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13;//�����������Եĵ���ģ�������ɱȣ�����ģ����
	double E1,Nu1,E2,Nu2,G2;	//��۸���ͬ�ԵĲ��ϵĵ���ģ�������ɱȣ�����ģ����1Ϊ����2Ϊ���
	double elas_matrix[6][6];		//Dij ���Ծ���
	double CTE_matrix[3][3];		//������ϵ������ 
	double aijhk[3][3][3][3];			//Dij ��Ӧ���Ľ�����
	//����20070103����
	double strength;						//ǿ�ȼ��ޣ�����ͬ�ԣ�
	double lon_str, tra_str;				//����ͺ���ǿ�ȼ��ޣ���۸���ͬ�ԣ�
	double str_11, str_22, str_33;	//���������ǿ�ȼ��ޣ��������ԣ�


	//���캯����
	MatPro(int itype=0);
	//---------------------------------------------
	//��Ա������
	//�����ĵ�����ȷ������ͬ�Բ��ϵ���ģ���Ͳ��ɱ�
	void set_ela_para(double E,double Mu);
	//�����ĵ�����ȷ��cte1,cte2;
	void set_expan_para(double a1,double a2);
	//�����ĵ�����ȷ�������������Բ��ϵ�E11,E22,E33,Nu12,Nu23,Nu13,G12,G23,G13��
	void set_ela_para(double iE11,double iE22,double iE33,double iNu12,double iNu23,double iNu13,double iG12,double iG23,double iG13);
	//�����ĵ�����ȷ����۸���ͬ�Բ��ϵ�
	void set_ela_para(double iE1,double iNu1,double iE2,double iNu2,double iG2);
	//����Dij���Ծ��󣬴˵��Ծ����Ӧ��Ӧ������Ϊ��e11,e22,e33,2*e12,2*e23,2*e31��ת��
	int Generate_elas_matrix();
	//����������ϵ������
	int Generate_CTE_matrix();
	//��Dij���Ծ���ӳ����֮��Ӧ��aijhk;
	int Generate_aijhk();
	//����aijhk���ĸ��±꣬�õ�������Ӧ���Ծ���Dij�е��Ǹ�ֵ��
	double Trans_aijhk2Dij(int ,int ,int ,int );
	//������Ծ���Dij��aijhk�Ľ�������
	//����20070103����
    void set_str_para(double s){ strength = s;}	//����ǿ�Ȳ���������ͬ�ԣ� 
    void set_str_para(double ls, double ts)			//����ǿ�Ȳ�������۸���ͬ�ԣ�
	{
		lon_str = ls;
		tra_str = ts;
    }                                         
    void set_str_para(double s1, double s2, double s3)	//����ǿ�Ȳ������������ԣ�
	{
		str_11 = s1 ;
		str_22 = s2 ;
		str_33 = s3 ;           
	}
	void print()const;
protected:
	//Dij��aijhk���±��Ӧ��ϵ��
	//�������Ϊaijhkǰ������������±꣬���������Dij�����±�����±ꣻ
	int Mapping(int i,int j);
};//------------------------------------------------

#endif
//==============================================

