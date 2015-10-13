//===========================================================================
// GloStiffMatrix.h
// ����ܸ�����ͷ�ļ�
// A Class of the Global Stiff Matrix
//===========================================================================
#ifndef GloStiffMatrixH
#define GloStiffMatrixH

#include"Fem.h"
#include"MatPro.h"
#include"Gauss.h"
#include"time.h"
//---------------------------------------------------------------------------
class GloStiffMatrix
{
public:
	//---------------------------------------------------------------------------
	//���캯��
	GloStiffMatrix( ){};
	//��������նȾ���
	int Gen_gsmatrix(double* &, int* &, int* &, const vector<Node> &, const vector<Element> &, const vector<MatPro> &);
	//���ɵ��գ�
	int Generate_elestiff(const Element &,const vector<Node> &,const vector<MatPro> &,const vector<Node> &,const vector<double> &);
	//������
	void Tetrahedron_line(const vector<Node> &,double ele_elas[][6]);
	//������
	void Threeprism_line(const vector<Node> &elenodes_vec,double ele_elas[][6],const vector<Node> &gauss,const vector<double> &wight);
	//��������ӵ��ܸ�
	int Add_to_gsmatrix(double* &, int* &, int* &, const Element &);
	vector<vector<double> > elestiff;
protected:
};

//---------------------------------------------------------------------------
#endif
//===========================================================================
