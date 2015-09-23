//===========================================================================
// PCMCell.h
// ������ǿ���ϲ��ϵ�����ͷ�ļ�
// A class of Cell of Particle reinforced Composite Material with random
// ��CMCell���������
//===========================================================================

#ifndef PCMCELL_H
#define PCMCELL_H

#include "CMCell.h"
#include "EllipseSurface.h"
#include "EllipseGen.h"
#include "EllipseMade.h"
//---------------------------------------------------------------------------
class PCMCell : public CMCell
{
	public:
		PCMCell();
		PCMCell(double len, double wid, double hei);
		PCMCell(double len, double wid, double hei, char* ellipsoidFile);

		int Cell_generate(ifstream &,string data_file);

		int cell_generrate();
	private:
		vector<EllipseSurface> elliSurfaces_vec;  //�洢���е��������;

		int import_ellipsoids(string file);
		//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
		char* Get_Line(ifstream& input_stream, char* read_line);
		//������Ϣһ�У�����ע���У���%��ͷ��
		string Get_Line(ifstream &infile)const;
}; 

//---------------------------------------------------------------------------
#endif
//===========================================================================
