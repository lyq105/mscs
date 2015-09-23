//===========================================================================
// PCMCell.h
// 颗粒增强复合材料单胞类头文件
// A class of Cell of Particle reinforced Composite Material with random
// 是CMCell类的派生类
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
		vector<EllipseSurface> elliSurfaces_vec;  //存储所有的椭球表面;

		int import_ellipsoids(string file);
		//读取输入文件一行，跳过注释行（以%开始）
		char* Get_Line(ifstream& input_stream, char* read_line);
		//读入信息一行，跳过注释行（以%开头）
		string Get_Line(ifstream &infile)const;
}; 

//---------------------------------------------------------------------------
#endif
//===========================================================================
