//===========================================================================
// PCMCell.cpp
// ������ǿ���ϲ��ϵ������Ա����
// Member Functions in a class of Cell of Particle reinforced Composite Material with random
//===========================================================================

#include "PCMCell.h"

//---------------------------------------------------------------------------
PCMCell::PCMCell()
        :CMCell()
{      
}
//---------------------------------------------------------------------------
PCMCell::PCMCell(double len, double wid, double hei)
        :CMCell(len, wid, hei)
{      
}
//---------------------------------------------------------------------------
PCMCell::PCMCell(double len, double wid, double hei, char* ellipsoidFile)
        :CMCell(len, wid, hei)
{
	import_ellipsoids(ellipsoidFile);
}
//---------------------------------------------------------------------------
//�����ɵ���
int PCMCell::Cell_generate(ifstream &infile,string data_file)
{
	int mod;
	istringstream istr(Get_Line(infile));		//��ȡ��ģ��Ϣ
	istr >> mod;
cout << mod << endl;
    size_t en = data_file.rfind('.');
    string title = data_file.substr(0,en);
	string ellips_file = title+"OutEllipDat.dat";
	
	if (mod==0)
	{
		EllipseGen *ellip = new EllipseGen;
		string ellips_para = title+"EllipPara.dat";
		if(ellip->uniell_generation(ellips_para,ellips_file,data_file)==0)
		{
			return 0;
		}
		delete ellip;
	}
	else if (mod==1)
	{
		EllipseMade *ellip = new EllipseMade;
		ellip->ellip_generation(infile,ellips_file,data_file,1);
		delete ellip;
	}
	else if (mod==2)
	{
		Get_Line(infile);	//����ȡ�����������������
		Get_Line(infile);	//������ȡ�����ж�����Ϣ
		Get_Line(infile);	//������ȡ����ٷֱ�
	}
	else
	{
		hout << "���󣡲�������ֲ���ģʽ���ڿ��ǵ�����ڡ�(mod>2)" << endl;
	}

	//���������ļ���������
	import_ellipsoids(ellips_file);

	return 1;
}

//---------------------------------------------------------------------------
//�Ӹ����ļ��ж�ȡ����
int PCMCell::import_ellipsoids(string file)
{
	ifstream input_stream;
	input_stream.open( file.c_str(), ios::in );
	if( !input_stream )
	{
		hout << "Can not open the file: " << file << endl;
		return 0;
	};

	int ellipsoid_count;
	double ellip_ratio;
	char read_line[250];
	istrstream istr(Get_Line(input_stream,read_line));
	istr >> ellipsoid_count >> ellip_ratio;

	//���뵥��ԭ��
	istrstream istr1(Get_Line(input_stream,read_line));
	istr1 >> origin_x >> origin_y >> origin_z;

	//���뵥�������
	istrstream istr2(Get_Line(input_stream,read_line));
	istr2 >> clength >> cwidth >> cheight;

	//���㵥�����
	Unitcell_V = clength*cwidth*cheight;

	elliSurfaces_vec.clear();
	elliSurfaces_vec.reserve(ellipsoid_count);

	surfaces_vec.clear();
	surfaces_vec.reserve(ellipsoid_count);

	for( int i=0; i<ellipsoid_count; i++ )
	{
		double x0,y0,z0,a,b,c,alpha1,alpha2,alpha3,beta1,beta2,beta3;
		double gama1,gama2,gama3;
		istrstream istr1(Get_Line(input_stream,read_line));
		istr1	>> x0 >> y0 >> z0
				>> a >> b >> c
				>> alpha1 >> alpha2 >> alpha3
				>> beta1 >> beta2 >> beta3
				>> gama1 >> gama2 >> gama3 ;

		EllipseSurface esur(	x0,y0,z0,a,b,c,alpha1,alpha2,alpha3,
										beta1,beta2,beta3,gama1,gama2,gama3	);
		esur.material_num = 1;
		elliSurfaces_vec.push_back(esur); 
        surfaces_vec.push_back(&elliSurfaces_vec[i]);
	};
	return 1;
}
//---------------------------------------------------------------------------
//��ȡ�����ļ�һ�У�����ע���У���%��ʼ��
char* PCMCell::Get_Line(ifstream& input_stream, char* read_line)
{
	//�����һ��char read_line[200]     
	input_stream.getline(read_line,200);
	//����ע����
	while(!input_stream.eof() && read_line[0] == '%')
	{
		input_stream.getline(read_line,200);
	}

	return read_line;
}
//����һ����Ϣ��������ע���У���"%"��ͷ����
string PCMCell::Get_Line(ifstream &infile)const
{
	string s;
	//������Ϣһ��
	getline(infile,s);
	//����ע����     
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
