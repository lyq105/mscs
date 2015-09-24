//===========================================================================
// PCMCell.cpp
// 颗粒增强复合材料单胞类成员函数
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
//从生成单胞
int PCMCell::Cell_generate(ifstream &infile,string data_file)
{
    int mod;
    istringstream istr(Get_Line(infile));		//读取建模信息
    istr >> mod;

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
        Get_Line(infile);	//不读取长方体区域的六条边
        Get_Line(infile);	//跳过读取椭球长中短轴信息
        Get_Line(infile);	//跳过读取体积百分比
    }
    else
    {
        hout << string("错误！产生椭球分布的模式不在考虑的情况内。(mod>2)") << endl;
    }

    //利用椭球文件生成椭球
    import_ellipsoids(ellips_file);

    return 1;
}

int PCMCell::Cell_generate(string input_file,string data_file)
{
    ifstream infile(input_file.c_str());
    Cell_generate(infile,data_file);
    infile.close();
}
//---------------------------------------------------------------------------
//从给定文件中读取椭球
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

    //读入单胞原点
    istrstream istr1(Get_Line(input_stream,read_line));
    istr1 >> origin_x >> origin_y >> origin_z;

    //读入单胞长宽高
    istrstream istr2(Get_Line(input_stream,read_line));
    istr2 >> clength >> cwidth >> cheight;

    //计算单胞体积
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
//读取输入文件一行，跳过注释行（以%开始）
char* PCMCell::Get_Line(ifstream& input_stream, char* read_line)
{
    //读入的一行char read_line[200]
    input_stream.getline(read_line,200);
    //跳过注释行
    while(!input_stream.eof() && read_line[0] == '%')
    {
        input_stream.getline(read_line,200);
    }

    return read_line;
}
//读入一行信息，并跳过注释行（以"%"开头）；
string PCMCell::Get_Line(ifstream &infile)const
{
    string s;
    //读入信息一行
    getline(infile,s);
    //跳过注释行
    while(!infile.eof() && s.substr(0,1)=="%")
        getline(infile,s);
    return s;
}
//===========================================================================
