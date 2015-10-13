#include "sotsinterface.h"

//namespace elastic{
#include "Mscm.h"
//}
using namespace std;
SOTSinterface::SOTSinterface()
{
    celldata = NULL;
    homodata = NULL;
}

SOTSinterface::~SOTSinterface()
{

}
void SOTSinterface::set_cell_mesh()
{
}
//void SOTSinterface::
void SOTSinterface::set_material(int mat_index)
{

}
void SOTSinterface::get_material(int mat_index)
{

}
void SOTSinterface::write_cell_solution(std::string filename)
{

}
void SOTSinterface::write_homo_solution(std::string filename)
{

}
void SOTSinterface::set_pmcell_info(PMCell_Info& mycell_info)
{
    //cell_info.coating = mycell_info.coating;

}
vtkSmartPointer<vtkUnstructuredGrid> SOTSinterface::get_cell_data()
{
    return celldata;
}

int SOTSinterface::Begin(string in_file, ifstream &infile, string data_file, int CNum)
{
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //将执行过程状况输出在文件
    //   clock_t ct_begin,ct_end;
    //   ct_begin = clock();
    cout << endl;
    cout << "*******************************************************************" << endl;
    cout << "-_- 开始计算......"<<endl;
    cout << "*******************************************************************" << endl;
    cout << endl;
    clock_t ct0,ct1;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //单胞建模；
lable_pcm: ct0 = clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始单胞建模......"<<endl;
    PCMCell *PCM = new PCMCell;
    if(PCM->Cell_generate(infile,data_file)==0)        //参数0：用于艳生成的椭球，1：用韩非生成的椭球，2表示读取数据
    {
        infile.close();
        infile.open(in_file.c_str());
        //跳过注释行（这样就可以跳过样本数信息）
        Get_Line(infile);
        delete PCM;
        goto lable_pcm;
    }
    ct1 = clock();
    cout << "    建模总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
    cout << "^_^ 单胞建模完毕！"<<endl<<endl;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //生成材料库信息
    ct0 = clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始读取材料信息......"<<endl;
    Matbase *Matrial = new Matbase;
    Matrial->Generate_matbase(infile);
    //	Matrial->mats_vec[0].print();
    ct1 = clock();
    cout << "    读取材料信息总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
    cout << "^_^ 读取材料信息完毕！"<<endl<<endl;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //建模并生成网格信息
    ct0 = clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始生成网格......"<<endl<<endl;
    Mesher *Mesh = new Mesher(PCM);
    if( Mesh->Mesh_generate(infile) == 0)
    {
        infile.close();
        infile.open(in_file.c_str());
        //跳过注释行（这样就可以跳过样本数信息）
        Get_Line(infile);
        delete PCM;
        delete Matrial;
        goto lable_pcm;
    }
    //	Mesh->Mesh_BinaryData(0,data_file, CNum);
    //	Mesh->Mesh_data(0);
    ct1 = clock();
    cout << "    生成网格总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
    cout << "^_^ 网格生成完毕！"<<endl<<endl;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //求解Na1和Na1a2
    ct0=clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始均匀化求解......"<<endl<<endl;
    //求解均匀化系数
    //读取数据，判断是否进行强度计算
    istringstream instr(Get_Line(infile));
    int elas_ana_only;
    instr >> elas_ana_only;

    HomoSolver *Hsolver = new HomoSolver(elas_ana_only,data_file);
    Hsolver->Solve(Mesh->nodes_vec, Mesh->bnodes_vec, Mesh->elements_vec, Matrial->mats_vec,PCM->Unitcell_V);
    //	Hsolver->Na_BinaryData(0,Mesh->nodes_vec,data_file, CNum);
    ct1=clock();
    cout << "    均匀化求解耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
    cout << "^_^ 均匀化求解过程完毕！"<<endl<<endl;



    delete Hsolver;
    delete PCM;
    delete Matrial;
    delete Mesh;
    //	delete Npara;

    return 1;
}
//读入一行信息，并跳过注释行（以"%"开头）；
string SOTSinterface::Get_Line(ifstream &infile)const
{
    string s;
    //读入信息一行
    getline(infile,s);
    //跳过注释行
    while(!infile.eof() && s.substr(0,1)=="%")
        getline(infile,s);
    return s;
}

int SOTSinterface::solve(string in_file,string out_file,string data_file)
{
    ofstream datafile(data_file.c_str());
    datafile.close();
    //开始计算

    //确定输出文件
    if( out_file.size() > 0 ) open_deffo_stream( (char*)out_file.c_str() );
    //打开文件
    ifstream infile;
    infile.open(in_file.c_str());
    if(!infile) { hout << "不能打开输入文件: "  << in_file << endl;  return 0; }

    string s;		//读入信息一行
    getline(infile,s);	//跳过注释行
    while(!infile.eof() && s.substr(0,1)=="%")	getline(infile,s);

    istringstream istr(s);		//读取计算样本数信息
    int sample_num;
    istr >> sample_num;
    infile.close();
    for(int i=1;i<=sample_num; i++)
    {
        cout <<"********************************************"<<endl;
        cout << "sample_num="<<i<<endl;
        cout <<"********************************************"<<endl;
        //打开样本结果数据文件
        ofstream out(data_file.c_str(),ios::app);
        out <<"%---------------------------------------------------------------" << endl;
        out << i <<" 个样本"<< endl;
        out.close();
        //重新打开文件
        ifstream infile;
        infile.open(in_file.c_str());
        //跳过注释行（这样就可以跳过样本数信息）
        getline(infile,s);
        while(!infile.eof() && s.substr(0,1)=="%")	getline(infile,s);
        //计算样本
        Mscm *Compute =  new Mscm;
        Compute->Begin(in_file, infile,data_file, i);
        delete Compute;
        //关闭文件
        infile.close();
    }
    //关闭输出文件
    close_deffo_stream();
    return 1;
}
