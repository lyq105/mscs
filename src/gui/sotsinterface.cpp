#include "sotsinterface.h"

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkFeatureEdges.h>
#include <vtkProperty.h>
#include <vtkPropPicker.h>
#include <vtkCallbackCommand.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCleanPolyData.h>
#include <vtkAxesActor.h>
#include <vtkSmartPointer.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkTransform.h>
#include <meshpara.h>
#include <vtkPolyData.h>
#include <vtkSTLReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkUnstructuredGridVolumeRayCastMapper.h>
#include <QTextEdit>
#include "analysistype.h"
#include "materiel.h"
#include <sstream>
#include <QScrollBar>
#include <vtkCellData.h>
#include "setcelldialog.h"
#include <vtkBorderWidget.h>

#include <vtkSmartPointer.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkTextMapper.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkProperty.h>
#include <vtkTextProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkParametricFunctionSource.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkParametricEllipsoid.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyData.h>
#include <vtkStripper.h>
#include <vtkFeatureEdges.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkSphereSource.h>
#include <vtkBox.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkClipPolyData.h>
#include <vtkCubeSource.h>
#include <vtkLookupTable.h>
#include <vtkAssignAttribute.h>
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
vtkUnstructuredGrid *SOTSinterface::get_cell_data()
{
    return celldata;
}

int SOTSinterface::build_pmcell(string in_file, string data_file, int CNum)
{
    //-----------------------------------------------------------------------------------------------------------------------------------------
    cout << endl;
    std::cout << "*******************************************************************" << endl;
    std::cout << "-_- 开始计算......"<<endl;
    cout << "*******************************************************************" << endl;
    cout << endl;
    clock_t ct0,ct1;

    ifstream infile;
    infile.open(in_file.c_str());
    //跳过注释行（这样就可以跳过样本数信息）
    Get_Line(infile);


    //-----------------------------------------------------------------------------------------------------------------------------------------
    //单胞建模；
lable_pcm: ct0 = clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始单胞建模......"<<endl;
    PCMCell *PCM = new PCMCell;
    if(PCM->Cell_generate(infile,data_file)==0)        //若单胞生成失败，则重新读取输入文件
    {
        infile.close();
        infile.open(in_file.c_str());
        //跳过注释行（这样就可以跳过样本数信息）
        Get_Line(infile);
        delete PCM;
        goto lable_pcm;
    }
    ct1 = clock();
//    cout << "    建模总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
//    cout << "^_^ 单胞建模完毕！"<<endl<<endl;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //生成材料库信息
    ct0 = clock();
    cout << "======================================================" << endl;
    cout << "-_- 开始读取材料信息......"<<endl;
    Matbase *Matrial = new Matbase;
    Matrial->Generate_matbase(infile);
    //	Matrial->mats_vec[0].print();
    ct1 = clock();
//    cout << "    读取材料信息总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
//    cout << "^_^ 读取材料信息完毕！"<<endl<<endl;

    //-----------------------------------------------------------------------------------------------------------------------------------------
    //建模并生成网格信息
    ct0 = clock();
//    cout << "======================================================" << endl;
//    cout << "-_- 开始生成网格......"<<endl<<endl;
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
//    cout << "    生成网格总耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
//    cout << "^_^ 网格生成完毕！"<<endl<<endl;

    process_pmcell(*Mesh);

    //-----------------------------------------------------------------------------------------------------------------------------------------
//    //求解Na1和Na1a2
//    ct0=clock();
//    cout << "======================================================" << endl;
//    cout << "-_- 开始均匀化求解......"<<endl<<endl;
//    //求解均匀化系数
//    //读取数据，判断是否进行强度计算
//    istringstream instr(Get_Line(infile));
//    int elas_ana_only;
//    instr >> elas_ana_only;

//    HomoSolver *Hsolver = new HomoSolver(elas_ana_only,data_file);
//    Hsolver->Solve(Mesh->nodes_vec, Mesh->bnodes_vec, Mesh->elements_vec, Matrial->mats_vec,PCM->Unitcell_V);
//    //	Hsolver->Na_BinaryData(0,Mesh->nodes_vec,data_file, CNum);
//    ct1=clock();
//    cout << "    均匀化求解耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
//    cout << "^_^ 均匀化求解过程完毕！"<<endl<<endl;



//    delete Hsolver;
    delete PCM;
    delete Matrial;
    delete Mesh;
    //	delete Npara;

    return 1;
}


int SOTSinterface::build_pmcell(string ellips_para, string ellips_file,string data_file)
{
    //-----------------------------------------------------------------------------------------------------------------------------------------
    //单胞建模；
    cout << "======================================================" << endl;
    cout << "-_- 开始颗粒增强单胞建模......"<<endl;


    EllipseGen *ellip = new EllipseGen;
    while(ellip->uniell_generation(ellips_para,ellips_file,data_file)==0)
    {
     //   cout << "======================================================" << endl;
    //    cout << "-_- 单胞建模结束......"<<endl;
    //    return 0;
    }
    delete ellip;
    cout << "======================================================" << endl;
    cout << "-_- 颗粒增强单胞建模成功......"<<endl;
    return 1;

//    PCMCell *PCM = new PCMCell();
//    if(PCM->Cell_generate(infile,data_file)==0)        //若单胞生成失败，则重新读取输入文件
//    {
//        infile.close();
//        infile.open(in_file.c_str());
//        //跳过注释行（这样就可以跳过样本数信息）
//        Get_Line(infile);
//        delete PCM;
//        goto lable_pcm;
//    }
//    ct1 = clock();

//    process_pmcell(*Mesh);

//-----------------------------------------------------------------------------------------------------------------------------------------
//    //求解Na1和Na1a2
//    ct0=clock();
//    cout << "======================================================" << endl;
//    cout << "-_- 开始均匀化求解......"<<endl<<endl;
//    //求解均匀化系数
//    //读取数据，判断是否进行强度计算
//    istringstream instr(Get_Line(infile));
//    int elas_ana_only;
//    instr >> elas_ana_only;

//    HomoSolver *Hsolver = new HomoSolver(elas_ana_only,data_file);
//    Hsolver->Solve(Mesh->nodes_vec, Mesh->bnodes_vec, Mesh->elements_vec, Matrial->mats_vec,PCM->Unitcell_V);
//    //	Hsolver->Na_BinaryData(0,Mesh->nodes_vec,data_file, CNum);
//    ct1=clock();
//    cout << "    均匀化求解耗时"<<(double)(ct1-ct0)/CLOCKS_PER_SEC<<"秒。"<<endl;
//    cout << "^_^ 均匀化求解过程完毕！"<<endl<<endl;



//    delete Hsolver;
//    delete PCM;
//    delete Matrial;
//    delete Mesh;
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

/// 从输入文件直接计算
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

int SOTSinterface::process_pmcell(const Mesher &mesh)
{
    vtkIdType number_of_points, number_of_tetra;
    vtkSmartPointer<vtkIntArray> intValue = vtkSmartPointer<vtkIntArray>::New();
    intValue->SetNumberOfComponents(1);
    intValue->SetName("materialId");

 //   infile >> number_of_points;
    number_of_points = mesh.nodes_vec.size();

    vtkSmartPointer<vtkPoints> points
            = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(number_of_points);
    for (vtkIdType i = 0; i < number_of_points; i++)
    {
        mesh.nodes_vec[i].x;
        double x, y, z;
        x = mesh.nodes_vec[i].x;
        y = mesh.nodes_vec[i].y;
        z = mesh.nodes_vec[i].z;
        points->SetPoint(i, x, y, z);
    }

    number_of_tetra = mesh.eles_vec.size();
    vtkSmartPointer<vtkCellArray> cellArray
            = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType i = 0; i < number_of_tetra; i++)
    {
        vtkIdType a, b, c, d;
//        vtkIdType subdomainId;
//        infile >> subdomainId >> a >> b >> c >> d;
        //       std::cout << a << b << c << d << std::endl;

        a = mesh.eles_vec[i].nodesId[0];
        b = mesh.eles_vec[i].nodesId[1];
        c = mesh.eles_vec[i].nodesId[2];
        d = mesh.eles_vec[i].nodesId[3];
        vtkSmartPointer<vtkTetra> tetra =
                vtkSmartPointer<vtkTetra>::New();
        if(mesh.eles_vec[i].materialId == 0) continue;
        tetra->GetPointIds()->SetId(0, a);
        tetra->GetPointIds()->SetId(1, b);
        tetra->GetPointIds()->SetId(2, c);
        tetra->GetPointIds()->SetId(3, d);
        intValue->InsertNextValue(mesh.eles_vec[i].materialId);
        cellArray->InsertNextCell(tetra);
    }
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TETRA, cellArray);
    unstructuredGrid->GetCellData()->AddArray(intValue);
    celldata = vtkUnstructuredGrid::New();
    celldata->SetPoints(points);
    celldata->SetCells(VTK_TETRA, cellArray);
    celldata->GetCellData()->SetScalars(intValue);

//    // Write file
//    vtkSmartPointer<vtkUnstructuredGridWriter> writer =
//            vtkSmartPointer<vtkUnstructuredGridWriter>::New();
//    writer->SetFileName("/home/yzh/mesher.vtu");
//#if VTK_MAJOR_VERSION <= 5
//    writer->SetInput(unstructuredGrid);
//#else
//    writer->SetInputData(unstructuredGrid);
//#endif
//    writer->Write();
    return 1;
}
