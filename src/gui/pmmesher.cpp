#include "pmmesher.h"

PMMesher::PMMesher(double ms,int c, int ln,double tr,CMCell* cell)
    :Mesher(cell)
{
    mesh_size = ms;
    coat= c;
    layer_num = ln;
    thick_ratio = tr;
}

int PMMesher::generate_mesh()
{
//    double mesh_size;
//    int coat, layer_num;
//    double thick_ratio;

    if(mesh_tet_map(mesh_size) == 0 )
    {
        cout << "~_~ 剖分失败！可能原因：增强项重叠或者最小体积单元与最大体积单元相差太大!" << endl;
        cout << "****************************************************************************************" <<endl;
        return 0;
    }

    //输出生成边界层前的tecplot网格数据
    //if (export_tecplot_data_before_coating() == 0)
    //{
    //	cout << "~_~ 输出生成边界层前的tecplot网格数据失败！" << endl;
    //	cout << "*************************************************" << endl;
    //	return 0;
    //}

    if (coat==1)
    {
        //经过确认materialId ==0 是基体单元，materialId ==1 是椭球单元
        clock_t ct0,ct1;
        ct0 = clock();
        cout << "-_- 开始生成界面层网格... " << endl;
        if(gen_coating_mesh(thick_ratio, layer_num) == 0)
        {
            cout << "~_~ 生成界面层网格失败！" << endl;
            cout << "*******************************************" <<endl;
            return 0;
        }

        cout << "    输出tecplot网格数据....." << endl;
        if (export_tecplot_data() == 0)
        {
            cout << "~_~ 输出tecplot网格数据失败！" << endl;
            cout << "*******************************************" << endl;
            return 0;
        }

        //cout<<"   输出三维空间可视化Ensight网格数据"
        //if (export_Ensight_data() == 0)
        //{
        //	cout << "~_~ 输出三维空间可视化Ensight网格数据失败！" << endl;
        //	cout << "*******************************************" << endl;
        //	return 0;
        //}

        ct1 = clock();
        cout << "    生成界面层耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒" << endl;
        cout << "^_^ 界面层网格生成完毕!" << endl<<endl;
    }
    else
    {
        //在不生成边界层的情况下，将四面体单元类型倒为混合单元类型
        clock_t ct0,ct1;
        ct0 = clock();
        cout << "-_- 将四面体网格转换...... " << endl;
        if (change_elements_to_blend() == 0)
        {
            cout << " 将四面体单元类型倒为混合单元类型操作失败！" << endl;
            cout << "*******************************************" << endl;
            return 0;
        }
        ct1 = clock();
        cout << "    转换四面体网格耗时：" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "秒" << endl;
        cout << "^_^ 转换四面体网格完毕!" << endl << endl;
    }

    //确定所有节点的相关单元信息
 //   deter_nodes_relative_eles();

    return 1;
}

