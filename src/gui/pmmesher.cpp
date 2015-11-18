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
//         cout << "~_~ �ʷ�ʧ�ܣ�����ԭ����ǿ���ص�������С�����Ԫ����������Ԫ���̫��!" << endl;
//         cout << "****************************************************************************************" <<endl;
        return 0;
    }

    //������ɱ߽��ǰ��tecplot��������
    //if (export_tecplot_data_before_coating() == 0)
    //{
    //	cout << "~_~ ������ɱ߽��ǰ��tecplot��������ʧ�ܣ�" << endl;
    //	cout << "*************************************************" << endl;
    //	return 0;
    //}

    if (coat==1)
    {
        //����ȷ��materialId ==0 �ǻ��嵥Ԫ��materialId ==1 ������Ԫ
        clock_t ct0,ct1;
        ct0 = clock();
//         cout << "-_- ��ʼ���ɽ��������... " << endl;
        if(gen_coating_mesh(thick_ratio, layer_num) == 0)
        {
//            cout << "~_~ ���ɽ��������ʧ�ܣ�" << endl;
//            cout << "*******************************************" <<endl;
            return 0;
        }

  
        //cout<<"   �����ά�ռ���ӻ�Ensight��������"
        //if (export_Ensight_data() == 0)
        //{
        //	cout << "~_~ �����ά�ռ���ӻ�Ensight��������ʧ�ܣ�" << endl;
        //	cout << "*******************************************" << endl;
        //	return 0;
        //}

        ct1 = clock();
 //       cout << "    ���ɽ�����ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��" << endl;
//        cout << "����������������!" << endl<<endl;
    }
    else
    {
        //�ڲ����ɱ߽�������£��������嵥Ԫ���͵�Ϊ��ϵ�Ԫ����
        clock_t ct0,ct1;
        ct0 = clock();
//        cout << "-_- ������������ת��...... " << endl;
        if (change_elements_to_blend() == 0)
        {
//            cout << " �������嵥Ԫ���͵�Ϊ��ϵ�Ԫ���Ͳ���ʧ�ܣ�" << endl;
//            cout << "*******************************************" << endl;
            return 0;
        }
        ct1 = clock();
//        cout << "    ת�������������ʱ��" << (double)(ct1 - ct0)/CLOCKS_PER_SEC << "��" << endl;
//        cout << "^_^ ת���������������!" << endl << endl;
    }

    //ȷ�����нڵ����ص�Ԫ��Ϣ
 //   deter_nodes_relative_eles();

    return 1;
}

