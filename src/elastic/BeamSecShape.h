//---------------------------------------------------------------------------

#ifndef BeamSecShapeH
#define BeamSecShapeH
//---------------------------------------------------------------------------
//����������
class BeamSecShape{
public:
        double rr;                                      //Բ�ν���İ뾶
        double hh,bb;                                   //���ν���ĸߺͿ�
        BeamSecShape( double r=1 );                       //Բ�ν��湹�캯��
        BeamSecShape( double h, double b );             //���ν��湹�캯��
        double Ix();                                    //������x��Ťת���Ծ�
        double Iy();                                    //������y���������Ծ�
        double Iz();                                    //������z���������Ծ�
        double area();                                  //�������
        int type(){ return shapeType; };                //������״����
protected:
        int shapeType ;                                 //0:Բ�Σ�1������
};
//---------------------------------------------------------------------------
#endif
