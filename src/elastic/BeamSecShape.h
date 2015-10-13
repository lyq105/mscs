//---------------------------------------------------------------------------

#ifndef BeamSecShapeH
#define BeamSecShapeH
//---------------------------------------------------------------------------
//定义梁截面
class BeamSecShape{
public:
        double rr;                                      //圆形截面的半径
        double hh,bb;                                   //方形截面的高和宽
        BeamSecShape( double r=1 );                       //圆形截面构造函数
        BeamSecShape( double h, double b );             //方形截面构造函数
        double Ix();                                    //截面绕x轴扭转惯性矩
        double Iy();                                    //截面绕y轴弯曲惯性矩
        double Iz();                                    //截面绕z轴弯曲惯性矩
        double area();                                  //截面面积
        int type(){ return shapeType; };                //截面形状类型
protected:
        int shapeType ;                                 //0:圆形，1：方形
};
//---------------------------------------------------------------------------
#endif
