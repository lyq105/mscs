//---------------------------------------------------------------------------
#include "BeamSecShape.h"
#define PI 3.141592654
//---------------------------------------------------------------------------
BeamSecShape::BeamSecShape(double r){
        rr        = r;
        shapeType = 0;
};
//---------------------------------------------------------------------------
BeamSecShape::BeamSecShape(double len, double wid){
        hh        = len;
        bb        = wid;
        shapeType = 1;
};
//---------------------------------------------------------------------------
double BeamSecShape::Ix(){
        double II_value = 0.0 ;
        switch( shapeType ){
        case 0:
                II_value = PI*rr*rr*rr*rr/4.0;
                break;
        case 1:
                II_value = bb*hh*hh*hh/12.0;
                break;
        default:
                break;
        };
        return II_value ;
};
//---------------------------------------------------------------------------
double BeamSecShape::Iy(){
        double II_value = 0.0 ;
        switch( shapeType ){
        case 0:
                II_value = PI*rr*rr*rr*rr/4.0;
                break;
        case 1:
                II_value = hh*bb*bb*bb/12.0;
                break;
        default:
                break;
        };
        return II_value ;
};
//---------------------------------------------------------------------------
double BeamSecShape::Iz(){
        double Iz_value = 0.0 ;
        return Iz_value ;
};
//---------------------------------------------------------------------------
double BeamSecShape::area(){
        double area_value = 0.0 ;
        switch( shapeType ){
        case 0:
                area_value = PI*rr*rr;
                break;
        case 1:
                area_value = hh*bb;
                break;
        default:
                break;
        };
        return area_value ;
};
//---------------------------------------------------------------------------
