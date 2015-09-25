
//==========================================
//二维向量类；
//Vector2D.h
//==========================================
#ifndef HEADER_VECTOR2D
#define HEADER_VECTOR2D
#include<iostream>
#include<cmath>
#include "Hns.h"
using namespace hns;
//==========================================

class Vector2D{
protected:
	
	
public:
      double x,y;

	  Vector2D(double ix=0,double iy=0);
//重载操作符；	
	 Vector2D operator+(const Vector2D& d);
	 Vector2D operator-(const Vector2D& d);
	 double   operator*(const Vector2D& d);
	 Vector2D operator*(double d);
	 Vector2D operator/(double d);
     bool     operator== (const Vector2D& d);
//二维向量单位化；
     Vector2D Normallize();
	 double   Length();
//重载输出流操作符；
	friend ostream& operator<<(ostream& o, const Vector2D& d);
	friend class Node;
};//========================================

inline ostream& operator<<(ostream& o,const Vector2D& d){
	return o<<d.x<<" "<<d.y<<endl;
}//=========================================
#endif  //HEADER_2DVECTOR  
