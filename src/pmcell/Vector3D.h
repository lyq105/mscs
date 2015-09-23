
//==========================================
//3ά�����ࣻ
//Vector3D.h
//==========================================
#ifndef HEADER_Vector3D
#define HEADER_Vector3D
#include<iostream>
#include<cmath>
#include "Hns.h"
using namespace hns;
//==========================================

class Vector3D{
protected:


public:
	double x,y,z;

	Vector3D(double ix=0,double iy=0,double iz=0);
	//���ز�������	
	Vector3D operator+(const Vector3D& d);
	Vector3D operator-(const Vector3D& d);
	double   operator*(const Vector3D& d);
	Vector3D operator*(double d);
	Vector3D operator/(double d);
	bool     operator== (const Vector3D& d);
	//��ά������λ����
	Vector3D Normallize();
	double   Length();
	//�����������������
	friend ostream& operator<<(ostream& o, const Vector3D& d);
	friend class Node;
};//========================================

inline ostream& operator<<(ostream& o,const Vector3D& d){
	return o<<d.x<<" "<<d.y<<" "<<d.z<<endl;
}//=========================================
#endif  //HEADER_2DVECTOR  
