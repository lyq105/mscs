//=========================================
//3Dvector.cpp
//=========================================
#include"Vector3D.h"
//-----------------------------------------
//¹¹Ôìº¯Êý
Vector3D::Vector3D(double ix,double iy,double iz){
	x=ix;
	y=iy;
	z=iz;
}//-------------------------------------------

Vector3D Vector3D::operator+(const Vector3D &d){
	return Vector3D(x+d.x,y+d.y,z+d.z);
}//-------------------------------------------


Vector3D Vector3D::operator-(const Vector3D &d){
	return Vector3D(x-d.x,y-d.y,z-d.z);
}//--------------------------------------------

double Vector3D::operator *(const Vector3D &d){
	return x*d.x+y*d.y+z*d.z;
}//---------------------------------------------

Vector3D Vector3D::operator *(double d){
	return Vector3D(x*d,y*d,z*d);
}//----------------------------------------------

Vector3D Vector3D::operator /(double d){
	return Vector3D(x/d,y/d,z/d);
}//----------------------------------------------

Vector3D Vector3D::Normallize(){
	return Vector3D(x/sqrt(x*x+y*y+z*z),y/sqrt(x*x+y*y+z*z),z/sqrt(x*x+y*y+z*z));
}//----------------------------------------------

bool Vector3D::operator ==(const Vector3D &d){
	if(x==d.x && y==d.y && z==d.z)
		return 1;
	else 
		return 0;
}//----------------------------------------------

double Vector3D::Length(){
	return sqrt(x*x+y*y+z*z);
}//------------------------------------------------

