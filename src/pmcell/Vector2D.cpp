//=========================================
//2Dvector.cpp
//=========================================
#include"Vector2D.h"
//-----------------------------------------
//¹¹Ôìº¯Êý
Vector2D::Vector2D(double ix,double iy){
	x=ix;
	y=iy;
}//-------------------------------------------

Vector2D Vector2D::operator+(const Vector2D &d){
	return Vector2D(x+d.x,y+d.y);
}//-------------------------------------------


Vector2D Vector2D::operator-(const Vector2D &d){
	return Vector2D(x-d.x,y-d.y);
}//--------------------------------------------

double Vector2D::operator *(const Vector2D &d){
	return x*d.x+y*d.y;
}//---------------------------------------------

Vector2D Vector2D::operator *(double d){
	return Vector2D(x*d,y*d);
}//----------------------------------------------

Vector2D Vector2D::operator /(double d){
	return Vector2D(x/d,y/d);
}//----------------------------------------------

Vector2D Vector2D::Normallize(){
	return Vector2D(x/sqrt(x*x+y*y),y/sqrt(x*x+y*y));
}//----------------------------------------------

bool Vector2D::operator ==(const Vector2D &d){
	if(x==d.x && y==d.y)
		return 1;
	else 
		return 0;
}//----------------------------------------------

double Vector2D::Length(){
	return sqrt(x*x+y*y);
}//------------------------------------------------

