#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>
#include <float.h>
//#include <queue>

using namespace std;

//***********************//
//		BASIC CLASSES
//***********************//
class Vector; //check
class Normal; 
class Point;
class Ray;
class LocalGeo;
class SmallMatrix;
class Matrix;
class Transformation;
class Color;
class BRDF;
class Sample;
class Intersection;
class Material;
class Color;
class Shape;
class Light;

#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; } 
//****************************Vector class**************************
class Vector {
       
public:
      float x;
      float y;
      float z;
      Vector();
      Vector(float a, float b, float c);
      void vector();
      void vector(float a, float b, float c);
      void reset(float x, float y, float z);
      Vector add(Vector added);
      Vector subtract(Vector subtractee);
      Vector scalarmultiply(float a);
      Vector scalardivide(float b);
      Vector normalize();
      Vector crossproduct(Vector other);
      float get_x();
      float get_y();
      float get_z();
      void set_x(float x);
      void set_y(float y);
      void set_z(float z);
      void printline();
      void print();
      float dot(Vector v0);
};

Vector::Vector() {
      this->x = 0;
      this->y = 0;
      this->z = 0;
}
Vector::Vector(float a, float b, float c){
      this->x = a;
      this->y = b;
      this->z = c;
}

void Vector::vector() {
      this->x = 0;
      this->y = 0;
      this->z = 0;
}
void Vector::vector(float a, float b, float c){
      this->x = a;
      this->y = b;
      this->z = c;
}
void Vector::reset(float x, float y, float z){
      this->x = x;
      this->y = y;
      this->z = z;
}
Vector Vector::add(Vector added){
      Vector sum;
      float sumx = this->x + added.get_x();
      float sumy = this->y + added.get_y();
      float sumz = this->z + added.get_z();
      sum.set_x(sumx);
      sum.set_y(sumy);
      sum.set_z(sumz);
      return sum;
}
Vector Vector::subtract(Vector subtractee){
      Vector subtracted;
      float subx = this->x - subtractee.get_x();
      float suby = this->y - subtractee.get_y();
      float subz = this->z - subtractee.get_z();
      subtracted.set_x(subx);
      subtracted.set_y(suby);
      subtracted.set_z(subz);
      return subtracted;
}
Vector Vector::scalarmultiply(float a){
      Vector temp;
      temp.set_x(this->x * a);
      temp.set_y(this->y * a);
      temp.set_z(this->z * a);
      return temp;
}
Vector Vector::scalardivide(float a){
      Vector temp;
      temp.set_x(this->x / a);
      temp.set_y(this->y / a);
      temp.set_z(this->z / a);
      return temp;
}
Vector Vector::normalize(){
      Vector result;
      result.vector();
      if (x != 0 || y != 0 || z != 0) {
            float temp = sqrt(sqrt(this->x)+sqrt(this->y)+sqrt(this->z));
            result.set_x(this->x / temp);
            result.set_y(this->y / temp);
            result.set_z(this->z / temp);
            return result;
      } else {
            return result;
      }
}
Vector Vector::crossproduct(Vector other){ //computes u x v with u->this and v->other
      Vector temp;
      temp.set_x((this->y * other.get_z()) - (other.get_y() * this->z));
      temp.set_y((this->z * other.get_x()) - (other.get_z() * this->x));
      temp.set_z((this->x * other.get_y()) - (other.get_x() * this->y));
      return temp;
}
float Vector::get_x(){
      return this->x;
}
float Vector::get_y(){
      return this->y;
}
float Vector::get_z(){
      return this->z;
}
void Vector::set_x(float x){
      this->x = x;
}
void Vector::set_y(float y){
      this->y = y;
}
void Vector::set_z(float z){
      this->z = z;
}
void Vector::printline(){           //for testing
      printf("Vector: x=%f, y=%f, z=%f\n", this->x, this->y, this->z);
}
void Vector::print(){              //for testing
      printf("x=%f, y=%f, z=%f", this->x, this->y, this->z);
}
float Vector::dot(Vector v0) {
      return (this->x * v0.get_x()) + (this->y * v0.get_y()) + (this->z * v0.get_z());
}

                        
//Point – Point = Vector 
 
//****************************Normal class**************************
//*********This is a child of Vector class, the results of operations of the Normal are automatically normalized****                       
class Normal : public Vector {
public:
      Normal();
      Normal(float a, float b, float c);
      void normal();
      void normal(float a, float b, float c);
      void reset(float x, float y, float z);
      Normal add(Vector addee);
      Normal subtract(Vector subtractee);
      void printline();
      void print();
};
Normal::Normal(){
      this->x = 0;
      this->y = 0;
      this->z = 0;
}
Normal::Normal(float a, float b, float c){
      if (x != 0 || y != 0 || z != 0) {
            Vector notnormalized;
            notnormalized.vector(x, y, z);
            this->x = notnormalized.normalize().get_x();
            this->y = notnormalized.normalize().get_y();
            this->z = notnormalized.normalize().get_z();
      } else {
            this->x = 0;
            this->y = 0;
            this->z = 0;
      }  
}
void Normal::normal(){
      this->x = 0;
      this->y = 0;
      this->z = 0;
      }
 
void Normal::normal(float x, float y, float z){
      if (x != 0 || y != 0 || z != 0) {
            Vector notnormalized;
            notnormalized.vector(x, y, z);
            this->x = notnormalized.normalize().get_x();
            this->y = notnormalized.normalize().get_y();
            this->z = notnormalized.normalize().get_z();
      } else {
            this->x = 0;
            this->y = 0;
            this->z = 0;
      }     
}
void Normal::reset(float x, float y, float z){
      if (x != 0 || y != 0 || z != 0) {
            Vector temp;
            temp.set_x(x);
            temp.set_y(y);
            temp.set_z(z);
            this->x = temp.normalize().get_x();
            this->y = temp.normalize().get_y();
            this->z = temp.normalize().get_z();
      } else {
            this->x = 0;
            this->y = 0;
            this->z = 0;
      }
}      
Normal Normal::add(Vector addee){
      float tempx, tempy, tempz;
      Normal temp;
      tempx = this->x + addee.get_x();
      tempy = this->y + addee.get_y();
      tempz = this->z + addee.get_z();
      temp.normal(tempx, tempy, tempz);
      return temp;
}
 
Normal Normal::subtract(Vector subtractee){
      float tempx, tempy, tempz;
      Normal temp;
      tempx = this->x - subtractee.get_x();
      tempy = this->y - subtractee.get_y();
      tempz = this->z - subtractee.get_z();
      temp.normal(tempx, tempy, tempz);
      return temp;
}
void Normal::printline(){          //for testing
      printf("Normal: x=%f, y=%f, z=%f\n", this->x, this->y, this->z);
}
void Normal::print(){              //for testing
      printf("x=%f, y=%f, z=%f", this->x, this->y, this->z);
}
 
//**************************Point class***************************
// Point - Point = Vector; Point - Vector = Point;
class Point {
public:
      float x;
      float y;
      float z;
      Point();
      Point(float x, float y, float z);
      void point();
      void point(float x, float y, float z);
      Point PaddvectorV(Vector addee);
      Point PsubtractV(Vector subtractee);
      Vector PsubtractP(Point subtractee);
      float get_x();
      float get_y();
      float get_z();
      void set_x(float x);
      void set_y(float y);
      void set_z(float z);
      void printline();
      void print();
 
};
Point::Point(){
      this->x = 0;
      this->y = 0;
      this->z = 0;
}
Point::Point(float x, float y, float z){
      this->x = x;
      this->y = y;
      this->z = z;
}
void Point::point(){
      this->x = 0;
      this->y = 0;
      this->z = 0;
}
void Point::point(float x, float y, float z){
      this->x = x;
      this->y = y;
      this->z = z;
}
Point Point::PaddvectorV(Vector addee){
      Point destination;
      destination.x = this->x + addee.x;
      destination.y = this->y + addee.y;
      destination.z = this->z + addee.z;
      return destination;
}
Point Point::PsubtractV(Vector subtractee){
      Point destination;
      destination.x = this->x - subtractee.x;
      destination.y = this->y - subtractee.y;
      destination.z = this->z - subtractee.z;
      return destination;
}
Vector Point::PsubtractP(Point subtractee){
      Vector result;
      result.x = this->x - subtractee.x;
      result.y = this->y - subtractee.y;
      result.z = this->z - subtractee.z;
      return result;
}
float Point::get_x(){
      return this->x;
}
float Point::get_y(){
      return this->y;
}
float Point::get_z(){
      return this->z;
}
void Point::set_x(float x){
      this->x = x;
}
void Point::set_y(float y){
      this->y = y;
}
void Point::set_z(float z){
      this->z = z;
}
void Point::printline(){                //for testing
      printf("Point: x=%f, y=%f, z=%f\n", this->x, this->y, this->z);
}
void Point::print(){              //for testing
      printf("x=%f, y=%f, z=%f", this->x, this->y, this->z);
}
 
 
//**************************Ray class***************************
//It represent the ray ray(t) = pos + t*dir, where t_min <= t <= t_max
class Ray {
 
public:
      Point pos;
      Vector dir;
      float t_min;
      float t_max;
      float t;

      Ray();
      Ray(Point a, Vector b, float c, float d);
      void ray();
      void ray(Point a, Vector b, float c, float d);
      Point currentposition(float t);
      Point get_pos();
      Vector get_dir();
      void set_pos(Point pos);
      void set_dir(Vector dir);
      void printline();
      void print();
};
Ray::Ray(){

}

Ray::Ray(Point a, Vector b, float c, float d){
      this->pos = a;
      this->dir = b;
      this->t_min = c;
      this->t_max = d;
}
void Ray::ray(){
      //do nothing
}
void Ray::ray(Point a, Vector b, float c, float d){

      this->pos = a;
      this->dir = b;
      this->t_min = c;
      this->t_max = d;
}
Point Ray::currentposition(float t){
      Point current;
      current = this->pos.PaddvectorV(this->dir.scalarmultiply(t));
      return current;
}
Point Ray::get_pos(){
      return this->pos;
}
Vector Ray::get_dir(){
      return this->dir;
}
void Ray::set_pos(Point pos){
      this->pos = pos;
}
void Ray::set_dir(Vector dir){
      this->dir = dir;
}
void Ray::printline(){
      printf("Ray: ");
      printf("Pos-> (");
      this->pos.print();
      printf(") Dir-> (");
      this->dir.print();
      printf(")\n");
}



/**********LocalGeo Class***********/
class LocalGeo {
public:
      Point pos;
      Normal normal;
      LocalGeo();
      LocalGeo(Point pos, Normal normal);
      void localGeo();
      void localGeo(Point pos, Normal normal);
      Point get_pos();
      Normal get_normal();
      void printline();
};
LocalGeo::LocalGeo(){
      //object is created
      this->pos = Point();
      this->normal = Normal();
}
LocalGeo::LocalGeo(Point pos, Normal normal){
      this->pos = pos;
      this->normal = normal;
}
void LocalGeo::localGeo(){
      this->pos = Point();
      this->normal = Normal();
}
void LocalGeo::localGeo(Point pos, Normal normal){
      this->pos = pos;
      this->normal = normal;
}
Point LocalGeo::get_pos(){
      return this->pos;
}
Normal LocalGeo::get_normal(){
      return this->normal;
}
void LocalGeo::printline(){
      printf("LocalGeo: ");
      printf("Pos-> (");
      this->pos.print();
      printf(") Normal-> (");
      this->normal.print();
      printf(")\n");
}




//for calculating intersections
class SmallMatrix {
public:
      float pos[3][3] = {
            {0,0,0},
            {0,0,0},
            {0,0,0}
      };
      SmallMatrix();
      SmallMatrix(float a, float b, float c, float d, float e, float f, float g, float h, float i);
      void smallMatrix();
      void smallMatrix(float a, float b, float c, float d, float e, float f, float g, float h, float i);
      float determinant();
      void print();
};
SmallMatrix::SmallMatrix(){
}
SmallMatrix::SmallMatrix(float a, float b, float c, float d, float e, float f, float g, float h, float i){
      this->pos[0][0] = a;
      this->pos[1][0] = b;
      this->pos[2][0] = c;
      this->pos[0][1] = d;
      this->pos[1][1] = e;
      this->pos[2][1] = f;
      this->pos[0][2] = g;
      this->pos[1][2] = h;
      this->pos[2][2] = i;
}
void SmallMatrix::smallMatrix() {
}
void SmallMatrix::smallMatrix(float a, float b, float c, float d, float e, float f, float g, float h, float i) {
      this->pos[0][0] = a;
      this->pos[1][0] = b;
      this->pos[2][0] = c;
      this->pos[0][1] = d;
      this->pos[1][1] = e;
      this->pos[2][1] = f;
      this->pos[0][2] = g;
      this->pos[1][2] = h;
      this->pos[2][2] = i;
}
float SmallMatrix::determinant() {
      return (this->pos[0][0] * this->pos[1][1] * this->pos[2][2])
      + (this->pos[1][0] * this->pos[2][1] * this->pos[0][2])
      + (this->pos[2][0] * this->pos[0][1] * this->pos[1][2])
      - (this->pos[2][0] * this->pos[1][1] * this->pos[0][2])
      - (this->pos[0][0] * this->pos[2][1] * this->pos[1][2])
      - (this->pos[1][0] * this->pos[0][1] * this->pos[2][2])
      ;
}

 
//**************************Matrix class***************************
//a matrix is 4x4
class Matrix{
//successfully debugged, now working: constructors, translation, scaling, m*m multiplication, identity, inverse, transpose,
       
public:
      float pos[4][4] = {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}};
      Matrix();
      Matrix(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p);
      void matrix();
      void matrix(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p);
      Matrix translation(float tx, float ty, float tz);
      Matrix scaling(float sx, float sy, float sz);
      Matrix rotation(float rx, float ry, float rz);
      Matrix arbitrary_rotation(float axisX, float axisY, float axisZ, float theta);
      Matrix multiplication(Matrix temp);
      Matrix operator*(Matrix temp);
      float determinant();
      Matrix inverse();
      Matrix transpose();
      Matrix identity();
      void print();
};
Matrix::Matrix(){
      pos[0][0] = 0;
}
Matrix::Matrix(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p){
      this->pos[0][0] = a;
      this->pos[1][0] = b;
      this->pos[2][0] = c;
      this->pos[3][0] = d;
      this->pos[0][1] = e;
      this->pos[1][1] = f;
      this->pos[2][1] = g;
      this->pos[3][1] = h;
      this->pos[0][2] = i;
      this->pos[1][2] = j;
      this->pos[2][2] = k;
      this->pos[3][2] = l;
      this->pos[0][3] = m;
      this->pos[1][3] = n;
      this->pos[2][3] = o;
      this->pos[3][3] = p;
 
}
void Matrix::matrix(){
      pos[0][0] = 0;
}
void Matrix::matrix(float a, float b, float c, float d, float e, float f, float g, float h, float i, float j, float k, float l, float m, float n, float o, float p){
      this->pos[0][0] = a;
      this->pos[1][0] = b;
      this->pos[2][0] = c;
      this->pos[3][0] = d;
      this->pos[0][1] = e;
      this->pos[1][1] = f;
      this->pos[2][1] = g;
      this->pos[3][1] = h;
      this->pos[0][2] = i;
      this->pos[1][2] = j;
      this->pos[2][2] = k;
      this->pos[3][2] = l;
      this->pos[0][3] = m;
      this->pos[1][3] = n;
      this->pos[2][3] = o;
      this->pos[3][3] = p;
 
}
Matrix translation(float tx, float ty, float tz){
      Matrix result;
      result.matrix(1,0,0,tx,0,1,0,ty,0,0,1,tz,0,0,0,1);
 
      return result;
}
Matrix scaling(float sx, float sy, float sz){
      Matrix result;
      result.pos[0][0] = sx;
      result.pos[1][1] = sy;
      result.pos[2][2] = sz;
      result.pos[3][3] = 1;
      return result;
}
Matrix rotation(float rx, float ry, float rz){
      //rx, ry, rz are agnles in degrees of rotation
      Matrix result;
      float radX = 3.14159265359*rx/180;
      float radY = 3.14159265359*ry/180;
      float radZ = 3.14159265359*rz/180;
      Matrix rxmatrix;
      Matrix rymatrix;
      Matrix rzmatrix;
 
      rxmatrix.matrix(1,0,0,0,      0,cos(radX),sin(radX),0,      0,-sin(radX),cos(radX),0,     0,0,0,1);
 
      rymatrix.matrix(cos(radY),0,-sin(radY),0,     0,1,0,0,    sin(radY),0,cos(radY),0,      0,0,0,1);
 
      rzmatrix.matrix(cos(radZ),sin(radZ),0,0,    -sin(radZ),cos(radZ),0,0,     0,0,1,0,    0,0,0,1);
 
      result = rxmatrix*rymatrix*rzmatrix; 
 
      return result;
}
Matrix arbitrary_rotation(float axisX, float axisY, float axisZ, float theta){
      //inputs are the x, y, z coordinates of the rotation axis and the angle of rotation about this axis
      Matrix result;
      Normal raxis;
      raxis.normal(axisX, axisY, axisZ); //automatically normalizes this vector
      float rad = 3.14159265359*theta/180;
      float cosine = cos(rad);
      float sine = sin(rad);
      float x = raxis.x;
      float y = raxis.y;
      float z = raxis.z;
      result.matrix(cosine+x*x*(1-cosine),y*x*(1-cosine)+z*sine,z*x*(1-cosine)-y*sine,0,   x*y*(1-cosine)-z*sine,cosine+y*y*(1-cosine),z*y*(1-cosine)+x*sine,0,      x*z*(1-cosine)+y*sine,y*z*(1-cosine)-x*sine,cosine+z*z*(1-cosine),0,      0,0,0,1);
      return result;
}
Matrix Matrix::multiplication(Matrix temp){
      Matrix result;
      for (int i = 0; i < 4; i = i + 1){
            for (int j = 0; j < 4; j = j + 1){
                  result.pos[i][j] = this->pos[0][i]*temp.pos[j][0] + this->pos[1][i]*temp.pos[j][1] + this->pos[2][i]*temp.pos[j][2] + this->pos[3][i]*temp.pos[j][3];
            }
      }
      return result;
}
Vector multiplicationV(Matrix m, Vector v){
      Vector result;
      float x = m.pos[0][0] * v.x + m.pos[0][1] * v.y + m.pos[0][2] * v.z + m.pos[0][3];
      float y = m.pos[1][0] * v.x + m.pos[1][1] * v.y + m.pos[1][2] * v.z + m.pos[1][3];
      float z = m.pos[2][0] * v.x + m.pos[2][1] * v.y + m.pos[2][2] * v.z + m.pos[2][3];
      return Vector(x, y, z);
}
Point multiplicationP(Matrix m, Point p){
      float x = m.pos[0][0] * p.x + m.pos[0][1] * p.y + m.pos[0][2] * p.z + m.pos[0][3];
      float y = m.pos[1][0] * p.x + m.pos[1][1] * p.y + m.pos[1][2] * p.z + m.pos[1][3];
      float z = m.pos[2][0] * p.x + m.pos[2][1] * p.y + m.pos[2][2] * p.z + m.pos[2][3];
      return Point(x, y, z);
}
Ray multiplicationR(Matrix m, Ray ray){
      Vector v = ray.get_dir();
      Point p = ray.get_pos();
      float x = m.pos[0][0] * v.x + m.pos[0][1] * v.y + m.pos[0][2] * v.z;
      float y = m.pos[1][0] * v.x + m.pos[1][1] * v.y + m.pos[1][2] * v.z;
      float z = m.pos[2][0] * v.x + m.pos[2][1] * v.y + m.pos[2][2] * v.z;
      Vector r = Vector(x, y, z);
      return Ray(p, r, 0, 0);
}
LocalGeo multiplicationL(Matrix m, LocalGeo localGeo){
      Vector v = localGeo.get_normal();
      Matrix minvt;
      minvt = m.inverse();
      float x = minvt.pos[0][0] * v.x + minvt.pos[0][1] * v.y + minvt.pos[0][2] * v.z;
      float y = minvt.pos[1][0] * v.x + minvt.pos[1][1] * v.y + minvt.pos[1][2] * v.z;
      float z = minvt.pos[2][0] * v.x + minvt.pos[2][1] * v.y + minvt.pos[2][2] * v.z;
      v = Normal(x, y, z);
      Point n;
      n = localGeo.pos;
      LocalGeo t;
      t = LocalGeo(n, v);
      return t;
}
Matrix Matrix::identity(){
      Matrix result;
      result.matrix(1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    0,0,0,1);
      return result;
}
Matrix Matrix::operator*(Matrix temp){
      Matrix result;
      result = this->multiplication(temp);
            return result;
}
float Matrix::determinant(){
      float det;
      det = this->pos[0][0] * this->pos[1][1] * this->pos[2][2] * this->pos[3][3] + 
      this->pos[0][0] * this->pos[1][2] * this->pos[2][3] * this->pos[3][1] + 
      this->pos[0][0] * this->pos[1][3] * this->pos[2][1] * this->pos[3][2] + 
      this->pos[0][1] * this->pos[1][0] * this->pos[2][3] * this->pos[3][2] + 
      this->pos[0][1] * this->pos[1][2] * this->pos[2][0] * this->pos[3][3] + 
      this->pos[0][1] * this->pos[1][3] * this->pos[2][2] * this->pos[3][0] +
      this->pos[0][2] * this->pos[1][0] * this->pos[2][1] * this->pos[3][3] +
      this->pos[0][2] * this->pos[1][1] * this->pos[2][3] * this->pos[3][0] +
      this->pos[0][2] * this->pos[1][3] * this->pos[2][0] * this->pos[3][1] +
      this->pos[0][3] * this->pos[1][0] * this->pos[2][2] * this->pos[3][1] +
      this->pos[0][3] * this->pos[1][1] * this->pos[2][0] * this->pos[3][2] +
      this->pos[0][3] * this->pos[1][2] * this->pos[2][1] * this->pos[3][0] - 
 
      this->pos[0][0] * this->pos[1][1] * this->pos[2][3] * this->pos[3][2] - 
      this->pos[0][0] * this->pos[1][2] * this->pos[2][1] * this->pos[3][3] - 
      this->pos[0][0] * this->pos[1][3] * this->pos[2][2] * this->pos[3][1] - 
      this->pos[0][1] * this->pos[1][0] * this->pos[2][2] * this->pos[3][3] - 
      this->pos[0][1] * this->pos[1][2] * this->pos[2][3] * this->pos[3][0] - 
      this->pos[0][1] * this->pos[1][3] * this->pos[2][0] * this->pos[3][2] - 
      this->pos[0][2] * this->pos[1][0] * this->pos[2][3] * this->pos[3][1] - 
      this->pos[0][2] * this->pos[1][1] * this->pos[2][0] * this->pos[3][3] - 
      this->pos[0][2] * this->pos[1][3] * this->pos[2][1] * this->pos[3][0] - 
      this->pos[0][3] * this->pos[1][0] * this->pos[2][1] * this->pos[3][2] - 
      this->pos[0][3] * this->pos[1][1] * this->pos[2][2] * this->pos[3][0] - 
      this->pos[0][3] * this->pos[1][2] * this->pos[2][0] * this->pos[3][1];
      return det;
}
Matrix Matrix::inverse(){
      float det;
      Matrix result;
 
      det = this->pos[0][0] * this->pos[1][1] * this->pos[2][2] * this->pos[3][3] + 
      this->pos[0][0] * this->pos[1][2] * this->pos[2][3] * this->pos[3][1] + 
      this->pos[0][0] * this->pos[1][3] * this->pos[2][1] * this->pos[3][2] + 
      this->pos[0][1] * this->pos[1][0] * this->pos[2][3] * this->pos[3][2] + 
      this->pos[0][1] * this->pos[1][2] * this->pos[2][0] * this->pos[3][3] + 
      this->pos[0][1] * this->pos[1][3] * this->pos[2][2] * this->pos[3][0] +
      this->pos[0][2] * this->pos[1][0] * this->pos[2][1] * this->pos[3][3] +
      this->pos[0][2] * this->pos[1][1] * this->pos[2][3] * this->pos[3][0] +
      this->pos[0][2] * this->pos[1][3] * this->pos[2][0] * this->pos[3][1] +
      this->pos[0][3] * this->pos[1][0] * this->pos[2][2] * this->pos[3][1] +
      this->pos[0][3] * this->pos[1][1] * this->pos[2][0] * this->pos[3][2] +
      this->pos[0][3] * this->pos[1][2] * this->pos[2][1] * this->pos[3][0] - 
 
      this->pos[0][0] * this->pos[1][1] * this->pos[2][3] * this->pos[3][2] - 
      this->pos[0][0] * this->pos[1][2] * this->pos[2][1] * this->pos[3][3] - 
      this->pos[0][0] * this->pos[1][3] * this->pos[2][2] * this->pos[3][1] - 
      this->pos[0][1] * this->pos[1][0] * this->pos[2][2] * this->pos[3][3] - 
      this->pos[0][1] * this->pos[1][2] * this->pos[2][3] * this->pos[3][0] - 
      this->pos[0][1] * this->pos[1][3] * this->pos[2][0] * this->pos[3][2] - 
      this->pos[0][2] * this->pos[1][0] * this->pos[2][3] * this->pos[3][1] - 
      this->pos[0][2] * this->pos[1][1] * this->pos[2][0] * this->pos[3][3] - 
      this->pos[0][2] * this->pos[1][3] * this->pos[2][1] * this->pos[3][0] - 
      this->pos[0][3] * this->pos[1][0] * this->pos[2][1] * this->pos[3][2] - 
      this->pos[0][3] * this->pos[1][1] * this->pos[2][2] * this->pos[3][0] - 
      this->pos[0][3] * this->pos[1][2] * this->pos[2][0] * this->pos[3][1];
 
      result.pos[0][0] = (this->pos[1][1] * this->pos[2][2] * this->pos[3][3] + 
            this->pos[1][2] * this->pos[2][3] * this->pos[3][1] + 
            this->pos[1][3] * this->pos[2][1] * this->pos[3][2] - 
            this->pos[1][1] * this->pos[2][3] * this->pos[3][2] - 
            this->pos[1][2] * this->pos[2][1] * this->pos[3][3] - 
            this->pos[1][3] * this->pos[2][2] * this->pos[3][1])/det;
      result.pos[0][1] = (this->pos[0][1] * this->pos[2][3] * this->pos[3][2] + 
            this->pos[0][2] * this->pos[2][1] * this->pos[3][3] + 
            this->pos[0][3] * this->pos[2][2] * this->pos[3][1] - 
            this->pos[0][1] * this->pos[2][2] * this->pos[3][3] - 
            this->pos[0][2] * this->pos[2][3] * this->pos[3][1] - 
            this->pos[0][3] * this->pos[2][1] * this->pos[3][2])/det;
      result.pos[0][2] = (this->pos[0][1] * this->pos[1][2] * this->pos[3][3] + 
            this->pos[0][2] * this->pos[1][3] * this->pos[3][1] + 
            this->pos[0][3] * this->pos[1][1] * this->pos[3][2] - 
            this->pos[0][1] * this->pos[1][3] * this->pos[3][2] - 
            this->pos[0][2] * this->pos[1][1] * this->pos[3][3] - 
            this->pos[0][3] * this->pos[1][2] * this->pos[3][1])/det;
      result.pos[0][3] = (this->pos[0][1] * this->pos[1][3] * this->pos[2][2] + 
            this->pos[0][2] * this->pos[1][1] * this->pos[2][3] + 
            this->pos[0][3] * this->pos[1][2] * this->pos[2][1] - 
            this->pos[0][1] * this->pos[1][2] * this->pos[2][3] - 
            this->pos[0][2] * this->pos[1][3] * this->pos[2][1] - 
            this->pos[0][3] * this->pos[1][1] * this->pos[2][2])/det;
      result.pos[1][0] = (this->pos[1][0] * this->pos[2][3] * this->pos[3][2] + 
            this->pos[1][2] * this->pos[2][0] * this->pos[3][3] + 
            this->pos[1][3] * this->pos[2][2] * this->pos[3][0] - 
            this->pos[1][0] * this->pos[2][2] * this->pos[3][3] - 
            this->pos[1][2] * this->pos[2][3] * this->pos[3][0] - 
            this->pos[1][3] * this->pos[2][0] * this->pos[3][2])/det;
      result.pos[1][1] = (this->pos[0][0] * this->pos[2][2] * this->pos[3][3] + 
            this->pos[0][2] * this->pos[2][3] * this->pos[3][0] + 
            this->pos[0][3] * this->pos[2][0] * this->pos[3][2] - 
            this->pos[0][0] * this->pos[2][3] * this->pos[3][2] - 
            this->pos[0][2] * this->pos[2][0] * this->pos[3][3] - 
            this->pos[0][3] * this->pos[2][2] * this->pos[3][0])/det;
      result.pos[1][2] = (this->pos[0][0] * this->pos[1][3] * this->pos[3][2] + 
            this->pos[0][2] * this->pos[1][0] * this->pos[3][3] + 
            this->pos[0][3] * this->pos[1][2] * this->pos[3][0] - 
            this->pos[0][0] * this->pos[1][2] * this->pos[3][3] - 
            this->pos[0][2] * this->pos[1][3] * this->pos[3][0] - 
            this->pos[0][3] * this->pos[1][0] * this->pos[3][2])/det;
      result.pos[1][3] = (this->pos[0][0] * this->pos[1][2] * this->pos[2][3] + 
            this->pos[0][2] * this->pos[1][3] * this->pos[2][0] + 
            this->pos[0][3] * this->pos[1][0] * this->pos[2][2] - 
            this->pos[0][0] * this->pos[1][3] * this->pos[2][2] - 
            this->pos[0][2] * this->pos[1][0] * this->pos[2][3] - 
            this->pos[0][3] * this->pos[1][2] * this->pos[2][0])/det;
      result.pos[2][0] = (this->pos[1][0] * this->pos[2][1] * this->pos[3][3] + 
            this->pos[1][1] * this->pos[2][3] * this->pos[3][0] + 
            this->pos[1][3] * this->pos[2][0] * this->pos[3][1] - 
            this->pos[1][0] * this->pos[2][3] * this->pos[3][1] - 
            this->pos[1][1] * this->pos[2][0] * this->pos[3][3] - 
            this->pos[1][3] * this->pos[2][1] * this->pos[3][0])/det;
      result.pos[2][1] = (this->pos[0][0] * this->pos[2][3] * this->pos[3][1] + 
            this->pos[0][1] * this->pos[2][0] * this->pos[3][3] + 
            this->pos[0][3] * this->pos[2][1] * this->pos[3][0] - 
            this->pos[0][0] * this->pos[2][1] * this->pos[3][3] - 
            this->pos[0][1] * this->pos[2][3] * this->pos[3][0] - 
            this->pos[0][3] * this->pos[2][0] * this->pos[3][1])/det;
      result.pos[2][2] = (this->pos[0][0] * this->pos[1][1] * this->pos[3][3] + 
            this->pos[0][1] * this->pos[1][3] * this->pos[3][0] + 
            this->pos[0][3] * this->pos[1][0] * this->pos[3][1] - 
            this->pos[0][0] * this->pos[1][3] * this->pos[3][1] - 
            this->pos[0][1] * this->pos[1][0] * this->pos[3][3] - 
            this->pos[0][3] * this->pos[1][1] * this->pos[3][0])/det;
      result.pos[2][3] = (this->pos[0][0] * this->pos[1][3] * this->pos[2][1] + 
            this->pos[0][1] * this->pos[1][0] * this->pos[2][3] + 
            this->pos[0][3] * this->pos[1][1] * this->pos[2][0] - 
            this->pos[0][0] * this->pos[1][1] * this->pos[2][3] - 
            this->pos[0][1] * this->pos[1][3] * this->pos[2][0] - 
            this->pos[0][3] * this->pos[1][0] * this->pos[2][1])/det;
      result.pos[3][0] = (this->pos[1][0] * this->pos[2][2] * this->pos[3][1] + 
            this->pos[1][1] * this->pos[2][0] * this->pos[3][2] + 
            this->pos[1][2] * this->pos[2][1] * this->pos[3][0] - 
            this->pos[1][0] * this->pos[2][1] * this->pos[3][2] - 
            this->pos[1][1] * this->pos[2][2] * this->pos[3][0] - 
            this->pos[1][2] * this->pos[2][0] * this->pos[3][1])/det;
      result.pos[3][1] = (this->pos[0][0] * this->pos[2][1] * this->pos[3][2] + 
            this->pos[0][1] * this->pos[2][2] * this->pos[3][0] + 
            this->pos[0][2] * this->pos[2][0] * this->pos[3][1] - 
            this->pos[0][0] * this->pos[2][2] * this->pos[3][1] - 
            this->pos[0][1] * this->pos[2][0] * this->pos[3][2] - 
            this->pos[0][2] * this->pos[2][1] * this->pos[3][0])/det;
      result.pos[3][2] = (this->pos[0][0] * this->pos[1][2] * this->pos[3][1] + 
            this->pos[0][1] * this->pos[1][0] * this->pos[3][2] + 
            this->pos[0][2] * this->pos[1][1] * this->pos[3][0] - 
            this->pos[0][0] * this->pos[1][1] * this->pos[3][2] - 
            this->pos[0][1] * this->pos[1][2] * this->pos[3][0] - 
            this->pos[0][2] * this->pos[1][0] * this->pos[3][1])/det;
      result.pos[3][3] = (this->pos[0][0] * this->pos[1][1] * this->pos[2][2] + 
            this->pos[0][1] * this->pos[1][2] * this->pos[2][0] + 
            this->pos[0][2] * this->pos[1][0] * this->pos[2][1] - 
            this->pos[0][0] * this->pos[1][2] * this->pos[2][1] - 
            this->pos[0][1] * this->pos[1][0] * this->pos[2][2] - 
            this->pos[0][2] * this->pos[1][1] * this->pos[2][0])/det;
 
      return result;
}
Matrix Matrix::transpose(){
      Matrix result;
      result.pos[0][0] = this->pos[0][0];
      result.pos[0][1] = this->pos[1][0];
      result.pos[0][2] = this->pos[2][0];
      result.pos[0][3] = this->pos[3][0];
      result.pos[1][0] = this->pos[0][1];
      result.pos[1][1] = this->pos[1][1];
      result.pos[1][2] = this->pos[2][1];
      result.pos[1][3] = this->pos[3][1];
      result.pos[2][0] = this->pos[0][2];
      result.pos[2][1] = this->pos[1][2];
      result.pos[2][2] = this->pos[2][2];
      result.pos[2][3] = this->pos[3][2];
      result.pos[3][0] = this->pos[0][3];
      result.pos[3][1] = this->pos[1][3];
      result.pos[3][2] = this->pos[2][3];
      result.pos[3][3] = this->pos[3][3];
      return result;
}
void Matrix::print(){
      int i,j;
 
      for (i = 0; i < 4; i++){
            for (j = 0; j < 4; j++){
                  printf("%f\t", pos[j][i]); 
            }
            printf("\n");
      }
      printf("\n");
}
 
/**************************Transformation Class*************************************/
 
class Transformation{

public:
      Matrix m;
      Matrix minvt;
      Transformation();
      Transformation(Matrix temp);
      void reset(Matrix new1);
      Transformation operator*(Transformation t);
      Point operator*(Point p);
      Ray operator*(Ray ray);
      LocalGeo operator*(LocalGeo localGeo);
};
Transformation::Transformation(){
      Matrix result;
      this->m = result;
      this->minvt = result.inverse();
}
Transformation::Transformation(Matrix temp){
      this->m = temp;
      this->minvt = temp.inverse();
}
void Transformation::reset(Matrix new1){
      this->m = new1;
      this->minvt = new1.inverse();
}
Transformation Transformation::operator*(Transformation t){
      Transformation result; 
      result.m = this->m * t.m;
      result.minvt = result.m.inverse();
      return result;
}
Point Transformation::operator*(Point p){
      Point result;
      result = multiplicationP(this->m, p);
      return result;
}
Ray Transformation::operator*(Ray ray){
      Ray result;
      result = multiplicationR(this->m, ray);
      return result;
}
LocalGeo Transformation::operator*(LocalGeo localgeo){
      LocalGeo result;
      result = multiplicationL(this->m, localgeo);
      return result;
}
 

  
 
/**********************BRDF***********************/
 class BRDF {

 
public:
      float kdr;
      float kdg;
      float kdb;
      float ksr;
      float ksg;
      float ksb;
      float kar;
      float kag;
      float kab;
      float krr;
      float krg;
      float krb;
      float p;
      BRDF();
      BRDF(float kdr,float kdg, float ksr, float ksg, float ksb, float kar, float kag, float kab, float krr, float krg, float krb, float p);
      void setKa(float r,float g, float b);
      void setKr(float r,float g, float b);
      void setKd(float r,float g, float b);
      void setKs(float r,float g, float b, float p);
      void setP(float p);
 
 
 };
 
BRDF::BRDF(){
      this->kdr = 0;
      this->kdg = 0;
      this->kdb = 0;
      this->ksr = 0;
      this->ksg = 0;
      this->ksb = 0;
      this->kar = 0;
      this->kag = 0;
      this->kab = 0;
      this->krr = 0;
      this->krg = 0;
      this->krb = 0;
      this->p = 0;
}
BRDF::BRDF(float kdr,float kdg, float ksr, float ksg, float ksb, float kar, float kag, float kab, float krr, float krg, float krb, float p){
      this->kdr = kdr;
      this->kdg = kdg;
      this->kdb = kdb;
      this->ksr = ksr;
      this->ksg = ksg;
      this->ksb = ksb;
      this->kar = kar;
      this->kab = kab;
      this->kag = kag;
      this->krr = krr;
      this->krg = krg;
      this->krb = krb;
      this->p = p;
 
} 
void BRDF::setKa(float r,float g, float b){
      this->kar = r;
      this->kag = g;
      this->kab = b;
}
void BRDF::setKd(float r,float g, float b){
      this->kdr = r;
      this->kdg = g;
      this->kdb = b;
}
void BRDF::setKs(float r,float g, float b, float p){
      this->ksr = r;
      this->ksg = g;
      this->ksb = b;
      this->p = p;
}
void BRDF::setKr(float r,float g, float b){
      this->krr = r;
      this->krg = g;
      this->krb = b;
}
 
 
 
 
/**************************Sample Class*******************/ 
/*    stores screen coordinates */
class Sample {
      float x;
      float y;
public:
      Sample();
      Sample(float x, float y);
};
Sample::Sample(){
      x = 0.0;
      y = 0.0;
}
Sample::Sample(float x, float y){
      this->x = x;
      this->y = y;
}




/*******************Primitive Class*********************/
class  Primitive{
public:
      virtual bool intersect(Ray& ray, float* thit, Intersection* in);
      virtual bool intersectP(Ray& ray);
      virtual void getBRDF(LocalGeo& local, BRDF* brdf);
 
};




 
/******************GeometricPrimitive Class**************/
class GeometricPrimitive : public Primitive {
      Transformation objToWorld;
      Transformation worldToObj;
      Shape* shape;
      BRDF* brdf;
 
public: 
      GeometricPrimitive(Shape* shape, Transformation transformation, BRDF* brdf);
      GeometricPrimitive(Shape *shape, float tx, float ty, float tz, float sx, float sy, float sz, float rotx, float roty, float rotz, float kar, float kag, float  kab, float kdr, float kdg, float kdb, float ksr, float ksg, float ksb);
      GeometricPrimitive(Shape *shape, float tx, float ty, float tz, float sx, float sy, float sz, float rx, float ry, float rz, float angle, float kar, float kag, float  kab, float kdr, float kdg, float kdb, float ksr, float ksg, float ksb);
      bool intersect(Ray& ray, float* thit, Intersection* in);
      bool intersectP(Ray& ray);
      void getBRDF(LocalGeo& local, BRDF* brdf);
}; 
GeometricPrimitive::GeometricPrimitive(Shape* shape, Transformation transformation, BRDF* brdf){
      this->objToWorld = transformation;
      Transformation temp;
      temp.m = transformation.m.inverse();
      temp.minvt = transformation.m;
      this->worldToObj = temp;
      this->shape = shape;
      this->brdf = brdf;
      }

GeometricPrimitive::GeometricPrimitive(Shape *shape, float tx, float ty, float tz, float sx, float sy, float sz, float rotx, float roty, float rotz, float kar, float kag, float  kab, float kdr, float kdg, float kdb, float ksr, float ksg, float ksb){

      Matrix rotate = rotation(rotx,roty,rotz);
      Matrix scale = scaling(sx,sy,sz);
      Matrix translate = translation(tx,ty,tz);
             
      BRDF brdf1;
      brdf1.kar = kar;
      brdf1.kag = kag;
      brdf1.kab = kab;
      brdf1.kdr = kdr;
      brdf1.kdg = kdg;
      brdf1.kdb = kdb;
      brdf1.ksr = ksr;
      brdf1.ksg = ksg;
      brdf1.ksb = ksb;
 
      this->objToWorld.m = translate*rotate*scale;
      this->worldToObj.m = this->objToWorld.m.inverse();
      this->shape = shape;
      this->brdf = brdf1;
      }
GeometricPrimitive::GeometricPrimitive(Shape *shape, float tx, float ty, float tz, float sx, float sy, float sz, float rx, float ry, float rz, float angle, float kar, float kag, float  kab, float kdr, float kdg, float kdb, float ksr, float ksg, float ksb){
      Matrix rotate = arbitrary_rotation(rx,ry,rz,angle);
      Matrix scale = scaling(sx,sy,sz);
      Matrix translate = translation(tx,ty,tz);
             
      BRDF brdf1;
      brdf1.kar = kar;
      brdf1.kag = kag;
      brdf1.kab = kab;
      brdf1.kdr = kdr;
      brdf1.kdg = kdg;
      brdf1.kdb = kdb;
      brdf1.ksr = ksr;
      brdf1.ksg = ksg;
      brdf1.ksb = ksb;
 
 
      this->objToWorld.m = translate*rotate*scale;
      this->worldToObj.m = this->objToWorld.m.inverse();
      this->shape = shape;
      this->brdf = brdf1;
}
bool GeometricPrimitive::intersect(Ray& ray, float* thit, Intersection* in)  {
      Ray oray = worldToObj*ray;
      LocalGeo olocal;                                 
      if (!shape->intersect(oray, thit, &olocal))  return false;
      in->primitive = this;
      in->local = objToWorld*olocal;
      return true;                               
}
bool GeometricPrimitive::intersectP(Ray& ray) {
      Ray oray = worldToObj*ray;
      return shape->intersectP(oray); 
                                                
}
void GeometricPrimitive::getBRDF(LocalGeo& local, BRDF* brdf) {
      this->brdf = brdf;
}
/************AggregatePrimitive********************/         
class AggregatePrimitive : public Primitive{
 
public:
    list<Primitive*> *primitives;
 
    AggregatePrimitive();
    AggregatePrimitive(vector<Primitive*> list);
    bool intersect(Ray& ray, float* thit, Intersection* in);
    bool intersectP(Ray& ray);
    void getBRDF(LocalGeo& local, BRDF* brdf){
        exit(1);
        //This should never get called, because in->primitve
        //will never be an aggregate primitive
    }
};
 
AggregatePrimitive::AggregatePrimitive(){
      this->primitives = NULL; 
}
 
AggregatePrimitive::AggregatePrimitive(vector<Primitive*> list){
     primitives = list;
}
void AggregatePrimitive::addPrimitive(Primitive temp){
     primitives.insert(temp);
}
 
 
bool AggregatePrimitive::intersect(Ray& ray, float* thit, Intersection* in){
    bool intersectobject = false;
    *thit = FLT_MAX;
    float newThit;
    for (auto primitive : primitives){
        Intersection newIn;
        if(primitive.shape.intersect(ray, &newThit, &newIn)){
            intersectobject = true;
            if (newThit < *thit){
                *thit = newThit;
                *in = newIn;
            }
        }
    }
    return intersectobject;
}
 
bool AggregatePrimitive::intersectP(Ray& ray){
    for (auto primitive: primitives){
        if (primitive.shape.intersectP(ray)) {
            return true;
        }
    }
    return false;
} 
/********************Material Class***********************/
class Material{
public:
      BRDF constantBRDF;
      BRDF getBRDF(LocalGeo& local, BRDF* brdf);
}; 

/***********************Color Class*****************/
//TODO: May support conversion from xyz
class Color {

public:
      float r;
      float g;
      float b;
	Color();
	Color(float r, float g, float b);
	void color();
	void color(float r, float g, float b);
	Color addColors(Color addeeColor);
	Color subColors(Color subtracteeColor);
	Color mulColorbyScalar(float scalar);
	Color divColorbyScalar(float scalar);
	float get_r();
	float get_g();
	float get_b();
	void set_r(float r);
	void set_g(float g);
    void set_b(float b);
 
};
Color::Color() {
	r = 0.0;
	g = 0.0;
	b = 0.0;
}
Color::Color(float r, float g, float b){
	this->r = r;
	this->g = g;
	this->b = b;
}
void Color::color() {
	r = 0.0;
	g = 0.0;
	b = 0.0;
}
 
void Color::color(float r, float g, float b){
	this->r = r;
	this->g = g;
	this->b = b;
}
 
void Color::addColors(Color addeeColor){
	this->r += addeeColor.get_r();
	this->g += addeeColor.get_g();
	this->b += addeeColor.get_b();
}
 
Color Color::subColors(Color subtracteeColor){
	this->r -= addeeColor.get_r();
	this->g -= addeeColor.get_g();
	this->b -= addeeColor.get_b();
}
 
Color Color::mulColorbyScalar(float scalar){
	this->r = scalar*(this->r);
	this->g = scalar*(this->g);
	this->b = scalar*(this->b);
}
 
Color Color::divColorbyScalar(float scalar){
	this->r = (this->r)/scalar;
	this->g = (this->g)/scalar;
	this->b = (this->b)/scalar;
}
float Color::get_r(){
	  return this->r;
}
float Color::get_g(){
	  return this->g;
}
float Color::get_b(){
	  return this->b;
}
void Color::set_r(float r){
	  this->r = r;
}
void Color::set_g(float g){
	  this->g = g;
}
void Color::set_b(float b){
	  this->b = b;
}


/*
      Implementations of triangle and sphere
      The intersection with the ray at t outside the range [t_min, t_max] should return false.
*/
class Shape {
      /* type=0 for sphere, type=1 for triangle, type=2 for not implemented */
      int type;

      /* implemented if sphere */
      float radius;
      Point center;

      /* if triangle */
      Point point0;
      Point point1;
      Point point2;

public:
      void shape();
      void makeSphere(float radius, Point center);
      void makeTriangle(Point point0, Point point1, Point point2);
      bool intersect(Ray& ray, float* thit, LocalGeo* local);
      bool intersectP(Ray& ray);

      float get_radius();
      Point get_center();
      Point get_point0();
      Point get_point1();
      Point get_point2();

      void set_radius(float radius);
      void set_center(Point center);
      void set_point0(Point point0);
      void set_point1(Point point1);
      void set_point2(Point point2);

      void printline();
};
void Shape::shape(){
      this->type = 2;
}

void Shape::makeSphere(float radius, Point center){
      this->type = 0;
      this->radius = radius;
      this->center = center;
}

void Shape::makeTriangle(Point point0, Point point1, Point point2){
      this->type = 1;
      this->point0 = point0;
      this->point1 = point1;
      this->point2 = point2;
}

float Shape::get_radius(){
      return this->radius;
}
Point Shape::get_center(){
      return this->center;
}
Point Shape::get_point0(){
      return this->point0;
}
Point Shape::get_point1(){
      return this->point1;
}
Point Shape::get_point2(){
      return this->point2;
}

void Shape::set_radius(float radius){
      this->radius = radius;
}
void Shape::set_center(Point center){
      this->center = center;
}
void Shape::set_point0(Point point0){
      this->point0 = point0;
}
void Shape::set_point1(Point point1){
      this->point1 = point1;
}
void Shape::set_point2(Point point2){
      this->point2 = point2;
}

void Shape::printline(){
      if (this->type==0){ //sphere
            printf("Sphere: ");
            printf("Center-> (");
            this->center.print();
            printf(") Radius-> (");
            printf("%f", this->radius);
            printf(")\n");
      }
      else if (this->type==1){ //triangular plane
            printf("Triangle: ");
            printf("Point0-> (");
            this->point0.print();
            printf(") Point1-> (");
            this->point1.print();
            printf(") Point2-> (");
            this->point2.print();
            printf(")\n");
      }
      else {
            printf("ERROR: Shape not implemented");
      }
}


// c = center, r = radius, e = starting position of ray, d = direction vector for ray
float find_discriminant(Point c, float R, Point e, Vector d) {
      float discriminant;
      discriminant = pow(d.scalarmultiply(2).dot(e.PsubtractP(c)), 2) - (4 * d.dot(d) * (e.PsubtractP(c).dot(e.PsubtractP(c)) - pow(R, 2)));
      return discriminant;
}


// c = center, r = radius, e = starting position of ray, d = direction vector for ray
float find_t_sphere(Point c, float R, Point e, Vector d) {
      return (d.scalarmultiply(-1).dot(e.PsubtractP(c)) - sqrt(pow(d.dot(e.PsubtractP(c)), 2) - (d.dot(d) * (e.PsubtractP(c).dot(e.PsubtractP(c)) - pow(R, 2)))))/(d.dot(d));  
}

// converts a vector to a normal
Normal vectorToNormal(Vector v) {
      Vector temp;
      temp = v.normalize();
      Normal n;
      n.normal(temp.get_x(), temp.get_y(), temp.get_z());
      return n;
}

// converts a pre-normalized vector to a normal
Normal normalizedVectorToNormal(Vector v) {
      Normal n;
      n.normal(v.get_x(), v.get_y(), v.get_z());
      return n;
}

//find distance between 2 points
float findDistance(Point a, Point b) {
      return sqrt(pow(b.get_x() - a.get_x(), 2) + pow(b.get_y() - a.get_y(), 2) + pow(b.get_z() - a.get_z(), 2));
}

//returns true if normal0 is the correct normal
bool isNormalCorrect(Vector normal0, Vector normal1, Point intersection, Point p) {
      if (findDistance(p, intersection.PaddvectorV(normal0)) < findDistance(p, intersection.PaddvectorV(normal1))) {
            return true;
      }
      return false;
}

/*
      Test if ray intersects with the shape or not (in object space), if so,
            return intersection point and normal
*/
bool Shape::intersect(Ray& ray, float* thit, LocalGeo* local){

      if (this->type == 0){ //sphere
            if (find_discriminant(this->center, this->radius, ray.get_pos(), ray.get_dir()) < 0) {
                  return false;
            }
            float myT;
            myT = find_t_sphere(this->center, this->radius, ray.get_pos(), ray.get_dir());
            if ((myT < ray.get_t_min()) || (myT > ray.get_t_max())) {
                  return false;
            }

            Point intersection;
            Vector normal;
            LocalGeo myLocal;

            *thit = myT;
            intersection = ray.get_pos().PaddvectorV(ray.get_dir().scalarmultiply(myT));
            normal = intersection.PsubtractP(this->center).scalardivide(this->radius);
            myLocal.localGeo(intersection, normalizedVectorToNormal(normal));
            *local = myLocal;

            myLocal.printline();
            return true;
      }
      else if (this->type == 1){ //triangle
            SmallMatrix A;
            A.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - this->point2.get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - this->point2.get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - this->point2.get_z(), ray.get_dir().get_z()
                  );

            float myT;
            SmallMatrix t_matrix;
            t_matrix.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - this->point2.get_x(), this->point0.get_x() - ray.get_pos().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - this->point2.get_y(), this->point0.get_y() - ray.get_pos().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - this->point2.get_z(), this->point0.get_z() - ray.get_pos().get_z()
                  );
            myT = t_matrix.determinant()/A.determinant();
            if (myT < ray.get_t_min() || myT > ray.get_t_max()) {
                  return false;
            }

            float myGamma;
            SmallMatrix gamma_matrix;
            gamma_matrix.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - ray.get_pos().get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - ray.get_pos().get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - ray.get_pos().get_z(), ray.get_dir().get_z()
                  );
            myGamma = gamma_matrix.determinant()/A.determinant();

            if (myGamma < 0 || myGamma > 1){
                  return false;
            }

            float myBeta;
            SmallMatrix beta_matrix;
            beta_matrix.smallMatrix(this->point0.get_x() - ray.get_pos().get_x(), this->point0.get_x() - this->point2.get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - ray.get_pos().get_y(), this->point0.get_y() - this->point2.get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - ray.get_pos().get_z(), this->point0.get_z() - this->point2.get_z(), ray.get_dir().get_z()
                  );
            myBeta = beta_matrix.determinant()/A.determinant();

            if (myBeta < 0 || myBeta > (1-myGamma)){
                  return false;
            }

            Point intersection;
            Vector normal0, normal1, v0, v1; //of normal0 and normal1, one will be the correct normal and one will be on the opposite side of plane
            LocalGeo myLocal;

            *thit = myT;
            intersection = ray.get_pos().PaddvectorV(ray.get_dir().scalarmultiply(myT));

            //cross product of 2 vectors is orthogonal to both these vectors
            v0 = intersection.PsubtractP(this->point0); 
            v1 = intersection.PsubtractP(this->point1);
            normal0 = v0.crossproduct(v1);
            normal1 = normal0.scalarmultiply(-1);

            if (isNormalCorrect(normal0, normal1, intersection, ray.get_pos())) {
                  myLocal.localGeo(intersection, vectorToNormal(normal0));
            }
            else {
                  myLocal.localGeo(intersection, vectorToNormal(normal1));
            }

            //printf("%f\n", myT);
            //myLocal.printline();
            *local = myLocal;
            //local->printline();
            return true;
      }
      else { //shape has not been set up
            return false;
      }
}

/*
      Same as intersect, but just return whether there is any intersection or not
*/
bool Shape::intersectP(Ray& ray){
      bool isIntersect = false;
      Point* intersection;
      Normal* normal;
      float* thit;
      if (this->type == 0){ //sphere
            if (find_discriminant(this->center, this->radius, ray.get_pos(), ray.get_dir()) < 0) {
                  return false;
            }
            float myT;
            myT = find_t_sphere(this->center, this->radius, ray.get_pos(), ray.get_dir());
            if ((myT < ray.get_t_min()) || (myT > ray.get_t_max())){
                  return false;
            }
            return true;
      }
      else if (this->type == 1){ //triangle
            SmallMatrix A;
            A.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - this->point2.get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - this->point2.get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - this->point2.get_z(), ray.get_dir().get_z()
                  );

            float myT;
            SmallMatrix t_matrix;
            t_matrix.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - this->point2.get_x(), this->point0.get_x() - ray.get_pos().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - this->point2.get_y(), this->point0.get_y() - ray.get_pos().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - this->point2.get_z(), this->point0.get_z() - ray.get_pos().get_z()
                  );
            myT = t_matrix.determinant()/A.determinant();
            if (myT < ray.get_t_min() || myT > ray.get_t_max()) {
                  return false;
            }

            float myGamma;
            SmallMatrix gamma_matrix;
            gamma_matrix.smallMatrix(this->point0.get_x() - this->point1.get_x(), this->point0.get_x() - ray.get_pos().get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - this->point1.get_y(), this->point0.get_y() - ray.get_pos().get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - this->point1.get_z(), this->point0.get_z() - ray.get_pos().get_z(), ray.get_dir().get_z()
                  );
            myGamma = gamma_matrix.determinant()/A.determinant();

            if (myGamma < 0 || myGamma > 1){
                  return false;
            }

            float myBeta;
            SmallMatrix beta_matrix;
            beta_matrix.smallMatrix(this->point0.get_x() - ray.get_pos().get_x(), this->point0.get_x() - this->point2.get_x(), ray.get_dir().get_x(),
                  this->point0.get_y() - ray.get_pos().get_y(), this->point0.get_y() - this->point2.get_y(), ray.get_dir().get_y(),
                  this->point0.get_z() - ray.get_pos().get_z(), this->point0.get_z() - this->point2.get_z(), ray.get_dir().get_z()
                  );
            myBeta = beta_matrix.determinant()/A.determinant();

            if (myBeta < 0 || myBeta > (1-myGamma)){
                  return false;
            }

            return true;
      }
      else { //shape has not been setup
            return false;
      }
}




/**************Light Class**************/
class Light {
      const static int POINTLIGHT = 0;
      const static int DIRECTIONALLIGHT = 1;
      //const static int AMBIENTLIGHT = 2;
public:
      int type;
      float x;
      float y;
      float z;
      Color color;
      int falloff;
      Light();
      makePointLight(float px, float py, float pz, Color c, int falloff);
      makeDirectionalLight(float dx, float dy, float dz, Color c);
      void generateLightRay(LocalGeo& local, Ray* lray, Color* lcolor);
};

Light::Light(){
      this->type = type;
      this->x = 0.0;
      this->y = 0.0;
      this->z = 0.0;
      this->color = Color();
      this->falloff = 0;
}

Light::makePointLight(float px, float py, float pz, Color c, int falloff){
      this->type = POINTLIGHT;
      this->x = px;
      this->y = py;
      this->z = pz;
      this->color = c;
      this->falloff = falloff;
}

Light::makeDirectionalLight(float dx, float dy, float dz, Color c){
      this->type = DIRECTIONALLIGHT;
      this->x = dx;
      this->y = dy;
      this->z = dz;
      this->color = c;
      this->falloff = 0;
}

void Light::generateLightRay(LocalGeo& local, Ray* lray, Color* lcolor){
      if (this->type==POINTLIGHT){
            // create point light ray
        Point p; 
		p = local.pos;
        Vector light_dir = this->ltp.PsubtractP(p);
        lray.ray(p, light_dir, 0.0001, FLT_MAX);
      } else if (this->type==DIRECTIONALLIGHT){
        Vector dirlight_dir = Vector(x,y,z);
        lray.ray(local.pos, dirlight_dir, 0.0001, FLT_MAX);
      }
}




