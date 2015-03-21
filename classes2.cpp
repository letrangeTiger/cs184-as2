#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <cstdlib>

//#include "classes.cpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>

#ifdef OSX


#include <time.h>
#include <math.h>

#endif

using namespace std;


#define PI 3.14159265  // Should be used from mathlib
inline float sqr(float x) { return x*x; }


#endif

class Vector;
class Normal;
class Point;
class Ray;
class SmallMatrix;
class Matrix;
class Color;
class Sample;
class LocalGeo;
class Shape;

float dot(Vector v0, Vector v1);
float length(Vector v0);
Point projectPointToLine(Point point, Ray ray);
float findScalarDivident(Vector v, Vector d);
void find_t(Ray ray, Point intersection, float thit);
int intersect3D_RayTriangle(Ray ray, Shape triangle, Point* intersection, Normal* normal, Point* thit);
int intersect3D_RaySphere(Ray ray, Shape sphere, Point* intersection, Normal* normal, Point* thit);
float dot(Vector v0);


class Vector {
      
public:
      float x;
      float y;
      float z;

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
      float subx = this->x + subtractee.get_x();
      float suby = this->y + subtractee.get_y();
      float subz = this->z + subtractee.get_z();
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
            float temp = sqrt(pow(this->x, 2)+pow(this->y, 2)+pow(this->z, 2));
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
      void normal();
      void normal(float a, float b, float c);
      void reset(float x, float y, float z);
      Normal add(Vector addee);
      Normal subtract(Vector subtractee);
      void printline();
      void print();
};
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
      float x;
      float y;
      float z;
public:
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
void Point::point(){
      x = 0;
      y = 0;
      z = 0;
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
      Point pos;
      Vector dir;
      float t_min;
      float t_max;
      float t;
public:
      void ray();
      void ray(Point a, Vector b, float c, float d);
      Point currentposition(float t);
      Point get_pos();
      Vector get_dir();
      float get_t_min();
      float get_t_max();
      void set_pos(Point pos);
      void set_dir(Vector dir);
      void printline();
      void print();
};
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
float Ray::get_t_min(){
      return this->t_min;
}
float Ray::get_t_max(){
      return this->t_max;
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


//for calculating intersections
class SmallMatrix {
public:
      float pos[3][3] = {
            {0,0,0},
            {0,0,0},
            {0,0,0}
      };
      void smallMatrix();
      void smallMatrix(float a, float b, float c, float d, float e, float f, float g, float h, float i);
      float determinant();
      void print();
};
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


class Matrix {
//successfully debugged, now working: constructors, translation, scaling, m*m multiplication, identity, inverse, transpose,
      
public:
      float pos[4][4] = {
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0},
            {0,0,0,0}};

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
Matrix Matrix::arbitrary_rotation(float axisX, float axisY, float axisZ, float theta){
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




/***********************Color Class*****************/
//TODO: May support conversion from xyz
class Color {
	float r;
	float g;
	float b;
public:
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

Color Color::addColors(Color addeeColor){
	Color newColor;
	newColor.r = this->r + addeeColor.get_r();
	newColor.g = this->g + addeeColor.get_g();
	newColor.b = this->b + addeeColor.get_b();
	return newColor;
}

Color Color::subColors(Color subtracteeColor){
	Color newColor;
	newColor.r = this->r + subtracteeColor.get_r();
	newColor.g = this->g + subtracteeColor.get_g();
	newColor.b = this->b + subtracteeColor.get_b();
	return newColor;
}

Color Color::mulColorbyScalar(float scalar){
	Color newColor;
	newColor.r = scalar*(this->r);
	newColor.g = scalar*(this->g);
	newColor.b = scalar*(this->b);
	return newColor;
}

Color Color::divColorbyScalar(float scalar){
	Color newColor;
	this->r = (this->r)/scalar;
	this->g = (this->g)/scalar;
	this->b = (this->b)/scalar;
	return newColor;
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

  

/*    stores screen coordinates */
class Sample {
	float x;
	float y;
public:
	void sample();
	void sample(float x, float y);
};
	void Sample::sample(){
		x = 0.0;
		y = 0.0;
	}
	void Sample::sample(float x, float y){
		this->x = x;
		this->y = y;
	}

	

class LocalGeo {
public:
      Point pos;
      Normal normal;
	LocalGeo();
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


class Intersection {
public:
      LocalGeo localGeo;
      Primitive primitive;
      void intersection();
      void intersection(LocalGeo localGeo, Primitive* primitive);
};
      void intersection(){
   }
      void intersection(LocalGeo localGeo, Primitive* primitive){
            this->localGeo = localGeo;
            this->primitive = primitive;
   }



//****************************************************
// TESTING SHAPE
//****************************************************
int main(int argc, char *argv[]) {
      // bool intersect_a, intersect_b, intersect_c, intersect_d, intersect_e;

      // //test a: triangle-ray intersection
      // printf("test a: triangle-ray intersection\n");

      // Point point0_a, point1_a, point2_a, ray_point_a;
      // Vector ray_dir_a, temp_ray_dir_a;
      // Shape triangle_a;
      // Ray ray_a;
      // float thit_a = 0;
      // LocalGeo localgeo_a = LocalGeo();

      // point0_a.point(0, 0, 0);
      // point1_a.point(0, 3, 0);
      // point2_a.point(3, 0, 0);


      // ray_point_a.point(1, 1, 3);
      // temp_ray_dir_a.vector(0, 0, -1);
      // ray_dir_a = temp_ray_dir_a.normalize();

      // ray_a.ray(ray_point_a, ray_dir_a, 0, 10000);

      // triangle_a.makeTriangle(point0_a, point1_a, point2_a);

      // intersect_a = triangle_a.intersect(ray_a, &thit_a, &localgeo_a);
      
      // localgeo_a.printline();
      // printf("%d\n", intersect_a);

      

      // //test b: sphere-ray intersection
      // printf("\n");
      // printf("test b: sphere-ray intersection\n");

      // Point center_b, ray_point_b;
      // float radius_b;
      // Vector ray_dir_b, temp_ray_dir_b;
      // Shape sphere_b;
      // Ray ray_b;
      // float thit_b;
      // LocalGeo localgeo_b;

      // center_b.point(0, 0, 0);
      // radius_b = 5;

      // ray_point_b.point(7, 2, 0);
      // temp_ray_dir_b.vector(-1, -1, 0);
      // ray_dir_b = temp_ray_dir_b.normalize();

      // ray_b.ray(ray_point_b, ray_dir_b, 0, 10000);

      // sphere_b.makeSphere(radius_b, center_b);

      // intersect_b = sphere_b.intersect(ray_b, &thit_b, &localgeo_b);

      // printf("%f\n", thit_b);

      // localgeo_b.printline();
      // printf("%d\n", intersect_b);
      

      // //test c: no intersection (triangle)
      // printf("\n");
      // printf("test c: no intersection (triangle)\n");

      // Point point0_c, point1_c, point2_c, ray_point_c;
      // Vector ray_dir_c, temp_ray_dir_c;
      // Shape triangle_c;
      // Ray ray_c;
      // float thit_c;
      // LocalGeo localgeo_c; 

      // point0_c.point(2, 0, 0);
      // point1_c.point(0, 0, 0);
      // point2_c.point(0, 2, 0);

      // ray_point_c.point(2, 2, 10);
      // temp_ray_dir_c.vector(0, 0, -1);
      // ray_dir_c = temp_ray_dir_c.normalize();

      // ray_c.ray(ray_point_c, ray_dir_c, 0, 10000);

      // triangle_c.makeTriangle(point0_c, point1_c, point2_c);

      // intersect_c = triangle_c.intersect(ray_c, &thit_c, &localgeo_c);

      // localgeo_c.printline();
      // printf("%d\n", intersect_c);


      // //test d: no intersection (sphere)
      // printf("\n");
      // printf("test d: no intersection (sphere)\n");

      // Point center_d, ray_point_d;
      // float radius_d;
      // Vector ray_dir_d, temp_ray_dir_d;
      // Shape sphere_d;
      // Ray ray_d;
      // float thit_d;
      // LocalGeo localgeo_d;

      // center_d.point(9, 9, 9);
      // radius_d = 1;

      // ray_point_d.point(11, 11, 0);
      // temp_ray_dir_d.vector(1, 1, 0);
      // ray_dir_d = temp_ray_dir_d.normalize();

      // ray_d.ray(ray_point_d, ray_dir_d, 0, 10000);

      // sphere_d.makeSphere(radius_d, center_d);

      // intersect_d = sphere_d.intersect(ray_d, &thit_d, &localgeo_d);

      // localgeo_d.printline();
      // printf("%d\n", intersect_d);
}

