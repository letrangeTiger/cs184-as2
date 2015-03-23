#include "class.cpp"
#include <iostream>

class Camera{
public:
	Point eye,ll,lr,ul,ur;
	float width, height;
	Camera(Point eye,Point ll,Point lr,Point ul,Point ur, float width, float height);
	void generateRay(Sample& sample, Ray* ray);
};

Camera::Camera(Point eye, Point ll, Point lr, Point ul, Point ur, float width, float height){
	this->eye = eye;
	this->ll = ll;
	this->lr = lr;
	this->ul = ul;
	this->ur = ur;
	this->width = width;
	this->height = height;
}

void Camera::generateRay(Sample& sample, Ray* ray){
	float u = ll.x+ (lr.x-ll.x)*sample.x/width;
	float v = lr.y+ (ur.y-lr.y)*sample.y/height;
	//float px = u*(v*ll.x+(1.0-v)*ul.x)+ (1.0-u)*(v*lr.x+(1.0-v)*ur.x);
	//float py = u*(v*ll.y+(1.0-v)*ul.y)+ (1.0-u)*(v*lr.y+(1.0-v)*ur.y);
	//float pz = u*(v*ll.z+(1.0-v)*ul.z)+ (1.0-u)*(v*lr.z+(1.0-v)*ur.z);
	Point p = Point(u,v,ll.z);
	//p.printline();
	Vector cam_dir = p.PsubtractP(this->eye);
	ray->ray(this->eye, cam_dir, 0.0, FLT_MAX);
}