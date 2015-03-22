#include "class.cpp"
#include <iostream>

class Camera{
	Point eye,ll,lr,ul,ur;
	Camera(Point eye,Point ll,Point lr,Point ul,Point ur);
	void generateRay(LocalGeo& local, Ray* ray)
};

Camera::Camera(Point eye, Point ll, Point lr, Point ul, Point ur){
	this->eye = eye;
	this->ll = ll;
	this->lr = lr;
	this->ul = ul;
	this->ur = ur;
}

void Camera::generateRay(Sample& sample, Ray* ray){
	float u = sample.get_x();
	float v = sample.get_y();
	float px = u*(v*ll.x+(1.0-v)*ul.x)+ (1.0-u)*(v*lr.x+(1.0-v)*ur.x);
	float py = u*(v*ll.y+(1.0-v)*ul.y)+ (1.0-u)*(v*lr.y+(1.0-v)*ur.y);
	float pz = u*(v*ll.z+(1.0-v)*ul.z)+ (1.0-u)*(v*lr.z+(1.0-v)*ur.z);
	Point p = Point(px,py,pz);
	Vector cam_dir = p.PsubtractP(this->eye);
	ray->ray(this->eye, cam_dir, 0.0, FLT_MAX);
}