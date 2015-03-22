#include "class.cpp"

class Sampler{
public:
	unsigned width, height;
	float current_u, current_v;
	Sampler();
	Sampler(float width, float height);
	bool generateSample(Sample* sample);
};

Sampler::Sampler(){
	this->width = 1000;
	this->height = 500;
	current_u = 0.5;
	current_v = 0.5;
}
Sampler::Sampler(float width, float height){
	this->width = width;
	this->height = height;
	current_u = 0.5;
	current_v = 0.5;
}
bool Sampler::generateSample(Sample* sample){
	cout << "generateSample called";
	if ((current_u > width) && (this->current_v > this->height)) {
		cout << "returning false";
		return false;
	} else if(this->current_u > this->width){
		sample->x = current_u;
		sample->y = current_v;
		this->current_u = 0.5;
		this->current_v += 1.0;
		cout << "returning the first";
		return true;
	} else {
		sample->x = current_u;
		sample->y = current_v;
		this->current_u += 1.0;
		cout << "returning the second";
		return true;
	}
}



