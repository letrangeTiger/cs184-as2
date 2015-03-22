#include <iostream>
#include <cmath>
#include <vector>
#include "lodepng.cpp" // from "lodev.org/lodepng/"
#include "class.cpp"

using namespace std;

class Film {
	private:
		std::vector<unsigned char> pixColors;
		unsigned width;
		unsigned height;	
	public:
		Film(unsigned width, unsigned height);
		void commit(Sample& sample, Color& color);
		void writeImage();
};

Film::Film(unsigned w, unsigned h){
	this->width = w;
	this->height = h;
	pixColors.resize(this->width*this->height*4);
	/* no need to initialize to black screen
	for(unsigned y = 0; y < height; y++)
		for(unsigned x = 0; x < width; x++) {
		image[4 * width * y + 4 * x + 0] = 0;
		image[4 * width * y + 4 * x + 1] = 0;
		image[4 * width * y + 4 * x + 2] = 0;
		image[4 * width * y + 4 * x + 3] = 0;
		}*/
}

void Film::commit(Sample& sample, Color& color){
	unsigned x = (unsigned) floor(sample.x);
	unsigned y = (unsigned) floor(sample.y);
	pixColors[4*width*y + 4*x + 0] = color.get_r()*255;
	pixColors[4*width*y + 4*x + 1] = color.get_g()*255;
	pixColors[4*width*y + 4*x + 2] = color.get_b()*255;
	pixColors[4*width*y + 4*x + 3] = 255; // transparency: 255=completely opaque
}

void Film::writeImage(){ 
	const char *filename = "result.png";
	lodepng::State state;
	std::vector<unsigned char> png;
	unsigned error = lodepng::encode(png, pixColors, width, height, state);
	lodepng::save_file(png, filename);
	if(error) std::cout << "png encode error " << error << std::endl;
}

//int main(int argc, char *argv[]){
	/*
	unsigned w=1000;
	unsigned h=500;
	Film film(w,h);
	for(unsigned i=0; i<w; i++)
		for(unsigned j=0; j<h; j++){
			film.commit(Sample(i,j),Color());
		}
	film.writeImage();*/
  
	//example encode from lodev.org/lodepng/
		/*
	unsigned width = 1000, height = 500;
	const char *filename = "result.png";
	std::vector<unsigned char> image;
	std::vector<unsigned char> png;
	image.resize(width * height * 4);
	for(unsigned y = 0; y < height; y++)
	for(unsigned x = 0; x < width; x++)
	{
	image[4 * width * y + 4 * x + 0] = 255*!(x&y);
	image[4 * width * y + 4 * x + 1] = x|y;
	image[4 * width * y + 4 * x + 2] = x^y;
	image[4 * width * y + 4 * x + 3] = 255;
	}
	lodepng::State state;
	lodepng::encode(png, image, width, height,state);
	lodepng::save_file(png, filename);
	*/

	//return 0;

//}
