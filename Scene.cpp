#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <sstream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "class.cpp"
#include "Film.cpp"
#include "Sampler.cpp"
#include "Camera.cpp"
#include "raytracer.cpp"
#include <algorithm>


using namespace std;

class Scene {
	public:		
		Point eye;
		Point ll;
		Point lr;
		Point ul;
		Point ur;
		int max_depth;
		AggregatePrimitive aggreprim;
		std::vector<Light> lights;
		Light amblight;
		unsigned width, height;
		Scene();
		Scene(unsigned width, unsigned height, int max_depth);
		void render();
};

Scene::Scene(){
	this->width = 1000;
	this->height = 500;
	this->max_depth = 5;
}

Scene::Scene(unsigned width, unsigned height, int max_depth) {
	this->width = width;
	this->height = height;
	this->max_depth = max_depth;
}

void Scene::render() {
	// implement this.
	Film film = Film(this->width,this->height);
	Sampler sampler = Sampler(this->width,this->height);
	Sample* sample;
	Camera camera = Camera(eye,ll,lr,ul,ur);
	RayTracer raytracer = RayTracer(max_depth, eye, aggreprim, lights);
	while (!sampler.generateSample(sample)){
		Ray *camray;
		Color *color;
		camera.generateRay(*sample, camray);
		raytracer.trace(*camray, 0, color);
		Color c = Color(color->r, color->g, color->b);
		film.commit(*sample, c);
	}

	film.writeImage();
}

int main(int argc, char *argv[]) {

	unsigned w = 1000;
	unsigned h = 1000;
	int maxdepth = 5;
	Scene scene = Scene(w,h,maxdepth);

	int counter = 1;
	Matrix trans_mat = Matrix();
	trans_mat = trans_mat.identity();
	BRDF *brdf;

	cout << "once";

  	while (counter<argc){
  		cout << "10000";
		std::string arg = argv[counter];
	  	if (!strncmp(argv[counter], "cam", 15)) {
	  		cout << "cam";
			float ex = atof(argv[counter+1]);
			float ey = atof(argv[counter+2]);
			float ez = atof(argv[counter+3]);
			float llx = atof(argv[counter+4]);
			float lly = atof(argv[counter+5]);
			float llz = atof(argv[counter+6]);
			float lrx = atof(argv[counter+7]);
			float lry = atof(argv[counter+8]);
			float lrz = atof(argv[counter+9]);
			float ulx = atof(argv[counter+10]);
			float uly = atof(argv[counter+11]);
			float ulz = atof(argv[counter+12]);
			float urx = atof(argv[counter+13]);
			float ury = atof(argv[counter+14]);
			float urz = atof(argv[counter+15]);
			counter = counter+16;	
			Point eye,ll,lr,ul,ur;
			eye.point(ex,ey,ez);
			ll.point(llx,lly,llz);
			lr.point(lrx,lry,lrz);
			ul.point(ulx,uly,ulz);
			ur.point(urx,ury,urz);
			scene.eye = eye;
			scene.ll = ll;
			scene.lr = lr;
			scene.ul = ul;
			scene.ur = ur;
	  	} else if (!strncmp(argv[counter], "sph", 4)){
	  		cout << "sph";
			float cx = atof(argv[counter+1]);
			float cy = atof(argv[counter+2]);
			float cz = atof(argv[counter+3]);
			float r = atof(argv[counter+4]);
			counter = counter+5;
			Shape *sphere;
			sphere->makeSphere(r, Point(cx,cy,cz));
			GeometricPrimitive geoprim = GeometricPrimitive(sphere,Transformation(trans_mat), brdf);
			scene.aggreprim.addPrimitive(&geoprim);
	    } else if (!strncmp(argv[counter], "tri", 9)){
	    	cout << "tri";
			float ax = atof(argv[counter+1]);
			float ay = atof(argv[counter+2]);
			float az = atof(argv[counter+3]);
			float bx = atof(argv[counter+4]);
			float by = atof(argv[counter+5]);
			float bz = atof(argv[counter+6]);
			float cx = atof(argv[counter+7]);
			float cy = atof(argv[counter+8]);
			float cz = atof(argv[counter+9]);	
			counter = counter+10;
			Shape *triangle;
			triangle->makeTriangle(Point(ax,ay,az), Point(bx,by,bz), Point(cx,cy,cz));
			GeometricPrimitive geoprim = GeometricPrimitive(triangle,Transformation(trans_mat), brdf);
			scene.aggreprim.addPrimitive(&geoprim);
	    } else if (arg=="obj"){
	    	/*
			NOTE: .obj file supporting:
			|	v [x_coord] [y_coord] [z_coord]
			|	f [v_num1] [v_num2] [v_num3]

				for now...
	    	*/
			std::string filename = argv[counter+1];

		    std::ifstream obj_file(filename);
		    string line;
		    string lineName;

		    std::vector<Point> points;

			while (std::getline(obj_file, line)){
		    	if (line.empty() || line[0] == '#'){	//skip empty line or comment and move onto next line
		    		continue;
		    	}
		    	std::istringstream buf(line);
		    	std::istream_iterator<std::string> beg(buf), end;
		    	std::vector<std::string> tokens(beg,end);
		    	lineName = tokens[0];

		    	if (lineName == "s"){
		    		continue; //smooth shading not enabled
		    	}
		    	if (lineName == "o"){
		    		continue; //object name - do nothing
		    	}
		    	if (lineName == "vn"){
		    		continue; // vertexnormals not enabled
		    	}
		    	if (lineName == "v"){	//vertex
		            float x,y,z;
		            string s;
		            s = tokens[1];
		            istringstream(s) >> x;
		            s = tokens[2];
		            istringstream(s) >> y;
		            s = tokens[3];
		            istringstream(s) >> z;
		            Point p = Point(x,y,z);
		            points.push_back(p);
      			}

      			if (lineName == "f"){	//face
      				int a,b,c;
      				string s;
		            s = tokens[1];
		            istringstream(s) >> a;
		            s = tokens[2];
		            istringstream(s) >> b;
		            s = tokens[3];
		            istringstream(s) >> c;
      				Shape* triangle;
      				triangle->makeTriangle(points[a-1],points[b-1],points[c-1]);
      				GeometricPrimitive geoprim = GeometricPrimitive(triangle, Transformation(trans_mat), brdf);
      				scene.aggreprim.addPrimitive(&geoprim);
      				//TODO: 3 items
      					// Shape triangle;
						// triangle.makeTriangle(Point(ax,ay,az), Point(bx,by,bz), Point(cx,cy,cz));
						// GeometricPrimitive geoprim(trans_mat, trans_mat.inverse(), triangle, brdf);
      				

      			} else continue;
      		}	
      			counter = counter+2;
	    } else if (!strncmp(argv[counter], "ltp", 7)){
	    	cout << "ltp";
	    	float px = atof(argv[counter+1]);
	    	float py = atof(argv[counter+2]);
	    	float pz = atof(argv[counter+3]);
	    	float r = atof(argv[counter+4]);
	    	float g = atof(argv[counter+5]);
	    	float b = atof(argv[counter+6]);
	    	int falloff = atof(argv[counter+7]);
			counter=counter+8;
			Light ptlight;
			ptlight.makePointLight(px,py,pz,Color(r,g,b),falloff); 
			scene.lights.push_back(ptlight);
	    } else if (!strncmp(argv[counter], "ltd", 6)){
	    	cout << "ltd";
	    	float dx = atof(argv[counter+1]);
	    	float dy = atof(argv[counter+2]);
	    	float dz = atof(argv[counter+3]);
	    	float r = atof(argv[counter+4]);
	    	float g = atof(argv[counter+5]);
	    	float b = atof(argv[counter+6]);	  		
	  		counter=counter+7;
	  		Light drlight;
	  		drlight.makeDirectionalLight(dx,dy,dz,Color(r,g,b));
	  		scene.lights.push_back(drlight);
	  	} else if (!strncmp(argv[counter], "lta", 3)){
	  		cout << "lta";
	  		float r = atof(argv[counter+1]);
	    	float g = atof(argv[counter+2]);
	    	float b = atof(argv[counter+3]);
	    	counter=counter+4;
	    	scene.amblight = Light(r,g,b);

	    } else if (!strncmp(argv[counter], "xft", 3)){
	    	cout << "xft";
	    	float tx = atof(argv[counter+1]);
	    	float ty = atof(argv[counter+2]);
	    	float tz = atof(argv[counter+3]);
	    	counter=counter+4;
	    	trans_mat = trans_mat * translation(tx,ty,tz);
	    } else if (!strncmp(argv[counter], "xfr", 3)){
	    	cout << "xfr";
	    	float rx = atof(argv[counter+1]);
	    	float ry = atof(argv[counter+2]);
	    	float rz = atof(argv[counter+3]);
	    	counter=counter+4;
	    	trans_mat = trans_mat * rotation(rx,ry,rz);
	    } else if (!strncmp(argv[counter], "xfs", 3)){
	    	cout << "xfs";
	    	float sx = atof(argv[counter+1]);
	    	float sy = atof(argv[counter+2]);
	    	float sz = atof(argv[counter+3]);
	    	counter=counter+4;
	    	trans_mat = trans_mat * scaling(sx,sy,sz);
	    } else if (!strncmp(argv[counter], "mat", 13)){
	    	float kar = atof(argv[counter+1]);
	    	float kag = atof(argv[counter+2]);
	    	float kab = atof(argv[counter+3]);
	    	float kdr = atof(argv[counter+4]);
	    	float kdg = atof(argv[counter+5]);
	    	float kdb = atof(argv[counter+6]);
	    	float ksr = atof(argv[counter+7]);
	    	float ksg = atof(argv[counter+8]);
	    	float ksb = atof(argv[counter+9]);
	    	float p = atof(argv[counter+10]);
	    	float krr = atof(argv[counter+11]);
	    	float krg = atof(argv[counter+12]);
	    	float krb = atof(argv[counter+13]);
	    	//brdf.setKa(kar,kag, kab);
           // brdf.setKr(krr,float g, float b);
    //  void setKd(float r,float g, float b);
      //void setKs(float r,float g, float b, float p);
      //void setP(float p);
	    	*brdf = BRDF(kar,kag,kab,
	    				kdr,kdg,kdb,
	    				ksr,ksg,ksb,p,
	    				krr,krg,krb);
	    	counter=counter+14;
	    } 
	    else if (!strncmp(argv[counter], "xfz", 0)){
	    	cout << "xfz";
			trans_mat = Matrix();
			trans_mat = trans_mat.identity();	    	
			counter+=1;
		}else {
	  		counter+=1;
		}
	
	}

	//scene.render();
	return 0;
}


