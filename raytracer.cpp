#include "class.cpp"
#include <iostream>

class RayTracer{

public:
      std::list<Light*> lights;
      AggregatePrimitive primitives;
      int maxrecursiondepth;
      Point eye;
      RayTracer();
      RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, std::list<Light*> lights);
      void trace(Ray& ray, int depth, Color* color);
};

RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, List<Lights*> lights){
      this->maxrecursiondepth = maxrecursiondepth;
      this->eye = eye;
      this->primitives = primitives;
      this->lights = lights;
}

void RayTracer::trace(Ray& ray, int depth, Color* color){
      BRDF brdf;
      float* thit;
      Intersection in = Intersection();
      /* Find the nearest shape that the ray intersects, if any */

      //no need to continue if reached maxrecursiondepth
      if (depth > this->maxrecursiondepth) {
         *color.color(0,0,0);
         return;
     }
      //if does not intersect anything return black
      if(!primitives.intersect(ray, &thit, &in)){
        *color = {0,0,0};
        return;
  }
      //find BRDF at intersection point
      in.primitive->getBRDF(in.localGeo, &brdf);

      //loop through lights
      for(int i = 0; i < lights.size(); i++) {
        Ray* currentray = Ray();
        Color* lcolor = Color();
        lights[i].generateLightRay(in.localGeo, &currentray, &lcolor);
        if (!primitives.intersectP(currentray)) {
            Vector n = in.localGeo.normal();
            Vector l = currentlight.get_dir().normalize(); //surface->light
            float NDotL = n.dot(l);
            Vector r = -l + 2*NDotL*n; //The reflection ray
            r = r.normalize();
            Vector v = (eye.PsubtractP(in.localGeo.get_pos())).normalize(); //The view vector

            Color diffuse_comp = Color(brdf.kdr * lcolor.get_r()*max<float>(NDotL, 0), brdf.kdg*lcolor.get_g()*max<float>(NDotL, 0), ), brdf.kdb*lcolor.get_b()*max<float>(NDotL, 0));

            Color spec_comp = Color(brdf.ksr*lcolor.get_r()*pow(max<float>(r.dot(v),0), brdf.p), brdf.ksg*lcolor.get_g()*pow(max<float>(r.dot(v),0), brdf.p), brdf.ksb*lcolor.get_b()*pow(max<float>(r.dot(v),0), brdf.p));

            Color ambient_comp = brdf.ka*lcolor; = Color(brdf.kar*lcolor.get_r(), brdf.kag*lcolor.get_g(), brdf.kab*lcolor.get_b());

            *color += diffuse_comp + spec_comp + ambient_comp;

        }else{
            *color += Color(brdf.kar*lcolor.get_r(), brdf.kag*lcolor.get_g(), brdf.kab*lcolor.get_b());
        }
    }
      if(brdf.krr > 0 || brdf.krg > 0 || brdf.krb > 0){
        

        Vector n = in.localGeo.normal();
        Vector l = -ray.direction.normalize();
        

        Ray reflectRay = Ray(in.localGeo.position(), -l + 2*n.dot(l)*n, 0.001f, FLT_MAX);
        
        //Make a recursive call to trace the reflected ray
        Color temp = Color(0,0,0);
        trace(reflectRay, depth+1, &temp);
        *color += Color(brdf.krr*temp.get_r, brdf.krg*temp.get_g, brdf.krb*temp.get_b);
    }

}