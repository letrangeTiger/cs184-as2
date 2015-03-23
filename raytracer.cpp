#include "class.cpp"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

class RayTracer{

public:
      std::vector<Light> lights;
      std::vector<Light>::iterator iter;
      AggregatePrimitive primitives;
      int maxrecursiondepth;
      Point eye;
      BRDF brdf;
      RayTracer();
      RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, std::vector<Light> lights);
      void trace(Ray& ray, int depth, Color* color);
};
RayTracer::RayTracer(){

}
RayTracer::RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, std::vector<Light> lights){
      this->maxrecursiondepth = maxrecursiondepth;
      this->eye = eye;
      this->primitives = primitives;
      this->lights = lights;
}

void RayTracer::trace(Ray& ray, int depth, Color* color){
      //cout << "inside raytracer loop";
      float thit= 0.0;
      Intersection in ;
      /* Find the nearest shape that the ray intersects, if any */

      //no need to continue if reached maxrecursiondepth
      if (depth > this->maxrecursiondepth) {
         *color = Color(0,0,0);
      }
      //if does not intersect anything return black
      else if(!primitives.intersect(ray, &thit, &in)){
         *color = Color(0,0,0);
      } else {  

      //find BRDF at intersection point
      in.primitive->getBRDF(in.localGeo, &brdf);
      //in.localGeo.printline();   
      

     // cout << "before lights loop gggggggggggggggggggggggggggggggggggggggg";
      for(auto light : lights) {
        //cout << "inside lights loop";
        Ray currentray;
        Color lcolor;
        //cout << "before generate LIGHT";
        light.generateLightRay(in.localGeo, &currentray, &lcolor);
        //currentray.get_dir().printline();
        if (!primitives.intersectP(currentray)) {

            Vector n = in.localGeo.normal;
            //n.printline();
            Vector l = currentray.get_dir().normalize();
            //l.printline();
            float NdotL = n.dot(l);
            printf("%f", NdotL);
            Vector r = l.reverse().add(n.scalarmultiply(2*NdotL));

            r = r.normalize();
            //r.printline();
            Vector v = (eye.PsubtractP(in.localGeo.get_pos())).normalize(); 
            //v.printline();
            //n.printline();
            //l.printline();
            //r.printline();
            float fmx = fmax(NdotL, 0.0f);

            float test = brdf.kdr;//*lcolor.get_r()*fmx;


            Color diffuse_comp = Color(brdf.kdr * lcolor.get_r()*fmx, brdf.kdg*lcolor.get_g()*fmx, brdf.kdb*lcolor.get_b()*fmx);
            
            float RdotV = r.dot(v);
            float rdv = fmax(RdotV,0.0f);
            Color spec_comp = Color(brdf.ksr*lcolor.get_r()*rdv, brdf.ksg*lcolor.get_g()*rdv, brdf.ksb*lcolor.get_b()*rdv);

            Color ambient_comp = Color(brdf.kar*lcolor.get_r(), brdf.kag*lcolor.get_g(), brdf.kab*lcolor.get_b());

            *color = *color + diffuse_comp + spec_comp + ambient_comp;
            

            //cout << "this has run";
        } else {
            //cout << "else";
            *color = *color + Color(brdf.kar*lcolor.get_r(), brdf.kag*lcolor.get_g(), brdf.kab*lcolor.get_b());
        }
    }}
      /*if(brdf->krr > 0 || brdf->krg > 0 || brdf->krb > 0){
        

        Vector n = in.localGeo.normal;
        Vector l = ray.get_dir().reverse().normalize();

        float NdotL = n.dot(l);
        Ray reflectRay = Ray(in.localGeo.pos, l.reverse().add(n.scalarmultiply(2*NdotL)), 0.001f, FLT_MAX);
        
        //Make a recursive call to trace the reflected ray
        Color temp = Color(0,0,0);
        trace(reflectRay, depth+1, &temp);
        *color = *color + Color(brdf->krr*temp.get_r(), brdf->krg*temp.get_g(), brdf->krb*temp.get_b());
    }*/
  
}