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
      BRDF *brdf;
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
      
      float thit= 0.0;
      Intersection in = Intersection();
      /* Find the nearest shape that the ray intersects, if any */

      //no need to continue if reached maxrecursiondepth
      if (depth > this->maxrecursiondepth) {
         *color = Color(0,0,0);
      }
      //if does not intersect anything return black
      if(!primitives.intersect(ray, &thit, &in)){
       // cout << "did not intersect!";
         *color = Color(0,0,0);

      }    
      //find BRDF at intersection point
      in.primitive->getBRDF(in.localGeo, brdf);
            
      

      
      //loop through lights
      int i = 0;
      for(iter = lights.begin(); iter < lights.end(); iter++, i++) {
        Ray* currentray;
        Color* lcolor;
        lights.at(i).generateLightRay(in.localGeo, currentray, lcolor);
        if (!primitives.intersectP(*currentray)) {
            Vector n = in.localGeo.normal;
            n.printline();
            Vector l = currentray->get_dir().normalize(); 
            float NdotL = n.dot(l);
            Vector r = l.reverse().add(n.scalarmultiply(2*NdotL));
            r = r.normalize();
            Vector v = (eye.PsubtractP(in.localGeo.get_pos())).normalize(); 
            float fmx = fmax(NdotL, 0.0f);
            
            //float test = brdf->kdr*lcolor->get_r()*fmx;
            Color diffuse_comp = Color(brdf->kdr * lcolor->get_r()*fmx, brdf->kdg*lcolor->get_g()*fmx, brdf->kdb*lcolor->get_b()*fmx);
            
            float RdotV = r.dot(v);
            float rdv = fmax(RdotV,0.0f);
            Color spec_comp = Color(brdf->ksr*lcolor->get_r()*rdv, brdf->ksg*lcolor->get_g()*rdv, brdf->ksb*lcolor->get_b()*rdv);

            Color ambient_comp = Color(brdf->kar*lcolor->get_r(), brdf->kag*lcolor->get_g(), brdf->kab*lcolor->get_b());

            *color = *color + diffuse_comp + spec_comp + ambient_comp;

        }else{
            *color = *color + Color(brdf->kar*lcolor->get_r(), brdf->kag*lcolor->get_g(), brdf->kab*lcolor->get_b());
        }
    }
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
