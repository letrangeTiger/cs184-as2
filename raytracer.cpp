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
      Light amblight;
      RayTracer();
      RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, std::vector<Light> lights, Light amblight);
      void trace(Ray& ray, int depth, Color* color);
      Ray createReflectRay(LocalGeo &local, Ray &ray);
};
RayTracer::RayTracer(){

}
RayTracer::RayTracer(int maxrecursiondepth, Point eye, AggregatePrimitive primitives, std::vector<Light> lights, Light amblight){
      this->maxrecursiondepth = maxrecursiondepth;
      this->eye = eye;
      this->primitives = primitives;
      this->lights = lights;
      this->amblight = amblight;
}

void RayTracer::trace(Ray& ray, int depth, Color* color){
      float thit;
      Intersection in;
      /* Find the nearest shape that the ray intersects, if any */

      //no need to continue if reached maxrecursiondepth
      if (depth > this->maxrecursiondepth) {
         *color = Color(0,0,0);
      }
      //if does not intersect anything return black
      else if(!primitives.intersect(ray, &thit, &in)){
         *color = Color(0,0,0);
      }else {

      //find BRDF at intersection point
      BRDF brdf;
      in.primitive->getBRDF(in.localGeo, &brdf);
      //in.localGeo.printline();   
     //brdf.printline();
      
      for(auto light : lights) {
        Ray currentray;
        Color lcolor;
        //cout << "before generate LIGHT";
        light.generateLightRay(in.localGeo, &currentray, &lcolor);
        //currentray.get_dir().printline();
        //lcolor.print();
        printf("%lu\n", primitives.primitives.size());
        cout << "Fdsfasdfsdfasdfds";
        if (!primitives.intersectP(currentray)) {
          //cout << "in raytracer loop";
            Vector n = in.localGeo.normal;
            Vector l =  currentray.get_dir().normalize();
            float NdotL = n.dot(l);
            //NdotL = - NdotL;
            //cout <<NdotL<<"\n";
            Vector r = l.reverse().add(n.scalarmultiply(2*NdotL));
            r = r.normalize();
            //r.printline();
            Vector v = (eye.PsubtractP(in.localGeo.get_pos())).normalize(); 
            float fmx = fmax(NdotL, 0.0f);
            //cout << fmx << "\n";

            Color diffuse_comp = Color(brdf.kdr * lcolor.get_r()*fmx, brdf.kdg*lcolor.get_g()*fmx, brdf.kdb*lcolor.get_b()*fmx);
            //diffuse_comp.print();
            float RdotV = r.dot(v);
            //RdotV = - RdotV;
            float rdv = pow(fmax(RdotV,0.0f),brdf.p);
            Color spec_comp = Color(brdf.ksr*lcolor.get_r()*rdv, brdf.ksg*lcolor.get_g()*rdv, brdf.ksb*lcolor.get_b()*rdv);

            Color ambient_comp = Color(brdf.kar*amblight.color.get_r(), brdf.kag*amblight.color.get_g(), brdf.kab*amblight.color.get_b());

            *color = *color + diffuse_comp + spec_comp + ambient_comp;

            //color->print();

            // printf("%f\n", brdf.kar);
            // printf("%f\n", brdf.kag);
            // printf("%f\n", brdf.kab);

            
            //cout << "this has run";
        }else{
            //cout << "else";

            *color = *color + Color(brdf.kar*amblight.color.get_r(), brdf.kag*amblight.color.get_g(), brdf.kab*amblight.color.get_b());

        }
    }
      if(brdf.krr > 0 || brdf.krg > 0 || brdf.krb > 0){
        

        Vector n = in.localGeo.normal;
        Vector l = ray.get_dir().reverse().normalize();

        float NdotL = n.dot(l);
        Ray reflectRay = Ray(in.localGeo.pos, l.reverse().add(n.scalarmultiply(2*NdotL)), 0.001f, FLT_MAX);
        
        //Make a recursive call to trace the reflected ray
        Color temp = Color(0,0,0);
        trace(reflectRay, depth+1, &temp);
        *color = *color + Color(brdf.krr*temp.get_r(), brdf.krg*temp.get_g(), brdf.krb*temp.get_b());
    }
  
}
}