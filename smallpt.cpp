// smallpt, a Path Tracer by Kevin Beason, 2008; Updated by Benjamin Stürz
// Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
//        Remove "-fopenmp" for g++ version < 4.2
// Usage: time ./smallpt 5000 && xv image.ppm
#include <algorithm>
#include <concepts>
#include <numbers>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <array>
#include "random.hpp"
#include "string.hpp"

// Define as 1 if you don't want runtime-support.
#define NORT 0

constexpr int width = 10, height = 10, num_samps = 1;


constexpr auto my_abs(auto n) noexcept {
   return n < 0.0 ? -n : n;
}
template<std::floating_point T>
constexpr T my_sin(T x) noexcept {
   if (std::is_constant_evaluated()) {
      const T t = x - (long)x;
      if (t >= 0.5) {
         return -16 * t * t + 8 * t;
      } else {
         return -16 * t * t - 24 * t + 8;
      }
      return t;
   } else {
      return std::sin(x);
   }
}
constexpr auto my_cos(auto x) noexcept {
   if (std::is_constant_evaluated()) {
      return my_sin(x + std::numbers::pi_v<decltype(x)> / 2);
   } else {
      return std::cos(x);
   }
}

constexpr auto my_sqrt(auto x) noexcept {
   if (std::is_constant_evaluated()) {
      using T = decltype(x);
      auto n = x / 2;
      auto lstX = T{};
      while (n != lstX) {
         lstX = n;
         n = (n + x / n) / 2;
      }
      return n;
   } else {
      return std::sqrt(x);
   }
}

struct Vec {
   // position, also color (r,g,b)
   double x, y, z;

   constexpr Vec(double x = 0, double y = 0, double z = 0) noexcept
      : x(x), y(y), z(z) {}

   [[nodiscard]]
   constexpr Vec operator+(const Vec &b) const noexcept {
      return Vec(x + b.x, y + b.y, z + b.z);
   }

   [[nodiscard]]
   constexpr Vec operator-(const Vec &b) const noexcept {
      return Vec(x - b.x, y - b.y, z - b.z);
   }

   [[nodiscard]]
   constexpr Vec operator*(double b) const noexcept {
      return Vec(x * b, y * b, z*b);
   }

   [[nodiscard]]
   constexpr Vec mult(const Vec &b) const noexcept {
      return Vec(x * b.x, y * b.y, z * b.z);
   }

   constexpr Vec& norm() noexcept {
      return *this = *this * (1 / my_sqrt(x*x + y*y + z*z));
   }

   [[nodiscard]]
   constexpr double dot(const Vec &b) const noexcept {
      return x*b.x+y*b.y+z*b.z;
   }

   // cross product
   [[nodiscard]]
   constexpr Vec operator%(const Vec&b) const noexcept {
      return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
   }
};

struct Ray {
   Vec o, d;
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()

struct Sphere {
   double rad;       // radius
   Vec p, e, c;      // position, emission, color
   Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

   // returns distance, 0 if nohit
   constexpr double intersect(const Ray &r) const noexcept {
      // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
      const Vec op = p - r.o;
      const double eps = 1e-4;
      const double b = op.dot(r.d);
      double det = b * b-op.dot(op) + rad * rad;
      double t;
      if (det < 0)
         return 0;
      det = my_sqrt(det);
      return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
   }
};
//Scene: radius, position, emission, color, material
constexpr Sphere spheres[] = {
  Sphere{1e5,  Vec{  1e5+1,40.8,81.6  }, Vec{        }, Vec{.75,.25,.25 },       DIFF}, // Left
  Sphere{1e5,  Vec{ -1e5+99,40.8,81.6 }, Vec{        }, Vec{.25,.25,.75 },       DIFF}, // Rght
  Sphere{1e5,  Vec{ 50,40.8, 1e5      }, Vec{        }, Vec{.75,.75,.75 },       DIFF}, // Back
  Sphere{1e5,  Vec{ 50,40.8,-1e5+170  }, Vec{        }, Vec{            },       DIFF}, // Front
  Sphere{1e5,  Vec{ 50, 1e5, 81.6     }, Vec{        }, Vec{.75,.75,.75 },       DIFF}, // Bottom
  Sphere{1e5,  Vec{ 50,-1e5+81.6,81.6 }, Vec{        }, Vec{.75,.75,.75 },       DIFF}, // Top
  Sphere{16.5, Vec{ 27,16.5,47        }, Vec{        }, Vec{1,1,1       } *.999, SPEC}, // Mirror
  Sphere{16.5, Vec{ 73,16.5,78        }, Vec{        }, Vec{1,1,1       } *.999, REFR}, // Glass
  Sphere{600,  Vec{ 50,681.6-.27,81.6 }, Vec{12,12,12}, Vec{            },       DIFF}  // Lite
};

[[nodiscard]]
constexpr auto clamp(auto x) noexcept {
   return x < 0 ? 0 : x > 1 ? 1 : x;
}

constexpr int toInt(double x) noexcept {
   return int(std::pow(clamp(x), 1 / 2.2) * 255 + 0.5);
}
constexpr bool intersect(const Ray &r, double &t, int &id){
   double n = sizeof(spheres) / sizeof(Sphere);
   double d;
   double inf=t=1e20;
   for(int i = int(n); i--;) {
      if ((d = spheres[i].intersect(r)) && d < t) {
         t = d;
         id = i;
      }
   }
   return t < inf;
}
constexpr Vec radiance(const Ray &r, int depth, auto& rnd) noexcept {
   double t;                               // distance to intersection
   int id=0;                               // id of intersected object

   // if miss, return black
   if (!intersect(r, t, id))
    return Vec{};

   // the hit object
   const Sphere &obj = spheres[id];

   const Vec x  = r.o + r.d * t;
   const Vec n  = (x - obj.p).norm();
   const Vec nl = n.dot(r.d) < 0 ? n : (n * -1);
   Vec f = obj.c;

   // max refl
   const double p = (f.x > f.y && f.x > f.z) ? f.x : (f.y > f.z ? f.y : f.z); 

   if (++depth > 5) {
      if (rnd() < p)
         f = f * (1 / p);
      else return obj.e; //R.R.
   }

   if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
      const double r1 = 2 * M_PI * rnd();
      const double r2 = rnd();
      const double r2s = my_sqrt(r2);
      const Vec w = nl;
      const Vec u = ((my_abs(w.x) > 0.1 ? Vec{0,1} : Vec{1}) % w).norm();
      const Vec v = w % u;
      const Vec d = (u * my_cos(r1) * r2s + v * my_sin(r1) * r2s + w * my_sqrt(1 - r2)).norm();
      return obj.e + f.mult(radiance(Ray{x, d}, depth, rnd));
   } else if (obj.refl == SPEC) {          // Ideal SPECULAR reflection
      return obj.e + f.mult(radiance(Ray{ x, r.d - n * 2 * n.dot(r.d) }, depth, rnd));
   }

   // Ideal dielectric REFRACTION
   const Ray reflRay{x, r.d - n * 2 * n.dot(r.d)};
   const bool into  = n.dot(nl)>0;                // Ray from outside going in?
   const double nc  = 1;
   const double nt  = 1.5;
   const double nnt = into ? nc / nt : (nt / nc);
   const double ddn = r.d.dot(nl);
   double cos2t     = 1 - nnt * nnt * (1 - ddn * ddn);

   // Total internal reflection
   if (cos2t < 0) {
      return obj.e + f.mult(radiance(reflRay, depth, rnd));
   }

   const Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + my_sqrt(cos2t)))).norm();
   const double a  = nt - nc;
   const double b  = nt + nc;
   const double R0 = a * a / (b * b);
   const double c  = 1 - (into ? -ddn : tdir.dot(n));
   const double Re = R0 + (1 - R0) * c * c * c * c * c;
   const double Tr = 1 - Re;
   const double P  = 0.25 + 0.5 * Re;
   const double RP = Re / P;
   const double TP = Tr / (1 - P);

   return obj.e + f.mult(depth > 2 ? (rnd() < P ?   // Russian roulette
      radiance(reflRay,depth,rnd)*RP:radiance(Ray{x,tdir},depth,rnd)*TP) :
      radiance(reflRay,depth,rnd)*Re+radiance(Ray{x,tdir},depth,rnd)*Tr);
}

constexpr auto do_compute(auto& c, int y, auto& r, auto& rnd, const int samps) {
   constexpr Ray cam{ Vec{50,52,295.6}, Vec{0, -0.042612, -1}.norm() };
   constexpr Vec cx{ width * 0.5135 / height };
   constexpr Vec cy = (cx % cam.d).norm() * 0.5135;

   if (!std::is_constant_evaluated()) {
      std::fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100.0 * y / (height - 1));
   }

   for (unsigned short x = 0; x < width; x++) {

      // 2x2 subpixel rows
      for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) {
         // 2x2 subpixel cols
         for (int sx = 0; sx < 2; sx++, r = Vec{}) {
            for (int s = 0; s < samps; s++) {
               const double r1 = 2 * rnd(), dx = r1 < 1 ? my_sqrt(r1) - 1 : (1 - my_sqrt(2 - r1));
               const double r2 = 2 * rnd(), dy = r2 < 1 ? my_sqrt(r2) - 1 : (1 - my_sqrt(2 - r2));
               Vec d = cx*( ( (sx+.5 + dx)/2 + x)/width - 0.5) +
                       cy*( ( (sy+.5 + dy)/2 + y)/height - 0.5) + cam.d;
               r = r + radiance(Ray{ cam.o + d * 140, d.norm() }, 0, rnd) * (1.0 / samps);
            } // Camera rays are pushed ^^^^^ forward to start in interior
            c[i] = c[i] + Vec{clamp(r.x), clamp(r.y), clamp(r.z)} * 0.25;
         }
      }
   }
}

#if !NORT
void compute_omp(auto& c, const int samps) {
   Vec r{};
   auto rnd = []{ return drand48(); };

#pragma omp parallel for schedule(dynamic, 1) private(r)
   for (int y = 0; y < height; y++) {
      do_compute(c, y, r, rnd, samps);
   }
}
#endif

constexpr auto compute_cexpr() {
   std::array<Vec, width * height> c{};
   Vec r{};
   F32Random rnd{};
   for (int y = 0; y < height; y++) {
      do_compute(c, y, r, rnd, num_samps);
   }
   return c;
}

constexpr auto generate_img(const auto& c) {
   String str{};
   str += "P3\n";
   str += String::parse_int(width);
   str += ' ';
   str += String::parse_int(height);
   str += "\n255\n";
   for (std::size_t i = 0; i < width*height; ++i) {
      str += String::parse_int(toInt(c[i].x));
      str += ' ';
      str += String::parse_int(toInt(c[i].y));
      str += ' ';
      str += String::parse_int(toInt(c[i].z));
      str += ' ';
   }
   return str;
}

static int write_file(const char* filename, const auto& c) {
   FILE* file = std::fopen(filename, "w");
   if (!file) {
      std::fprintf(stderr, "smallpt: failed to open '%s'\n", filename);
      return 1;
   }
   std::fputs(c.c_str(), file);
   /*std::fprintf(file, "P3\n%d %d\n%d\n", width, height, 255);
   for (int i = 0; i < width*height; i++) {
     fprintf(file, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
   }*/
   std::fclose(file);
   return 0;
}

int main(int argc, char *argv[]){
   if (argc == 1) {
      static constinit auto ppm = generate_img(compute_cexpr());
      return write_file("image_cexpr.ppm", ppm);
   } else if (argc == 2) {
#if NORT
      std::fputs("smallpt: runtime execution not supporte.\n", stderr);
#else
      // static, because of stack overflow
      static std::array<Vec, width * height> c{};
      compute_omp(c, std::atoi(argv[1]));
      auto ppm = generate_img(c);
      return write_file("image.ppm", ppm);
#endif
   } else {
      std::fputs("Usage: smallpt [num_samps]\n", stderr);
      return 1;
   }
}
