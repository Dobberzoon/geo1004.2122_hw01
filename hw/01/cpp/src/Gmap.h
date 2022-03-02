#pragma once

#include "Point.h"

struct Point;
struct Dart;
struct Vertex;
struct Edge;
struct Face;
struct Volume;

struct Dart {
  // involutions:
  Dart* a0 = nullptr;
  Dart* a1 = nullptr;
  Dart* a2 = nullptr;
  Dart* a3 = nullptr;

  // cells:
  Vertex* v = nullptr;
  Edge* e = nullptr;
  Face* f = nullptr;
  Volume* vo = nullptr;

  // constructor without arguments
  Dart(){}

  void invol_a0(Dart *dart_a, Dart *dart_b) {
      if ((dart_a->e==dart_b->e) && (dart_a->f==dart_b->f) && (dart_a->v!=dart_b->v)) {
          //std::cout << "involution a0 found\n";
          dart_a->a0 = dart_b;
          //return a0;
      }
  }

  void invol_a1(Dart *dart_a, Dart *dart_b) {
      if (dart_a->v==dart_b->v && dart_a->f==dart_b->f && dart_a->e!=dart_b->e) {
          //std::cout << "involution a1 found\n";
          dart_a->a1 = dart_b;
      }
  }

  void invol_a2(Dart *dart_a, Dart *dart_b) {
      if (dart_a->v==dart_b->v && dart_a->e==dart_b->e && dart_a->f!=dart_b->f) {
          //std::cout << "involution a2 found\n";
          dart_a->a2 = dart_b;
      }
  }

};

struct Vertex {
  // the coordinates of this vertex:
  Point point;

  // also in string format for making keys
  std::string xyz;

  // constructor without arguments
  Vertex() : point(Point()) {}

  // constructor with x,y,z arguments to immediately initialise the point member on this Vertex.
  Vertex(const double &x, const double &y, const double &z) : point(Point(x,y,z)) {}

  // a dart incident to this Vertex:
  Dart* dart = nullptr;

  // function to convert point.x/y/z into concatenated string
  std::string xyz_tostring(const double &x, const double &y, const double &z) {
      std::string xS, yS, zS;
      xS = std::to_string(x);
      yS = std::to_string(y);
      zS = std::to_string(z);
      xyz = xS + yS + zS;
      return xyz;
  }

};

struct Edge {
  // a dart incident to this Edge:
  Dart* dart = nullptr;

  // begin and end vertex of edge
  int origin_v, end_v;

  // also in string format for keys
  std::string edgeS;

  // barycenter of edge
  Point barCE;

  // constructor without arguments
  Edge(){}

  // constructor
  Edge(const int &origin_v, const int &end_v) {
      this->origin_v = origin_v;
      this->end_v = end_v;
  }

  // function to convert point.x/y/z into concatenated string
  std::string edge_tostring(const int &origin_v, const int &end_v) {

        std::string origin_vS, end_vS;
        origin_vS = std::to_string(origin_v);
        end_vS = std::to_string(end_v);
        edgeS = origin_vS + end_vS;
        return edgeS;
    }

  // function to compute the barycenter for this Edge (needed for triangulation output):
  // Point barycenter() {}
  // function to compute the barycenter for this Edge (needed for triangulation output):
  Point barycenter(std::vector<Vertex>&vertices) {
      Vertex start, end;
      start = vertices[origin_v];
      end = vertices[end_v];
      barCE = (start.point+end.point)/2;
      return barCE;
  }
};

struct Face {
  // a dart incident to this Face:
  Dart* dart = nullptr;

  // face vertices
  std::vector<int> face_vertices;

  // also in string format for keys
  std::string faceS;

  // barycenter of face
  Point barCF;

  // function to convert point.x/y/z into concatenated string
  void face_tostring() {
      std::vector<std::string> vfaceS;
      for (auto i : face_vertices) {
          std::string str;
          str = std::to_string(i);
          faceS += str;
      }
  }
  // function to compute the barycenter for this Face (needed for triangulation output):
  // Point barycenter() {}
  // function to compute the barycenter for this Face (needed for triangulation output):
  // Point barycenter() {}

  Point barycenter(std::vector<Vertex>&vertices) {
      Vertex a, b, c, d;
      a = vertices[face_vertices[0]];
      b = vertices[face_vertices[1]];
      c = vertices[face_vertices[2]];
      d = vertices[face_vertices[3]];
      barCF=(a.point+b.point+c.point+d.point)/4;
      return barCF;}

};

struct Volume {
  // a dart incident to this Volume:
  Dart* dart = nullptr;

};

