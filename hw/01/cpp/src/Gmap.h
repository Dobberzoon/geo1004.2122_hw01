#pragma once

#include "Point.h"

struct Point;
struct Dart;
struct Vertex;
struct Edge;
struct Face;
struct Volume;

/*
Below you find the basic elements that you need to build the generalised map.
The main thing you need to fill out are the links between the elements:
  * the involutions and cells on the Dart
  * the darts on the cells

One way to do this is by using pointers. eg. define a member on the dart struct like

  Struct Dart {
    // involutions:
    Dart* a0 = nullptr;
    // ...

    // cells:
    // ...
  
  };

Then you could create and link Darts like:
  
  Dart* dart_a = new Dart();
  Dart* dart_b = new Dart();

  dart_a->a0 = dart_b;
*/

struct Dart {
  // involutions:
  Dart* a0 = nullptr;
  Dart* a1 = nullptr;
  Dart* a2 = nullptr;
  Dart* a3 = nullptr;

  // cells:
  Dart* v = nullptr;
  Dart* e = nullptr;
  Dart* f = nullptr;
  Dart* vo = nullptr;

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

  /*
  bool operator==(const Vertex& v) const {
      return point.x == v.point.x && point.y == v.point.y && point.z == v.point.z;
  }
   */

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
  int origin_v, end_v, sum;

  // also in string format for keys
  std::string edgeS, sumS;

  // constructor without arguments
  Edge(){}

  // constructor
  Edge(const int &origin_v, const int &end_v) {
      this->origin_v = origin_v;
      this->end_v = end_v;
  }

  // function to convert point.x/y/z into concatenated string
  std::string edge_tostring(const int &origin_v, const int &end_v) {

        //sum = origin_v + end_v + 1;

        std::string origin_vS, end_vS;

        //sumS = std::to_string(sum);
        origin_vS = std::to_string(origin_v);
        end_vS = std::to_string(end_v);
        edgeS = origin_vS + end_vS;
        return edgeS;
    }

  // function to compute the barycenter for this Edge (needed for triangulation output):
  // Point barycenter() {}
};

struct Face {
  // a dart incident to this Face:
  Dart* dart = nullptr;

  // face vertices
  int v0, v1, v2, v3;

  // constructor
  // input vertices should be given in CCW
  Face(const int &v0, const int &v1, const int &v2, const int &v3) {
      this->v0 = v0;
      this->v1 = v1;
      this->v2 = v2;
      this->v3 = v3;
  }

  // function to compute the barycenter for this Face (needed for triangulation output):
  // Point barycenter() {}

};

struct Volume {
  // a dart incident to this Volume:
  Dart* dart = nullptr;

};

/*
class HashVertexFunction {
public:
    size_t operator()(const Vertex& v) const {
        return (std::hash<float>()(v.point.x)) ^ (std::hash<float>()(v.point.y)) ^
               (std::hash<float>()(v.point.z));
    }
};
 */