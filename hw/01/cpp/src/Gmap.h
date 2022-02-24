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
  int id;

  // constructor without arguments
  Vertex() : point(Point()) 
  {}

  // constructor with x,y,z arguments to immediately initialise the point member on this Vertex.
  Vertex(const double &x, const double &y, const double &z) : point(Point(x,y,z))
  {}

  // a dart incident to this Vertex:
  // ...

};

struct Edge {
  // a dart incident to this Edge:
  // ...

  // function to compute the barycenter for this Edge (needed for triangulation output):
  // Point barycenter() {}
  Point bcentre;

  // function to compute the barycenter for this Edge (needed for triangulation output):
  Point barycenter(Point a, Point b){
    bcentre = (a + b) / 2;
    return bcentre;
  }
};

struct Face {
    int id;
  // a dart incident to this Face:
  // ...

  // function to compute the barycenter for this Face (needed for triangulation output):
  // Point barycenter() {}

  int e, f2, f3, d;
  int f_id;

};

struct Volume {
  // a dart incident to this Volume:
  // ...


};