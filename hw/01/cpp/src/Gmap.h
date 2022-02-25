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
  // functions to determine the darts 
  int invols_a0(Dart da, Dart db)
  {if (da.e_id==db.e_id && da.f_id==db.f_id && da.v_id!=db.v_id) {da.a0=db.id;}
      return da.a0;}

      int invols_a1(Dart da, Dart db)
      {if (da.v_id==db.v_id && da.f_id==db.f_id && da.e_id!=db.e_id) {da.a1=db.id;}
        return da.a1;}

    int invols_a2(Dart da, Dart db)
    {if (da.v_id==db.v_id && da.e_id==db.e_id && da.f_id!=db.f_id) {da.a2=db.id;}
        return da.a2;}

  // cells:
  // ...

};

struct Vertex {
  // the coordinates of this vertex:
  Point point;

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
};

struct Face {
  // a dart incident to this Face:
  // ...

  // function to compute the barycenter for this Face (needed for triangulation output):
  // Point barycenter() {}

};

struct Volume {
  // a dart incident to this Volume:
  // ...

};
