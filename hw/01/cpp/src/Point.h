#pragma once

/*
Structure to represent a 3D point. Various operations are defined. eg creation 
using `Point p(x,y,z)`, accessing a coordinate as `p[0]` (equivalent to `p.x`), vector addition/substraction 
as `p+q` or `p-q`, multiplication/division with a scalar as `p*s` or `p/s`, dot/cross 
product as `p.dot(q)` or `p.cross(q)`.
*/
struct Point {
  float x, y, z;
  
  // constructor without arguments; initialises all coordinates on 0.
  Point() : x(0), y(0), z(0)
  {}
  
  // constructor with x,y,z arguments;
  Point(const float &x, const float &y, const float &z) : x(x), y(y), z(z)
  {}
  
  // non-const coordinate accessor with `[]`
  float &operator[](const int &coordinate) {
    if (coordinate == 0) return x;
    else if (coordinate == 1) return y;
    else if (coordinate == 2) return z;
  }
  
  // const coordinate accessor with `[]`
  float operator[](const int &coordinate) const {
    if (coordinate == 0) return x;
    else if (coordinate == 1) return y;
    else if (coordinate == 2) return z;
  }

  // assignment operator. Assigns values from other point to this point
  void operator=(const Point &other) {
    x=other.x;
    y=other.y;
    z=other.z;
  }
  
  // addition operator (yields a new point)
  const Point operator+(const Point &other) const {
    return Point(x+other.x, y+other.y, z+other.z);
  }
  
  // assignment/addition operator. Adds other point to this point.
  void operator+=(const Point &other) {
    x+=other.x;
    y+=other.y;
    z+=other.z;
  }
  
  // subtraction operator (yields a new point)
  const Point operator-(const Point &other) const {
    return Point(x-other.x, y-other.y, z-other.z);
  }
  
  // element-wise multiplication operator (yields a new point)
  const Point operator*(const float &other) const {
    return Point(x*other, y*other, z*other);
  }
  
  // division operator (yields a new point)
  const Point operator/(const float &other) const {
    return Point(x/other, y/other, z/other);
  }
  
  // dot product operator
  float dot(const Point &other) const {
    return x*other.x + y*other.y + z*other.z;
  }
  
  // cross product operator
  const Point cross(const Point &other) const {
    return Point(y*other.z-z*other.y, -(x*other.z-z*other.x), x*other.y-y*other.x);
  }
};

// function to print the coordinates of a point using eg `std::cout << p`, where p is a Point struct
std::ostream& operator<<(std::ostream& os, const Point p) {
  os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
  return os;
}