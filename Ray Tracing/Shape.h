#ifndef _SHAPE_H_
#define _SHAPE_H_

#include <vector>

#include "Ray.h"
#include "defs.h"

// Base class for any shape object
class Shape {
public:
    int id;             // Id of the shape
    int matIndex;       // Material index of the shape

    virtual ReturnVal intersect(const Ray & ray) const = 0; // Pure virtual method for intersection test. You must implement this for sphere, triangle, and mesh.

    Shape(void);
    Shape(int id, int matIndex); // Constructor

private:
    // Write any other stuff here
};

// Class for sphere
class Sphere: public Shape {
public:
    Sphere(void);       // Constructor
    Sphere(int id, int matIndex, int cIndex, float R, std::vector<Vector3f> *vertices);      // Constructor
    ReturnVal intersect(const Ray & ray) const;     // Will take a ray and return a structure related to the intersection information. You will implement this.

private:
    // Write any other stuff here

    Vector3f center;
    float radius;
    float R_square;
    Vector3f get_normal(const Vector3f &point) const;
};

// Class for triangle
class Triangle: public Shape {
public:
    Triangle(void);     // Constructor
    Triangle(int id, int matIndex, int p1Index, int p2Index, int p3Index, std::vector<Vector3f> *vertices);  // Constructor
    ReturnVal intersect(const Ray & ray) const;     // Will take a ray and return a structure related to the intersection information. You will implement this.

private:
    // Write any other stuff here

    Vector3f point1;
    Vector3f point2;
    Vector3f point3;
};

// Class for mesh
class Mesh: public Shape {
public:
    Mesh(void);     // Constructor
    Mesh(int id, int matIndex, const std::vector<Triangle>& faces, std::vector<int> *pIndices, std::vector<Vector3f> *vertices);       // Constructor
    ReturnVal intersect(const Ray & ray) const;     // Will take a ray and return a structure related to the intersection information. You will implement this.

private:
    // Write any other stuff here

    std::vector<Triangle> faces;
};

#endif
