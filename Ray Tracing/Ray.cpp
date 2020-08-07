#include "Ray.h"
#include "VectorUtility.hpp"

#include <limits>

Ray::Ray() {}

Ray::Ray(const Vector3f& origin, const Vector3f& direction)
    : origin(origin), direction(direction) {}

/* Takes a parameter t and returns the point accoring to t. t is the parametric variable in the ray equation o+t*d.*/
Vector3f Ray::getPoint(float t) const {
    return origin + t * this->direction;
}

/*
Here are possible causes for concern about the following function: 
We're checking for the equality of floating point numbers, a small error interval may be introduced for robustness.  
Checking whether if the point is on the ray is much more costly than assuming it is.
*/

float Ray::gett(const Vector3f & p) const {
    Vector3f temp = (p - origin);
    float t = std::numeric_limits<float>::infinity();

    if (this->direction.x == 0) {
        if (temp.x != 0) {
            return std::numeric_limits<float>::infinity();
        }
    }

    else {
        t = temp.x / this->direction.x;
    }

    if (this->direction.y == 0){
        if (temp.y != 0) {
            return std::numeric_limits<float>::infinity();
        }
    }

    else {
        if (t == -1) {
            t = temp.y / this->direction.y;
        }

        else if (t * this->direction.y != temp.y) {
            return std::numeric_limits<float>::infinity();
        }
    }

    if (this->direction.z == 0){
        if (temp.z != 0) {
            return std::numeric_limits<float>::infinity();
        }
    }

    else {
        if (t == -1) {
            t = temp.z / this->direction.z;
        }
   
        else if (t * this->direction.z != temp.z) {
            return std::numeric_limits<float>::infinity();
        }
    }

    // If indeed the point is on the ray, all 3 equations of the following form must be satisfied: (p - o) / d = t 
    return t;
}
