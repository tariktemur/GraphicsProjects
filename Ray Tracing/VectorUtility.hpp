#ifndef VECTORUTILITY_H
#define VECTORUTILITY_H

#include "defs.h"

#include <ostream>

//
//	Vector addition
//
Vector3f operator+(const Vector3f &vec1, const Vector3f &vec2);

Vector3f operator+=(Vector3f &vector1, const Vector3f &vector2);

//
// Vector subtraction
//
Vector3f operator-(const Vector3f &vec1, const Vector3f &vec2);

//
// Vector element-wise multiplication
//
Vector3f operator*(const Vector3f &vec1, const Vector3f &vec2);

//
//	Element-wise division
//
Vector3f operator/(const Vector3f &vec1, const Vector3f &vec2);

// vector scalar division
Vector3f operator/(const Vector3f &vector, float scalar);

//
// Equality check
//
bool operator==(const Vector3f &vec1, const Vector3f &vec2);

bool operator!=(const Vector3f &vec1, const Vector3f &vec2);

//
// Print vector:
//
std::ostream& operator<<(std::ostream& os, const Vector3f &vec);

//
// Vector scalar multiplication
//
Vector3f operator*(const float &scalar, const Vector3f &vec);

Vector3f operator*(const Vector3f &vec, const float &scalar);

namespace vector_utility {
    //
    // Vector dot product
    //
    float dot(const Vector3f &vec1, const Vector3f &vec2);

    //
    // Vector cross product
    //
    Vector3f cross(const Vector3f &vec1, const Vector3f &vec2);

    //
    // Magnitude of the vector
    //
    float magnitude(const Vector3f &vec);

    //
    // Normalize vector
    //
    Vector3f normalize(const Vector3f &vec);
}

#endif
