#include "VectorUtility.hpp"

#include <cmath>

//
//	Vector addition
//
Vector3f operator+(const Vector3f &vec1, const Vector3f &vec2) {
    return {{vec1.x + vec2.x}, {vec1.y + vec2.y}, {vec1.z + vec2.z}};
}

Vector3f operator+=(Vector3f &vector1, const Vector3f &vec2) {
    vector1.r += vec2.r;
    vector1.g += vec2.g;
    vector1.b += vec2.b;

    return vector1;
}

//
// Vector subtraction
//
Vector3f operator-(const Vector3f &vec1, const Vector3f &vec2) {
    return {{vec1.x - vec2.x}, {vec1.y - vec2.y}, {vec1.z - vec2.z}};
}

//
// Vector element-wise multiplication
//
Vector3f operator*(const Vector3f &vec1, const Vector3f &vec2) {
    return {
        {vec1.x * vec2.x},
        {vec1.y * vec2.y},
        {vec1.z * vec2.z}
    };
}

//
//	Element-wise division
//
Vector3f operator/(const Vector3f &vec1, const Vector3f &vec2) {
    return {
        {vec1.x / vec2.x},
        {vec1.y / vec2.y},
        {vec1.z / vec2.z}
    };
}

// vector scalar division
Vector3f operator/(const Vector3f &vector, float scalar) {
    return {
        {vector.x / scalar},
        {vector.y / scalar},
        {vector.z / scalar}
    };
}

//
// Equality check
//
bool operator==(const Vector3f &vec1, const Vector3f &vec2) {
	return (vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z);
}

bool operator!=(const Vector3f &vec1, const Vector3f &vec2) {
	return !(vec1.x == vec2.x && vec1.y == vec2.y && vec1.z == vec2.z);
}


//
// Print vector:
//
std::ostream& operator<<(std::ostream& os, const Vector3f &vec) {
	os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
	return os;
}

//
// Vector scalar multiplication
//
Vector3f operator*(const float &scalar, const Vector3f &vec) {
    return {{scalar * vec.x}, {scalar * vec.y}, {scalar * vec.z}};
}

Vector3f operator*(const Vector3f &vec, const float &scalar) {
    return {{scalar * vec.x}, {scalar * vec.y}, {scalar * vec.z}};
}

namespace vector_utility {
    //
    // Vector dot product
    //
    float dot(const Vector3f &vec1, const Vector3f &vec2) {
        return (vec1.x * vec2.x) + (vec1.y * vec2.y) + (vec1.z * vec2.z);
    }

    //
    // Vector cross product
    //
    Vector3f cross(const Vector3f &vec1, const Vector3f &vec2) {
        return {
            {(vec1.y * vec2.z) - (vec1.z * vec2.y)},
            {(vec1.z * vec2.x) - (vec1.x * vec2.z)},
            {(vec1.x * vec2.y) - (vec1.y * vec2.x)}
        };
    }

    //
    // Magnitude of the vector
    //
    float magnitude(const Vector3f &vec) {
        return std::sqrt(vector_utility::dot(vec, vec));
    }

    //
    // Normalize vector
    //
    Vector3f normalize(const Vector3f &vec){
        float mag = vector_utility::magnitude(vec);

        return vec / mag;
    }
}
