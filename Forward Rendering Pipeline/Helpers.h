#ifndef __HELPERS_H__
#define __HELPERS_H__

#define EPSILON 0.000000001

#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b);

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b);

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v);

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v);

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v);

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b);

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b);

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c);

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v);

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b);
/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 getIdentityMatrix();

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v);

// ----------------------------------------------------------------------------

//
//	Vector addition
//
Vec3 operator+(const Vec3 &vec1, const Vec3 &vec2);

Vec3& operator+=(Vec3 &vector1, const Vec3 &vec2);

//
// Vector subtraction
//
Vec3 operator-(const Vec3 &vec1, const Vec3 &vec2);

//
// Vector element-wise multiplication
//
Vec3 operator*(const Vec3 &vec1, const Vec3 &vec2);

//
//	Element-wise division
//
Vec3 operator/(const Vec3 &vec1, const Vec3 &vec2);

// vector scalar division
Vec3 operator/(const Vec3 &vector, double scalar);

//
// Equality check
//
bool operator==(const Vec3 &vec1, const Vec3 &vec2);

bool operator!=(const Vec3 &vec1, const Vec3 &vec2);

//
// Print vector:
//
std::ostream& operator<<(std::ostream& os, const Vec3 &vec);

//
// Vector scalar multiplication
//
Vec3 operator*(double scalar, const Vec3 &vec);

Vec3 operator*(const Vec3 &vec, double scalar);

namespace vector_utility {
    //
    // Vector dot product
    //
    double dot(const Vec3 &vec1, const Vec3 &vec2);

    //
    // Vector cross product
    //
    Vec3 cross(const Vec3 &vec1, const Vec3 &vec2);

    //
    // Magnitude of the vector
    //
    double magnitude(const Vec3 &vec);

    //
    // Normalize vector
    //
    void normalize(Vec3 &vec);
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 operator*(const Matrix4 & m1, const Matrix4 & m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 operator*(const Matrix4 & m, const Vec4 & v);

Vec3 operator*(const Matrix4 & m, const Vec3 & v);

namespace matrix_utility {
    Matrix4 get_identity_matrix();
    Matrix4 transpose(const Matrix4 & matrix);
}
#endif
