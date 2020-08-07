#include <iostream>
#include <cmath>        // std::abs(float)
#include "Helpers.h"
#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"

/*
 * Calculate cross product of vec3 a, vec3 b and return resulting vec3.
 */
Vec3 crossProductVec3(Vec3 a, Vec3 b) {
    Vec3 result;

    result.x = a.y * b.z - b.y * a.z;
    result.y = b.x * a.z - a.x * b.z;
    result.z = a.x * b.y - b.x * a.y;

    return result;
}

/*
 * Calculate dot product of vec3 a, vec3 b and return resulting value.
 */
double dotProductVec3(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Find length (|v|) of vec3 v.
 */
double magnitudeOfVec3(Vec3 v) {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

/*
 * Normalize the vec3 to make it unit vec3.
 */
Vec3 normalizeVec3(Vec3 v) {
    Vec3 result;
    double d;

    d = magnitudeOfVec3(v);
    result.x = v.x / d;
    result.y = v.y / d;
    result.z = v.z / d;

    return result;
}

/*
 * Return -v (inverse of vec3 v)
 */
Vec3 inverseVec3(Vec3 v) {
    Vec3 result;
    result.x = -v.x;
    result.y = -v.y;
    result.z = -v.z;

    return result;
}

/*
 * Add vec3 a to vec3 b and return resulting vec3 (a+b).
 */
Vec3 addVec3(Vec3 a, Vec3 b) {
    Vec3 result;
    result.x = a.x + b.x;
    result.y = a.y + b.y;
    result.z = a.z + b.z;

    return result;
}

/*
 * Subtract vec3 b from vec3 a and return resulting vec3 (a-b).
 */
Vec3 subtractVec3(Vec3 a, Vec3 b) {
    Vec3 result;
    result.x = a.x - b.x;
    result.y = a.y - b.y;
    result.z = a.z - b.z;

    return result;
}

/*
 * Multiply each element of vec3 with scalar.
 */
Vec3 multiplyVec3WithScalar(Vec3 v, double c) {
    Vec3 result;
    result.x = v.x * c;
    result.y = v.y * c;
    result.z = v.z * c;

    return result;
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void printVec3(Vec3 v) {
    std::cout << "(" << v.x << "," << v.y << "," << v.z << ")" << '\n';
}

/*
 * Check whether vec3 a and vec3 b are equal.
 * In case of equality, returns 1.
 * Otherwise, returns 0.
 */
int areEqualVec3(Vec3 a, Vec3 b) {
    /* if x difference, y difference and z difference is smaller than threshold, then they are equal */
    if ((std::abs(a.x - b.x) < EPSILON) && (std::abs(a.y - b.y) < EPSILON) && (std::abs(a.z - b.z) < EPSILON)) {
        return 1;
    }

    else {
        return 0;
    }
}

/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 getIdentityMatrix() {
    Matrix4 result;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) {
                result.val[i][j] = 1.0;
            }

            else {
                result.val[i][j] = 0.0;
            }
        }
    }

    return result;
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 multiplyMatrixWithMatrix(Matrix4 m1, Matrix4 m2) {
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            total = 0;

            for (int k = 0; k < 4; k++){
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 multiplyMatrixWithVec4(Matrix4 m, Vec4 v) {
    double values[4];
    double total;

    for (int i = 0; i < 4; i++) {
        total = 0;

        for (int j = 0; j < 4; j++){
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

// ----------------------------------------------------------------------------

//
//	Vector addition
//
Vec3 operator+(const Vec3 &vec1, const Vec3 &vec2) {
    return {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z, -1};
}

Vec3& operator+=(Vec3 &vector1, const Vec3 &vec2) {
    vector1.x += vec2.x;
    vector1.y += vec2.y;
    vector1.z += vec2.z;

    return vector1;
}

//
// Vector subtraction
//
Vec3 operator-(const Vec3 &vec1, const Vec3 &vec2) {
    return {vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z, -1};
}

//
// Vector element-wise multiplication
//
Vec3 operator*(const Vec3 &vec1, const Vec3 &vec2) {
    return {
        vec1.x * vec2.x,
        vec1.y * vec2.y,
        vec1.z * vec2.z,
        -1,
    };
}

//
//	Element-wise division
//
Vec3 operator/(const Vec3 &vec1, const Vec3 &vec2) {
    return {
        vec1.x / vec2.x,
        vec1.y / vec2.y,
        vec1.z / vec2.z,
        -1,
    };
}

// vector scalar division
Vec3 operator/(const Vec3 &vector, double scalar) {
    return {
        vector.x / scalar,
        vector.y / scalar,
        vector.z / scalar,
        -1,
    };
}

Vec3& operator/=(Vec3 &vector, double scalar) {
    vector.x /= scalar;
    vector.y /= scalar;
    vector.z /= scalar;

    return vector;
}

//
// Equality check
//
bool operator==(const Vec3 &vec1, const Vec3 &vec2) {
	return (std::abs(vec1.x - vec2.x) < EPSILON && std::abs(vec1.y - vec2.y) < EPSILON && std::abs(vec1.z - vec2.z) < EPSILON);
}

bool operator!=(const Vec3 &vec1, const Vec3 &vec2) {
	return !(vec1 == vec2);
}


//
// Vector scalar multiplication
//
Vec3 operator*(double scalar, const Vec3 &vec) {
    return {scalar * vec.x, scalar * vec.y, scalar * vec.z, -1};
}

Vec3 operator*(const Vec3 &vec, double scalar) {
    return {scalar * vec.x, scalar * vec.y, scalar * vec.z, -1};
}

namespace vector_utility {
    //
    // Vector dot product
    //
    double dot(const Vec3 &vec1, const Vec3 &vec2) {
        return (vec1.x * vec2.x) + (vec1.y * vec2.y) + (vec1.z * vec2.z);
    }

    //
    // Vector cross product
    //
    Vec3 cross(const Vec3 &vec1, const Vec3 &vec2) {
        return {
            (vec1.y * vec2.z) - (vec1.z * vec2.y),
            (vec1.z * vec2.x) - (vec1.x * vec2.z),
            (vec1.x * vec2.y) - (vec1.y * vec2.x),
            -1,
        };
    }

    //
    // Magnitude of the vector
    //
    double magnitude(const Vec3 &vec) {
        return std::sqrt(vector_utility::dot(vec, vec));
    }

    //
    // Normalize vector
    //
    void normalize(Vec3 &vec){
        double mag = vector_utility::magnitude(vec);

        vec /= mag;
    }
}

/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 operator*(const Matrix4 & m1, const Matrix4 & m2) {
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            total = 0;
            for (int k = 0; k < 4; k++) {
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
        }
    }

    return result;
}

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 operator*(const Matrix4 & m, const Vec4 & v) {
    double values[4];
    double total;

    for (int i = 0; i < 4; i++) {
        total = 0;
        for (int j = 0; j < 4; j++) {
            total += m.val[i][j] * v.getElementAt(j);
        }
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}

Vec3 operator*(const Matrix4 & m, const Vec3 & v) {
    Vec4 v4(v);

    v4 = m * v4;
    v4 /= v4.t;

    return Vec3(v4);
}

namespace matrix_utility {
    Matrix4 get_identity_matrix() {
        Matrix4 result;

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    result.val[i][j] = 1.0;
                }

                else {
                    result.val[i][j] = 0.0;
                }
            }
        }

        return result;
    }

    Matrix4 transpose(const Matrix4 & matrix) {
        Matrix4 transposed_matrix;

        for (unsigned int i = 0; i < 4; ++i) {
            for (unsigned int j = 0; j < 4; ++j) {
                transposed_matrix.val[j][i] = matrix.val[i][j];
            }
        }

        return transposed_matrix;
    }
}
