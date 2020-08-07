#include "Light.h"
#include "VectorUtility.hpp"

/* Constructor. Implemented for you. */
PointLight::PointLight(const Vector3f & position, const Vector3f & intensity)
    : position(position), intensity(intensity) {}

// Compute the contribution of light at point p using the
// inverse square law formula
Vector3f PointLight::computeLightContribution(const Vector3f& p) {
    Vector3f distance = this->position - p;

    return this->intensity / vector_utility::dot(distance, distance);
}
