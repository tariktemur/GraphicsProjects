#include "Vec3.h"
#include <iomanip>

#include "Vec4.h"

Vec3::Vec3() :
    x(0.0), y(0.0), z(0.0), colorId(-1) {}

Vec3::Vec3(double x, double y, double z) :
    x(x), y(y), z(z), colorId(-1) {}

Vec3::Vec3(double x, double y, double z, int colorId) :
    x(x), y(y), z(z), colorId(colorId) {}

Vec3::Vec3(const Vec3 &other) :
    x(other.x), y(other.y), z(other.z), colorId(other.colorId) {}

Vec3::Vec3(const Vec4 &other) :
    x(other.x), y(other.y), z(other.z), colorId(other.colorId) {}

Vec3& Vec3::operator=(const Vec3 &other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->colorId = other.colorId;

    return *this;
}

double Vec3::getElementAt(int index) {
    switch (index) {
    case 0:
        return this->x;

    case 1:
        return this->y;

    case 2:
        return this->z;

    default:
        return this->z;
    }
}

std::ostream& operator<<(std::ostream& os, const Vec3& v) {
    os << std::fixed << std::setprecision(6) << "[" << v.x << ", " << v.y << ", " << v.z << "]";

    return os;
}
