#include "Color.h"
#include <iostream>
#include <iomanip>

Color::Color() {}

Color::Color(double r, double g, double b) {
    this->r = r;
    this->g = g;
    this->b = b;
}

Color::Color(const Color &other) {
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
}

Color Color::operator*(double scalar) {
    return Color{
        this->r * scalar,
        this->g * scalar,
        this->b * scalar,
    };
}

Color Color::operator/(double scalar) {
    return Color{
        this->r / scalar,
        this->g / scalar,
        this->b / scalar,
    };
}

Color Color::operator+(const Color & other) {
    return Color{
        this->r + other.r,
        this->g + other.g,
        this->b + other.b,
    };
}

Color Color::operator-(const Color & other) {
    return Color{
        this->r - other.r,
        this->g - other.g,
        this->b - other.b,
    };
}

Color & Color::operator=(const Color & other) {
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;

    return *this;
}

Color & Color::operator+=(const Color & other) {
    *this = *this + other;

    return *this;
}

std::ostream& operator<<(std::ostream& os, const Color& c) {
    os << std::fixed << std::setprecision(0) << "rgb(" << c.r << ", " << c.g << ", " << c.b << ")";
    return os;
}
