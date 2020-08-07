#ifndef __VEC4_H__
#define __VEC4_H__

#include <iostream>

class Vec3;

class Vec4 {
public:
    double x, y, z, t;
    int colorId;

    Vec4();
    Vec4(double x, double y, double z, double t, int colorId);
    Vec4(const Vec3 &other);
    Vec4(const Vec4 &other);
    Vec4 operator/(double scalar);
    Vec4& operator=(const Vec4 &other);
    Vec4& operator/=(double scalar);
    
    double getElementAt(int index) const;

    friend std::ostream& operator<<(std::ostream& os, const Vec4& v);
};

#endif
