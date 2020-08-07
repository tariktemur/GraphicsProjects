#ifndef __VEC3_H__
#define __VEC3_H__

#include <iostream>

class Vec4;

class Vec3 {
public:
    double x, y, z;
    int colorId;

    Vec3();
    Vec3(double x, double y, double z);
    Vec3(double x, double y, double z, int colorId);
    Vec3(const Vec3 &other);
    Vec3(const Vec4 &other);
    Vec3& operator=(const Vec3 &other);

    double getElementAt(int index);
    
    friend std::ostream& operator<<(std::ostream& os, const Vec3& v);
};

#endif
