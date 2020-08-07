#include "Shape.h"
#include "VectorUtility.hpp"

#include <cmath>
#include <limits>

Shape::Shape(void) {}

Shape::Shape(int id, int matIndex)
    : id(id), matIndex(matIndex) {}

Sphere::Sphere(void) {}

/* Constructor for sphere. You will implement this. */
Sphere::Sphere(int id, int matIndex, int cIndex, float R, std::vector<Vector3f> *pVertices)
    : Shape(id, matIndex), center((*pVertices)[cIndex - 1]), radius(R), R_square(R * R) {}

/* Sphere-ray intersection routine. You will implement this.
Note that ReturnVal structure should hold the information related to the intersection point, e.g., coordinate of that point, normal at that point etc.
You should to declare the variables in ReturnVal structure you think you will need. It is in defs.h file. */
ReturnVal Sphere::intersect(const Ray & ray) const {
    Vector3f e_minus_c = ray.origin - this->center;
    float A = vector_utility::dot(ray.direction, ray.direction);
    float B_over_2 = vector_utility::dot(ray.direction, e_minus_c);
    float C = vector_utility::dot(e_minus_c, e_minus_c) - this->R_square;
    float discriminant_over_4 = (B_over_2 * B_over_2 - A * C);     // B * B / 4 - * A * C;

    if (discriminant_over_4 < 0) {

        return ReturnVal{
            Vector3f{},
            Vector3f{},
            Vector3f{},
            -1,
            std::numeric_limits<float>::infinity()
        };
    }

    else {
        Vector3f p;     // intersection point
        float t;
        float discriminant_sqrt_over_2 = 0;

        if (discriminant_over_4 > 0) {
            discriminant_sqrt_over_2 = std::sqrt(discriminant_over_4);
        }   

        t = (-B_over_2 - discriminant_sqrt_over_2) / A;
    
        if (t < 0) {
            t = (-B_over_2 + discriminant_sqrt_over_2) / A;

            if (t < 0) {
                return ReturnVal{
                    Vector3f{},
                    Vector3f{},
                    Vector3f{},
                    -1,
                    std::numeric_limits<float>::infinity()
                };
            }
        }

        p = ray.getPoint(t);

        return ReturnVal{
            ray.direction,
            p,
            this->get_normal(p),
            this->matIndex,
            t
        };
    }
}

// return the unit normal vector on the given point
Vector3f Sphere::get_normal(const Vector3f &point) const {
    return (point - this->center) / this->radius;
}

Triangle::Triangle(void) {}

/* Constructor for triangle. You will implement this. */
Triangle::Triangle(int id, int matIndex, int p1Index, int p2Index, int p3Index, std::vector<Vector3f> *pVertices)
    : Shape(id, matIndex), point1((*pVertices)[p1Index - 1]), point2((*pVertices)[p2Index - 1]), point3((*pVertices)[p3Index - 1]) {}

/* Triangle-ray intersection routine. You will implement this.
Note that ReturnVal structure should hold the information related to the intersection point, e.g., coordinate of that point, normal at that point etc.
You should to declare the variables in ReturnVal structure you think you will need. It is in defs.h file. */
ReturnVal Triangle::intersect(const Ray & ray) const {
    float a = this->point1.x - this->point2.x;
    float b = this->point1.y - this->point2.y;
    float c = this->point1.z - this->point2.z;
    float d = this->point1.x - this->point3.x;
    float e = this->point1.y - this->point3.y;
    float f = this->point1.z - this->point3.z;
    float g = ray.direction.x;
    float h = ray.direction.y;
    float i = ray.direction.z;
    float j = this->point1.x - ray.origin.x;
    float k = this->point1.y - ray.origin.y;
    float l = this->point1.z - ray.origin.z;
    float ei_minus_hf = e * i - h * f;
    float gf_minus_di = g * f - d * i;
    float dh_minus_eg = d * h - e * g;
    float M = a * ei_minus_hf + b * gf_minus_di + c * dh_minus_eg;
    float beta = (j * ei_minus_hf + k * gf_minus_di + l * dh_minus_eg) / M;


    if (beta < 0 || beta > 1) {
        return ReturnVal{
            Vector3f{},
            Vector3f{},
            Vector3f{},
            -1,
            std::numeric_limits<float>::infinity()
        };
    }

    else {
        float ak_minus_jb = a * k - j * b;
        float jc_minus_al = j * c - a * l;
        float bl_minus_kc = b * l - k * c;
        float gamma = (i * ak_minus_jb + h * jc_minus_al + g * bl_minus_kc) / M;

        if (gamma < 0) {
            return ReturnVal{
                Vector3f{},
                Vector3f{},
                Vector3f{},
                -1,
                std::numeric_limits<float>::infinity()
            };
        }

        float gamma_beta_sum = gamma + beta; 

        if (gamma_beta_sum > 1) {
            return ReturnVal{
                Vector3f{},
                Vector3f{},
                Vector3f{},
                -1,
                std::numeric_limits<float>::infinity()
            };
        }

        float t = -1 * (f * ak_minus_jb + e * jc_minus_al + d * bl_minus_kc) / M;

        if (t < 0) {
            return ReturnVal{
                Vector3f{},
                Vector3f{},
                Vector3f{},
                -1,
                std::numeric_limits<float>::infinity()
            };
        }

        else {
            Vector3f p = ray.getPoint(t);   // intersection point
            
            Vector3f cross;

            if(gamma_beta_sum != 1){
                cross = vector_utility::cross(p - this->point2, p - this->point3);
            }
            else if(beta != 1){
                cross = vector_utility::cross(p - this->point1, p - this->point2);                
            }
            else{
                cross = vector_utility::cross(p - this->point3, p - this->point1);
            }
            
            return ReturnVal{
                ray.direction,
                p,
                vector_utility::normalize(cross),
                this->matIndex,
                t
            };
        }
    }
}

Mesh::Mesh() {}

/* Constructor for mesh. You will implement this. */
Mesh::Mesh(int id, int matIndex, const std::vector<Triangle>& faces, std::vector<int> *pIndices, std::vector<Vector3f> *pVertices)
    : Shape(id, matIndex), faces(faces) {}

/* Mesh-ray intersection routine. You will implement this.
Note that ReturnVal structure should hold the information related to the intersection point, e.g., coordinate of that point, normal at that point etc.
You should to declare the variables in ReturnVal structure you think you will need. It is in defs.h file. */
ReturnVal Mesh::intersect(const Ray & ray) const {
    ReturnVal first_hit_record = {
        Vector3f{},
        Vector3f{},
        Vector3f{},
        -1,
        std::numeric_limits<float>::infinity()
    };

    for (const auto & triangle : this->faces) {    
        ReturnVal hit_record = triangle.intersect(ray);

        if (hit_record.t < first_hit_record.t) {
            first_hit_record = hit_record;
        }
    }

    return first_hit_record;
}
