#ifndef _SCENE_H_
#define _SCENE_H_

#include <vector>

#include "Ray.h"
#include "Image.h"
#include "Light.h"
#include "defs.h"

// Forward declarations to avoid cyclic references
class Camera;
class PointLight;
class Material;
class Shape;

// Class to hold everything related to a scene.
class Scene {
public:
    int maxRecursionDepth;          // Maximum recursion depth
    float intTestEps;               // IntersectionTestEpsilon. You will need this one while implementing intersect routines in Shape class
    float shadowRayEps;             // ShadowRayEpsilon. You will need this one while generating shadow rays.
    Vector3f backgroundColor;       // Background color
    Vector3f ambientLight;          // Ambient light radiance

    std::vector<Camera *> cameras;       // Vector holding all cameras
    std::vector<PointLight *> lights;    // Vector holding all point lights
    std::vector<Material *> materials;   // Vector holding all materials
    std::vector<Vector3f> vertices;      // Vector holding all vertices (vertex data)
    std::vector<Shape *> objects;        // Vector holding all shapes

    Scene(const char *xmlPath);     // Constructor. Parses XML file and initializes vectors above. Implemented for you.

    void renderScene(void);         // Method to render scene, an image is created for each camera in the scene. You will implement this.

private:
    // Write any other stuff here

    unsigned char number_of_cpus;

    void render_image(const Camera * camera, Image & image, unsigned char cpu_index, int number_of_rows_per_cpu) const;
    ReturnVal send_ray(const Ray & ray) const;
    Vector3f evaluate_shading(const ReturnVal & hit_record, unsigned char bounces_remaining) const;
    bool point_in_shadow(const Ray & ray) const;
};

#endif
