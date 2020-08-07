#include "Scene.h"
#include "Camera.h"
#include "Material.h"
#include "Shape.h"
#include "tinyxml2.h"
#include "VectorUtility.hpp"

#include <algorithm>
#include <cmath>
#include <thread>   // for std::thread
#include <functional>   // for std::ref

using namespace tinyxml2;

/* 
 * Must render the scene from each camera's viewpoint and create an image.
 * You can use the methods of the Image class to save the image as a PPM file. 
 */
void Scene::renderScene(void) {
    for (const Camera * camera : this->cameras) {
        Image image(camera->imgPlane.nx, camera->imgPlane.ny);
        int number_of_rows_per_cpu = image.height / this->number_of_cpus;
        std::vector<std::thread> threads;
        threads.reserve(this->number_of_cpus);

        for (unsigned char i = 0; i < this->number_of_cpus; ++i) {
            threads.push_back(std::thread(&Scene::render_image, this, camera, std::ref(image), i, number_of_rows_per_cpu));
        }

        for (auto & thread : threads) {
            thread.join();
        }

        image.saveImage(camera->imageName);
    }
}

void Scene::render_image(const Camera * camera, Image & image, unsigned char cpu_index, int number_of_rows_per_cpu) const {
    int row = cpu_index * number_of_rows_per_cpu;
    int row_index_end;

    if (cpu_index != this->number_of_cpus) {
        row_index_end = row + number_of_rows_per_cpu;
    }

    else {
        row_index_end = image.height;
    }

    for ( ; row < row_index_end; ++row) {
        for (int col = 0; col < image.width; ++col) {
            ReturnVal hit_record = this->send_ray(camera->getPrimaryRay(col, row));

            Vector3f color_intensity = evaluate_shading(hit_record, this->maxRecursionDepth); 
            image.setPixelValue(col, row, Color{{
                static_cast<unsigned char>(std::min(255.0f, color_intensity.r)),
                static_cast<unsigned char>(std::min(255.0f, color_intensity.g)),
                static_cast<unsigned char>(std::min(255.0f, color_intensity.b))
            }});
        }
    }
}

// returns the first hit record if there is, otherwise return a dummy one with 
// background color as the color intensity and -1 as the t value
ReturnVal Scene::send_ray(const Ray & ray) const {
    ReturnVal first_hit_record = {
        Vector3f{},
        Vector3f{},
        Vector3f{},
        -1,
        std::numeric_limits<float>::infinity()
    };

    for (const Shape * shape : this->objects) {
        ReturnVal hit_record = shape->intersect(ray);
        if (hit_record.t < first_hit_record.t) {
            first_hit_record = hit_record;
        }
    }

    return first_hit_record;
}

Vector3f Scene::evaluate_shading(const ReturnVal & hit_record, unsigned char bounces_remaining) const {
    if (hit_record.MaterialIndex == -1) {
        return this->backgroundColor;
    }

    Material * material = this->materials[hit_record.MaterialIndex - 1];
    Vector3f color_intensity = material->ambientRef * this->ambientLight;

    for (PointLight * light : this->lights) {
        Vector3f light_distance = light->position - hit_record.intersection_point;

        if (point_in_shadow(Ray{hit_record.intersection_point + shadowRayEps * light_distance, light_distance})) {
            continue;
        }

        Vector3f surface_color = material->diffuseRef;
        Vector3f specular_color = material->specularRef;
        Vector3f light_intensity = light->computeLightContribution(hit_record.intersection_point);
        Vector3f light_direction = vector_utility::normalize(light_distance);
        Vector3f viewer_direction = -1 * vector_utility::normalize(hit_record.ray_direction);
        Vector3f half_vector = vector_utility::normalize(light_direction + viewer_direction);


        Vector3f diffuse_reflection = surface_color * light_intensity * std::max(0.0f, vector_utility::dot(hit_record.surface_normal, light_direction));
        color_intensity += diffuse_reflection;

        Vector3f specular_reflection = specular_color * light_intensity * std::pow(std::max(0.0f, vector_utility::dot(hit_record.surface_normal, half_vector)), material->phongExp);
        color_intensity += specular_reflection;


    }

    if (material->mirrorRef != Vector3f{{0}, {0}, {0}} && bounces_remaining > 0) {
        Vector3f d = vector_utility::normalize(hit_record.ray_direction);
        Ray ray;

        ray.direction =  d - 2 * vector_utility::dot(d, hit_record.surface_normal) * hit_record.surface_normal; 
        ray.origin = hit_record.intersection_point + this->shadowRayEps * ray.direction;

        ReturnVal hit_record = this->send_ray(ray);
        color_intensity += material->mirrorRef * evaluate_shading(hit_record, bounces_remaining - 1);
    }

    return color_intensity;
}

bool Scene::point_in_shadow(const Ray & ray) const {
    for (const Shape * shape : this->objects) {
        ReturnVal hit_record = shape->intersect(ray);

        if (hit_record.t >= 0 && hit_record.t < 1) {
            return true;
        }
    }

    return false;
}

// Parses XML file. 
Scene::Scene(const char *xmlPath) : number_of_cpus(std::thread::hardware_concurrency()) {
        const char *str;
        XMLDocument xmlDoc;
        XMLError eResult;
        XMLElement *pElement;

        maxRecursionDepth = 1;
        shadowRayEps = 0.001;

        eResult = xmlDoc.LoadFile(xmlPath);

        XMLNode *pRoot = xmlDoc.FirstChild();

        pElement = pRoot->FirstChildElement("MaxRecursionDepth");
        if(pElement != nullptr)
                pElement->QueryIntText(&maxRecursionDepth);

        pElement = pRoot->FirstChildElement("BackgroundColor");
        str = pElement->GetText();
        sscanf(str, "%f %f %f", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

        pElement = pRoot->FirstChildElement("ShadowRayEpsilon");
        if(pElement != nullptr)
                pElement->QueryFloatText(&shadowRayEps);

        pElement = pRoot->FirstChildElement("IntersectionTestEpsilon");
        if(pElement != nullptr)
                eResult = pElement->QueryFloatText(&intTestEps);

        // Parse cameras
        pElement = pRoot->FirstChildElement("Cameras");
        XMLElement *pCamera = pElement->FirstChildElement("Camera");
        XMLElement *camElement;
        while(pCamera != nullptr)
        {
        int id;
        char imageName[64];
        Vector3f pos, gaze, up;
        ImagePlane imgPlane;

                eResult = pCamera->QueryIntAttribute("id", &id);
                camElement = pCamera->FirstChildElement("Position");
                str = camElement->GetText();
                sscanf(str, "%f %f %f", &pos.x, &pos.y, &pos.z);
                camElement = pCamera->FirstChildElement("Gaze");
                str = camElement->GetText();
                sscanf(str, "%f %f %f", &gaze.x, &gaze.y, &gaze.z);
                camElement = pCamera->FirstChildElement("Up");
                str = camElement->GetText();
                sscanf(str, "%f %f %f", &up.x, &up.y, &up.z);
                camElement = pCamera->FirstChildElement("NearPlane");
                str = camElement->GetText();
                sscanf(str, "%f %f %f %f", &imgPlane.left, &imgPlane.right, &imgPlane.bottom, &imgPlane.top);
                camElement = pCamera->FirstChildElement("NearDistance");
                eResult = camElement->QueryFloatText(&imgPlane.distance);
                camElement = pCamera->FirstChildElement("ImageResolution");     
                str = camElement->GetText();
                sscanf(str, "%d %d", &imgPlane.nx, &imgPlane.ny);
                camElement = pCamera->FirstChildElement("ImageName");
                str = camElement->GetText();
                strcpy(imageName, str);

                cameras.push_back(new Camera(id, imageName, pos, gaze, up, imgPlane));

                pCamera = pCamera->NextSiblingElement("Camera");
        }

        // Parse materals
        pElement = pRoot->FirstChildElement("Materials");
        XMLElement *pMaterial = pElement->FirstChildElement("Material");
        XMLElement *materialElement;
        while(pMaterial != nullptr)
        {
                materials.push_back(new Material());

                int curr = materials.size() - 1;
        
                eResult = pMaterial->QueryIntAttribute("id", &materials[curr]->id);
                materialElement = pMaterial->FirstChildElement("AmbientReflectance");
                str = materialElement->GetText();
                sscanf(str, "%f %f %f", &materials[curr]->ambientRef.r, &materials[curr]->ambientRef.g, &materials[curr]->ambientRef.b);
                materialElement = pMaterial->FirstChildElement("DiffuseReflectance");
                str = materialElement->GetText();
                sscanf(str, "%f %f %f", &materials[curr]->diffuseRef.r, &materials[curr]->diffuseRef.g, &materials[curr]->diffuseRef.b);
                materialElement = pMaterial->FirstChildElement("SpecularReflectance");
                str = materialElement->GetText();
                sscanf(str, "%f %f %f", &materials[curr]->specularRef.r, &materials[curr]->specularRef.g, &materials[curr]->specularRef.b);
                materialElement = pMaterial->FirstChildElement("MirrorReflectance");
                if(materialElement != nullptr)
                {
                        str = materialElement->GetText();
                        sscanf(str, "%f %f %f", &materials[curr]->mirrorRef.r, &materials[curr]->mirrorRef.g, &materials[curr]->mirrorRef.b);
                }
                                else
                {
                        materials[curr]->mirrorRef.r = 0.0;
                        materials[curr]->mirrorRef.g = 0.0;
                        materials[curr]->mirrorRef.b = 0.0;
                }
                materialElement = pMaterial->FirstChildElement("PhongExponent");
                if(materialElement != nullptr)
                        materialElement->QueryIntText(&materials[curr]->phongExp);

                pMaterial = pMaterial->NextSiblingElement("Material");
        }

        // Parse vertex data
        pElement = pRoot->FirstChildElement("VertexData");
        int cursor = 0;
        Vector3f tmpPoint;
        str = pElement->GetText();
        while(str[cursor] == ' ' || str[cursor] == '\t' || str[cursor] == '\n')
                cursor++;
        while(str[cursor] != '\0')
        {
                for(int cnt = 0 ; cnt < 3 ; cnt++)
                {
                        if(cnt == 0)
                                tmpPoint.x = atof(str + cursor);
                        else if(cnt == 1)
                                tmpPoint.y = atof(str + cursor);
                        else
                                tmpPoint.z = atof(str + cursor);
                        while(str[cursor] != ' ' && str[cursor] != '\t' && str[cursor] != '\n')
                                cursor++; 
                        while(str[cursor] == ' ' || str[cursor] == '\t' || str[cursor] == '\n')
                                cursor++;
                }
                vertices.push_back(tmpPoint);
        }

        // Parse objects
        pElement = pRoot->FirstChildElement("Objects");
        
        // Parse spheres
        XMLElement *pObject = pElement->FirstChildElement("Sphere");
        XMLElement *objElement;
        while(pObject != nullptr)
        {
                int id;
                int matIndex;
                int cIndex;
                float R;

                eResult = pObject->QueryIntAttribute("id", &id);
                objElement = pObject->FirstChildElement("Material");
                eResult = objElement->QueryIntText(&matIndex);
                objElement = pObject->FirstChildElement("Center");
                eResult = objElement->QueryIntText(&cIndex);
                objElement = pObject->FirstChildElement("Radius");
                eResult = objElement->QueryFloatText(&R);

                objects.push_back(new Sphere(id, matIndex, cIndex, R, &vertices));

                pObject = pObject->NextSiblingElement("Sphere");
        }

        // Parse triangles
        pObject = pElement->FirstChildElement("Triangle");
        while(pObject != nullptr)
        {
                int id;
                int matIndex;
                int p1Index;
                int p2Index;
                int p3Index;

                eResult = pObject->QueryIntAttribute("id", &id);
                objElement = pObject->FirstChildElement("Material");
                eResult = objElement->QueryIntText(&matIndex);
                objElement = pObject->FirstChildElement("Indices");
                str = objElement->GetText();
                sscanf(str, "%d %d %d", &p1Index, &p2Index, &p3Index);

                objects.push_back(new Triangle(id, matIndex, p1Index, p2Index, p3Index, &vertices));

                pObject = pObject->NextSiblingElement("Triangle");
        }

        // Parse meshes
        pObject = pElement->FirstChildElement("Mesh");
        while(pObject != nullptr)
        {
                int id;
                int matIndex;
                int p1Index;
                int p2Index;
                int p3Index;
                int cursor = 0;
                int vertexOffset = 0;
                std::vector<Triangle> faces;
                std::vector<int> *meshIndices = new std::vector<int>;

                eResult = pObject->QueryIntAttribute("id", &id);
                objElement = pObject->FirstChildElement("Material");
                eResult = objElement->QueryIntText(&matIndex);
                objElement = pObject->FirstChildElement("Faces");
                objElement->QueryIntAttribute("vertexOffset", &vertexOffset);
                str = objElement->GetText();
                while(str[cursor] == ' ' || str[cursor] == '\t' || str[cursor] == '\n')
                        cursor++;
                while(str[cursor] != '\0')
                {
                        for(int cnt = 0 ; cnt < 3 ; cnt++)
                        {
                                if(cnt == 0)
                                        p1Index = atoi(str + cursor) + vertexOffset;
                                else if(cnt == 1)
                                        p2Index = atoi(str + cursor) + vertexOffset;
                                else
                                        p3Index = atoi(str + cursor) + vertexOffset;
                                while(str[cursor] != ' ' && str[cursor] != '\t' && str[cursor] != '\n')
                                        cursor++; 
                                while(str[cursor] == ' ' || str[cursor] == '\t' || str[cursor] == '\n')
                                        cursor++;
                        }
                        faces.push_back(*(new Triangle(-1, matIndex, p1Index, p2Index, p3Index, &vertices)));
                        meshIndices->push_back(p1Index);
                        meshIndices->push_back(p2Index);
                        meshIndices->push_back(p3Index);
                }

                objects.push_back(new Mesh(id, matIndex, faces, meshIndices, &vertices));

                pObject = pObject->NextSiblingElement("Mesh");
        }

        // Parse lights
        int id;
        Vector3f position;
        Vector3f intensity;
        pElement = pRoot->FirstChildElement("Lights");

        XMLElement *pLight = pElement->FirstChildElement("AmbientLight");
        XMLElement *lightElement;
        str = pLight->GetText();
        sscanf(str, "%f %f %f", &ambientLight.r, &ambientLight.g, &ambientLight.b);

        pLight = pElement->FirstChildElement("PointLight");
        while(pLight != nullptr)
        {
                eResult = pLight->QueryIntAttribute("id", &id);
                lightElement = pLight->FirstChildElement("Position");
                str = lightElement->GetText();
                sscanf(str, "%f %f %f", &position.x, &position.y, &position.z);
                lightElement = pLight->FirstChildElement("Intensity");
                str = lightElement->GetText();
                sscanf(str, "%f %f %f", &intensity.r, &intensity.g, &intensity.b);

                lights.push_back(new PointLight(position, intensity));

                pLight = pLight->NextSiblingElement("PointLight");
        }
}

