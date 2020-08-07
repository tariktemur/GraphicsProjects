#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <utility>
#include <algorithm>
#include <limits>
#include <functional>

#include "Scene.h"
#include "Triangle.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;

/*
    Transformations, clipping, culling, rasterization are done here.
    You can define helper functions inside Scene class implementation.
*/
void Scene::forwardRenderingPipeline(Camera *camera) {
    this->initializeImage(camera);

    if (!this->modeling_transformation_happened) {
        this->perform_modeling_transformation();
    }

    auto M_cam = this->get_camera_transformation_matrix(*camera);
    auto M_pro = this->get_projection_transformation_matrix(*camera);
    auto M_vp = this->get_viewport_transformation_matrix(*camera);
    auto M = M_pro * M_cam;

    this->rasterize(*camera, M, M_vp);
}

void Scene::perform_modeling_transformation() {
    this->image_width = static_cast<int>(this->image.size());
    this->image_height = static_cast<int>(this->image[0].size());
    this->vertices.insert(this->vertices.begin(), nullptr); // placeholder value to make the index 1-based
                                                            // ideally this shouldn't be here but in the constructor
    this->colorsOfVertices.insert(this->colorsOfVertices.begin(), nullptr); // placeholder value to make the index 1-based
                                                                            // ideally this shouldn't be here but in the constructor

    this->create_modeling_tranformation_matrices();

    for (Model * model : this->models) {
        this->transform_model(*model);
    }

    this->modeling_transformation_happened = true;
}

void Scene::create_modeling_tranformation_matrices() {
    this->create_rotation_matrices();
    this->create_scaling_matrices();
    this->create_translation_matrices();
}

void Scene::create_rotation_matrices() {
    this->rotation_matrices.push_back({});      // placeholder value to make the index 1-based

    for (const Rotation * rotation : this->rotations) {
        this->rotation_matrices.push_back(this->create_rotation_matrix(*rotation));
    }
}

Matrix4 Scene::create_rotation_matrix(const Rotation& rotation) {
    Vec3 u = {rotation.ux, rotation.uy, rotation.uz};
    Vec3 v = {-rotation.uy, rotation.ux, 0};

    if (u.z != 0) {
        v = {0, -rotation.uz, rotation.uy};
    }

    Vec3 w = vector_utility::cross(u, v);

    vector_utility::normalize(u);
    vector_utility::normalize(v);
    vector_utility::normalize(w);

    Matrix4 reverse_rotation_matrix = matrix_utility::get_identity_matrix();
    reverse_rotation_matrix.val[0][0] = u.x;
    reverse_rotation_matrix.val[0][1] = v.x;
    reverse_rotation_matrix.val[0][2] = w.x;
    reverse_rotation_matrix.val[1][0] = u.y;
    reverse_rotation_matrix.val[1][1] = v.y;
    reverse_rotation_matrix.val[1][2] = w.y;
    reverse_rotation_matrix.val[2][0] = u.z;
    reverse_rotation_matrix.val[2][1] = v.z;
    reverse_rotation_matrix.val[2][2] = w.z;

    Matrix4 rotation_matrix = this->get_rotation_x_matrix(rotation.angle);

    rotation_matrix = rotation_matrix * matrix_utility::transpose(reverse_rotation_matrix);
    rotation_matrix = reverse_rotation_matrix * rotation_matrix;

    return rotation_matrix;
}

void Scene::create_scaling_matrices() {
    this->scaling_matrices.push_back({});       // placeholder value to make the index 1-based

    for (const Scaling * scaling : this->scalings) {
        this->scaling_matrices.push_back(this->get_scaling_matrix(scaling->sx, scaling->sy, scaling->sz));
    }
}

void Scene::create_translation_matrices() {
    this->translation_matrices.push_back({});   // placeholder value to make the index 1-based

    for (const Translation * translation : this->translations) {
        this->translation_matrices.push_back(this->get_translation_matrix(translation->tx, translation->ty, translation->tz));
    }
}

void Scene::transform_model(Model & model) {
    auto vertices_in_the_model = std::unordered_set<int>{};

    for (const auto & triangle : model.triangles) {
        vertices_in_the_model.insert(triangle.vertexIds[0]);
        vertices_in_the_model.insert(triangle.vertexIds[1]);
        vertices_in_the_model.insert(triangle.vertexIds[2]);
    }

    for (int i = 0; i < model.numberOfTransformations; ++i) {
        switch (model.transformationTypes[i]) {
            case 'r':
                transform_vertices(vertices_in_the_model, this->rotation_matrices[model.transformationIds[i]]);
                break;

            case 's':
                transform_vertices(vertices_in_the_model, this->scaling_matrices[model.transformationIds[i]]);
                break;

            case 't':
                transform_vertices(vertices_in_the_model, this->translation_matrices[model.transformationIds[i]]);
                break;

            default:
                // Error
                /* std::cout << "Error: Undefined Transformation\n"; */
                break;
        }
    }
}

void Scene::transform_vertices(std::unordered_set<int> & vertices_in_the_model, const Matrix4 & transformation) {
    for (auto vertex_id : vertices_in_the_model) {
        *this->vertices[vertex_id] = transformation * *this->vertices[vertex_id];
    }
}

Matrix4 Scene::get_camera_transformation_matrix(Camera& camera) {
    Matrix4 camera_translation_matrix = this->get_translation_matrix(-camera.pos.x, -camera.pos.y, -camera.pos.z);
    Matrix4 camera_rotation_matrix = matrix_utility::get_identity_matrix();

    camera_rotation_matrix.val[0][0] = camera.u.x;
    camera_rotation_matrix.val[0][1] = camera.u.y;
    camera_rotation_matrix.val[0][2] = camera.u.z;
    camera_rotation_matrix.val[1][0] = camera.v.x;
    camera_rotation_matrix.val[1][1] = camera.v.y;
    camera_rotation_matrix.val[1][2] = camera.v.z;
    camera_rotation_matrix.val[2][0] = camera.w.x;
    camera_rotation_matrix.val[2][1] = camera.w.y;
    camera_rotation_matrix.val[2][2] = camera.w.z;

    return camera_rotation_matrix * camera_translation_matrix;
}

Matrix4 Scene::get_rotation_x_matrix(double theta) {
    auto matrix = matrix_utility::get_identity_matrix();
    auto theta_as_radian = theta * this->degrees_to_radian;
    auto cos_theta = std::cos(theta_as_radian);
    auto sin_theta = std::sin(theta_as_radian);

    matrix.val[1][1] = cos_theta;
    matrix.val[1][2] = -sin_theta;
    matrix.val[2][1] = sin_theta;
    matrix.val[2][2] = cos_theta;

    return matrix;
}

Matrix4 Scene::get_rotation_y_matrix(double theta) {
    auto matrix = matrix_utility::get_identity_matrix();
    auto theta_as_radian = theta * this->degrees_to_radian;
    auto cos_theta = std::cos(theta_as_radian);
    auto sin_theta = std::sin(theta_as_radian);

    matrix.val[0][0] = cos_theta;
    matrix.val[0][2] = sin_theta;
    matrix.val[2][0] = -sin_theta;
    matrix.val[2][2] = cos_theta;

    return matrix;
}

Matrix4 Scene::get_rotation_z_matrix(double theta) {
    auto matrix = matrix_utility::get_identity_matrix();
    auto theta_as_radian = theta * this->degrees_to_radian;
    auto cos_theta = std::cos(theta_as_radian);
    auto sin_theta = std::sin(theta_as_radian);

    matrix.val[0][0] = cos_theta;
    matrix.val[0][1] = -sin_theta;
    matrix.val[1][0] = sin_theta;
    matrix.val[1][1] = cos_theta;

    return matrix;
}

Matrix4 Scene::get_scaling_matrix(double s_x, double s_y, double s_z) {
    auto matrix = matrix_utility::get_identity_matrix();

    matrix.val[0][0] = s_x;
    matrix.val[1][1] = s_y;
    matrix.val[2][2] = s_z;

    return matrix;
}

Matrix4 Scene::get_translation_matrix(double t_x, double t_y, double t_z) {
    auto matrix = matrix_utility::get_identity_matrix();

    matrix.val[0][3] = t_x;
    matrix.val[1][3] = t_y;
    matrix.val[2][3] = t_z;

    return matrix;
}

Matrix4 Scene::get_projection_transformation_matrix(Camera & camera) {
    Matrix4 projection_transformation_matrix = matrix_utility::get_identity_matrix();

    if (this->projectionType == 1) {
        projection_transformation_matrix = this->get_perspective_projection_transform_matrix(camera);
    }

    projection_transformation_matrix = this->get_orthographic_projection_transform_matrix(camera) * projection_transformation_matrix;

    return projection_transformation_matrix;
}

Matrix4 Scene::get_orthographic_projection_transform_matrix(Camera& camera) {
    auto r_minus_l = camera.right - camera.left;
    auto t_minus_b = camera.top - camera.bottom;
    auto f_minus_n = camera.far - camera.near;

    auto orthographic_projection_transform_matrix = this->get_translation_matrix(-camera.left, -camera.bottom, -camera.near);
    orthographic_projection_transform_matrix = this->get_scaling_matrix(2 / r_minus_l, 2 / t_minus_b, 2 / f_minus_n) * orthographic_projection_transform_matrix;
    orthographic_projection_transform_matrix = this->get_translation_matrix(-1, -1, -1) * orthographic_projection_transform_matrix;

    return orthographic_projection_transform_matrix;
}

Matrix4 Scene::get_perspective_projection_transform_matrix(Camera& camera) {
    auto perspective_projection_transform_matrix = Matrix4{};

    perspective_projection_transform_matrix.val[0][0] = camera.near;
    perspective_projection_transform_matrix.val[1][1] = camera.near;
    perspective_projection_transform_matrix.val[2][2] = camera.far + camera.near;
    perspective_projection_transform_matrix.val[2][3] = camera.far * camera.near;
    perspective_projection_transform_matrix.val[3][2] = -1;

    return perspective_projection_transform_matrix;
}

Matrix4 Scene::get_viewport_transformation_matrix(Camera& camera) {
    auto viewport_transformation_matrix = matrix_utility::get_identity_matrix();

    viewport_transformation_matrix = this->get_translation_matrix(1, 1, 0) * viewport_transformation_matrix;
    viewport_transformation_matrix = this->get_scaling_matrix(camera.horRes / 2.0, camera.verRes / 2.0, 0) * viewport_transformation_matrix;
    viewport_transformation_matrix = this->get_translation_matrix(-0.5, -0.5, 0) * viewport_transformation_matrix;

    return viewport_transformation_matrix;
}

void Scene::rasterize(const Camera & camera, const Matrix4 & M, const Matrix4 & M_vp) {
    for (const Model * model : this->models) {
        if (model->type == 0) {
            this->render_wireframe(camera, *model, M, M_vp);
        }

        else {
            this->render_solid(camera, *model, M_vp * M);
        }
    }
}

void Scene::render_wireframe(const Camera & camera, const Model & model, const Matrix4 & M, const Matrix4 & M_vp) {
    for (const auto & triangle : model.triangles) {
        if (this->cullingEnabled && this->is_back_face(camera, triangle)) {
            continue;
        }

        auto p0_id = triangle.vertexIds[0];
        auto p1_id = triangle.vertexIds[1];
        auto p2_id = triangle.vertexIds[2];

        std::vector<Vec3> l01 =  clip(M * *this->vertices[p0_id], M * *this->vertices[p1_id]);
        std::vector<Vec3> l12 =  clip(M * *this->vertices[p1_id], M * *this->vertices[p2_id]);
        std::vector<Vec3> l20 =  clip(M * *this->vertices[p2_id], M * *this->vertices[p0_id]);

        this->draw_line(M_vp * l01[0], M_vp * l01[1]);
        this->draw_line(M_vp * l12[0], M_vp * l12[1]);
        this->draw_line(M_vp * l20[0], M_vp * l20[1]);
      
        // this->draw_line(M * *this->vertices[p0_id], M * *this->vertices[p1_id]);
        // this->draw_line(M * *this->vertices[p1_id], M * *this->vertices[p2_id]);
        // this->draw_line(M * *this->vertices[p2_id], M * *this->vertices[p0_id]);
    }
}

std::set<std::pair<int, int>> Scene::create_line_segments(const Model & model) {
    auto line_segments = std::set<std::pair<int, int>>{};

    for (const auto & triangle : model.triangles) {
        line_segments.insert(this->create_line_segment(triangle.vertexIds[0], triangle.vertexIds[1]));
    }

    return line_segments;
}

std::pair<int, int> Scene::create_line_segment(int vertex_id0, int vertex_id1) {
    if (vertex_id0 > vertex_id1) {
        return {vertex_id1, vertex_id0};
    }

    else {
        return {vertex_id0, vertex_id1};
    }
}

void Scene::set_pixel(int x, int y, Color & color) {
    if (x >= 0 && x < this->image_width && y >= 0 && y < this->image_height) {
        this->image[x][y] = color;
    }
}

bool Scene::visible(double den, double num, double& t_E, double& t_L){
    if(den == 0) {  
        /* std::cout << "The line is parrallel to one of the sides!!! This case is not yet implemented\n"; */ 
        return den <= 0; 
    }

    double t = num / den;

    /* std::cout << "The intersection point is at t = " << t << "\n"; */
    
    if(den > 0){ // Potentially entering...
        if(t > t_L) // We have left before entering, this means that the line is outside of the box!
            return false; 
        if(t > t_E)// We have entered before, but this is a later entry. 
            t_E = t;
    }
    else{ // Potentially leaving...

        if(t < t_E) // We are leaving before we enter!! This means that the line is outside of the box!
            return false;
        if(t < t_L) // This is a valid exit, if its the most recent exit, keep it!
            t_L = t; 
    }

    /* std::cout << "Done!!\n"; */
    return true;
}

std::vector<Vec3> Scene::clip(Vec3 p0, Vec3 p1) {
    // TODO

    // auto c0 = *this->colorsOfVertices[v0.colorId];
    // auto c1 = *this->colorsOfVertices[v1.colorId];

    Vec3 v0 = p0; 
    Vec3 v1 = p1; 

    double x_min = -1;
    double y_min = -1;
    double z_min = -1;

    double x_max = 1; 
    double y_max = 1; 
    double z_max = 1; 


    double dx = v1.x - v0.x; // x_1 - x_0
    double dy = v1.y - v0.y; // y_1 - y_0
    double dz = v1.z - v0.z; // z_1 - z_0


    double t_E = 0; // Originally, the entry point is set to be the beginning of the line (v0,v1) 
    double t_L = 1; // and the exit point is set to be the end of line.


    if(this->visible(dx, x_min - v0.x, t_E, t_L)){ // Does it have a part this->visible to the left plane?
        if(this->visible(-dx, v0.x - x_max, t_E, t_L)){ // Does it have a part this->visible to the right plane?
            if(this->visible(dy, y_min - v0.y, t_E, t_L)){ // Does it have a part this->visible to the bottom plane? 
                if(this->visible(-dy, v0.y - y_max, t_E, t_L)){ // Does it have a part this->visible to the top plane?
                    if(this->visible(dz, z_min - v0.z, t_E, t_L)){ // Does it have a part this->visible to the front plane?
                        if(this->visible(-dz, v0.z - z_max, t_E, t_L)){ // Does it have a part this->visible to the back plane?
                            // If the code is here, some part of the line is inside the viewing volume. 
                            /* std::cout << "Intersection detected, the entry and exit times are as follows\nEntry: " << t_E << "\nExit: " << t_L << "\n"; */
                            if(t_L < 1){
                                v1 = Vec3(
                                    v0.x + t_L * dx, // x_1 = x_0 + t_L * (x_1 - x_0)
                                    v0.y + t_L * dy, // y_1 = y_0 + t_L * (y_1 - y_0)
                                    v0.z + t_L * dz, // z_1 = z_0 + t_L * (z_1 - z_0)
                                    v1.colorId
                                );
                            } // update v1 

                            if(t_E > 0){
                                v0 = Vec3(
                                    v0.x + t_E * dx, // x_0 = x_0 + t_E * (x_1 - x_0)
                                    v0.y + t_E * dy, // y_0 = y_0 + t_E * (y_1 - y_0)
                                    v0.z + t_E * dz, // z_0 = z_0 + t_E * (z_1 - z_0)
                                    v0.colorId
                                ); // update v0 
                            }
                        }
                    }
                }
            }
        }
    }
    return std::vector<Vec3>({v0, v1});
}

void Scene::draw_line(Vec3 p0, Vec3 p1) {
    if (p0.x > p1.x) {
        auto temp = p0;
        p0 = p1;
        p1 = temp;
    }

    int x0 = std::round(p0.x);
    int y0 = std::round(p0.y);
    int x1 = std::round(p1.x);
    int y1 = std::round(p1.y);

    auto m = std::numeric_limits<double>::infinity();

    if (x0 != x1) {
        m = (y1 - y0) / static_cast<double>((x1 - x0));
    }

    else if (y1 < y0) {
        m = -std::numeric_limits<double>::infinity();
    }

    // f(x, y) = (y0 − y1) * x + (x1 − x0) * y + x0 * y1 − x1 * y0 = 0
    auto x_delta = y0 - y1;
    auto y_delta = x1 - x0;

    if (m > 1) {        // (1, \infty)
        auto x = x0;
        auto d = 2 * y_delta + x_delta;
        auto color = *this->colorsOfVertices[p0.colorId];
        auto color_delta = (*this->colorsOfVertices[p1.colorId] - color) / (y1 - y0);

        for (auto y = y0; y <= y1; ++y) {
            this->set_pixel(x, y, color);

            if (d > 0) {    // NE
                ++x;
                d += 2 * (y_delta + x_delta);
            }

            else {          // N
                d += 2 * y_delta;
            }

            color += color_delta;
        }
    }

    else if (m > 0) {   // (0, 1]
        auto y = y0;
        auto d = 2 * x_delta + y_delta;
        auto color = *this->colorsOfVertices[p0.colorId];
        auto color_delta = (*this->colorsOfVertices[p1.colorId] - color) / (x1 - x0);

        for (auto x = x0; x <= x1; ++x) {
            this->set_pixel(x, y, color);

            if (d < 0) {    // NE
                ++y;
                d += 2 * (x_delta + y_delta);
            }

            else {          // E
                d += 2 * x_delta;
            }

            color += color_delta;
        }
    }

    else if (m > -1) {  // (-1, 0]
        auto y = y0;
        auto d = 2 * x_delta - y_delta;
        auto color = *this->colorsOfVertices[p0.colorId];
        auto color_delta = (*this->colorsOfVertices[p1.colorId] - color) / (x1 - x0);

        for (auto x = x0; x <= x1; ++x) {
            this->set_pixel(x, y, color);

            if (d > 0) {    // SE
                --y;
                d += 2 * (x_delta - y_delta);
            }

            else {          // E
                d += 2 * x_delta;
            }

            color += color_delta;
        }
    }

    else {              // (-\infty, -1]
        auto x = x0;
        auto d = x_delta - 2 * y_delta;
        auto color = *this->colorsOfVertices[p0.colorId];
        auto color_delta = (*this->colorsOfVertices[p1.colorId] - color) / (y0 - y1);

        for (auto y = y0; y >= y1; --y) {
            this->set_pixel(x, y, color);

            if (d < 0) {    // SE
                ++x;
                d += 2 * (x_delta - y_delta);
            }

            else {          // S
                d -= 2 * y_delta;
            }

            color += color_delta;
        }
    }
}

void Scene::render_solid(const Camera & camera, const Model & model, const Matrix4 & M) {
    for (const auto & triangle : model.triangles) {
        if (this->cullingEnabled && this->is_back_face(camera, triangle)) {
            continue;
        }

        auto p0_id = triangle.vertexIds[0];
        auto p1_id = triangle.vertexIds[1];
        auto p2_id = triangle.vertexIds[2];

        this->draw_triangle(M * *this->vertices[p0_id], M * *this->vertices[p1_id], M * *this->vertices[p2_id]);
    }
}

bool Scene::is_back_face(const Camera & camera, const Triangle & triangle) {
    auto p0 = *this->vertices[triangle.vertexIds[0]];
    auto p1 = *this->vertices[triangle.vertexIds[1]];
    auto p2 = *this->vertices[triangle.vertexIds[2]];

    auto v = p0 - camera.pos;
    auto n = vector_utility::cross(p1 - p0, p2 - p0);

    return vector_utility::dot(v, n) > 0;
}

void Scene::draw_triangle(const Vec3 & p0, const Vec3 & p1, const Vec3 & p2) {
    int x0 = std::round(p0.x);
    int y0 = std::round(p0.y);
    int x1 = std::round(p1.x);
    int y1 = std::round(p1.y);
    int x2 = std::round(p2.x);
    int y2 = std::round(p2.y);

    auto x_min = std::max(0, std::min({this->image_width - 1, x0, x1, x2}));
    auto x_max = std::min(this->image_width - 1, std::max({0, x0, x1, x2}));
    auto y_min = std::max(0, std::min({this->image_height - 1, y0, y1, y2}));
    auto y_max = std::min(this->image_height - 1, std::max({0, y0, y1, y2}));

    for (int x = x_min; x <= x_max; ++x) {
        for (int y = y_min; y <= y_max; ++y) {
            auto alpha = this->f_12(x, y, x1, y1, x2, y2) / this->f_12(x0, y0, x1, y1, x2, y2);
            auto beta = this->f_20(x, y, x2, y2, x0, y0) / this->f_20(x1, y1, x2, y2, x0, y0);
            auto gamma = this->f_01(x, y, x0, y0, x1, y1) / this->f_01(x2, y2, x0, y0, x1, y1);

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                this->image[x][y] =
                    *this->colorsOfVertices[p0.colorId] * alpha +
                    *this->colorsOfVertices[p1.colorId] * beta +
                    *this->colorsOfVertices[p2.colorId] * gamma;
            }
        }
    }
}

double Scene::f_01(int x, int y, int x0, int y0, int x1, int y1) {
    return static_cast<double>((y0 - y1) * x + (x1 - x0) * y + x0 * y1 - x1 * y0);
}

double Scene::f_12(int x, int y, int x1, int y1, int x2, int y2) {
    return static_cast<double>((y1 - y2) * x + (x2 - x1) * y + x1 * y2 - x2 * y1);
}

double Scene::f_20(int x, int y, int x2, int y2, int x0, int y0) {
    return static_cast<double>((y2 - y0) * x + (x0 - x2) * y + x2 * y0 - x0 * y2);
}

/*
    Parses XML file
*/
Scene::Scene(const char *xmlPath) {
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL) {
        pElement->QueryBoolText(&cullingEnabled);
    }

    // read projection type
    pElement = pRoot->FirstChildElement("ProjectionType");
    if (pElement != NULL) {
        pElement->QueryIntText(&projectionType);
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL) {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze = normalizeVec3(cam->gaze);
        cam->u = crossProductVec3(cam->gaze, cam->v);
        cam->u = normalizeVec3(cam->u);

        cam->w = inverseVec3(cam->gaze);
        cam->v = crossProductVec3(cam->u, cam->gaze);
        cam->v = normalizeVec3(cam->v);

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = std::string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL) {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL) {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL) {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL) {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read models
    pElement = pRoot->FirstChildElement("Models");

    XMLElement *pModel = pElement->FirstChildElement("Model");
    while (pModel != NULL) {
        Model *model = new Model();

        pModel->QueryIntAttribute("id", &model->modelId);
        pModel->QueryIntAttribute("type", &model->type);

        // read model transformations
        XMLElement *pTransformations = pModel->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        pTransformations->QueryIntAttribute("count", &model->numberOfTransformations);

        while (pTransformation != NULL) {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            model->transformationTypes.push_back(transformationType);
            model->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        // read model triangles
        XMLElement *pTriangles = pModel->FirstChildElement("Triangles");
        XMLElement *pTriangle = pTriangles->FirstChildElement("Triangle");

        pTriangles->QueryIntAttribute("count", &model->numberOfTriangles);

        while (pTriangle != NULL) {
            int v1, v2, v3;

            str = pTriangle->GetText();
            sscanf(str, "%d %d %d", &v1, &v2, &v3);

            model->triangles.push_back(Triangle(v1, v2, v3));

            pTriangle = pTriangle->NextSiblingElement("Triangle");
        }

        models.push_back(model);

        pModel = pModel->NextSiblingElement("Model");
    }
}

/*
    Initializes image with background color
*/
void Scene::initializeImage(Camera *camera) {
    if (this->image.empty()) {
        for (int i = 0; i < camera->horRes; i++) {
            std::vector<Color> rowOfColors;

            for (int j = 0; j < camera->verRes; j++) {
                rowOfColors.push_back(this->backgroundColor);
            }

            this->image.push_back(rowOfColors);
        }
    }
    // if image is filled before, just change color rgb values with the background color
    else {
        for (int i = 0; i < camera->horRes; i++) {
            for (int j = 0; j < camera->verRes; j++) {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
            }
        }
    }
}

/*
    If given value is less than 0, converts value to 0.
    If given value is more than 255, converts value to 255.
    Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value) {
    if (value >= 255.0) {
        return 255;
    }

    if (value <= 0.0) {
        return 0;
    }

    return (int)(value);
}

/*
    Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera) {
    std::ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << '\n';
    fout << "# " << camera->outputFileName << '\n';
    fout << camera->horRes << " " << camera->verRes << '\n';
    fout << "255" << '\n';

    for (int j = camera->verRes - 1; j >= 0; j--) {
        for (int i = 0; i < camera->horRes; i++) {
            fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
                 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
        }
        fout << '\n';
    }
    fout.close();
}

/*
    Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
    os_type == 1        -> Ubuntu
    os_type == 2        -> Windows
    os_type == other    -> No conversion
*/
void Scene::convertPPMToPNG(std::string ppmFileName, int osType) {
    std::string command;

    // call command on Ubuntu
    if (osType == 1) {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // call command on Windows
    else if (osType == 2) {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // default action - don't do conversion
    else {
    }
}
