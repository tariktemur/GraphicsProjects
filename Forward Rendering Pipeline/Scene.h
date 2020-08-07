#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <set>

#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Vec3.h"
#include "Matrix4.h"

class Scene {
public:
	Color backgroundColor;
	bool cullingEnabled;
	int projectionType;

	std::vector< std::vector<Color> > image;
	std::vector< Camera* > cameras;
	std::vector< Vec3* > vertices;
	std::vector< Color* > colorsOfVertices;
	std::vector< Scaling* > scalings;
	std::vector< Rotation* > rotations;
	std::vector< Translation* > translations;
	std::vector< Model* > models;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(std::string ppmFileName, int osType);

private:
    bool modeling_transformation_happened = false;
	std::vector<Matrix4> rotation_matrices;
	std::vector<Matrix4> scaling_matrices;
	std::vector<Matrix4> translation_matrices;
    const double degrees_to_radian = std::acos(-1) / 180;
    int image_width;
    int image_height;

    Matrix4 get_rotation_x_matrix(double theta);
    Matrix4 get_rotation_y_matrix(double theta);
    Matrix4 get_rotation_z_matrix(double theta);
    Matrix4 get_scaling_matrix(double s_x, double s_y, double s_z);
    Matrix4 get_translation_matrix(double t_x, double t_y, double t_z);

    void perform_modeling_transformation();
    void create_modeling_tranformation_matrices();
    void create_rotation_matrices();
    Matrix4 create_rotation_matrix(const Rotation & rotation);
    void create_scaling_matrices();
    void create_translation_matrices();
    void transform_model(Model & model);
    void transform_vertices(std::unordered_set<int> & vertices_in_the_model, const Matrix4 & transformation);

    Matrix4 get_camera_transformation_matrix(Camera& camera);

    Matrix4 get_projection_transformation_matrix(Camera& camera);
    Matrix4 get_orthographic_projection_transform_matrix(Camera& camera);
    Matrix4 get_perspective_projection_transform_matrix(Camera& camera);

    Matrix4 get_viewport_transformation_matrix(Camera& camera);

    void rasterize(const Camera & camera, const Matrix4 & M, const Matrix4 & M_vp);
    void render_wireframe(const Camera & camera, const Model & model, const Matrix4 & M, const Matrix4 & M_vp);

    std::set<std::pair<int, int>> create_line_segments(const Model & model);
    std::pair<int, int> create_line_segment(int vertex_id0, int vertex_id1);

    void set_pixel(int x, int y, Color & color);
    std::vector<Vec3> clip(Vec3, Vec3);
    bool visible(double den, double num, double& t_E, double& t_L);
    void draw_line(Vec3 p0, Vec3 p1);

    void render_solid(const Camera & camera, const Model & model, const Matrix4 & M);

    void draw_triangle(const Vec3 & p0, const Vec3 & p1, const Vec3 & p2);
    double f_01(int x, int y, int x0, int y0, int x1, int y1);
    double f_12(int x, int y, int x1, int y1, int x2, int y2);
    double f_20(int x, int y, int x2, int y2, int x0, int y0);

    bool is_back_face(const Camera & camera, const Triangle & triangle);
};

#endif
