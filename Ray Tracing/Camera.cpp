#include "Camera.h"
#include "VectorUtility.hpp"

Camera::Camera(
    int id,                      // Id of the camera
    const char* imageName,       // Name of the output PPM file
    const Vector3f& pos,         // Camera position
    const Vector3f& gaze,        // Camera gaze direction
    const Vector3f& up,          // Camera up direction
    const ImagePlane& imgPlane)  // Image plane parameters
    :
    id(id),
    imgPlane(imgPlane),
    position(pos),
    gaze(gaze),
    up(up),
    right(vector_utility::cross(gaze, up)),     // -w X v = u
    pixel_width((imgPlane.right - imgPlane.left) / imgPlane.nx * right),
    pixel_height((imgPlane.bottom - imgPlane.top) / imgPlane.ny * up)
{
    unsigned int i = 0;

    for ( ; *imageName != '\0' && i < 31; ) {
        this->imageName[i++] = *imageName++;
    }
    this->imageName[i] = '\0';  
}

/* Takes coordinate of an image pixel as row and col, and
 * returns the ray going through that pixel. 
 */
Ray Camera::getPrimaryRay(int col, int row) const {
    Vector3f m_minus_e = this->imgPlane.distance * this->gaze;
    Vector3f q = m_minus_e + this->imgPlane.left * right + this->imgPlane.top * up;
    Vector3f s = q + (col + 0.5f) * pixel_width + (row + 0.5f) * pixel_height;

    return Ray{this->position, s};
}
