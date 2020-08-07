#include "Ray.h"
#include "VectorUtility.hpp"

#include <iostream>
#include <vector>




int main(void){

    Vector3f orig1 = {2, 2, 2};
    Vector3f dir1 = {1,0,0};
    Vector3f orig2 = {0, 0, 0};
    Vector3f dir2 = {0,1,0};
    Vector3f orig3 = {0, 0, 0};
    Vector3f dir3 = {0,0,1};


    Vector3f point_true = {9,2,2}; // for t = 7
    Vector3f point_false = {9,2,3}; // wrong direction. 
    Vector3f point_false2 = {9,1,1}; // wrong alignment. 
    Vector3f point4 = {0,8,0}; // t = 8
    Vector3f point5 = {0,0,5}; // t = 5

    Ray Ray1(orig1, dir1);
    Ray Ray2(orig2, dir2);
    Ray Ray3(orig3, dir3);

    std::cout << "----Checking get point for different values of t, we should get (2+t, 2, 2)----\n";
    for(int i = 0; i < 10; i++){
        std::cout << "Get point for t = " << ((float)i) / 10 << " is --> " << Ray1.getPoint(((float)i) / 10) << std::endl;
    }
    

    std::cout << "Checking t value for a true point (we should get t = 7) --> " << Ray1.gett(point_true) << std::endl;
    std::cout << "Checking t value for a false point (direction) --> " << Ray1.gett(point_false) << std::endl;
    std::cout << "Checking t value for a false point (alignment) --> " << Ray1.gett(point_false2) << std::endl;
    std::cout << "Checking t value for a true point (we should get t = 8) --> " << Ray2.gett(point4) << std::endl;
    std::cout << "Checking t value for a true point (we should get t = 5) --> " << Ray3.gett(point5) << std::endl;

   


    return 0;
}
