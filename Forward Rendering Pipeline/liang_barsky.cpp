#include <iostream>
#include <vector>


#define x_min 0.0
#define x_max 1.0
#define y_min 0.0
#define y_max 1.0
#define z_min 0.0 // The near plane
#define z_max 1.0 // The far plane


bool visible(float den, float num, float& t_E, float& t_L){
    if(den == 0) {  
        std::cout << "The line is parrallel to one of the sides!!! This case is not yet implemented\n"; 
        return false; 
    }

    float t = num / den;

    std::cout << "The intersection point is at t = " << t << "\n";
    
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

    std::cout << "Done!!\n";
}


bool liang_barsky_tester(std::vector<float>& v0, std::vector<float>& v1){

    float dx = v1[0] - v0[0]; // x_1 - x_0
    float dy = v1[1] - v0[1]; // y_1 - y_0
    float dz = v1[2] - v0[2]; // z_1 - z_0  

    float t_E = 0; // Originally, the entry point is set to be the beginning of the line (v0,v1) 
    float t_L = 1; // and the exit point is set to be the end of line.

    if(visible(dx, x_min - v0[0], t_E, t_L)){ // Does it have a part visible to the left plane?
        if(visible(-dx, v0[0] - x_max, t_E, t_L)){ // Does it have a part visible to the right plane?
            if(visible(dy, y_min - v0[1], t_E, t_L)){ // Does it have a part visible to the bottom plane? 
                if(visible(-dy, v0[1] - y_max, t_E, t_L)){ // Does it have a part visible to the top plane?
                    if(visible(dz, z_min - v0[2], t_E, t_L)){ // Does it have a part visible to the front plane?
                        if(visible(-dz, v0[2] - z_max, t_E, t_L)){ // Does it have a part visible to the back plane?
                            // If the code is here, some part of the line is inside the viewing volume. 
                            std::cout << "Intersection detected, the entry and exit times are as follows\nEntry: " << t_E << "\nExit: " << t_L << "\n";
                            if(t_L < 1){
                                v1 = std::vector<float>{
                                    v0[0] + t_L * dx, // x_1 = x_0 + t_L * (x_1 - x_0)
                                    v0[1] + t_L * dy, // y_1 = y_0 + t_L * (y_1 - y_0)
                                    v0[2] + t_L * dz // z_1 = z_0 + t_L * (z_1 - z_0)
                                    };
                            } // update v1 

                            if(t_E > 0){
                                v0 = std::vector<float>{
                                    v0[0] + t_E * dx, // x_0 = x_0 + t_E * (x_1 - x_0)
                                    v0[1] + t_E * dy, // y_0 = y_0 + t_E * (y_1 - y_0)
                                    v0[2] + t_E * dz // z_0 = z_0 + t_E * (z_1 - z_0)
                                }; // update v0 
                            }
                            return true; 
                        }
                    }
                }
            }
        }
    }


    return false;     
}


int main(){


    std::cout << "The 3D box has the following intervals: \n[x_min, x_max] = [" << x_min << ", " << x_max << 
                    "]\n[y_min, y_max] = [" << y_min << ", " << y_max << "]\n[z_min, z_max] = [" << z_min << ", " << z_max << "]\n";


    std::vector<float> v0{.25, .25, .25};
    std::vector<float> v1{7, 3, 7}; 

    if(liang_barsky_tester(v0, v1))
        std::cout << "The intersection points are :\n[" << v0[0] << ", " << v0[1] << "]\n[" << v1[0] << ", " << v1[1] << "]\n";
    
    else 
        std::cout << "No intersection detected, the line is should be culled!\n";


    return 0; 
}