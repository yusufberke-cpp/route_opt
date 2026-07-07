#include <math.h> 

const float pi = 3.14159f;
const float min_Vreal = 0.5f;    
const float lambda_turn = 0.02f;

const float temp_factor = 0.5f;

const float k_head = 0.05f;     
const float k_cross = 0.04f;    


inline float deg2rad(float deg) {
    return deg * (pi / 180.0f);
}

inline float angleDiff(float a, float b) {
    float d = fabs(a - b);
    if (d > 180.0f) d = 360.0f - d;
    return d;
}



float calcCostAF(float inv_ship_speed, float distance, float temp, float prev_ship_dir, float new_ship_dir) {

    float time_cost = distance * inv_ship_speed;

    float cost = time_cost + (temp * temp_factor);

    if (prev_ship_dir >= 0.0f) {
        float turn_angle = angleDiff(new_ship_dir, prev_ship_dir);

    
        if (turn_angle > 20.0f) {
            cost += lambda_turn * (turn_angle - 20.0f);
        }
    }

    return cost;
}