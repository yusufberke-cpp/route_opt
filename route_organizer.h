#pragma once
#include <vector>

using namespace std;

struct Node
{
    float lat, lon, temp;
    float total_cost;
    int back_pointer = -1;
    float incoming_heading;

    Node()
        : lat(0), lon(0), temp(0),
        total_cost(numeric_limits<float>::infinity()),
        back_pointer(-1), incoming_heading(-1) {
    }

    Node(float la, float lo, float t)
        : lat(la), lon(lo), temp(t),
        total_cost(numeric_limits<float>::infinity()),
        back_pointer(-1),
        incoming_heading(-1) {
    }
};

struct Stage { //her noktanın node struct'ını tutcak array

    vector<Node> nodes;

};



struct Point {

    double lat, lon;

};
extern std::vector<uint8_t> landMask;
extern vector<uint8_t> rawMask;
void loadLandMask();
bool isLand(double lat, double lon);
bool isSafePath(float lat1, float lon1, float lat2, float lon2);
bool isCriticalZone(float lat, float lon);
float haversine(float lat1, float lon1, float lat2, float lon2);
float calculateHeading(float lat1, float lon1, float lat2, float lon2);
int checkPathSafety(float lat1, float lon1, float lat2, float lon2);
vector<Stage> loadDPMap(const string& filename);
void saveDPMap(const vector<Stage>& dp_map, const string& filename);
vector<vector<float>> readCSV_fast(const string& filename, int colNum, bool skipHeader = true);

vector<Point> adaptiveInterpolate(const vector<vector<float>>& original_waypoints, float max_distance_km);
vector<Point> interpolateRoute(const vector<vector<float>>& original_waypoints, float interval_km);
vector<Stage> buildDPMap(const vector<Point>& centers, const vector<vector<float>>& weather_data, float radius_km);