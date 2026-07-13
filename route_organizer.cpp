#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <algorithm>
#include <cstdint>
#include <nanoflann.hpp>

using namespace std;

bool isLand(double lat, double lon);

vector<uint8_t> landMask;
vector<uint8_t> rawMask;
const float pi = 3.14159f;

struct WeatherPoint {
    float lat, lon, temp;
};


struct WeatherCloud {
    std::vector<WeatherPoint> pts;

    inline size_t kdtree_get_point_count() const { return pts.size(); }


    inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].lat;
        else return pts[idx].lon;

    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }

};

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



static inline bool parseLine(const char* start, const char* end,
    vector<float>& row, int colNum) {
    row.clear();
    row.reserve(colNum);
    const char* ptr = start;

    while (ptr < end && (int)row.size() < colNum) {
        while (ptr < end && (*ptr == ' ' || *ptr == ',' ||
            *ptr == '\t' || *ptr == '\r'))
            ptr++;

        if (ptr >= end || *ptr == '\n') break;

        char* next;
        float val = strtof(ptr, &next);

        if (next == ptr) { ptr++; continue; }

        ptr = next;

        row.push_back(val);
    }
    return (int)row.size() == colNum;

}





static vector<size_t> findLineStarts(const char* buf, size_t size, bool skipHeader) {
    vector<size_t> starts;
    starts.reserve(size / 16);

    size_t i = 0;

    if (skipHeader) {
        while (i < size && buf[i] != '\n') i++;
        if (i < size) i++;
    }

    if (i < size) starts.push_back(i);

    for (; i < size; i++) {
        if (buf[i] == '\n' && i + 1 < size)
            starts.push_back(i + 1);
    }
    return starts;

}





vector<vector<float>> readCSV_fast(const string& filename, int colNum, bool skipHeader = true) {

    ifstream file(filename, ios::binary);
    if (!file.is_open()) {
        cerr << "Dosya acilamadi: " << filename << "\n";
        return {};
    }

    file.rdbuf()->pubsetbuf(nullptr, 1 << 23);
    file.seekg(0, ios::end);
    size_t fileSize = (size_t)file.tellg();
    file.seekg(0, ios::beg);
    string buf(fileSize, '\0');
    file.read(&buf[0], fileSize);
    file.close();

    vector<size_t> starts = findLineStarts(buf.data(), fileSize, skipHeader);
    size_t nLines = starts.size();

    vector<vector<float>> data;

    data.reserve(nLines);

    vector<float> row;

    for (size_t i = 0; i < nLines; i++) {
        const char* lineStart = buf.data() + starts[i];
        const char* lineEnd = (i + 1 < nLines) ? buf.data() + starts[i + 1] : buf.data() + fileSize;
        if (parseLine(lineStart, lineEnd, row, colNum))
            data.push_back(move(row));
    }
    return data;

}



vector<vector<float>> readCSV(const string& filename, int colNum) {

    vector<vector<float>> data;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Dosya acilamadi: " << filename << endl;
        return data;
    }

    string line;
    int line_count = 0;

    while (getline(file, line)) {
        line_count++;
        if (line.empty()) continue;
        replace(line.begin(), line.end(), ',', ' ');
        vector<float> row;
        stringstream ss(line);
        ss.imbue(locale("C"));
        float value;

        while (ss >> value) {
            row.push_back(value);
            if (row.size() == colNum) break; // rowun eleman sayisi coldur
        }

        if (row.size() == colNum) {
            data.push_back(row);
        }

        else if (data.empty() && line_count == 1) {
            cout << "Ilk satir dogru okunamadi. Sadece " << row.size() << " sayi bulundu.\n";
            cout << "Okunmaya calisilan satir: " << line << "\n";
        }
    }

    file.close();
    return data;

}



void saveDPMap(const vector<Stage>& dp_map, const string& filename) {

    ofstream out(filename, ios::binary);
    size_t stage_count = dp_map.size();
    out.write((char*)&stage_count, sizeof(stage_count));

    for (const auto& stage : dp_map) {
        size_t node_count = stage.nodes.size();
        out.write((char*)&node_count, sizeof(node_count));

        for (const auto& n : stage.nodes) {
            out.write((char*)&n.lat, sizeof(float));
            out.write((char*)&n.lon, sizeof(float));
            out.write((char*)&n.temp, sizeof(float));
        }
    }

    out.close();

}

vector<Stage> loadDPMap(const string& filename) {

    vector<Stage> dp_map;
    ifstream in(filename, ios::binary);

    if (!in.is_open()) return dp_map;

    size_t stage_count;

    in.read((char*)&stage_count, sizeof(stage_count));

    dp_map.resize(stage_count);

    for (size_t i = 0; i < stage_count; ++i) {
        size_t node_count;
        in.read((char*)&node_count, sizeof(node_count));
        dp_map[i].nodes.resize(node_count);

        for (size_t j = 0; j < node_count; ++j) {
            Node& n = dp_map[i].nodes[j];

            in.read((char*)&n.lat, sizeof(float));
            in.read((char*)&n.lon, sizeof(float));
            in.read((char*)&n.temp, sizeof(float));

            n.total_cost = numeric_limits<float>::infinity();
            n.back_pointer = -1;
            n.incoming_heading = -1;
        }
    }

    in.close();

    return dp_map;

}





// iki koordinat arasındaki gidiş açısını hesaplamak içim

float calculateHeading(float lat1, float lon1, float lat2, float lon2) {

    float dLon = (lon2 - lon1) * pi / 180.0f;
    lat1 = lat1 * pi / 180.0f;
    lat2 = lat2 * pi / 180.0f;

    float y = sin(dLon) * cos(lat2);
    float x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon);
    float heading = atan2(y, x) * 180.0f / pi;
    return fmod((heading + 360.0f), 360.0f);
}



float haversine(float lat1, float lon1, float lat2, float lon2) {
    float R = 6371.0f;
    float dLat = (lat2 - lat1) * (pi / 180.0f);
    float dLon = (lon2 - lon1) * (pi / 180.0f);

    float a = sinf(dLat * 0.5f) * sinf(dLat * 0.5f) +
        cosf(lat1 * (pi / 180.0f)) * cosf(lat2 * (pi / 180.0f)) *
        sinf(dLon * 0.5f) * sinf(dLon * 0.5f);

    float c = 2.0f * atan2f(sqrtf(a), sqrtf(1.0f - a));
    return R * c;
}





vector<Point> interpolateRoute(const vector<vector<float>>& original_waypoints, float interval_km) {

    vector<Point> centers;

    for (size_t i = 0; i < original_waypoints.size() - 1; ++i) {

        float start_lat = original_waypoints[i][0];
        float start_lon = original_waypoints[i][1];
        float end_lat = original_waypoints[i + 1][0];
        float end_lon = original_waypoints[i + 1][1];
        float segment_dist = haversine(start_lat, start_lon, end_lat, end_lon);

        // segmentler arası 5km lik kaç nokta olur

        int num_points = max(1, (int)floor(segment_dist / interval_km));

        // noktaları üret ve listeye at

        for (int j = 0; j < num_points; ++j) {

            float fraction = (float)j / num_points;
            float inter_lat = start_lat + (end_lat - start_lat) * fraction;
            float inter_lon = start_lon + (end_lon - start_lon) * fraction;
            centers.push_back({ inter_lat, inter_lon });
        }
    }

    // varış noktasını en son merkez olarak ekle

    centers.push_back({ original_waypoints.back()[0], original_waypoints.back()[1] });
    return centers;

}

vector<Point> adaptiveInterpolate(const vector<vector<float>>& original_waypoints, float max_distance_km) {

    vector<Point> centers;
    if (original_waypoints.empty()) return centers;

    for (size_t i = 0; i < original_waypoints.size() - 1; ++i) {

        float start_lat = original_waypoints[i][0];
        float start_lon = original_waypoints[i][1];
        float end_lat = original_waypoints[i + 1][0];
        float end_lon = original_waypoints[i + 1][1];
        centers.push_back({ start_lat, start_lon });
        float segment_dist = haversine(start_lat, start_lon, end_lat, end_lon);

        if (segment_dist > max_distance_km) {
            int num_splits = max(1, (int)ceil(segment_dist / max_distance_km));
            for (int j = 1; j < num_splits; ++j) {
                float fraction = (float)j / num_splits;
                float inter_lat = start_lat + (end_lat - start_lat) * fraction;
                float inter_lon = start_lon + (end_lon - start_lon) * fraction;
                centers.push_back({ inter_lat, inter_lon });
            }
        }
    }
    centers.push_back({ original_waypoints.back()[0], original_waypoints.back()[1] });

    return centers;

}



bool isCriticalZone(float lat, float lon) {

    // Süveyş Kanalı (algoritma sadece burayı handle'layamıyor saç teli gibi biyer...) 

    if (lat >= 30.55f && lat <= 31.35f && lon >= 32.20f && lon <= 32.40f) return true; //ilk kısmı
    if (lat >= 29.90f && lat <= 30.55f && lon >= 32.30f && lon <= 32.65f) return true; //ikinci kısmı
    //if (lat >= 29.90f && lat <= 31.30f && lon >= 32.25f && lon <= 32.65f) return true;
    //if (lat >= 27.00f && lat <= 31.50f && lon >= 32.00f && lon <= 34.50f) return true;
    // Panama Kanalı ve girişleri
    if (lat >= 8.85f && lat <= 9.40f && lon >= -80.00f && lon <= -79.50f) return true;

    // Cebelitarık Boğazı
    if (lat >= 35.70f && lat <= 36.20f && lon >= -6.00f && lon <= -5.20f) return true;

    // Singapur Boğazı geçişi
    if (lat >= 1.10f && lat <= 1.45f && lon >= 103.50f && lon <= 104.10f) return true;

    // Babülmendep Boğazı
    if (lat >= 12.30f && lat <= 12.80f && lon >= 43.20f && lon <= 43.50f) return true;

    // İstanbul Boğazı, Marmara Denizi ve Çanakkale Boğazı koridoru
    if (lat >= 40.00f && lat <= 41.30f && lon >= 26.10f && lon <= 29.20f) return true;

    // Oresund (Danimarka) ve Kiel Kanalı (Almanya) bölgesi
    if (lat >= 53.30f && lat <= 56.30f && lon >= 9.80f && lon <= 13.00f) return true;

    return false;

}





vector<Stage> buildDPMap(const vector<Point>& centers, const vector<vector<float>>& weather_data, float radius_km) {

    vector<Stage> dp_map;

    int emptyNodes = 0;

    WeatherCloud cloud;

    cloud.pts.reserve(weather_data.size());

    for (const auto& row : weather_data) {
        if (row.size() >= 3) {
            cloud.pts.push_back({ row[0], row[1], row[2] });

        }

    }



    typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, WeatherCloud>, WeatherCloud, 2> my_kd_tree_t;

    my_kd_tree_t index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));

    index.buildIndex();

    float search_radius_deg = radius_km / 111.0f;

    float search_radius_sq = search_radius_deg * search_radius_deg;

    for (const auto& center : centers) {

        Stage current_stage;

        if (!isCriticalZone(center.lat, center.lon)) {
            float query_pt[2] = { (float)center.lat, (float)center.lon };
            vector<nanoflann::ResultItem<uint32_t, float>> ret_matches;
            nanoflann::SearchParameters params;
            index.radiusSearch(&query_pt[0], search_radius_sq, ret_matches, params);
            for (const auto& match : ret_matches) {

                const WeatherPoint& wp = cloud.pts[match.first];
                float dist = haversine(center.lat, center.lon, wp.lat, wp.lon);

                if (dist > radius_km) continue;

                if (isLand(wp.lat, wp.lon)) continue;

                current_stage.nodes.push_back(Node(wp.lat, wp.lon, wp.temp));
            }
        }


        if (current_stage.nodes.empty()) {
            float query_pt[2] = { (float)center.lat, (float)center.lon };
            size_t ret_index = 0;
            float out_dist_sqr = 0.0f;
            nanoflann::KNNResultSet<float> resultSet(1);
            resultSet.init(&ret_index, &out_dist_sqr);
            index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParameters(10));
            float real_temp = cloud.pts[ret_index].temp;
            current_stage.nodes.push_back(Node(center.lat, center.lon, real_temp));

        }
        dp_map.push_back(current_stage);

    }

    ofstream json_out("points.json");

    if (json_out.is_open()) {
        json_out << "[\n";
        bool first_node = true;
        int total_survivors = 0;

        for (size_t i = 0; i < dp_map.size(); ++i) {
            for (size_t j = 0; j < dp_map[i].nodes.size(); ++j) {
                if (!first_node) {
                    json_out << ",\n";
                }
                json_out << "  [" << dp_map[i].nodes[j].lat << ", " << dp_map[i].nodes[j].lon << "]";
                first_node = false;
                total_survivors++;
            }
        }

        json_out << "\n]\n";
        json_out.close();
        cout << "Toplam hayatta kalan nokta: " << total_survivors << endl;

    }

    else {
        cerr << "surviving_nodes.json dosyasi olusturulamadi" << endl;

    }
    
    return dp_map;

}


void loadLandMask() {
    std::ifstream is("dilated_landmask12km.bin", std::ios::binary);
    if (is) {
        is.seekg(0, std::ios::end);
        size_t size = is.tellg();
        is.seekg(0, std::ios::beg);
        landMask.resize(size);
        is.read((char*)landMask.data(), size);
    }
    else {

        exit(1);
    }

    std::ifstream is_raw("land_mask_grid.bin", std::ios::binary);
    if (is_raw) {
        is_raw.seekg(0, std::ios::end);
        size_t size = is_raw.tellg();
        is_raw.seekg(0, std::ios::beg);
        rawMask.resize(size);
        is_raw.read((char*)rawMask.data(), size);
    }
    else {
        
        exit(1);
    }

    cout << "maskler yuklendi" << endl;
}



bool isLand(double lat, double lon) {
    int64_t x = static_cast<int64_t>(std::floor((lon + 180.0) / 0.005));
    int64_t y = static_cast<int64_t>(std::floor((90.0 - lat) / 0.005));

    if (x < 0) x = 0; if (x >= 72000) x = 71999;

    if (y < 0) y = 0; if (y >= 36000) y = 35999;

    uint64_t index = static_cast<uint64_t>(y) * 72000 + x;

    uint64_t byteIndex = index / 8;
    uint8_t bitOffset = index % 8;

    if (byteIndex >= landMask.size()) return false;

    return (landMask[byteIndex] & (1U << bitOffset)) != 0;

}

bool isSafePath(float lat1, float lon1, float lat2, float lon2) {
    int x1 = static_cast<int>((lon1 + 180.0) * 200.0);
    int y1 = static_cast<int>((90.0 - lat1) * 200.0);
    int x2 = static_cast<int>((lon2 + 180.0) * 200.0);
    int y2 = static_cast<int>((90.0 - lat2) * 200.0);
    int dx = abs(x2 - x1), sx = x1 < x2 ? 1 : -1;
    int dy = -abs(y2 - y1), sy = y1 < y2 ? 1 : -1;
    int err = dx + dy, e2;

    while (true) {
        if (x1 >= 0 && x1 < 72000 && y1 >= 0 && y1 < 36000) {
            uint64_t index = static_cast<uint64_t>(y1) * 72000 + x1;

            uint64_t byteIndex = index / 8;
            uint8_t bitOffset = index % 8;

            if (byteIndex < landMask.size()) {
                if ((landMask[byteIndex] & (1U << bitOffset)) != 0) return false; // kara
            }
            else {
				return false; // harita dışı alanlar
            }
        }

        if (x1 == x2 && y1 == y2) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x1 += sx; }
        if (e2 <= dx) { err += dx; y1 += sy; }
    }
    return true;
}



int checkPathSafety(float lat1, float lon1, float lat2, float lon2) {
    int x1 = static_cast<int>((lon1 + 180.0) * 200.0);
    int y1 = static_cast<int>((90.0 - lat1) * 200.0);
    int x2 = static_cast<int>((lon2 + 180.0) * 200.0);
    int y2 = static_cast<int>((90.0 - lat2) * 200.0);
    int dx = abs(x2 - x1), sx = x1 < x2 ? 1 : -1;
    int dy = -abs(y2 - y1), sy = y1 < y2 ? 1 : -1;
    int err = dx + dy, e2;

    int penalty_score = 0;

    // Bu çizginin kritik bir kanalda (Süveyş vb.) olup olmadığını kontrol et
    bool in_critical = isCriticalZone(lat1, lon1) || isCriticalZone(lat2, lon2);

    while (true) {
        if (x1 >= 0 && x1 < 72000 && y1 >= 0 && y1 < 36000) {
            uint64_t index = static_cast<uint64_t>(y1) * 72000 + x1;
            uint64_t byteIndex = index / 8;
            uint8_t bitOffset = index % 8;

            // 1. GERÇEK KARA KONTROLÜ (RAW MASK)
            if ((rawMask[byteIndex] & (1U << bitOffset)) != 0) {
                if (in_critical) {
                    // KANALDAYIZ: Zinciri koparma! Sadece devasa ceza ver.
                    // DP mecburen en az kumdan (en çok sudan) geçen çizgiyi seçecek.
                    penalty_score += 100;
                }
                else {
                    // AÇIK DENİZ: Gerçek karaya çarptı, affetme, yolu kapat.
                    return -1;
                }
            }
            // 2. GÜVENLİK BÖLGESİ KONTROLÜ (LAND MASK 12KM)
            else if ((landMask[byteIndex] & (1U << bitOffset)) != 0) {
                penalty_score += 1; // Normal dolgu cezası
            }
        }
        else {
            return -1; // Harita dışı
        }

        if (x1 == x2 && y1 == y2) break;
        e2 = 2 * err;
        if (e2 >= dy) { err += dy; x1 += sx; }
        if (e2 <= dx) { err += dx; y1 += sy; }
    }

    return penalty_score;
}