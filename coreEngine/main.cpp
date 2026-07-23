#define NOMINMAX 


#include <stdio.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <chrono>
#include <functional>
#include <iomanip>
#include <fstream>
#include <utility> 
#include <cmath>   
#include <cstdlib>

#include "dp_opt.h"
#include "route_organizer.h"
#include <nanoflann.hpp>

using namespace std;

template <typename Func, typename... Args>
auto time_it(Func&& func, Args&&... args) {
    auto start = std::chrono::high_resolution_clock::now();

    auto result = std::invoke(std::forward<Func>(func), std::forward<Args>(args)...);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    std::cout << " Gecen Sure: " << elapsed.count() << " ms\n";
    return result;
}

struct WeatherPoint {
    float lat, lon, temp;
};

struct WeatherCloud {
    std::vector<WeatherPoint> pts;

    // nanoflann arayüz fonksiyonları
    inline size_t kdtree_get_point_count() const { return pts.size(); }

    inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
        if (dim == 0) return pts[idx].lat;
        else return pts[idx].lon;
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};


void clearPortArea(float lat, float lon, int radius_pixels = 10) {
    int center_x = (int)((lon + 180.0) / 0.005);
    int center_y = (int)((90.0 - lat) / 0.005);

    for (int dy = -radius_pixels; dy <= radius_pixels; ++dy) {
        for (int dx = -radius_pixels; dx <= radius_pixels; ++dx) {
            if (dx * dx + dy * dy <= radius_pixels * radius_pixels) {
                int x = center_x + dx;
                int y = center_y + dy;

                if (x >= 0 && x < 72000 && y >= 0 && y < 36000) {
                    uint64_t index = static_cast<uint64_t>(y) * 72000 + x;
                    uint64_t byteIndex = index / 8;
                    uint8_t bitOffset = index % 8;

                    landMask[byteIndex] &= static_cast<uint8_t>(~(1U << bitOffset));

                    rawMask[byteIndex] &= static_cast<uint8_t>(~(1U << bitOffset));
                }
            }
        }
    }
}

// her çemberde 69 nokta var
// lat-lon / enlem-boylam
// hızlar knot cinsinden 1 knot = 0.514444444 m / s
int main() {

    loadLandMask();

    vector<vector<float>> original_waypoints = readCSV_fast("original_wp.csv", 2);
    if (original_waypoints.empty()) {
        cerr << "Hata: original_wp.csv dosyasi bulunamadi veya bos" << endl;
    }
    else {
        cout << "Waypointler okundu: " << original_waypoints.size() << " adet." << endl;
        clearPortArea(original_waypoints.front()[0], original_waypoints.front()[1], 15); // Kalkış
        clearPortArea(original_waypoints.back()[0], original_waypoints.back()[1], 15);   // Varış
    }

    // orijinal rotadaki noktalar arasinda cizgi cekip etrafindaki pikselleri deniz (0) yaptim ki kanallarda rota bulunmazligi yasanmasin direk ana rotayı deniz gorsun çünkü buffer boğazları birleştiriyor
    for (size_t i = 0; i < original_waypoints.size() - 1; ++i) {

        double lat1 = original_waypoints[i][0];
        double lon1 = original_waypoints[i][1];
        double lat2 = original_waypoints[i + 1][0];
        double lon2 = original_waypoints[i + 1][1];

        if (true) {
            double dLat = lat2 - lat1;
            double dLon = lon2 - lon1;
            double dist = std::sqrt(dLat * dLat + dLon * dLon);

            int steps = static_cast<int>(dist / 0.001) + 1;

            int tunnel_radius = 6;

            for (int j = 0; j <= steps; ++j) {
                double t = (double)j / steps;
                double currLat = lat1 + t * dLat;
                double currLon = lon1 + t * dLon;
                bool in_canal = isCriticalZone(currLat, currLon);

                for (int dy = -tunnel_radius; dy <= tunnel_radius; ++dy) {
                    for (int dx = -tunnel_radius; dx <= tunnel_radius; ++dx) {


                        if (dx * dx + dy * dy <= tunnel_radius * tunnel_radius) {

                            int x = (int)((currLon + 180.0) / 0.005) + dx;
                            int y = (int)((90.0 - currLat) / 0.005) + dy;

                            if (x >= 0 && x < 72000 && y >= 0 && y < 36000) {
                                uint64_t index = static_cast<uint64_t>(y) * 72000 + x;
                                uint64_t byteIndex = index / 8;
                                uint8_t bitOffset = index % 8;

                                //if (in_canal || (rawMask[byteIndex] & (1U << bitOffset)) == 0) {

                                    landMask[byteIndex] &= static_cast<uint8_t>(~(1U << bitOffset));

                                //}

                            }
                        }
                    }
                }
            }
        }

    }

	float max_gap_km = 150.0f; // kaç km'den fazla boşluk varsa araya sanal nokta ekleyecek
    vector<Point> centers = adaptiveInterpolate(original_waypoints, max_gap_km);

    cout << "Orijinal wp sayisi: " << original_waypoints.size() << endl;
    cout << "Uretilen toplam dp merkezi (Sanal noktalarla): " << centers.size() << endl;

    cout << "hava durumu verisi okunuyor" << endl;

    vector<vector<float>> weather_data = readCSV_fast("veri.txt", 3);
    if (weather_data.size() == 0) {
        cerr << "Hata: weather.csv okunamadi buildDPMap bos node uretecektir" << endl;
    }

    cout << "okunan toplam veri satiri: " << weather_data.size() << endl;
    cout << "veriler merkezlere yerlestiriliyor" << endl;

    float search_radius_km = 50.0f; // her merkez etrafindaki cemberin yaricapi (rotadan ne kadar sapcağınıda bu belirliyor)

    vector<Stage> dp_map;
    ifstream test("dp_cache.bin", ios::binary);

    if (test.good()) {
        cout << "cache yukleniyor" << endl;
        dp_map = loadDPMap("dp_cache.bin");
    }
    else {
        cout << "cache yok, hesaplaniyor" << endl;
        dp_map = buildDPMap(centers, weather_data, search_radius_km);
    }

    cout << "dp haritasi hazir: " << dp_map.size() << endl;
    for (size_t s = 1; s < dp_map.size() - 1; ++s) {
        bool center_exists = false;
        for (auto& n : dp_map[s].nodes) {
            // Toleransı biraz yüksek tutuyoruz ki interpolasyon sapmalarını yakalasın
            if (std::abs(n.lat - centers[s].lat) < 0.02f && std::abs(n.lon - centers[s].lon) < 0.02f) {
                center_exists = true;
                break;
            }
        }

        // Eğer merkez nokta karada sanılıp listeye eklenmemişse, biz zorla ekliyoruz
        if (!center_exists) {
            Node n;
            n.lat = centers[s].lat;
            n.lon = centers[s].lon;
            n.total_cost = numeric_limits<float>::infinity();
            n.back_pointer = -1;
            n.temp = 15.0f; // Güvenlik değeri

            // Sıcaklık ataması
            float min_d = 999999.0f;
            for (const auto& w : weather_data) {
                if (w.size() >= 3) {
                    float d = haversine(n.lat, n.lon, w[0], w[1]);
                    if (d < min_d) {
                        min_d = d;
                        n.temp = w[2];
                    }
                }
            }
            dp_map[s].nodes.push_back(n);
        }
    }
    auto start = std::chrono::steady_clock::now();
    int empty_stages = 0;
    int single_stages = 0;
    for (size_t i = 0; i < dp_map.size(); i++) {
        if (dp_map[i].nodes.size() == 0) {
            cout << "BOS STAGE: " << i << " (" << centers[i].lat << ", " << centers[i].lon << ")" << endl;
            empty_stages++;
        }
        else if (dp_map[i].nodes.size() == 1) {
            single_stages++;
        }
    }
    cout << "Bos stage sayisi: " << empty_stages << endl;
    cout << "Tek nodlu stage sayisi: " << single_stages << endl;

    // -------------------------------DP baslangici--------------------------
    cout << "dp baslatiliyor" << endl;

    size_t max_nodes_in_any_stage = 0;
    for (const auto& stage : dp_map) {
        if (stage.nodes.size() > max_nodes_in_any_stage) {
            max_nodes_in_any_stage = stage.nodes.size();
        }
    }

    float ship_speed = 10.0f; // geminin ortalama hızı (10 m/s 19.4 knot) bu çok önemli değil statik değer olması lazım
    float inv_ship_speed = 1.0f / ship_speed;

    for (int i = 0; i < 2; i++)
    {
        dp_map[i].nodes.clear();

        Node node;
        node.lat = original_waypoints[i][0];
        node.lon = original_waypoints[i][1];

        node.total_cost = (i == 0) ? 0.0f : std::numeric_limits<float>::infinity();

        float min_d = std::numeric_limits<float>::max();

        for (const auto& w : weather_data)
        {
            float d = haversine(node.lat, node.lon, w[0], w[1]);
            if (d < min_d)
            {
                min_d = d;
                node.temp = w[2];
            }
        }

        dp_map[i].nodes.push_back(node);
    }

    int last_idx = dp_map.size() - 1;
    dp_map[last_idx].nodes.clear();

    Node final_node;
    final_node.lat = original_waypoints.back()[0];
    final_node.lon = original_waypoints.back()[1];

    float min_dist = 999999.0f;
    float real_temp = 15.0f; // eğer haritada hiç veri yoksa güvenlik amaçlı

    for (const auto& w : weather_data) {
        if (w.size() >= 3) {
            float d = haversine(final_node.lat, final_node.lon, w[0], w[1]);
            if (d < min_dist) {
                min_dist = d;
                real_temp = w[2];
            }
        }
    }
    final_node.temp = real_temp;
    dp_map[last_idx].nodes.push_back(final_node);

    // ---forward pass---
#pragma omp parallel
    {
        // her thread kendi lokal struct'ını ve vektörünü 1 kez oluşturur ki çakışma fln olmasın
        struct Offer {
            float cost;
            int index;
            float heading;
        };

        vector<Offer> offers;
        // vektör kapasitesini kilitler hafıza tahsis etme maliyetini kaldırması şçin
        offers.reserve(max_nodes_in_any_stage);

        // her thread tüm stage'leri sırayla gezcek
        for (size_t s = 0; s < dp_map.size() - 1; ++s) {

            // ekrana yazdırma işlemi tek threadin
#pragma omp single
            {
                if (s % 10 == 0 || s == dp_map.size() - 2) {
                    float progress = (float)s / (dp_map.size() - 2) * 100.0f;
                    printf("\rilerleme: [%%%.2f] - Stage: %zu / %zu", progress, s, dp_map.size() - 2);
                    fflush(stdout);
                }
            }

            Stage& current_stage = dp_map[s];
            Stage& next_stage = dp_map[s + 1];

#pragma omp for schedule(dynamic)
            for (int j = 0; j < (int)next_stage.nodes.size(); ++j) {
                offers.clear();

                Node& nodeB = next_stage.nodes[j];

                for (size_t i = 0; i < current_stage.nodes.size(); ++i) {
                    Node& nodeA = current_stage.nodes[i];
                    if (nodeA.total_cost == numeric_limits<float>::infinity()) continue;

                    // Trigonometrik hesaplamalar
                    float dist = haversine(nodeA.lat, nodeA.lon, nodeB.lat, nodeB.lon);
                    float ship_dir = calculateHeading(nodeA.lat, nodeA.lon, nodeB.lat, nodeB.lon);
                    //----sarmal engelleme (rotanın kendi üstünden tekrar dolanmaması için, yani çember gibi işte)
                    float ref_heading = calculateHeading(centers[s].lat, centers[s].lon, centers[s + 1].lat, centers[s + 1].lon);

                    // Geminin gitmek istediği yön ile rotanın ana yönü arasındaki fark
                    float angle_diff = std::abs(ref_heading - ship_dir);
                    if (angle_diff > 180.0f) angle_diff = 360.0f - angle_diff;

                    // Eğer gemi ana ilerleme yönünden 85 dereceden fazla sapıyorsa (yani geriye doğru akıyorsa)
                    // direk o adımı geç
                    if (angle_diff > 85.0f) {
                        continue;
                    }
                    //----

                    float step_cost = calcCostAF(inv_ship_speed, dist, nodeB.temp, nodeA.incoming_heading, ship_dir);
                    float tentative_cost = nodeA.total_cost + step_cost;

                    offers.push_back({ tentative_cost, (int)i, ship_dir });
                }

                sort(offers.begin(), offers.end(), [](const Offer& a, const Offer& b) {
                    if (a.cost != b.cost) return a.cost < b.cost;
                    return a.index < b.index;
                    });

                float best_total_cost = numeric_limits<float>::infinity();
                int best_parent_index = -1;
                float best_heading = -1;

                for (const auto& offer : offers) {
                    // Eğer önceden bulduğumuz en iyi maliyet, bu offer'ın saf maliyetinden bile küçükse,
                    // listenin geri kalanına bakmamıza gerek yok (çünkü offers cost'a göre sıralı), raporda belirttim
                    if (offer.cost >= best_total_cost) {
                        break;
                    }

                    int a_index = offer.index;
                    Node& best_nodeA = current_stage.nodes[a_index];

                    int path_status = checkPathSafety(best_nodeA.lat, best_nodeA.lon, nodeB.lat, nodeB.lon);
                    if (path_status == -1) {
                        bool is_center_A = (std::abs(best_nodeA.lat - centers[s].lat) < 0.02f && std::abs(best_nodeA.lon - centers[s].lon) < 0.02f);
                        bool is_center_B = (std::abs(nodeB.lat - centers[s + 1].lat) < 0.02f && std::abs(nodeB.lon - centers[s + 1].lon) < 0.02f);

                        if (is_center_A && is_center_B) {
                            path_status = 5; 
                        }
                    }
                  
                    if (path_status != -1) {
                        float buffer_penalty = path_status * 100.0f;
                        float final_cost = offer.cost + buffer_penalty;

                        if (final_cost < best_total_cost) {
                            best_total_cost = final_cost;
                            best_parent_index = a_index;
                            best_heading = offer.heading;
                        }
                    }
                }

                if (best_parent_index != -1) {
                    nodeB.total_cost = best_total_cost;
                    nodeB.back_pointer = best_parent_index;
                    nodeB.incoming_heading = best_heading;
                }
            }
        }
    }
    bool total_failure = false;
    for (size_t s = 0; s < dp_map.size(); ++s) {
        bool stage_alive = false;
        for (const auto& node : dp_map[s].nodes) {
            if (node.total_cost != numeric_limits<float>::infinity()) {
                stage_alive = true;
                break;
            }
        }

        if (!stage_alive) {
            cout << ">>> hata:" << s << "'den itibaren zincir koptu" << endl;
            cout << ">>> kopan koordinat: " << centers[s].lat << ", " << centers[s].lon << endl;
            cout.flush();
            total_failure = true;
            break;
        }
    }

    if (total_failure) {
        return -1;
    }
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "\n gecen sure: " << duration.count() << " ms" << std::endl;

    // --- backward pass ---
    cout << "\nbackward pass baslatiliyor" << endl;

    vector<Point> optimized_route;
    float opt_total_temp = 0.0f, opt_total_dist = 0.0f;
    float orig_total_temp = 0.0f, orig_total_dist = 0.0f;
    int last_stage_idx = dp_map.size() - 1;

    int best_final_node_idx = -1;
    float min_final_cost = numeric_limits<float>::infinity();

    for (size_t j = 0; j < dp_map[last_stage_idx].nodes.size(); ++j) {
        if (dp_map[last_stage_idx].nodes[j].total_cost < min_final_cost) {
            min_final_cost = dp_map[last_stage_idx].nodes[j].total_cost;
            best_final_node_idx = j;
        }
    }

    if (best_final_node_idx == -1) {
        cerr << "bir seyler ters gitti (final node yok)" << endl;
        return -1;
    }

    int current_node_idx = best_final_node_idx;

    for (int s = last_stage_idx; s >= 0; --s) {
        if (current_node_idx < 0 || current_node_idx >= (int)dp_map[s].nodes.size()) {
            cerr << "hata: Stage " << s << " icin gecersiz node index: " << current_node_idx << endl;
            return -1;
        }

        Node& n = dp_map[s].nodes[current_node_idx];

        if (s > 0 && (n.back_pointer < 0 || n.back_pointer >= (int)dp_map[s - 1].nodes.size())) {
            cerr << "hata: Stage " << s << " back_pointer gecersiz: " << n.back_pointer << endl;
            return -1;
        }

        optimized_route.push_back({ n.lat, n.lon });
        opt_total_temp += n.temp;


        if (s > 0) {
            current_node_idx = n.back_pointer;
        }
    }

    reverse(optimized_route.begin(), optimized_route.end());

    vector<float> opt_node_temps;
    opt_node_temps.reserve(optimized_route.size());

    current_node_idx = best_final_node_idx;
    for (int s = last_stage_idx; s >= 0; --s) {
        opt_node_temps.push_back(dp_map[s].nodes[current_node_idx].temp);
        if (s > 0) {
            current_node_idx = dp_map[s].nodes[current_node_idx].back_pointer;
        }
    }
    reverse(opt_node_temps.begin(), opt_node_temps.end());

    opt_total_dist = 0.0f;
    float opt_weighted_temp_sum = 0.0f;

    for (size_t i = 0; i < optimized_route.size() - 1; ++i) {
        float seg_dist = haversine(optimized_route[i].lat, optimized_route[i].lon,
            optimized_route[i + 1].lat, optimized_route[i + 1].lon);

        // İki düğüm arasındaki segmentin ortalama sıcaklığı
        float seg_temp = (opt_node_temps[i] + opt_node_temps[i + 1]) * 0.5f;

        opt_total_dist += seg_dist;
        opt_weighted_temp_sum += (seg_temp * seg_dist);
    }

    // Mesafe ağırlıklı ortalama sıcaklık: Sum(T * d) / Toplam Mesafe
    float opt_avg_temp = (opt_total_dist > 0.0f) ? (opt_weighted_temp_sum / opt_total_dist) : 0.0f;


    WeatherCloud cloud;
    cloud.pts.reserve(weather_data.size());
    for (const auto& row : weather_data) {
        if (row.size() >= 3) cloud.pts.push_back({ row[0], row[1], row[2] });
    }

    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<float, WeatherCloud>,
        WeatherCloud, 2> my_kd_tree_t;

    my_kd_tree_t index(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    index.buildIndex();


    vector<float> orig_node_temps(centers.size(), 15.0f);
    for (size_t s = 0; s < centers.size(); ++s) {
        float query_pt[2] = { (float)centers[s].lat, (float)centers[s].lon };
        size_t ret_index = 0;
        float out_dist_sqr = 0.0f;
        nanoflann::KNNResultSet<float> resultSet(1);
        resultSet.init(&ret_index, &out_dist_sqr);

        index.findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParameters(10));
        orig_node_temps[s] = cloud.pts[ret_index].temp;
    }

    orig_total_dist = 0.0f;
    for (size_t i = 0; i < original_waypoints.size() - 1; ++i) {
        orig_total_dist += haversine(original_waypoints[i][0], original_waypoints[i][1],
            original_waypoints[i + 1][0], original_waypoints[i + 1][1]);
    }


    float orig_temp_dist_sum = 0.0f;
    float orig_weighted_temp_sum = 0.0f;

    for (size_t i = 0; i < centers.size() - 1; ++i) {
        float seg_dist = haversine(centers[i].lat, centers[i].lon,
            centers[i + 1].lat, centers[i + 1].lon);

        float seg_temp = (orig_node_temps[i] + orig_node_temps[i + 1]) * 0.5f;

        orig_temp_dist_sum += seg_dist;
        orig_weighted_temp_sum += (seg_temp * seg_dist);
    }

    float orig_avg_temp = (orig_temp_dist_sum > 0.0f) ? (orig_weighted_temp_sum / orig_temp_dist_sum) : 0.0f;



    cout << "\n--- ROTA OPTIMIZASYON SONUCLARI ---" << endl;
    cout << "Orijinal merkez sayisi     : " << centers.size() << endl;
    cout << "Optimize rota nokta sayisi : " << optimized_route.size() << endl;

    cout << "----------------------------------------" << endl;
    cout << "ORIJINAL ROTA:" << endl;
    cout << "Gercek Toplam Mesafe       : " << orig_total_dist << " km" << endl;
    cout << "Agırlıklı Ort. Sicaklik    : " << orig_avg_temp << " C" << endl;
    cout << "----------------------------------------" << endl;
    cout << "OPTIMIZE EDILMIS ROTA:" << endl;
    cout << "Yeni Toplam Mesafe         : " << opt_total_dist << " km" << endl;
    cout << "Agırlıklı Ort. Sicaklik    : " << opt_avg_temp << " C" << endl;
    cout << "----------------------------------------" << endl;

    float dist_diff_percent = ((opt_total_dist - orig_total_dist) / orig_total_dist) * 100.0f;
    float temp_diff = opt_avg_temp - orig_avg_temp;

    cout << "Mesafe Degisimi            : %" << (dist_diff_percent > 0 ? "+" : "") << dist_diff_percent << endl;
    cout << "Sicaklik Farki             : " << (temp_diff > 0 ? "+" : "") << temp_diff << " C" << endl;
    cout << "----------------------------------------" << endl;


    cout << "Disa aktariliyor..." << endl;

    ofstream json_file("route_statistics.json");
    if (json_file.is_open()) {
        json_file << std::fixed << std::setprecision(6);
        json_file << "{\n";
        json_file << "  \"original_route\": {\n";
        json_file << "    \"center_count\": " << centers.size() << ",\n";
        json_file << "    \"total_distance_km\": " << orig_total_dist << ",\n";
        json_file << "    \"average_temperature_c\": " << orig_avg_temp << "\n";
        json_file << "  },\n";

        json_file << "  \"optimized_route\": {\n";
        json_file << "    \"point_count\": " << optimized_route.size() << ",\n";
        json_file << "    \"total_distance_km\": " << opt_total_dist << ",\n";
        json_file << "    \"average_temperature_c\": " << opt_avg_temp << "\n";
        json_file << "  },\n";

        json_file << "  \"comparison\": {\n";
        json_file << "    \"distance_change_percent\": " << dist_diff_percent << ",\n";
        json_file << "    \"temperature_difference_c\": " << temp_diff << "\n";
        json_file << "  }\n";
        json_file << "}\n";

        json_file.close();
        cout << "Istatistikler route_statistics.json dosyasina kaydedildi." << endl;
    }
    else {
        cerr << "Hata: JSON dosyasi olusturulamadi" << endl;
    }

    ofstream out_file("optimized_route_output.js");
    if (out_file.is_open()) {
        out_file << "var dp_optimized_route = [\n";
        out_file << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < optimized_route.size(); ++i) {
            out_file << "  [" << optimized_route[i].lat << ", " << optimized_route[i].lon << "]";
            if (i < optimized_route.size() - 1) {
                out_file << ",\n";
            }
            else {
                out_file << "\n";
            }
        }
        out_file << "];\n";
        out_file.close();
        cout << "optimized_route_output.js dosyasi olusturuldu." << endl;
    }
    else {
        cerr << "Hata: JS dosyasi olusturulamadi" << endl;
    }

    ofstream orig_js_file("original_wp.js");
    if (orig_js_file.is_open()) {
        orig_js_file << "const original_csv_data = `";
        orig_js_file << std::fixed << std::setprecision(8);

        for (size_t i = 0; i < original_waypoints.size(); ++i) {
            orig_js_file << original_waypoints[i][0] << "," << original_waypoints[i][1];
            if (i < original_waypoints.size() - 1) {
                orig_js_file << "\n";
            }
        }

        orig_js_file << "`;\n";
        orig_js_file.close();
        cout << "original_wp.js dosyasi olusturuldu." << endl;
    }
    else {
        cerr << "Hata: original_wp.js dosyasi olusturulamadi" << endl;
    }
    return 0;
}
