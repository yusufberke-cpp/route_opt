import http.server
import json
import urllib.request
import math
import time
import os
import concurrent.futures
import subprocess


PORT = 5000
land_geojson = None

def load_land_data():
    global land_geojson
    if land_geojson is None:
        print("Dünya sınırları verisi indiriliyor (İlk kurulum için)...")
        url = "https://raw.githubusercontent.com/johan/world.geo.json/master/countries.geo.json"
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode('utf-8'))
                
                for feature in data.get('features', []):
                    geom = feature.get('geometry')
                    if not geom: continue
                    minX, maxX, minY, maxY = 180, -180, 90, -90
                    
                    def update_bbox(ring):
                        nonlocal minX, maxX, minY, maxY
                        for pt in ring:
                            if pt[0] < minX: minX = pt[0]
                            if pt[0] > maxX: maxX = pt[0]
                            if pt[1] < minY: minY = pt[1]
                            if pt[1] > maxY: maxY = pt[1]
                    
                    g_type = geom.get('type')
                    coords = geom.get('coordinates', [])
                    if g_type == "Polygon":
                        update_bbox(coords[0])
                    elif g_type == "MultiPolygon":
                        for poly in coords:
                            update_bbox(poly[0])
                            
                    feature['bbox'] = {'minX': minX, 'maxX': maxX, 'minY': minY, 'maxY': maxY}
                    
                land_geojson = data
                print("Dünya sınırları başarıyla yüklendi! Kara/Deniz filtresi hazır.")
        except Exception as e:
            print(f"HATA: Sınır verisi indirilemedi: {e}")

def in_ring(point, ring):
    x, y = point[0], point[1]
    inside = False
    j = len(ring) - 1
    for i in range(len(ring)):
        xi, yi = ring[i][0], ring[i][1]
        xj, yj = ring[j][0], ring[j][1]
        
        intersect = ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        if intersect:
            inside = not inside
        j = i
    return inside

def in_polygon(point, polygon):
    outer_ring = polygon[0]
    if not in_ring(point, outer_ring):
        return False
    for i in range(1, len(polygon)):
        if in_ring(point, polygon[i]):
            return False
    return True

def is_point_on_land(lat, lon):
    if not land_geojson: return False
    
    x, y = lon, lat
    point = [x, y]
    
    for feature in land_geojson.get('features', []):
        bbox = feature.get('bbox')
        if not bbox: continue
        
        if x < bbox['minX'] or x > bbox['maxX'] or y < bbox['minY'] or y > bbox['maxY']:
            continue
            
        geom = feature.get('geometry')
        g_type = geom.get('type')
        coords = geom.get('coordinates', [])
        
        if g_type == "Polygon":
            if in_polygon(point, coords): return True
        elif g_type == "MultiPolygon":
            for poly in coords:
                if in_polygon(point, poly): return True
    return False

def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    dLat = math.radians(lat2 - lat1)
    dLon = math.radians(lon2 - lon1)
    a = math.sin(dLat / 2) ** 2 + \
        math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * \
        math.sin(dLon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

def safe_cos(deg):
    val = math.cos(math.radians(deg))
    return val if abs(val) > 0.0001 else 0.0001

def generate_grid_points(waypoints, interval_km, max_width_km, corridor_count):
    points_map = {}
    global_step_index = 0

    segment_angles = []
    for i in range(len(waypoints) - 1):
        p1 = waypoints[i]
        p2 = waypoints[i + 1]
        dLat = p2[0] - p1[0]
        latMid = (p1[0] + p2[0]) / 2.0
        cosLatMid = safe_cos(latMid)
        dLon = (p2[1] - p1[1]) * cosLatMid
        segment_angles.append(math.atan2(dLat, dLon))

    for i in range(len(waypoints) - 1):
        p1 = waypoints[i]
        p2 = waypoints[i + 1]
        dist = haversine(p1[0], p1[1], p2[0], p2[1])

        stepSize = interval_km / 2.0
        steps = max(1, math.ceil(dist / stepSize))

        dLat = p2[0] - p1[0]
        dLon = p2[1] - p1[1]
        current_angle = segment_angles[i]

        for j in range(steps + 1):
            t = j / steps
            cLat = p1[0] + dLat * t
            cLon = p1[1] + dLon * t

            if j == 0 and i > 0:
                prev_angle = segment_angles[i - 1]
                theta_start = prev_angle - math.pi / 2
                theta_end = current_angle - math.pi / 2

                diff = theta_end - theta_start
                while diff < -math.pi: diff += 2 * math.pi
                while diff > math.pi: diff -= 2 * math.pi

                M = max(1, math.ceil(abs(diff) / 0.15))

                for m in range(M + 1):
                    phi = theta_start + (m / M) * diff
                    dy = math.sin(phi)
                    dx = math.cos(phi)

                    rng = corridor_count
                    for k in range(-rng, rng + 1):
                        dist_offset = (k * (max_width_km / rng)) if rng > 0 else 0
                        offset_lat = cLat + dy * (dist_offset / 111.0)
                        offset_lon = cLon + (dx * (dist_offset / 111.0)) / safe_cos(cLat)

                        pt_type = 'center' if k == 0 else 'side'
                        key = f"{offset_lat:.4f},{offset_lon:.4f}"

                        if key not in points_map:
                            points_map[key] = {'lat': offset_lat, 'lon': offset_lon, 'type': pt_type, 'step': global_step_index}
            elif j > 0 or i == 0:
                phi = current_angle - math.pi / 2
                dy = math.sin(phi)
                dx = math.cos(phi)

                rng = corridor_count
                for k in range(-rng, rng + 1):
                    dist_offset = (k * (max_width_km / rng)) if rng > 0 else 0
                    offset_lat = cLat + dy * (dist_offset / 111.0)
                    offset_lon = cLon + (dx * (dist_offset / 111.0)) / safe_cos(cLat)

                    pt_type = 'center' if k == 0 else 'side'
                    key = f"{offset_lat:.4f},{offset_lon:.4f}"

                    if key not in points_map:
                        points_map[key] = {'lat': offset_lat, 'lon': offset_lon, 'type': pt_type, 'step': global_step_index}
            
            global_step_index += 1

    return list(points_map.values())

def process_points(waypoints, interval_km=20, max_width=150, corridor_count=20, filter_land=False):
    raw_points = generate_grid_points(waypoints, interval_km, max_width, corridor_count)
    filtered = []
    grid = {}
    CELL_SIZE = 0.01

    for p in raw_points:
        if filter_land and p['type'] != 'center' and is_point_on_land(p['lat'], p['lon']):
            continue
        
        cellY = round(p['lat'] / CELL_SIZE)
        cellX = round(p['lon'] / CELL_SIZE)
        
        cosLat = safe_cos(p['lat'])
        maxDegreeLon = 1.0 / (111.0 * max(0.01, cosLat))
        
        limitY = 1
        limitX = max(1, math.ceil(maxDegreeLon / CELL_SIZE))
        
        isTooClose = False
        
        for dy in range(-limitY, limitY + 1):
            for dx in range(-limitX, limitX + 1):
                key = f"{cellX + dx},{cellY + dy}"
                cell_points = grid.get(key)
                if cell_points:
                    for fp in cell_points:
                        if haversine(p['lat'], p['lon'], fp['lat'], fp['lon']) < 1.0:
                            isTooClose = True
                            break
                if isTooClose: break
            if isTooClose: break
            
        if not isTooClose:
            filtered.append(p)
            key = f"{cellX},{cellY}"
            if key not in grid:
                grid[key] = []
            grid[key].append(p)
            
    return filtered

def fetch_weather(filtered_points):
    if not filtered_points:
        print("HATA: Hava durumu çekilecek koordinat bulunamadı (Tümü elenmiş olabilir).")
        return ""
        
    chunk_size = 100
    final_data = []
    api_key = "ypQTP99k2XyDdIgE"
    
    chunks = [filtered_points[i:i + chunk_size] for i in range(0, len(filtered_points), chunk_size)]
    print(f"Toplam {len(filtered_points)} nokta, {len(chunks)} grup halinde çekilecek...")

    if api_key:
        concurrency_limit = 15
        print(f"API anahtarı algılandı. {len(chunks)} grup isteği paralel seriler halinde (eşzamanlı limit: {concurrency_limit}) gönderiliyor...")
        for i in range(0, len(chunks), concurrency_limit):
            batch_chunks = chunks[i:i + concurrency_limit]
            
            def fetch_chunk(args):
                idx_in_all, chunk = args
                lats = ",".join([f"{p['lat']:.4f}" for p in chunk])
                lons = ",".join([f"{p['lon']:.4f}" for p in chunk])
                url = f"https://customer-api.open-meteo.com/v1/forecast?latitude={lats}&longitude={lons}&current=temperature_2m&apikey={api_key}"
                
                req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
                # HTML mantığına uygun olarak: HTTP hatası oluşursa Exception fırlatılır, devam edilmez
                with urllib.request.urlopen(req, timeout=15) as response:
                    data = json.loads(response.read().decode('utf-8'))
                    results = data if isinstance(data, list) else [data]
                    
                    chunk_data = []
                    for idx_res, res in enumerate(results):
                        temp = 15.0
                        if 'current' in res and res['current'].get('temperature_2m') is not None:
                            temp = res['current']['temperature_2m']
                        chunk_data.append(f"{chunk[idx_res]['lat']:.4f}\t{chunk[idx_res]['lon']:.4f}\t{temp}\t{chunk[idx_res]['step']}")
                    return (idx_in_all, chunk_data)

            with concurrent.futures.ThreadPoolExecutor(max_workers=concurrency_limit) as executor:
                batch_args = [(i + j, batch_chunks[j]) for j in range(len(batch_chunks))]
                
                batch_results = list(executor.map(fetch_chunk, batch_args))
                
                batch_results.sort(key=lambda x: x[0])
                for _, chunk_data in batch_results:
                    final_data.extend(chunk_data)
                    
            progress = min(100, round(((i + len(batch_chunks)) / len(chunks)) * 100))
            print(f"İlerleme: %{progress}")
    else:
        print("API anahtarı bulunamadı. Ücretsiz API üzerinden kısıtlamalara uyarak sırayla çekiliyor...")
        for i, chunk in enumerate(chunks):
            lats = ",".join([f"{p['lat']:.4f}" for p in chunk])
            lons = ",".join([f"{p['lon']:.4f}" for p in chunk])
            url = f"https://api.open-meteo.com/v1/forecast?latitude={lats}&longitude={lons}&current=temperature_2m"
            
            print(f"İstek atılıyor: Grup {i + 1} / {len(chunks)}")
            
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        
            with urllib.request.urlopen(req, timeout=15) as response:
                data = json.loads(response.read().decode('utf-8'))
                results = data if isinstance(data, list) else [data]
                
                for idx_res, res in enumerate(results):
                    temp = 15.0
                    if 'current' in res and res['current'].get('temperature_2m') is not None:
                        temp = res['current']['temperature_2m']
                    final_data.append(f"{chunk[idx_res]['lat']:.4f}\t{chunk[idx_res]['lon']:.4f}\t{temp}\t{chunk[idx_res]['step']}")
                    
            time.sleep(1)
            
    return "\n".join(final_data)

class RequestHandler(http.server.BaseHTTPRequestHandler):
    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'POST, GET, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        self.end_headers()

    def do_POST(self):
        content_length = int(self.headers.get('Content-Length', 0))
        post_data = self.rfile.read(content_length)
        
        self.send_response(200)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        
        try:
            data = json.loads(post_data.decode('utf-8'))
            
            # Gelen veriden GeoJSON FeatureCollection rotayı ayıklama
            if data.get("type") != "FeatureCollection" or "features" not in data:
                raise ValueError("Geçersiz JSON formatı. Lütfen bir GeoJSON 'FeatureCollection' gönderin.")
                
            waypoints = []
            for feature in data["features"]:
                geom = feature.get("geometry", {})
                if geom.get("type") == "Point" and "coordinates" in geom:
                    coords = geom["coordinates"]
                    if len(coords) >= 2:
                
                        lon, lat = float(coords[0]), float(coords[1])
                        waypoints.append([lat, lon])
                        
            if not waypoints:
                raise ValueError("GeoJSON verisinde geçerli 'Point' (nokta) geometrisi bulunamadı.")
            
            # original_wp.txt Olarak Kaydetme
            script_dir = os.path.dirname(os.path.abspath(__file__))
            original_wp_path = os.path.join(script_dir, 'original_wp.csv')
            
            with open(original_wp_path, 'w', encoding='utf-8') as f:
                for wp in waypoints:
                    
                    f.write(f"{wp[0]},{wp[1]}\n")
                    
            print(f"Rota alındı: {len(waypoints)} nokta. 'original_wp.txt' olarak kaydedildi. İşleniyor...")
            
            interval_km = 20
            max_width_km = 150
            corridor_count = 20
            filter_land = False
            
            # 1. Koordinat Ağı Hesaplama ve Filtreleme
            filtered_points = process_points(waypoints, interval_km, max_width_km, corridor_count, filter_land)
            print(f"Hesaplandı! Toplam nokta sayısı: {len(filtered_points)}. Hava durumu API'den çekiliyor...")
            
            # 2. Open-Meteo'dan Hava Durumunu Çekme
            final_data_text = fetch_weather(filtered_points)
            
            # 3. veri.txt Olarak Kaydetme
            veri_txt_path = os.path.join(script_dir, 'veri.txt')
            with open(veri_txt_path, 'w', encoding='utf-8') as f:
                f.write(final_data_text)
                
            print("İşlem tamamlandı, 'veri.txt' başarıyla oluşturuldu.")
            
            # 4. Binary (sr.exe, ./sr veya RouteOptimizer.exe) Çalıştırma
            possible_bins = [
                os.path.join(script_dir, "sr.exe"),
                os.path.join(script_dir, "RouteOptimizer.exe"),
                os.path.join(script_dir, "sr"),
                "sr.exe",
                "RouteOptimizer.exe",
                "./sr"
            ]

            optimized_route = []
            executed_bin = None
            for candidate in possible_bins:
                if os.path.exists(candidate):
                    try:
                        print(f"Binary dosyası çalıştırılıyor: {candidate}...")
                        subprocess.run([candidate], input="\n", text=True, cwd=script_dir, timeout=300, check=True)
                        print(f"'{candidate}' başarıyla çalıştırıldı.")
                        executed_bin = candidate
                        break
                    except Exception as bin_e:
                        print(f"Uyarı: '{candidate}' çalıştırılırken hata: {bin_e}")

            opt_path = os.path.join(script_dir, 'optimized_route_output.js')
            if os.path.exists(opt_path):
                try:
                    import re
                    with open(opt_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                       
                        matches = re.findall(r'\[\s*(-?\d+(?:\.\d+)?)\s*,\s*(-?\d+(?:\.\d+)?)\s*\]', content)
                        for m in matches:
                            optimized_route.append([float(m[0]), float(m[1])])
                except Exception as read_e:
                    print(f"optimized_route_output.js okuma hatası: {read_e}")

            statistics = {}
            stats_path = os.path.join(script_dir, 'route_statistics.json')
            if os.path.exists(stats_path):
                try:
                    with open(stats_path, 'r', encoding='utf-8') as f:
                        statistics = json.load(f)
                except Exception as read_e:
                    print(f"route_statistics.json okuma hatası: {read_e}")

            response = {
                "status": "success",
                "message": f"Veri işlendi, ./sr çalıştırıldı ({len(optimized_route)} nokta bulundu).",
                "optimized_route": optimized_route,
                "statistics": statistics,
            }
            self.wfile.write(json.dumps(response).encode('utf-8'))
            
        except Exception as e:
            print(f"HATA OLUŞTU: {e}")
            response = {"status": "error", "message": str(e)}
            self.wfile.write(json.dumps(response).encode('utf-8'))

def run():
    load_land_data()
    server_address = ('', PORT)
    httpd = http.server.HTTPServer(server_address, RequestHandler)
    print(f"Python Endpoint Sunucusu başlatıldı! (Port: {PORT})")
    print("Siteden POST isteklerini bekliyor...")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nSunucu kapatılıyor.")
        httpd.server_close()

if __name__ == '__main__':
    run()

