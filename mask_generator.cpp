#define NOMINMAX
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <omp.h>
#include <chrono>
#include "gdal_priv.h"
#include "gdal_alg.h"
#include "ogrsf_frmts.h"

using namespace std;

// Harita boyutları ve piksel sayısı
const int WIDTH = 72000;
const int HEIGHT = 36000;
const uint64_t TOTAL_PIXELS = static_cast<uint64_t>(WIDTH) * HEIGHT;
const uint64_t PACKED_SIZE = TOTAL_PIXELS / 8;
const int PADDING_RADIUS = 22; // 0.005 derece için 12 km

// Piksel okuyucu
inline bool getBit(const vector<uint8_t>& mask, int x, int y) {
    if (x < 0 || x >= WIDTH || y < 0 || y >= HEIGHT) return false;
    uint64_t index = static_cast<uint64_t>(y) * WIDTH + x;
    return (mask[index / 8] & (1U << (index % 8))) != 0;
}

int main() {
   
	//shapefile okuma ve rasterizasyon
    GDALAllRegister();
    OGRRegisterAll();

    const char* shp_path = "GSHHS_f_L1.shp";
    GDALDataset* poShapeDS = (GDALDataset*)GDALOpenEx(shp_path, GDAL_OF_VECTOR, NULL, NULL, NULL);
    if (!poShapeDS) {
        cerr << "Shapefile bulunamadi dosya yolunu kontrol et." << endl;
        return 1;
    }

    // RAM üzerinde 72000x36000 boyutunda sanal harita 
    GDALDriver* poMemDriver = GetGDALDriverManager()->GetDriverByName("MEM");
    GDALDataset* poMemDS = poMemDriver->Create("", WIDTH, HEIGHT, 1, GDT_Byte, NULL);

    // Koordinat sistemi [Sol_Ust_X, X_Cozunurluk, 0, Sol_Ust_Y, 0, Y_Cozunurluk]
    double adfGeoTransform[6] = { -180.0, 0.005, 0.0, 90.0, 0.0, -0.005 };
    poMemDS->SetGeoTransform(adfGeoTransform);

    OGRLayer* poLayer = poShapeDS->GetLayer(0);
    int panBandList[1] = { 1 };
    double padfBurnValues[1] = { 1.0 }; // Karalar 1

    CPLErr err = GDALRasterizeLayers(poMemDS, 1, panBandList, 1, (OGRLayerH*)&poLayer, NULL, NULL, padfBurnValues, NULL, NULL, NULL);
    if (err != CE_None) {
        cerr << "Rasterizasyon basarisiz oldu" << endl;
        return 1;
    }

	// 1 bitlik maske oluşturmak için raster veriyi paketleme

    vector<uint8_t> rawMask(PACKED_SIZE, 0);
    GDALRasterBand* poBand = poMemDS->GetRasterBand(1);
    vector<GByte> scanline(WIDTH);

    // RAM haritasını satır satır okuyup bitlere (0-1) çevirme
    for (int y = 0; y < HEIGHT; y++) {
        poBand->RasterIO(GF_Read, 0, y, WIDTH, 1, scanline.data(), WIDTH, 1, GDT_Byte, 0, 0);
        for (int x = 0; x < WIDTH; x++) {
            if (scanline[x] == 1) { // Eğer piksel karaysa
                uint64_t index = static_cast<uint64_t>(y) * WIDTH + x;
                rawMask[index / 8] |= (1U << (index % 8));
            }
        }
    }


    GDALClose(poMemDS);
    GDALClose(poShapeDS);

    // Orijinal haritayı diske kaydetme

    ofstream outFileRaw("land_mask_grid.bin", ios::binary);
    outFileRaw.write(reinterpret_cast<const char*>(rawMask.data()), PACKED_SIZE);
    outFileRaw.close();

	// 12 km genişletme işlemi (dilation)
    vector<uint8_t> dilatedMask = rawMask; // Kopyasını al

#pragma omp parallel for schedule(dynamic, 100)
    for (int y = 0; y < HEIGHT; ++y) {
        if (omp_get_thread_num() == 0 && y % 500 == 0) {
            printf("\r Ilerleme: %d / %d Satir islendi...", y, HEIGHT);
            fflush(stdout);
        }

        for (int x = 0; x < WIDTH; ++x) {
            if (getBit(rawMask, x, y)) { 
                for (int dy = -PADDING_RADIUS; dy <= PADDING_RADIUS; ++dy) {
                    for (int dx = -PADDING_RADIUS; dx <= PADDING_RADIUS; ++dx) {

                        // Dairesel şişirme
                        if (dx * dx + dy * dy <= PADDING_RADIUS * PADDING_RADIUS) {
                            int ny = y + dy;
                            int nx = x + dx;

                            if (ny >= 0 && ny < HEIGHT && nx >= 0 && nx < WIDTH) {
                                uint64_t targetIndex = static_cast<uint64_t>(ny) * WIDTH + nx;
                                uint64_t targetByte = targetIndex / 8;
                                uint8_t targetBit = targetIndex % 8;

#pragma omp atomic
                                dilatedMask[targetByte] |= (1U << targetBit);
                            }
                        }
                    }
                }
            }
        }
    }


    ofstream outFileDilated("dilated_landmask_12km.bin", ios::binary);
    outFileDilated.write(reinterpret_cast<const char*>(dilatedMask.data()), PACKED_SIZE);
    outFileDilated.close();




    return 0;
}