#include <omp.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;

void read_image(vector<vector<int>>& v_image, const char* file_name, int& width, int& height, int& bpp){

    int channel_num = 1;
    uint8_t* rgb_image = stbi_load(file_name, &width, &height, &bpp, channel_num);

    v_image.resize(height);

    #pragma omp parallel for
    for(int i = 0; i < height; i++){

        v_image[i].resize(width);

        #pragma omp parallel for
        for(int j = 0; j < width; j++){
            v_image[i][j] = (int)rgb_image[i * width + j];
        }
    }

    stbi_image_free(rgb_image);
}

void write_image(vector<vector<int>>& v_image, int width, int height){
    
    int channel_num = 1;
    uint8_t* rgb_image = new uint8_t[width * height * channel_num];

    #pragma omp parallel for
    for (int i = 0; i < height; ++i){
        #pragma omp parallel for
        for (int j = 0; j < width; ++j){
            rgb_image[i * width + j] = (uint8_t)v_image[i][j];
        }
    }

    stbi_write_png("image_new.png", width, height, channel_num, (void*)rgb_image, width*channel_num);
}

vector<vector<int>> median_filter(vector<vector<int>>& v_image, const int width, const int height, const int core_size){
    
    vector<vector<int>> new_image(v_image.size());
    copy(v_image.begin(), v_image.end(), new_image.begin());

    int edge_x = floor(core_size/2);
    int edge_y = floor(core_size/2);

    #pragma omp parallel for
    for (int i = edge_x; i < height - edge_x; ++i){
        #pragma omp parallel for
        for (int j = edge_y; j < width - edge_y;  ++j){

            vector<int> arr(edge_x*edge_y);

            for (int di = 0; di < edge_x; ++di){
                for (int dj = 0; dj < edge_y; ++dj){
                    arr[di * edge_y + dj] = v_image[i + di - edge_x][j + dj - edge_y];
                }
            }

            sort(begin(arr), end(arr));
            new_image[i][j] = arr[(edge_x * edge_y)/2];
        }
    }

    return new_image;
}

int main(int argc, char *argv[]) {

    if (argc < 3){
        cout << "Wrong parameter!" << endl;
        exit(-1);
    }

    int width, height, bpp;
    vector<vector<int>> v_image;
    const char* file_name = argv[1];
    const int core_size = atoi(argv[2]);

    cout << "Input file: " << file_name << ", core size = " <<  core_size << ';' << endl;

    //Время считывания
    auto t1 = chrono::high_resolution_clock::now();
    read_image(v_image, file_name, width, height, bpp);
    auto t2 = chrono::high_resolution_clock::now();
    double time = chrono::duration<double>(t2 - t1).count();

    cout << "Image read for --- " << time << " sec ---"<< endl;

    //Время фильтрации
    t1 = chrono::high_resolution_clock::now();
    //Фильтрация
    vector<vector<int>> result = median_filter(v_image, width, height, core_size);
    t2 = chrono::high_resolution_clock::now();
    time = chrono::duration<double>(t2 - t1).count();

    write_image(result, width, height);

    cout << "Image filtered for --- " << time << " sec ---"<< endl;

    return 0;
}