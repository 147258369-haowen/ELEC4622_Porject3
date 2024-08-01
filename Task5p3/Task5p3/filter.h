#ifndef FILTER_H
#define FILTER_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "image_comps.h"
#include <cmath>
#include <vector>
#include <omp.h>
#define PI 3.14159265358979323846
#define SINC 1
class Filter {

public:
    int length;
    int D;
    double* sinc_buffer;
    double* hann;
    double* t;

    Filter(int H, int D) {
        this->length = H;
        this->D = D;
        this->sinc_buffer = new double[2 * H + 1];
        this->hann = new double[2 * H + 1];
        this->t = new double[2 * H + 1];

    }
    ~Filter() {
        delete[] this->sinc_buffer;
        delete[] this->hann;
        delete[] this->t;
    }
    inline double sinc(double x) {
        return (fabs(x) <0.01) ? 1.0 : sin(PI * x) / (PI * x);
    }

    void hannwindow(double* hann) {
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            hann[i] = 0.5 * (1 - cos(2 * PI * i / (2 * this->length)));
        }
    }
    inline float hann_(double x) {
        return  (0.5f * (1.0f + cos(2.0f * PI * x / (2.0f * 7.0f))));
    }
    void SincKernelGenerate(double* sinc_buffer) {
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            this->t[i] = -(this->length) + i;
        }
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            sinc_buffer[i] = sinc(this->t[i] / (3.5));/// 
        }
        hannwindow(this->hann);
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            sinc_buffer[i] *= this->hann[i];
        }
    }
    float sinc_interpolation(int* image, int width, int height, float x, float y, int stride,int printfflag) {
        int x_int = static_cast<int>(x);
        int y_int = static_cast<int>(y);
        float result = 0.0f;
        float norm_factor = 0.0f;
        static float prev_y = 0;
        static int update_x_flag = 0;
        const int window_size = 7;
        const int half_window_size = (window_size - 1) >> 1;

        // Cache sinc values to avoid recomputation
        std::vector<float> sinc_x(2 * half_window_size + 1);
        std::vector<float> sinc_y(2 * half_window_size + 1);
        #pragma omp parallel for
        for (int i = -half_window_size; i <= half_window_size; ++i) {
            sinc_x[i + half_window_size] = sinc(x - (x_int + i));
            sinc_y[i + half_window_size] = sinc(y - (y_int + i));
        }
        //if (printfflag) {
        //    printf("%f,%f\r\n", x, y);
        //}
        #pragma omp parallel for collapse(2) reduction(+:result, norm_factor)
        for (int m = -half_window_size; m <= half_window_size; ++m) {
            for (int n = -half_window_size; n <= half_window_size; ++n) {
                int x_move = x_int + n;
                int y_move = y_int + m;
                if (x_move >= 0 && x_move < width && y_move >= 0 && y_move < height) {
                    float weight = sinc_x[n + half_window_size] * sinc_y[m + half_window_size];
                    result += image[y_move * stride + x_move] * weight;
                    norm_factor += weight;
                }
            }
        }

        if (norm_factor > 0) {
            result /= norm_factor;
        }
        //prev_y = y;
        return result;
    }
    int sinc_interpolation_(int* image, int width, int height, float x, float y, int stride) {
        int x_int = x;
        int y_int = y;
        float result = 0.0;
        float norm_factor = 0.0;
        int window_size = 7;
        int half_window_size = (window_size + 1) >> 1;
        for (int m = -half_window_size; m <= half_window_size; m++) {
            for (int n = -half_window_size; n <= half_window_size; n++) {
                int x_move = x_int + n;
                int y_move = y_int + m;
                if (x_move >= 0 && x_move < width && y_move >= 0 && y_move < height) {
                    float weight = this->sinc(x - x_move) * this->sinc(y - y_move);
                    result += image[y_move * stride + x_move] * weight;
                    norm_factor += weight;
                }
            }
        }
        if (norm_factor > 0) {
            result /= norm_factor;
        }
        return (int)result;
    }
    void horizontal_filter(my_image_comp* in, my_image_comp* out, double* inputfilter, int width, int G_MF_flag) {

        int filter_extent = (width - 1) / 2;
        //double* mirror_buffer = new double[width];
        double* mirror_psf = inputfilter + filter_extent;//中间点
        /* for (int j = -filter_extent; j <= filter_extent; j++) {
             mirror_psf[j] = inputfilter[j + filter_extent];
         }*/
        float temp = 0.0f;
        for (int j = -filter_extent; j <= filter_extent; j++) {
            temp += mirror_psf[j];
        }
        if (temp > 1.03f || temp < 0.98) {
            for (int j = -filter_extent; j <= filter_extent; j++) {
                mirror_psf[j] = (mirror_psf[j] * (1.0f / temp));
            }
        }
        int inital_flag = 0;
        float sum = 0.0F;
        float* buffer = new float[width];
        float* bufptr = buffer;
        for (int r = 0; r < out->height; r++)//进行卷积操作

            for (int c = 0; c < out->width; c++)
            {
                float* ip = in->buf + r * in->stride + c;
                float* op = out->buf + r * out->stride + c;
                //mirror_psf = mirror_psf + r * out->stride + c;
                if (!G_MF_flag) {//gaussian
                    sum = 0.0F;
                }
                // for (int y = -filter_extent; y <= filter_extent; y++)//列
                if (inital_flag == 0) {
                    for (int x = -filter_extent; x <= filter_extent; x++)//行
                    {
                        sum += ip[x] * mirror_psf[x];
                        //if (G_MF_flag) {//MV
                        //    *bufptr = ip[x] * mirror_psf[x];
                        //    listshift(buffer, &bufptr, width);
                        //}

                    }

                }
                else {
                    sum += (ip[filter_extent] * mirror_psf[filter_extent]);
                    sum -= *bufptr;
                    *bufptr = (ip[filter_extent] * mirror_psf[filter_extent]);
                    listshift(buffer, &bufptr, width);
                }
                CLAMP_TO_BYTE(sum);
                inital_flag = G_MF_flag ? 1 : 0;
                *op = sum;
            }
    }
    void vertical_filter(my_image_comp* in, my_image_comp* out, double* inputfilter, int width, int G_MF_flag) {

        int filter_extent = (width - 1) / 2;
        //double* mirror_buffer = new double[width];
        double* mirror_psf = inputfilter + filter_extent;//中间点
        /* for (int j = -filter_extent; j <= filter_extent; j++) {
             mirror_psf[j] = inputfilter[j + filter_extent];
         }*/

        float temp = 0.0f;
        for (int j = -filter_extent; j <= filter_extent; j++) {
            temp += mirror_psf[j];
        }
        if (temp > 1.03f || temp < 0.98) {
            for (int j = -filter_extent; j <= filter_extent; j++) {
                mirror_psf[j] = (mirror_psf[j] * (1.0f / temp));
            }
        }
        float* buffer = new float[width];
        float* bufptr = buffer;
        int inital_flag = 0;
        float sum = 0.0F;
        for (int r = 0; r < out->width; r++)//进行卷积操作
            for (int c = 0; c < out->height; c++)
            {
                float* ip = in->buf + c * in->stride + r;
                float* op = out->buf + c * out->stride + r;

                if (!G_MF_flag) {
                    sum = 0.0F;
                }
                if (!inital_flag) {
                    for (int y = -filter_extent; y <= filter_extent; y++)//
                    {

                        sum += ip[y * in->stride] * mirror_psf[y];
                        //if (G_MF_flag) {
                        //    *bufptr = ip[y * in->stride] * mirror_psf[y];
                        //    listshift(buffer, &bufptr, width);
                        //}

                    }
                }
                else {
                    sum += ip[filter_extent * in->stride] * mirror_psf[filter_extent];
                    sum -= *bufptr;
                    *bufptr = (ip[filter_extent * in->stride] * mirror_psf[filter_extent]);
                    listshift(buffer, &bufptr, width);
                }
                CLAMP_TO_BYTE(sum);
                inital_flag = G_MF_flag ? 1 : 0;
                *op = sum;

            }
    }

};
#endif