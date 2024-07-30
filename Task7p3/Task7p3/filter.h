#ifndef FILTER_H
#define FILTER_H
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "image_comps.h"
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
    double sinc(double x) {
        return (x == 0) ? 1.0 : sin(PI * x) / (PI * x);
    }

    void hannwindow(double* hann) {
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            hann[i] = 0.5 * (1 - cos(2 * PI * i / (2 * this->length)));
        }
    }
    void SincKernelGenerate(double* sinc_buffer) {
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            this->t[i] = -(this->length) + i;
        }
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            sinc_buffer[i] = sinc(this->t[i] / (2.5));/// 
        }
        hannwindow(this->hann);
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            sinc_buffer[i] *= this->hann[i];
        }
    }
    int sinc_interpolation(int* image, int width, int height ,float x,float y,int stride) {
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
                int* ip = in->buf_ + r * in->stride + c;
                int* op = out->buf_ + r * out->stride + c;
                //mirror_psf = mirror_psf + r * out->stride + c;
                if (!G_MF_flag) {//gaussian
                    sum = 0.0F;
                }
                // for (int y = -filter_extent; y <= filter_extent; y++)//列
                if (inital_flag == 0) {
                    for (int x = -filter_extent; x <= filter_extent; x++)//行
                    {
                        sum += ip[x] * mirror_psf[x];
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
                int* ip = in->buf_ + 2*c * in->stride + 2*r;
                int* op = out->buf_ + c * out->stride + r;

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
    void sinc_2d_filter(my_image_comp* in, my_image_comp* out, double* filter, int filter_width) {
        int filter_extent = (filter_width - 1) / 2;
        double* mirror_psf = filter + filter_extent;  // 中间点

        // 归一化滤波器
        float temp = 0.0f;
        for (int j = -filter_extent; j <= filter_extent; j++) {
            temp += mirror_psf[j];
        }
        if (temp > 1.03f || temp < 0.98f) {
            for (int j = -filter_extent; j <= filter_extent; j++) {
                mirror_psf[j] = (mirror_psf[j] * (1.0f / temp));
            }
        }

        for (int r = 0; r < out->height; r++) {
            for (int c = 0; c < out->width; c++) {
                double sum = 0.0;
                for (int y = -filter_extent; y <= filter_extent; y++) {
                    for (int x = -filter_extent; x <= filter_extent; x++) {
                        int in_r = 2 * r + y;  // 上采样因子为2
                        int in_c = 2 * c + x;  // 上采样因子为2
                        if (in_r < 0) in_r = 0;
                        if (in_c < 0) in_c = 0;
                        if (in_r >= in->height) in_r = in->height - 1;
                        if (in_c >= in->width) in_c = in->width - 1;
                        int* ip = in->buf_ + in_r * in->stride + in_c;
                        sum += (*ip) * mirror_psf[y] * mirror_psf[x];
                    }
                }
                int* op = out->buf_ + r * out->stride + c;
                CLAMP_TO_BYTE(sum);
                *op = (int)sum;
            }
        }
    }


};
#endif