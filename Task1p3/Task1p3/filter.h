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
            sinc_buffer[i] = sinc(this->t[i] / (3.5));/// 
        }
        hannwindow(this->hann);
        for (int i = 0; i < (this->length * 2 + 1); i++) {
            sinc_buffer[i] *= this->hann[i];
        }
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