﻿// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "filter.h"
#include "image_comps.h"
#include <math.h>
#include "motion.h"
#include <time.h>
extern int height_offset;
extern float laplacianKernel[3][3];

int main(int argc, char* argv[]) {
    STATE state = CHECK_INPUT;
    int filterChoose = 0;
    bmp_in in;
    bmp_in in2;
    my_image_comp mono[2];
    bmp_out out;
    my_image_comp* input_comps = nullptr;
    my_image_comp* input_comps2 = nullptr;
    my_image_comp* output_comps = nullptr;
    io_byte* line = nullptr;
    float** matrix = nullptr;
    my_image_comp* temp_comps = nullptr;
    my_image_comp* temp_comps_out = nullptr;
    my_image_comp* overlaid = nullptr;
    my_image_comp* nextinput = nullptr;
    my_image_comp** lapalacin = nullptr;
    my_image_comp** lapalacin2 = nullptr;
    my_image_comp** lapalacin3 = nullptr;
    my_image_comp* temp_image_2 = nullptr;
    my_image_comp* store_out_image_1 = new my_image_comp[3];
    my_image_comp* store_out_image_2 = new my_image_comp[3];
    ImageParam imageParam;
    float sigma = 0;
    float** lapkernel = allocateMatrix(3);
    Filter* sinc = nullptr;
    int newheight = 0;
    int tempheight = 0;
    bool change_picture = false;
    int first_width,first_height,first_initHeight;
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    while (true) {
        switch (state) {
        case CHECK_INPUT:
            CheckInput(argc, argv, &sigma, &filterChoose, &imageParam);
            state = LOAD_PICTURE;
            break;
        case LOAD_PICTURE:
            //printf("%f\r\n", FilterInit(matrix, HIGHT, WIDTH));//归一化
            LoadImage(&in,&in2, &input_comps, &input_comps2 ,&output_comps, &line, &imageParam, &filterChoose, argv);
            state = MOTION_ESTIMATION_COMPENSATATION;
            break;
        case MOTION_ESTIMATION_COMPENSATATION: {
            // Now perform simple motion estimation and compensation
            int nominal_block_width = imageParam.B;
            int nominal_block_height = imageParam.B;
            int block_width, block_height;
            int height = imageParam.height;
            int width = imageParam.width;
            int S = imageParam.S;
            int mid;
            if (imageParam.B % 2 == 0) { mid = imageParam.B >> 1; }
            else { mid = (imageParam.B + 1) >> 1; }
            my_image_comp* input_upsample = new  my_image_comp[imageParam.num_comp];
          /*  Image_copy_no_offset(input_comps2, output_comps, &imageParam);
            Image_copy_no_offset(input_comps2, output_comps+1, &imageParam);
            Image_copy_no_offset(input_comps2, output_comps + 2, &imageParam);*/
            //Image_upsample(&input_comps, &input_upsample, &imageParam);
            for (int n = 0; n < imageParam.num_comp; n++) {
                for (int r = 0; r < height; r += block_height)//height is the image hight
                {
                    block_height = nominal_block_height;
                    if ((r + block_height) > height)
                        block_height = height - r;
                    for (int c = 0; c < width; c += block_width)//width is th image width
                    {
                        block_width = nominal_block_width;
                        if ((c + block_width) > width)
                            block_width = width - c;
                        mvector vec = Coarse_find_motion(input_comps +n, input_comps2 + n,//input_comps reference frame
                        r, c, block_width, block_height, S);
                        //printf("vec.x_:%d,vec.y_:%d\r\n", vec.x_, vec.y_);
                        mvector vec_pixel = Increment_find_motion(input_comps + n, input_comps2 + n,//input_comps reference frame
                        r,c, block_width, block_height, S, vec.y_, vec.x_);
                        //printf("vec_pixel.x:%d,vec_pixel.y:%d\r\n", vec_pixel.x_, vec_pixel.y_);
                        mvector vec_half_pixel = Half_pixel_find_motion(input_comps + n, input_comps2 + n,//input_comps reference frame
                        r,c, block_width, block_height, S, vec_pixel.y_, vec_pixel.x_);

                        mvector vec_total;
                        //vec_total.x = (float)vec.x_ + (float)vec_pixel.x_ + vec_half_pixel.x;
                        //vec_total.y = (float)vec.y_ + (float)vec_pixel.y_ + vec_half_pixel.y;
                        vec_total.x = vec_half_pixel.x;
                        vec_total.y = vec_half_pixel.y;
                        //printf("vec_half_pixel.x:%f,vec_half_pixel.y:%f\r\n", vec_half_pixel.x, vec_half_pixel.y);
                        motion_comp_float(input_comps + n, output_comps + n, vec_total,
                            r, c, block_width, block_height);
                        if (imageParam.num_comp == 1) {
                            motion_copy(output_comps + n, output_comps + 1, vec_total,
                                r, c, block_width, block_height);
                            motion_copy(output_comps + n, output_comps + 2, vec_total,
                                r, c, block_width, block_height);
                        }
                        int y_start = r +mid;
                        int x_start = c +mid;
                        int x_end = x_start - vec_total.x;
                        int y_end = y_start - vec_total.y;
                        /*int x_end = c - (vec.x + vec_pixel.x + vec_half_pixel.x);
                        int y_end = r - (vec.y + vec_pixel.y + vec_half_pixel.y);*/
                        draw_vector_(&output_comps[1], y_start, x_start, vec_total.y, vec_total.x, n);
                        //draw_vector(&output_comps[1], y_start, x_start, y_end, x_end, n);
                        
                    }
                }
            }
            //printf("total mse %d\r\n", get_global_mse());
            Calculate_mse(&input_comps2[0], &output_comps[0]);
            end = clock();
            cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
            printf("Run time: %f second\n", cpu_time_used);
            state = OUTPUT_PICTURE;
        }   
        
            break;
        case OUTPUT_PICTURE:
            delete line;
            line = new  io_byte[output_comps[0].width * 3];
            OutputImage(&out, input_comps, &output_comps, &line, &imageParam, argv);
            bmp_in__close(&in);
            delete[] line;
            delete[] input_comps;
            delete[] output_comps;

            return 1;
        default:
            break;
        }
    }
}
