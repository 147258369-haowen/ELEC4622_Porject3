// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "filter.h"
#include "image_comps.h"
#include <math.h>
#include "motion.h"
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
            int mid;
            if (imageParam.B % 2 == 0) { mid = imageParam.B >> 2; } else { mid = (imageParam.B + 1) >> 2; }
            int block_width, block_height;
            int height = imageParam.height;
            int width = imageParam.width;
            int S = imageParam.S;
            Image_copy_no_offset(input_comps, &output_comps[0], &imageParam);
            Image_copy_no_offset(input_comps, &(output_comps[1]), &imageParam);
            Image_copy_no_offset(input_comps, &(output_comps[2]), &imageParam);
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
                        mvector vec = find_motion(input_comps+n, input_comps2 + n,
                            r, c, block_width, block_height, S);
                        motion_comp(input_comps + n, output_comps + n, vec,
                            r, c, block_width, block_height);

                        int y_start = r + mid;
                        int x_start = c + mid;
                        int x_end = x_start - vec.x;
                        int y_end = y_start - vec.y;
                        
                        //draw_vector_(&output_comps[1], y_start, x_start, vec.y, vec.x, n);
                        draw_vector(&output_comps[1], y_start, x_start, y_end, x_end,n);
                        //draw_vectors2(&output_comps[1], vec, y_start, x_start, nominal_block_width, nominal_block_width);
                    }
                }
            }
           
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
