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
void motion_comp(my_image_comp* ref, my_image_comp* tgt, mvector vec,
    int start_row, int start_col, int block_width, int block_height);
mvector
find_motion(my_image_comp* ref, my_image_comp* tgt,
    int start_row, int start_col, int block_width, int block_height, int S);
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
            int block_width, block_height;
            int height = imageParam.height;
            int width = imageParam.width;
            int S = imageParam.S;
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
                    }
                }
            }
            state = OUTPUT_PICTURE;
        } 
        
            break;
        case OUTPUT_PICTURE:

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
void
motion_comp(my_image_comp* ref, my_image_comp* tgt, mvector vec,
    int start_row, int start_col, int block_width, int block_height)
    /* This function transfers data from the `ref' frame to a block within the
       `tgt' frame, thereby realizing motion compensation.  The motion in
       question has already been found by `find_motion' and is captured by
       the `vec' argument.  The block in the `tgt' frame commences
       at the coordinates given by `start_row' and `start_col' and extends
       for `block_width' columns and `block_height' rows. */
{
    int r, c;
    int ref_row = start_row - vec.y;
    int ref_col = start_col - vec.x;
    int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
    int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
    for (r = block_height; r > 0; r--,
        rp += ref->stride, tp += tgt->stride)
        for (c = 0; c < block_width; c++)
            tp[c] = rp[c];
}



/*****************************************************************************/

mvector
find_motion(my_image_comp* ref, my_image_comp* tgt,
    int start_row, int start_col, int block_width, int block_height,int S)
    /* This function finds the motion vector which best describes the motion
       between the `ref' and `tgt' frames, over a specified block in the
       `tgt' frame.  Specifically, the block in the `tgt' frame commences
       at the coordinates given by `start_row' and `start_col' and extends
       over `block_width' columns and `block_height' rows.  The function finds
       the translational offset (the returned vector) which describes the
       best matching block of the same size in the `ref' frame, where
       the "best match" is interpreted as the one which minimizes the sum of
       absolute differences (SAD) metric. */
{
    mvector vec, best_vec;//vec就是运动矢量
    int sad, best_sad = 256 * block_width * block_height;
    for (vec.y = -S; vec.y <= S; vec.y++)
        for (vec.x = -S; vec.x <= S; vec.x++)
        {
            int ref_row = start_row - vec.y;
            int ref_col = start_col - vec.x;
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
            int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
            for (sad = 0, r = block_height; r > 0; r--,
                rp += ref->stride, tp += tgt->stride)//换行
                for (c = 0; c < block_width; c++)
                {
                    int diff = tp[c] - rp[c];
                    sad += (diff < 0) ? (-diff) : diff;//abs
                }
            if (sad < best_sad)
            {
                best_sad = sad;
                best_vec = vec;
            }
        }

    return best_vec;
}