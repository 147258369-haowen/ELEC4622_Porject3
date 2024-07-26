#include "io_bmp.h"
#include "filter.h"
#include "image_comps.h"
#include <math.h>
#include "motion.h"
static int global_mse = 0;

int get_global_mse() {
    return global_mse;
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
float bilinear_interpolate(int* image, int width, int height, float x, float y, int stride) {
    int x1 = (int)x;
    int y1 = (int)y;
    int x2 = x1 + 1;
    int y2 = y1 + 1;

    if (x2 >= width) x2 = width - 1;
    if (y2 >= height) y2 = height - 1;

    float R1 = (x2 - x) * image[y1 * stride + x1] + (x - x1) * image[y1 * stride + x2];
    float R2 = (x2 - x) * image[y2 * stride + x1] + (x - x1) * image[y2 * stride + x2];
    float P = (y2 - y) * R1 + (y - y1) * R2;

    return P;
}


/*****************************************************************************/

mvector
find_motion(my_image_comp* ref, my_image_comp* tgt,
    int start_row, int start_col, int block_width, int block_height, int S)
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
    float N = (float)block_height * (float)block_width;
    mvector vec, best_vec;//vec就是运动矢量
    int mse, best_mse = 256 * block_width * block_height;
    for (vec.y = -S; vec.y <= S; vec.y +=0.5)
        for (vec.x = -S; vec.x <= S; vec.x +=0.5)
        {

            float ref_row = (float)start_row - vec.y;
            float ref_col = (float)start_col - vec.x;
           // float rvalue = bilinear_interpolate(ref->buf_, ref->width, ref->height, ref_row, ref_col, ref->stride);
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            //int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
           
            int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
            for (mse = 0, r = 0; r < block_height; r++,
                 tp += tgt->stride)//换行 rp += 2 * ref->stride,
                for (c = 0; c < block_width; c++)
                {
                    float  rvalue = bilinear_interpolate(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
                    int diff = (tp[c] - rvalue)* (tp[c] - rvalue);
                   
                    mse += diff;
                }
            mse = (1.0 / N) * mse;
           // printf("count:%d\r\n", count);
            if (mse < best_mse)
            {
                best_mse = mse;
                best_vec = vec;
                
            }
        }
    global_mse += best_mse;
    printf("%d\r\n", best_mse);
    return best_vec;
}


void draw_vector(my_image_comp* tgt, int y_start, int x_start, int y_end, int x_end,int n) {
    int dx = abs(x_end - x_start), dy = abs(y_end - y_start);
    int sx = x_start < x_end ? 1 : -1, sy = y_start < y_end ? 1 : -1;
    int err = dx - dy, e2;
    int offset = n ? 0 : 2;
    while (1) {
        int*pixel = (tgt + offset)->buf_ + y_start * (tgt + offset)->stride + x_start;
        *pixel = 254;
        if (x_start == x_end && y_start == y_end) break;
        e2 = err << 1;
        if (e2 > -dy) { err -= dy; x_start += sx; }
        if (e2 < dx) { err += dx; y_start += sy; }
    }
}
