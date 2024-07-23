#include "io_bmp.h"
#include "filter.h"
#include "image_comps.h"
#include <math.h>
#include "motion.h"
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
                    int diff = (tp[c] - rp[c])* (tp[c] - rp[c]);
                   
                    sad += diff;
                }
            sad = (1.0 / N) * sad;
           // printf("count:%d\r\n", count);
            if (sad < best_sad)
            {
                best_sad = sad;
                best_vec = vec;
            }
        }

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
