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
    int ref_row = start_row - vec.y_;
    int ref_col = start_col - vec.x_;
    int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
    int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
    for (r = block_height; r > 0; r--,
        rp += ref->stride, tp += tgt->stride)
        for (c = 0; c < block_width; c++)
            tp[c] = ((rp[c] >> 1) + 128);
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
void
motion_comp_float(my_image_comp* ref, my_image_comp* tgt, mvector vec,
    int start_row, int start_col, int block_width, int block_height)
    /* This function transfers data from the `ref' frame to a block within the
       `tgt' frame, thereby realizing motion compensation.  The motion in
       question has already been found by `find_motion' and is captured by
       the `vec' argument.  The block in the `tgt' frame commences
       at the coordinates given by `start_row' and `start_col' and extends
       for `block_width' columns and `block_height' rows. */
{
    int r, c;
    float ref_row = start_row - vec.y;
    float ref_col = start_col - vec.x;
    //Filter sinc(7, 7);
    //int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
    int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
    for (r = 0; r < block_height; r++,
        tp += tgt->stride)
        for (c = 0; c < block_width; c++) {
            float  rvalue = bilinear_interpolate(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
            //float  rvalue = sinc.sinc_interpolation(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
            tp[c] = ((rvalue / 2.0) + 128);
        }
}

mvector Coarse_find_motion(my_image_comp* ref, my_image_comp* tgt, int start_row, int start_col, int block_width, int block_height, int S) {
    float N = (float)block_height * (float)block_width;
    mvector vec, best_vec;//vec就是运动矢量
    float mse, best_mse = 256 * block_width * block_height;
    for (vec.y_ = -S; vec.y_ <= S; vec.y_ += 4)
        for (vec.x_ = -S; vec.x_ <= S; vec.x_ += 4)
        {
            int ref_row = start_row - vec.y_;
            int ref_col = start_col - vec.x_;
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
            int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
            for (mse = 0, r = block_height; r > 0; r--,
                rp += ref->stride, tp += tgt->stride)//换行
                for (c = 0; c < block_width; c++)
                {
                    float diff = (tp[c] - rp[c]) * (tp[c] - rp[c]);
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
    //printf("%d\r\n", best_mse);
    return best_vec;

}
mvector Increment_find_motion(my_image_comp* ref, my_image_comp* tgt, int start_origion_row,
    int star_origion_col, int block_width, int block_height, int S, int vec_y, int vec_x) {
    float N = (float)block_height * (float)block_width;
    mvector vec, best_vec;//vec就是运动矢量
    float mse, best_mse = 256 * block_width * block_height;
    /*printf("vecpre.x_:%d,vecpre.y_:%d\r\n", vecpre.x_, vecpre.y_);*/
    for (vec.y_ = (vec_y - 2); vec.y_ <= (vec_y + 2); vec.y_ += 1)
        for (vec.x_ = (vec_x - 2); vec.x_ <= (vec_x + 2); vec.x_ += 1)
        {
            int ref_row = start_origion_row - vec.y_;
            int ref_col = star_origion_col - vec.x_;
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
            int* tp = tgt->buf_ + start_origion_row * tgt->stride + star_origion_col;
            for (mse = 0, r = block_height; r > 0; r--,
                rp += ref->stride, tp += tgt->stride)//换行
                for (c = 0; c < block_width; c++)
                {
                    float diff = (tp[c] - rp[c]) * (tp[c] - rp[c]);
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
    // printf("%d\r\n", best_mse);
    return best_vec;
}
mvector Half_pixel_find_motion(my_image_comp* ref, my_image_comp* tgt,
    int start_origion_row, int star_origion_col, int block_width, int block_height, int S, int vec_y, int vec_x)

{
    float N = (float)block_height * (float)block_width;
    mvector vec, best_vec;//vec就是运动矢量
    float mse, best_mse = 256 * block_width * block_height;
    //Filter sinc(7,7);
    //printf("vecpre.x_:%d,vecpre.y_:%d\r\n", vecpre.x_, vecpre.y_);
    for (vec.y = ((float)vec_y - 0.5f); vec.y <= ((float)vec_y + 0.5f); vec.y += 0.5f)
        for (vec.x = ((float)vec_x - 0.5f); vec.x <= ((float)vec_x + 0.5f); vec.x += 0.5f)
        {

            float ref_row = (float)start_origion_row - vec.y;
            float ref_col = (float)star_origion_col - vec.x;
            // float rvalue = bilinear_interpolate(ref->buf_, ref->width, ref->height, ref_row, ref_col, ref->stride);
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            //int* rp = ref->buf_ + ref_row * ref->stride + ref_col;

            int* tp = tgt->buf_ + start_origion_row * tgt->stride + star_origion_col;
            for (mse = 0, r = 0; r < block_height; r++,
                tp += tgt->stride)//换行 rp += 2 * ref->stride,
                for (c = 0; c < block_width; c++)
                {
                    float  rvalue = bilinear_interpolate(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
                    //float rvalue = sinc.sinc_interpolation(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
                    float diff = (tp[c] - rvalue) * (tp[c] - rvalue);
                    mse += diff;
                }
            mse = (1.0 / N) * mse;
            //printf("count:%d\r\n", count);
            if (mse < best_mse)
            {
                best_mse = mse;
                best_vec = vec;
            }
        }
    global_mse += best_mse;
    //printf("%d\r\n", best_mse);
    return best_vec;
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
    //Filter sinc(7,7);
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
                    //float rvalue = sinc.sinc_interpolation(ref->buf_, ref->width, ref->height, ref_col + c, ref_row + r, ref->stride);
                    int diff = (tp[c] - rvalue)* (tp[c] - rvalue);
                   
                    mse += diff;
                }
            mse = (1.0 / N) * mse;
            //printf("count:%d\r\n", count);
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


void draw_vector(my_image_comp* tgt, int y_start, int x_start, int y_end, int x_end, int n) {
    int dx = abs(x_end - x_start), dy = abs(y_end - y_start);
    int sx = x_start < x_end ? 1 : -1, sy = y_start < y_end ? 1 : -1;
    int err = dx - dy, e2;
    //int offset = n ? 0 : 2;
    while (1) {
        int* pixel = (tgt)->buf_ + y_start * (tgt)->stride + x_start;
        *pixel = 0;
        if (x_start == x_end && y_start == y_end) break;
        e2 = err << 1;
        if (e2 > -dy) { err -= dy; x_start += sx; }
        if (e2 < dx) { err += dx; y_start += sy; }
    }
}
void Calculate_mse(my_image_comp* tgt, my_image_comp* out) {
    float sum = 0;
    float temp = 0;
    for (int r = 0; r < tgt->height; r++) {
        for (int c = 0; c < tgt->width; c++) {
            int* tg = tgt->buf_ + r * tgt->stride + c;
            int* ref = out->buf_ + r * out->stride + c;
            temp = (*tg - ((*ref - 128) * 2));
            sum += temp * temp;
        }
    }
    float mse = sum / (tgt->height * tgt->width);
    printf("Total MSE : %f\n", mse);
}

void
motion_copy(my_image_comp* ref, my_image_comp* tgt, mvector vec,
    int start_row, int start_col, int block_width, int block_height) {
    int r, c;
    int* rp = ref->buf_ + start_row * ref->stride + start_col;
    int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
    for (r = block_height; r > 0; r--,
        rp += ref->stride, tp += tgt->stride)
        for (c = 0; c < block_width; c++)
            tp[c] = rp[c];

}

mvector
find_motion_origion(my_image_comp* ref, my_image_comp* tgt,
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
    for (vec.y_ = -S; vec.y_ <= S; vec.y_++)
        for (vec.x = -S; vec.x <= S; vec.x++)
        {
            int ref_row = start_row - vec.y_;
            int ref_col = start_col - vec.x_;
            if ((ref_row < 0) || (ref_col < 0) ||
                ((ref_row + block_height) > ref->height) ||
                ((ref_col + block_width) > ref->width))//如果计算出的块超出了参考帧的边界，则跳过这个运动矢量
                continue; // Translated block not containe within reference frame
            int r, c;
            int* rp = ref->buf_ + ref_row * ref->stride + ref_col;
            int* tp = tgt->buf_ + start_row * tgt->stride + start_col;
            for (mse = 0, r = block_height; r > 0; r--,
                rp += ref->stride, tp += tgt->stride)//换行
                for (c = 0; c < block_width; c++)
                {
                    int diff = (tp[c] - rp[c]) * (tp[c] - rp[c]);

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
    //printf("%d\r\n", best_mse);
    return best_vec;
}

void draw_vector_(my_image_comp* img, int start_row, int start_col, int vec_y, int vec_x, int color_plane) {
    int abs_v_2 = abs(vec_y);
    int abs_v_1 = abs(vec_x);
    int c1 = start_col;
    int c2 = start_row;

    if (abs_v_2 > abs_v_1) {
        if (vec_x != 0) {
            float k = (float)vec_y / vec_x;
            if (vec_x >= 0) {
                for (int n1 = c1; n1 <= (c1 + vec_x); n1++) {
                    int n2 = c2 + (n1 - c1) * k;
                    if (n1 >= 0 && n1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
            else {
                for (int n1 = (c1 + vec_x); n1 <= c1; n1++) {
                    int n2 = c2 + (n1 - c1) * k;
                    if (n1 >= 0 && n1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
        }
        else {
            // vec_x == 0, 只绘制垂直线
            if (vec_y >= 0) {
                for (int n2 = c2; n2 <= (c2 + vec_y); n2++) {
                    if (c1 >= 0 && c1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + c1;
                        *pixel = 0;
                    }
                }
            }
            else {
                for (int n2 = (c2 + vec_y); n2 <= c2; n2++) {
                    if (c1 >= 0 && c1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + c1;
                        *pixel = 0;
                    }
                }
            }
        }
    }
    else {
        if (vec_y != 0) {
            float k = (float)vec_x / vec_y;
            if (vec_y >= 0) {
                for (int n2 = c2; n2 <= (c2 + vec_y); n2++) {
                    int n1 = c1 + (n2 - c2) * k;
                    if (n1 >= 0 && n1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
            else {
                for (int n2 = (c2 + vec_y); n2 <= c2; n2++) {
                    int n1 = c1 + (n2 - c2) * k;
                    if (n1 >= 0 && n1 < img->width && n2 >= 0 && n2 < img->height) {
                        int* pixel = img->buf_ + n2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
        }
        else {
            // vec_y == 0, 只绘制水平线
            if (vec_x >= 0) {
                for (int n1 = c1; n1 <= (c1 + vec_x); n1++) {
                    if (n1 >= 0 && n1 < img->width && c2 >= 0 && c2 < img->height) {
                        int* pixel = img->buf_ + c2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
            else {
                for (int n1 = (c1 + vec_x); n1 <= c1; n1++) {
                    if (n1 >= 0 && n1 < img->width && c2 >= 0 && c2 < img->height) {
                        int* pixel = img->buf_ + c2 * img->stride + n1;
                        *pixel = 0;
                    }
                }
            }
        }
    }
}