// File: motion.h
// Author: David Taubman
// Last Revised: 30 September, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

#ifndef MOTION_H
#define MOTION_H
struct my_image_comp;
struct mvector {
    mvector() { x = y = 0; }
    float x, y;
};
void
motion_comp(my_image_comp* ref, my_image_comp* tgt, mvector vec,
    int start_row, int start_col, int block_width, int block_height);
mvector
find_motion(my_image_comp* ref, my_image_comp* tgt,
    int start_row, int start_col, int block_width, int block_height, int S);
void draw_vector(my_image_comp* tgt, int y_start, int x_start, int y_end, int x_end, int n);
int get_global_mse();
#endif // MOTION_H
