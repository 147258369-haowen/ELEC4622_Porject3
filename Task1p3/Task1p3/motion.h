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
    int x, y;
};

#endif // MOTION_H
