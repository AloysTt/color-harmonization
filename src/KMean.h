#ifndef KMEAN_H
#define KMEAN_H

#include "HSV.h"

//
// indexImage and colors have to be initialized
void KMean(const unsigned char* imageIn, unsigned int * classes, ColorRGB * colors, int nbColors, int h, int w);

#endif // KMEAN_H
