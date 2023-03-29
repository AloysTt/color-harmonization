#ifndef SEGMENTATION_H
#define SEGMENTATION_H

float getPxStdDev(const float * img, int h, int w, int row, int col);
void computeSimilarityImage(const float * imageYCbCr, float * imageSimilarity, int h, int w);

#endif // SEGMENTATION_H
