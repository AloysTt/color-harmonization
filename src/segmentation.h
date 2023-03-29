#ifndef SEGMENTATION_H
#define SEGMENTATION_H

enum class OtsuCriterion
{
	BCW_WCV,
	TV_WCV,
	BCV_TV
};

float getPxStdDev(const float * img, int h, int w, int row, int col);
void computeSimilarityImage(const float * imageYCbCr, float * imageSimilarity, int h, int w);
double otsu(int k, double *ddp, int grayLevelCount, OtsuCriterion criterion);

#endif // SEGMENTATION_H
