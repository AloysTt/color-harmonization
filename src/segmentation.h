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
void computeSeedsFromSimilarity(const float *imageYCbCr, const float *imageSimilarity, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance);
void filterOtsuSeedsEuclidean(const float *imageYCbCr, const unsigned char *imageSeedOtsu, unsigned char * out, int h, int w, float maxEuclideanDistance);
void computeSeeds(const float *imageYCbCr, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance);

#endif // SEGMENTATION_H
