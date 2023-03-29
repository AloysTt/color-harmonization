#include "segmentation.h"
#include "util.h"
#include "ColorSpaces.h"
#include <cmath>
#include <cstring>
#include <algorithm>

float getPxStdDev(const float * img, int h, int w, int row, int col)
{
	int w3 = 3*w;
	float avgY, avgCb, avgCr;
	avgY = avgCb = avgCr = 0.0f;
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			avgY += img[(row + i) * w3 + 3*(col + j)];
			avgCb += img[(row + i) * w3 + 3*(col + j) + 1];
			avgCr += img[(row + i) * w3 + 3*(col + j) + 2];
		}
	}
	avgY /= 9.0f;
	avgCb /= 9.0f;
	avgCr /= 9.0f;

	float stdY, stdCb, stdCr;
	stdY = stdCb = stdCr = 0.0f;
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			stdY += (img[(row + i) * w3 + 3*(col + j)] - avgY) * (img[(row + i) * w3 + 3*(col + j)] - avgY);
			stdCb += (img[(row + i) * w3 + 3*(col + j) + 1] - avgCb) * (img[(row + i) * w3 + 3*(col + j) + 1] - avgCb);
			stdCr += (img[(row + i) * w3 + 3*(col + j) + 2] - avgCr) * (img[(row + i) * w3 + 3*(col + j) + 2] - avgCr);
		}
	}
	stdY /= 9.0f;
	stdCb /= 9.0f;
	stdCr /= 9.0f;
	stdY = std::sqrt(stdY);
	stdCb = std::sqrt(stdCb);
	stdCr = std::sqrt(stdCr);

	return stdY + stdCb + stdCr;
}

void computeSimilarityImage(const float * imageYCbCr, float * imageSimilarity, int h, int w)
{
	int size = h*w;
	std::memset(imageSimilarity, 0, size * sizeof(float));
	for (int row=1; row<h-1; ++row)
	{
		for (int col=1; col<w-1; ++col)
		{
			imageSimilarity[row * w + col] = getPxStdDev(imageYCbCr, h, w, row, col);
		}
	}
	float maxStdDev = *std::max_element(imageSimilarity, imageSimilarity + size);
	for (int i=0; i<size; ++i)
		imageSimilarity[i] = 1.0f - imageSimilarity[i] / maxStdDev;
}

double otsu(int k, double *ddp, int grayLevelCount, OtsuCriterion criterion)
{
	double w0, u;
	w0 = u = 0.0;
	for(int i=0; i<=k; ++i)
	{
		w0+=ddp[i];
		u+=i*ddp[i];
	}
	double w1 = 1.0-w0;

	double ut = u;
	for (int i=k+1; i < grayLevelCount; i++)
		ut += i * ddp[i];

	double u0 = u/w0;
	double u1 = (ut - u) / w1;

	double v0 = 0;
	double v1 = 0;
	for(int i=0; i<=k; i++)
		v0 += ((double)i-u0)*((double)i-u0)
			  * (ddp[i] / w0);
	for(int i=k+1; i < grayLevelCount; i++)
		v1 += ((double)i-u1)*((double)i-u1)
			  * (ddp[i] / w1);

	switch (criterion)
	{
		case OtsuCriterion::BCW_WCV:
		{
			double WCV = w0 * v0 + w1 * v1;
			double BCV = w0 * w1 * (u1 - u0) * (u1 - u0);
			return BCV / WCV;
		}
		case OtsuCriterion::TV_WCV:
		{
			double TV = 0;
			for (int i=0; i<grayLevelCount; ++i)
				TV += ((double)i-ut)*((double)i-ut)*ddp[i];
			double WCV = w0 * v0 + w1 * v1;
			return TV/WCV;
		}
		case OtsuCriterion::BCV_TV:
		{
			double TV = 0;
			for (int i=0; i<grayLevelCount; ++i)
				TV += ((double)i-ut)*((double)i-ut)*ddp[i];
			double BCV = w0*w1*(u1-u0)*(u1-u0);
			return BCV/TV;
		}
		default:
			return 0;
	}
}

void computeSeedsFromSimilarity(const float *imageYCbCr, const float *imageSimilarity, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance)
{
	int size = h*w;
	// split similarity image into buckets
	uint * imageOtsu = new uint[size];
	for (int i=0; i<size; ++i)
		imageOtsu[i] = static_cast<float>(otsuBuckets-1)*imageSimilarity[i];

	// compute ddp
	double * ddp = new double[otsuBuckets];
	create_ddp(imageOtsu, h, w, ddp, otsuBuckets);

	// get threshold using Otsu's method
	int threshold = 0;
	double thresholdVal = 0;
	for (int i=0; i<otsuBuckets; ++i)
	{
		double otsuVal = otsu(i, ddp, otsuBuckets, OtsuCriterion::BCW_WCV);
		if (otsuVal > thresholdVal)
		{
			threshold = i;
			thresholdVal = otsuVal;
		}
	}

	// create binary seed image
	unsigned char * imageSeedOtsu = new unsigned char[size];
	for (int i=0; i<size; ++i)
		imageSeedOtsu[i] = imageOtsu[i] >= threshold ? 255 : 0;

	// 2nd condition : euclidean distance
	filterOtsuSeedsEuclidean(imageYCbCr, imageSeedOtsu, imageSeeds, h, w, maxEuclideanDistance);

	delete [] imageSeedOtsu;
	delete [] imageOtsu;
	delete [] ddp;
}

void filterOtsuSeedsEuclidean(const float *imageYCbCr, const unsigned char *imageSeedOtsu, unsigned char * out, int h, int w, float maxEuclideanDistance)
{
	for (int row=1; row<h-1; ++row)
	{
		for (int col=1; col<w-1; ++col)
		{
			if (imageSeedOtsu[row*w+col] == 0)
			{
				out[row*w+col] = 0;
				continue;
			}
			float max = 0.0f;
			for (int i=-1; i<=1; ++i)
			{
				for (int j=-1; j<=1; ++j)
				{
					float dist = ycbcr_distance_relative_euclidean(
						imageYCbCr[row*w*3+col*3],
						imageYCbCr[row*w*3+col*3+1],
						imageYCbCr[row*w*3+col*3+2],
						imageYCbCr[(row+i)*w*3+(col+j)*3],
						imageYCbCr[(row+i)*w*3+(col+j)*3+1],
						imageYCbCr[(row+i)*w*3+(col+j)*3+2]
					);
					if (dist > max)
						max = dist;
				}
			}
			if (max > maxEuclideanDistance)
				out[row*w+col] = 0;
			else
				out[row*w+col] = 255;
		}
	}
}

void computeSeeds(const float *imageYCbCr, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance)
{
	int size = h*w;
	float * imageSimilarity = new float[size];
	computeSimilarityImage(imageYCbCr, imageSimilarity, h, w);

	computeSeedsFromSimilarity(imageYCbCr, imageSimilarity, imageSeeds, h, w, otsuBuckets, maxEuclideanDistance);

	delete [] imageSimilarity;
}