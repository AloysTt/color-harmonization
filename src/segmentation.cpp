#include "segmentation.h"
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