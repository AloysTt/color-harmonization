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
