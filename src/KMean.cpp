#include <random>
#include <cstring>
#include "KMean.h"

#define EPSILON 1e-8f

inline float squareLength(const ColorRGB & px1, const ColorRGB & px2)
{
	return powf(px1.r - px2.r, 2)
		   + powf(px1.g - px2.g, 2)
		   + powf(px1.b - px2.b, 2);
}

inline bool checkClasses(ColorRGB * colors, ColorRGB * means, int size)
{
	for (int i=0; i<size; ++i)
	{
		float sql = squareLength(colors[i], means[i]);
		if (sql > EPSILON)
		{
			return false;
		}
	}
	return true;
}

void KMean(const unsigned char* imageIn, unsigned int * classes, ColorRGB * colors, int nbColors, int h, int w)
{
	// init colors
	ColorRGB means[nbColors];
	for (int i=0; i<nbColors; ++i)
	{
		ColorHSV hsv{(float)i/((float)nbColors-1.0f), 1.0f, 1.0f};
		ColorRGB rgb{0.0f, 0.0f, 0.0f};
		hsv_to_rgb(hsv, rgb);
		colors[i] = rgb;
		means[i] = rgb;
	}

	// store classes for each pixel
	int size = h*w;
	bool done = false;
	int iterations = 0;
	do
	{
		std::memcpy(colors, means, nbColors*sizeof(ColorRGB)); // old means become the colors

		// update classes
		for (int i=0; i<size; ++i)
		{
			ColorRGB px{imageIn[3*i]/255.0f, imageIn[3*i+1]/255.0f, imageIn[3*i+2]/255.0f};
			float minDist = MAXFLOAT;
			int minDistColor = -1;
			for (int col=0; col<nbColors; ++col)
			{
				float dist = squareLength(colors[col], px);
				if (dist < minDist)
				{
					minDist = dist;
					minDistColor = col;
				}
			}
			classes[i] = minDistColor;
		}

		// compute means
		unsigned int count[nbColors];
		for (int col=0; col<nbColors; ++col)
		{
			means[col] = {0.0f, 0.0f, 0.0f};
			count[col] = 0;
		}
		for (int i = 0; i < size; ++i)
		{
			ColorRGB px{imageIn[3*i]/255.0f, imageIn[3*i+1]/255.0f, imageIn[3*i+2]/255.0f};
			means[classes[i]].r += px.r;
			means[classes[i]].g += px.g;
			means[classes[i]].b += px.b;
			count[classes[i]]++;
		}

		for (int cl=0; cl < nbColors; ++cl)
		{
			means[cl].r /= static_cast<float>(count[cl]);
			means[cl].g /= static_cast<float>(count[cl]);
			means[cl].b /= static_cast<float>(count[cl]);
		}
		done = checkClasses(colors, means, nbColors);
	}
	while (!done);
}