#ifndef UTIL_H
#define UTIL_H

#include <cmath>

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923

inline float gaussian(float x, float avg, float sd)
{
	return (1.0f / (sd * std::sqrt(2.0f * PI)))
		   * std::exp(
		-0.5*pow((x-avg)/sd, 2.0f)
	);
}

inline float rad(float angle)
{
	return angle*PI/180.0f;
}

// tab contents must be initialized to 0
inline void create_histo(const uint *image, int h, int w, uint *tab)
{
	int size = h * w;
	for (int i=0; i < size; i++)
		++tab[image[i]];
}

inline void create_ddp(const uint *image, int h, int w, double *tab, int tabLength)
{
	uint * histogram = new uint[tabLength];
	std::memset(histogram, 0u, tabLength*sizeof(uint));
	create_histo(image, h, w, histogram);

	int size = h*w;
	for (int i=0; i<tabLength; ++i)
		tab[i] = (double)histogram[i]/(double)size;

	delete [] histogram;
}

#endif // UTIL_H
