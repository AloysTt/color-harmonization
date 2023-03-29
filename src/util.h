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

inline void histo(const uint *ImgIn, uint* T, int h, int w)
{
	int size = h * w;
	for (int i=0; i < size; i++)
		++T[ImgIn[i]];
}

#endif // UTIL_H
