#include "ColorSpaces.h"

#include <cmath>

#define PI_6 0.5235987755982988730771072305465838140328615665625176368291574320


void rgb_to_hsv(const ColorRGB & rgb, ColorHSV & hsv)
{
	float X_max = rgb.r;
	float X_min = rgb.r;
	size_t maxComponent = 0;
	size_t minComponent = 0;
	for (int c=1; c<3; ++c)
	{
		if (rgb.val[c] > X_max)
		{
			X_max = rgb.val[c];
			maxComponent = c;
		}
		if (rgb.val[c] < X_min)
		{
			X_min = rgb.val[c];
			minComponent = c;
		}
	}

	float chroma = X_max - X_min;

	hsv.v = (X_max+X_min+chroma)/2.0f;

	if (chroma == 0)
		hsv.h = 0.0f;
	else
	{
		switch (maxComponent)
		{
			case 0:
				hsv.h = 60.0f*std::fmod((rgb.g - rgb.b)/chroma, 6.0);
				break;
			case 1:
				hsv.h = 60.0f*(((rgb.b - rgb.r)/chroma)+2.0f);
				break;
			case 2:
				hsv.h = 60.0f*(((rgb.r - rgb.g)/chroma)+4.0f);
				break;
			default:
				hsv.h = 0.0f;
				break;
		}
	}

	if (hsv.v == 0.0f)
		hsv.s = 0.0f;
	else
		hsv.s = chroma/hsv.v;
}

void hsv_to_rgb(const ColorHSV & hsv, ColorRGB & rgb)
{
	float chroma = hsv.v * hsv.s;
	auto f = [&hsv](float n){
		float k = std::fmod(n+(hsv.h/60.0f), 6.0f);
		return hsv.v - hsv.v*hsv.s*std::max(0.0f, std::min(std::min(k, 4.0f-k), 1.0f));
	};
	rgb.r = f(5.0f);
	rgb.g = f(3.0f);
	rgb.b = f(1.0f);
}


void rgb_to_ycbcr(const ColorRGB & rgb, ColorYCbCr & ycbcr)
{
	ycbcr.y = 65.481*rgb.r + 128.553*rgb.g + 24.966*rgb.b + 16;
	ycbcr.cb = -39.797*rgb.r -74.203*rgb.g + 112.0*rgb.b + 128;
	ycbcr.cr = 112*rgb.r -93.786*rgb.g - 18.214 * rgb.b + 128;
}

void rgb_to_ycbcr(const unsigned char * image, float * imageYCbCr, int h, int w)
{
	int size = h*w;
	int size3 = size*3;
	for (int i=0; i<size; ++i)
	{
		ColorYCbCr * ycbcr = new (imageYCbCr+3*i) ColorYCbCr;
		ColorRGB rgb{image[3*i]/255.0f, image[3*i+1]/255.0f, image[3*i+2]/255.0f};
		rgb_to_ycbcr(rgb, *ycbcr);
	}
}

void ycbcr_to_rgb(const ColorYCbCr & ycbcr, ColorRGB & rgb)
{
	float y = ycbcr.y - 16;
	float cb = ycbcr.cb - 128;
	float cr = ycbcr.cr - 128;

	rgb.r = 0.00456621 * y + 0.00625893*cr;
	rgb.g = 0.00456621 * y -0.00153632*cb -0.00318811*cr;
	rgb.b = 0.00456621 * y + 0.00791071*cb;
}

float ycbcr_distance_relative_euclidean(float y1, float cb1, float cr1, float y2, float cb2, float cr2)
{
	return std::sqrt(
			(y1-y2)*(y1-y2)
			+ (cb1-cb2)*(cb1-cb2)
			+ (cr1-cr2)*(cr1-cr2))
		/ std::sqrt(y1*y1 + cb1*cb1 + cr1*cr1);
}

double ycbcr_distance_relative_euclidean_squared(float y1, float cb1, float cr1, float y2, float cb2, float cr2)
{
	return (
		((double)y1-(double)y2)*((double)y1-(double)y2)
		+ ((double)cb1-(double)cb2)*((double)cb1-(double)cb2)
		+ ((double)cr1-(double)cr2)*((double)cr1-(double)cr2))
		   / ((double)y1*(double)y1 + (double)cb1*(double)cb1 + (double)cr1*(double)cr1);
}