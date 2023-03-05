#include "HSV.h"

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
