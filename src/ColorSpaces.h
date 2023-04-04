#ifndef COLORSPACES_H
#define COLORSPACES_H

struct ColorRGB
{
	union
	{
		float val[3];
		struct
		{
			float r, g, b;
		};
	};
};

struct ColorHSV
{
	union
	{
		float val[3];
		struct
		{
			float h, s, v;
		};
	};
};

struct ColorYCbCr
{
	union
	{
		float val[3];
		struct
		{
			float y, cb, cr;
		};
	};
};

void rgb_to_hsv(const ColorRGB & rgb, ColorHSV & hsv);
void hsv_to_rgb(const ColorHSV & hsv, ColorRGB & rgb);

void rgb_to_ycbcr(const ColorRGB & rgb, ColorYCbCr & ycbcr);
void rgb_to_ycbcr(const unsigned char * image, float * imageYCbCr, int h, int w);
void ycbcr_to_rgb(const ColorYCbCr & ycbcr, ColorRGB & rgb);

float ycbcr_distance_relative_euclidean(float y1, float cb1, float cr1, float y2, float cb2, float cr2);
double ycbcr_distance_relative_euclidean_squared(float y1, float cb1, float cr1, float y2, float cb2, float cr2);
inline float ycbcr_distance_relative_euclidean(const float *px1, const float *px2)
{
	return ycbcr_distance_relative_euclidean(px1[0], px1[1], px1[2], px2[0], px2[1], px2[2]);
}
inline double ycbcr_distance_relative_euclidean_squared(const float *px1, const float *px2)
{
	return ycbcr_distance_relative_euclidean_squared(px1[0], px1[1], px1[2], px2[0], px2[1], px2[2]);
}

#endif // COLORSPACES_H
