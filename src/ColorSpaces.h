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

#endif // COLORSPACES_H
