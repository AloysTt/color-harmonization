#ifndef HSV_H
#define HSV_H

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

void rgb_to_hsv(const ColorRGB & rgb, ColorHSV & hsv);
void hsv_to_rgb(const ColorHSV & hsv, ColorRGB & rgb);

#endif // HSV_H
