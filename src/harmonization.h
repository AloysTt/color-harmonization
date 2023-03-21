#ifndef HARMONIZATION_H
#define HARMONIZATION_H

#include <vector>

// hue of 0 degrees = red

// Distributions, from:
// Daniel Cohen-Or et al. « Color Harmonization ». In : ACM Transactions on Graphics
// (Proceedings of ACM SIGGRAPH) 25.3 (2006), p. 624-630.

enum class HarmonyType
{
	COMPLEMENTARY = 0,
	TRIADIC,
	TETRADIC_RECTANGLE,
	TETRADIC_SQUARE,
	SPLIT_COMPLEMENTARY,
	ANALOGOUS,
};

struct DistributionColor
{
	float hue;
	float arcWidth;
};

class Distribution
{
public:
	Distribution();

	unsigned int getColorCount() const;
	// for every color, we have 1 sector and 2 sector borders
	float getSectorBorder1(unsigned int colorIndex) const;
	float getSectorBorder2(unsigned int colorIndex) const;

	const DistributionColor & getColor(unsigned int index) const;
	float getRotation() const;

	void addColor(float hue, float arcWidth);

	static Distribution create_I_template(float angle, float arcWidth);
	static Distribution create_i_template(float angle, float arcWidth);
private:
	std::vector<DistributionColor> colors;
	float rotationAngle;
};

// shifting
void KMeanShift(HarmonyType type, int h, int w, const unsigned char *image, unsigned char *imageOut);
void shift(const Distribution & distrib, int h, int w, const unsigned char *image, unsigned char *imageOut);
void
find_sector_borders(const Distribution & distrib, const unsigned char *image, unsigned char *imageSectorBorders, int h,
					int w);

#endif // HARMONIZATION_H
