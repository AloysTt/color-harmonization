#ifndef HARMONIZATION_H
#define HARMONIZATION_H

#include <vector>

// hue of 0 degrees = red

// Distributions, from:
// Daniel Cohen-Or et al. « Color Harmonization ». In : ACM Transactions on Graphics
// (Proceedings of ACM SIGGRAPH) 25.3 (2006), p. 624-630.


struct DistributionColor
{
	float hue;
	float arcWidth;
};

class Distribution
{
public:

	unsigned int getColorCount() const;
	// for every color, we have 1 sector and 2 sector borders
	float getSectorBorder1(unsigned int colorIndex) const;
	float getSectorBorder2(unsigned int colorIndex) const;

	const DistributionColor & getColor(unsigned int index) const;
	float getRotation() const;

	static Distribution create_I_template(float angle, float arcWidth);
	static Distribution create_i_template(float angle, float arcWidth);
private:
	std::vector<DistributionColor> colors;
	float rotationAngle;
};

#endif // HARMONIZATION_H
