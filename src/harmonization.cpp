#include "harmonization.h"

#include <cmath>


Distribution Distribution::create_I_template(float angle, float arcWidth) {
	Distribution d;
	d.colors.push_back(DistributionColor{0, arcWidth});
	d.colors.push_back(DistributionColor{180.0f, arcWidth});
	d.rotationAngle = std::fmod(angle, 360.0f);
	return d;
}

unsigned int Distribution::getColorCount() const {
	return colors.size();
}

float Distribution::getSectorBorder1(unsigned int colorIndex) const {
	return std::fmod(rotationAngle+colors[colorIndex].hue-colors[colorIndex].arcWidth/2.0f, 360.0f);
}

float Distribution::getSectorBorder2(unsigned int colorIndex) const {
	return std::fmod(rotationAngle+colors[colorIndex].hue+colors[colorIndex].arcWidth/2.0f, 360.0f);
}
