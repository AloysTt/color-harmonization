#include "harmonization.h"
#include "HSV.h"
#include <cmath>
#include <cfloat>
#include "util.h"

Distribution Distribution::create_I_template(float angle, float arcWidth) {
	Distribution d;
	d.colors.push_back(DistributionColor{0, arcWidth});
	d.colors.push_back(DistributionColor{180.0f, arcWidth});
	d.rotationAngle = std::fmod(angle, 360.0f);
	return d;
}

Distribution Distribution::create_i_template(float angle, float arcWidth) {
	Distribution d;
	d.colors.push_back(DistributionColor{0, arcWidth});
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

const DistributionColor & Distribution::getColor(unsigned int index) const
{
	return colors[index];
}

float Distribution::getRotation() const
{
	return rotationAngle;
}

void shift(const Distribution & distrib, int h, int w, const unsigned char *image, unsigned char *imageOut)
{
	int size = h*w;
	unsigned char *imageSectorBorders = new unsigned char[size];
	find_sector_borders(distrib, image, size, imageSectorBorders);

	for (int i=0; i < size; ++i)
	{
		ColorRGB rgbIn{(float) image[3 * i] / 255.0f, (float) image[3 * i + 1] / 255.0f,
					 (float) image[3 * i + 2] / 255.0f};
		ColorHSV hsv{};
		rgb_to_hsv(rgbIn, hsv);

		DistributionColor closestDistribColor = distrib.getColor(imageSectorBorders[i]/2);
		closestDistribColor.hue+= distrib.getRotation();

		// compute difference to find out whether we should add or subtract from the sector's hue
		float diff = std::atan2(
			std::sin(rad(hsv.h) - rad(closestDistribColor.hue)),
			std::cos(rad(hsv.h) - rad(closestDistribColor.hue))
		);
		diff = diff*180.0f/PI;
		float absdiff = std::abs(diff);
		float invmax = 1.0f/gaussian(0.0f, 0.0f, closestDistribColor.arcWidth/2.0f);

		// apply shift
		float newHue;
		if (diff < 0.0f)
			newHue = closestDistribColor.hue
			- (closestDistribColor.arcWidth/2.0f)
			* (1.0f -
				invmax*gaussian(absdiff, 0.0f, closestDistribColor.arcWidth/2.0f));
		else
			newHue = closestDistribColor.hue
					 + (closestDistribColor.arcWidth/2.0f)
					   * (1.0f -
				invmax*gaussian(absdiff, 0.0f, closestDistribColor.arcWidth/2.0f));
		hsv.h = newHue;

		// convert back to rgb and write to output image
		ColorRGB rgbOut{};
		hsv_to_rgb(hsv, rgbOut);
		imageOut[3 * i] = rgbOut.r * 255.0f;
		imageOut[3 * i + 1] = rgbOut.g  * 255.0f;
		imageOut[3 * i + 2] = rgbOut.b * 255.0f;
	}
	delete [] imageSectorBorders;
}

void find_sector_borders(const Distribution & distrib, const unsigned char *image, int size,
						 unsigned char *imageSectorBorders)
{
	int colorCount = distrib.getColorCount();
	// assign each pixel to a sector border
	// there are 2 borders for each sector (and 1 sector per color)
	// ex: 0 and 1 are borders for the color 0; 2 and 3 for the color 1
	for (int i=0; i < size; ++i)
	{
		ColorRGB rgb{(float)image[3*i] / 255.0f, (float)image[3*i + 1] / 255.0f, (float)image[3*i + 2] / 255.0f};
		ColorHSV hsv{};
		rgb_to_hsv(rgb, hsv);

		float Hrad = hsv.h*PI/180.0f;
		float minDist = FLT_MAX;
		int closestSector = 0;
		for (int color=0; color<colorCount; ++color)
		{
			// check border 1
			float sectorBorder1 = distrib.getSectorBorder1(color)*PI/180.0f;
			float dist = atan2(
					sinf(Hrad-sectorBorder1),
					cosf(Hrad-sectorBorder1)
					);
			dist = std::abs(dist)*180.0f/PI;
			if (dist < minDist)
			{
				minDist = dist;
				closestSector = color*2;
			}

			// check border 2
			float sectorBorder2 = distrib.getSectorBorder2(color)*PI/180.0f;
			dist = atan2(
					sinf(Hrad-sectorBorder2),
					cosf(Hrad-sectorBorder2)
			);
			dist = std::abs(dist)*180.0f/PI;
			if (dist < minDist)
			{
				minDist = dist;
				closestSector = color*2+1;
			}
		}
		imageSectorBorders[i] = closestSector;
	}
}
