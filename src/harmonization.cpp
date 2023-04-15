#include "harmonization.h"
#include "ColorSpaces.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <cstring>
#include "segmentation.h"
#include "util.h"
#include "KMean.h"


Distribution Distribution::create_I_template(float angle, float arcWidth) {
	Distribution d;
	d.colors.push_back(DistributionColor{0, arcWidth});
	d.colors.push_back(DistributionColor{180.0f, arcWidth});
	d.rotationAngle = std::fmod(360.0f+angle, 360.0f);
	return d;
}

Distribution Distribution::create_i_template(float angle, float arcWidth) {
	Distribution d;
	d.colors.push_back(DistributionColor{0, arcWidth});
	d.rotationAngle = std::fmod(360.0f+angle, 360.0f);
	return d;
}

unsigned int Distribution::getColorCount() const {
	return colors.size();
}

float Distribution::getSectorBorder1(unsigned int colorIndex) const {
	return std::fmod(360.0f+rotationAngle+colors[colorIndex].hue-colors[colorIndex].arcWidth/2.0f, 360.0f);
}

float Distribution::getSectorBorder2(unsigned int colorIndex) const {
	return std::fmod(360.0f+rotationAngle+colors[colorIndex].hue+colors[colorIndex].arcWidth/2.0f, 360.0f);
}

const DistributionColor & Distribution::getColor(unsigned int index) const
{
	return colors[index];
}

float Distribution::getRotation() const
{
	return rotationAngle;
}

void Distribution::addColor(float hue, float arcWidth)
{
	colors.push_back(DistributionColor{hue, arcWidth});
}

Distribution::Distribution()
: colors()
, rotationAngle(0)
{
}

void shift(const Distribution & distrib, int h, int w, const unsigned char *image, unsigned char *imageOut)
{
	int size = h*w;
	unsigned char *imageSectorBorders = new unsigned char[size];
	find_sector_borders(distrib, image, imageSectorBorders, h, w);

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

void shiftRegions(const Distribution & distrib, int h, int w, const unsigned char *image, unsigned char *imageOut, const Region * regions, const int * imageRegions)
{
	int size = h*w;
	unsigned char * imageColorRegions = new unsigned char[size*3];
	for (int i=0; i<size; ++i)
	{
		ColorRGB rgb{};
		ycbcr_to_rgb(regions[imageRegions[i]].avg, rgb);
		imageColorRegions[3*i] = 255.0f*rgb.r;
		imageColorRegions[3*i+1] = 255.0f*rgb.g;
		imageColorRegions[3*i+2] = 255.0f*rgb.b;
	}
	unsigned char *imageSectorBorders = new unsigned char[size];
	find_sector_borders(distrib, imageColorRegions, imageSectorBorders, h, w);
	delete [] imageColorRegions;

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

void
find_sector_borders(const Distribution & distrib, const unsigned char *image, unsigned char *imageSectorBorders, int h,
					int w)
{
	int colorCount = distrib.getColorCount();
	int size = h*w;
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

void KMeanShift(HarmonyType type, int h, int w, const unsigned char *image, unsigned char *imageOut)
{
	int nbColors;
	switch (type)
	{
		case HarmonyType::COMPLEMENTARY:
			nbColors = 2;
			break;
		case HarmonyType::TRIADIC:
			nbColors = 3;
			break;
		case HarmonyType::TETRADIC_RECTANGLE:
			nbColors = 4;
			break;
		case HarmonyType::TETRADIC_SQUARE:
			nbColors = 4;
			break;
		case HarmonyType::SPLIT_COMPLEMENTARY:
			nbColors = 3;
			break;
		case HarmonyType::ANALOGOUS:
			nbColors = 3;
			break;
		default:
			break;
	}

	int size = h*w;
	uint * classes = new uint[size];
	ColorRGB * colors = new ColorRGB[nbColors];
	KMean(image, classes, colors, nbColors, h, w);

	// find which color has more pixels
	uint * pxCount = new uint[nbColors];
	std::memset(pxCount, 0, nbColors);
	for (int i=0; i<size; ++i)
		++pxCount[classes[i]];
	int max = -1;
	int maxVal = -1;
	for (int i=0; i<nbColors; ++i)
	{
		if (pxCount[i] > maxVal)
		{
			max = i;
			maxVal = pxCount[i];
		}
	}
	delete [] pxCount;

	ColorHSV * colorsHSV = new ColorHSV[nbColors];
	for (int i=0; i<2; ++i)
		rgb_to_hsv(colors[i], colorsHSV[i]);

	switch (type)
	{
		case HarmonyType::COMPLEMENTARY:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		case HarmonyType::TRIADIC:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+120.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+240.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		case HarmonyType::TETRADIC_RECTANGLE:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+60.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+240.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		case HarmonyType::TETRADIC_SQUARE:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+90.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+270.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		case HarmonyType::SPLIT_COMPLEMENTARY:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+150.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+210.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		case HarmonyType::ANALOGOUS:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h-30.0f+360.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+30.0f, 360.0f), 18.0f);
			shift(distrib, h, w, image, imageOut);
		}
			break;
		default:
			break;
	}
	delete [] classes;
	delete [] colorsHSV;
	delete [] colors;
}

void KMeanShiftRegions(HarmonyType type, int h, int w, const unsigned char *image, unsigned char *imageOut, const Region * regions, const int * imageRegions)
{
	int nbColors;
	switch (type)
	{
		case HarmonyType::COMPLEMENTARY:
			nbColors = 2;
			break;
		case HarmonyType::TRIADIC:
			nbColors = 3;
			break;
		case HarmonyType::TETRADIC_RECTANGLE:
			nbColors = 4;
			break;
		case HarmonyType::TETRADIC_SQUARE:
			nbColors = 4;
			break;
		case HarmonyType::SPLIT_COMPLEMENTARY:
			nbColors = 3;
			break;
		case HarmonyType::ANALOGOUS:
			nbColors = 3;
			break;
		default:
			break;
	}

	int size = h*w;
	uint * classes = new uint[size];
	ColorRGB * colors = new ColorRGB[nbColors];
	KMean(image, classes, colors, nbColors, h, w);

	// find which color has more pixels
	uint * pxCount = new uint[nbColors];
	std::memset(pxCount, 0, nbColors);
	for (int i=0; i<size; ++i)
		++pxCount[classes[i]];
	int max = -1;
	int maxVal = -1;
	for (int i=0; i<nbColors; ++i)
	{
		if (pxCount[i] > maxVal)
		{
			max = i;
			maxVal = pxCount[i];
		}
	}
	delete [] pxCount;

	ColorHSV * colorsHSV = new ColorHSV[nbColors];
	for (int i=0; i<2; ++i)
		rgb_to_hsv(colors[i], colorsHSV[i]);

	switch (type)
	{
		case HarmonyType::COMPLEMENTARY:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		case HarmonyType::TRIADIC:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+120.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+240.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		case HarmonyType::TETRADIC_RECTANGLE:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+60.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+240.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		case HarmonyType::TETRADIC_SQUARE:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+90.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+180.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+270.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		case HarmonyType::SPLIT_COMPLEMENTARY:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+150.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+210.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		case HarmonyType::ANALOGOUS:
		{
			Distribution distrib;
			distrib.addColor(colorsHSV[max].h, 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h-30.0f+360.0f, 360.0f), 18.0f);
			distrib.addColor(std::fmod(colorsHSV[max].h+30.0f, 360.0f), 18.0f);
			shiftRegions(distrib, h, w, image, imageOut, regions, imageRegions);
		}
			break;
		default:
			break;
	}
	delete [] classes;
	delete [] colorsHSV;
	delete [] colors;
}
