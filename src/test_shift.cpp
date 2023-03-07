#include <cmath>
#include <cfloat>
#include <random>
#include <algorithm>
#include "image_ppm.h"

#include "harmonization.h"
#include "HSV.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


typedef unsigned char uchar;

inline float gaussian(float x, float avg, float sd)
{
	return (1.0f / (sd * std::sqrt(2.0f * PI)))
	* std::exp(
		-0.5*pow((x-avg)/sd, 2.0f)
		);
}

inline float rad(float angle)
{
	return angle*PI/180.0f;
}

int main(int argc, char **argv)
{

	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

	Distribution distrib = Distribution::create_I_template(60.0f, 60.0f);
	int h, w, size3, size;
	lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h*w;
	uchar * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	int colorCount = distrib.getColorCount();
	// assign each pixel to a sector border
	// there are 2 borders for each sector (and 1 sector per color)
	// ex: 0 and 1 are borders for the color 0; 2 and 3 for the color 1
	uchar * imageSectorBorders = new uchar[size];
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
					sin(Hrad-sectorBorder1),
					cos(Hrad-sectorBorder1)
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
					sin(Hrad-sectorBorder2),
					cos(Hrad-sectorBorder2)
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

	uchar * imageOut = new uchar[size3];
	for (int i=0; i < size; ++i)
	{
		ColorRGB rgbIn{(float) image[3 * i] / 255.0f, (float) image[3 * i + 1] / 255.0f,
					 (float) image[3 * i + 2] / 255.0f};
		ColorHSV hsv{};
		rgb_to_hsv(rgbIn, hsv);

		DistributionColor closestDistribColor = distrib.getColor(imageSectorBorders[i]/2);
		closestDistribColor.hue+= distrib.getRotation();

		float diff = atan2(
			sin(rad(hsv.h)-rad(closestDistribColor.hue)),
			cos(rad(hsv.h)-rad(closestDistribColor.hue))
		);
		diff = diff*180.0f/PI;
		float absdiff = std::abs(diff);
		float invmax = 1.0f/gaussian(0.0f, 0.0f, closestDistribColor.arcWidth/2.0f);

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

		ColorRGB rgbOut{};
		hsv_to_rgb(hsv, rgbOut);

		imageOut[3 * i] = rgbOut.r * 255.0f;
		imageOut[3 * i + 1] = rgbOut.g  * 255.0f;
		imageOut[3 * i + 2] = rgbOut.b * 255.0f;
	}
	ecrire_image_ppm((baseName+"_shift.ppm").data(), imageOut, h, w);


	delete [] image;
	delete [] imageSectorBorders;
	delete [] imageOut;
	return 0;
}