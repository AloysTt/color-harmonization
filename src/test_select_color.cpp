#include <cmath>
#include <cfloat>
#include "image_ppm.h"

#include "harmonization.h"
#include "ColorSpaces.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


typedef unsigned char uchar;

int main(int argc, char **argv)
{
	Distribution distrib = Distribution::create_i_template(0, 90);
	int h, w, size3, size;
	lire_nb_lignes_colonnes_image_ppm("lena.ppm", &h, &w);
	size3 = h * w * 3;
	size = h*w;
	uchar * image = new uchar[size3];
	lire_image_ppm("lena.ppm", image, size);

	int colorCount = distrib.getColorCount();
	// assign each pixel to a sector border
	// there are 2 borders for each sector (and 1 sector per color)
	// ex: 0 and 1 are borders for the color 0; 2 and 3 for the color 1
	uchar * imageSectorBorders = new uchar[size];
	for (int i=0; i<size; ++i)
		imageSectorBorders[i] = 0;
	for (int i=0; i < size; ++i)
	{
		ColorRGB rgb{(float)image[3*i] / 255.0f, (float)image[3*i + 1] / 255.0f, (float)image[3*i + 2] / 255.0f};
		ColorHSV hsv{};
		rgb_to_hsv(rgb, hsv);

		for (int color=0; color<colorCount; ++color)
		{
			float sectorBorder1 = distrib.getSectorBorder1(color);
			float sectorBorder2 = distrib.getSectorBorder2(color);
			if (hsv.h > sectorBorder1 && hsv.h < sectorBorder2)
				imageSectorBorders[i] = color+1;
		}
	}

	// debug image
	for (int i=0; i<size; ++i)
		imageSectorBorders[i] *= 127;
	ecrire_image_pgm("lena_test2.pgm", imageSectorBorders, h, w);


	delete [] image;
	delete [] imageSectorBorders;
	return 0;
}