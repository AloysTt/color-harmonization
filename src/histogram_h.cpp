#include <cmath>
#include <cfloat>
#include <string>
#include <fstream>
#include "image_ppm.h"

#include "ColorSpaces.h"

#define PI 3.1415926535897932384626433832795028841971693993751058209749445923


typedef unsigned char uchar;
typedef unsigned int uint;

inline float gaussian(float x, float avg, float sd)
{
	return (1.0f / (sd * std::sqrt(2.0f * PI))) * std::exp(-(std::pow(x - avg, 2.0f) / (2.0f * sd * sd)));
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

	int h, w, size3, size;
	lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h*w;
	uchar * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	uint histo[360];
	for (int i=0; i < 360; ++i)
		histo[i] = 0;
	for (int i=0; i < size; ++i)
	{
		ColorRGB rgb{(float)image[3*i] / 255.0f, (float)image[3*i + 1] / 255.0f, (float)image[3*i + 2] / 255.0f};
		ColorHSV hsv{};
		rgb_to_hsv(rgb, hsv);
		++histo[(int)hsv.h];
	}

	std::ofstream f(baseName+"_histo.dat", std::ios::out | std::ios::trunc);
	for (int i=0; i < 360; ++i)
		f << i << "\t" << histo[i] << "\n";
	f.close();

	delete [] image;
	return 0;
}