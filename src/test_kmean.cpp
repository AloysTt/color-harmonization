#include <random>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <iostream>
#include "image_ppm.h"
#include "harmonization.h"
#include "KMean.h"

typedef unsigned char uchar;

int main(int argc, char **argv)
{

	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

	int size3, size, h, w;
	lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h*w;
	unsigned char * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	unsigned char * imageOut = new uchar[size3];
	KMeanShift(HarmonyType::COMPLEMENTARY, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftComplementary.ppm").data(), imageOut, h, w);
	KMeanShift(HarmonyType::TRIADIC, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftTriadic.ppm").data(), imageOut, h, w);
	KMeanShift(HarmonyType::TETRADIC_RECTANGLE, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftTetradicRectangle.ppm").data(), imageOut, h, w);
	KMeanShift(HarmonyType::TETRADIC_SQUARE, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftTetradicSquare.ppm").data(), imageOut, h, w);
	KMeanShift(HarmonyType::SPLIT_COMPLEMENTARY, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftSplitComplementary.ppm").data(), imageOut, h, w);
	KMeanShift(HarmonyType::ANALOGOUS, h, w, image, imageOut);
	ecrire_image_ppm((baseName+"_kmeanShiftAnalogous.ppm").data(), imageOut, h, w);

	delete [] image;
	delete [] imageOut;
	return 0;
}

