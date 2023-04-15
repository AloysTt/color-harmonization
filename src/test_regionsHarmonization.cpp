#include <iostream>
#include <vector>
#include <cstring>
#include <cfloat>
#include <algorithm>
#include <queue>
#include <set>
#include "image_ppm.h"
#include "ColorSpaces.h"
#include "segmentation.h"
#include "util.h"
#include "harmonization.h"

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		printf("Usage: ImageIn.ppm\n");
		return 1;
	}
	std::string imageName{argv[1]};
	std::string baseName = {imageName.substr(0, imageName.size()-4)};

    int h, w, size, size3;
    lire_nb_lignes_colonnes_image_ppm(imageName.data(), &h, &w);
	size3 = h * w * 3;
	size = h * w;
	uchar * image = new uchar[size3];
	lire_image_ppm(imageName.data(), image, size);

	// convert to ycbcr float
	float * imageYCbCr = new float[size3];
	rgb_to_ycbcr(image, imageYCbCr, h, w);

	// compute seeds
	std::cout << "computing region seeds" << std::endl;
	uchar * imageSeed = new uchar[size];
	std::memset(imageSeed, 0u, sizeof(uchar)*size);
	computeSeeds(imageYCbCr, imageSeed, h, w, 100, 0.02f);

	// create regions based solely on seeds
	std::cout << "creating regions" << std::endl;
	int * regionIds = new int[size];
	int regionCount = createRegionIDFromSeeds(h, w, imageSeed, regionIds);

	// compute region size and average
	std::cout << "computing region sizes and averages" << std::endl;
	Region * regions = new Region[regionCount];
	computeRegionSizeAndAvg(regions, regionCount, regionIds, imageYCbCr, h, w);

	// compute frontiers
	std::cout << "computing frontier" << std::endl;
	std::multiset<PxDist, CmpPxDist> frontier;
	computeFrontier(h, w, imageYCbCr, regionIds, regions, frontier);

	// grow regions
	std::cout << "growing regions" << std::endl;
	growRegions(h, w, imageYCbCr, regionIds, regions, frontier);

	// merge
	std::cout << "merging regions" << std::endl;
	// compute neighbors of each regions
	std::set<int> * regionNeighbours = new std::set<int>[regionCount];
	computeRegionNeighbors(h, w, regionIds, regionNeighbours);
	int * regionsAssociations = new int[regionCount];
	std::memset(regionsAssociations, -1, sizeof(int)*regionCount);

	// merging by average
	std::cout << "\t1. by average color" << std::endl;
	constexpr float THRESHOLD_AVG = 0.05f;
	int regionMin;
	mergeByAverage(regionCount, regions, regionNeighbours, THRESHOLD_AVG, regionsAssociations);

	std::cout << "\t2. by size" << std::endl;
	const int THRESHOLD_SIZE = (int)(size/1000.0f);
	mergeBySize(regionCount, regions, regionNeighbours, regionsAssociations, regionMin, THRESHOLD_SIZE);


	// update regionsIds
	updateRegionsID(size, regionIds, regionsAssociations);



	unsigned char * imageOut = new uchar[size3];
	std::cout << "Shift using K-mean for complementary colors\n";
	KMeanShiftRegions(HarmonyType::COMPLEMENTARY, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftComplementary.ppm").data(), imageOut, h, w);
	std::cout << "Shift using K-mean for triadic colors\n";
	KMeanShiftRegions(HarmonyType::TRIADIC, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftTriadic.ppm").data(), imageOut, h, w);
	std::cout << "Shift using K-mean for tetradic rectangle colors\n";
	KMeanShiftRegions(HarmonyType::TETRADIC_RECTANGLE, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftTetradicRectangle.ppm").data(), imageOut, h, w);
	std::cout << "Shift using K-mean for tetradic square colors\n";
	KMeanShiftRegions(HarmonyType::TETRADIC_SQUARE, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftTetradicSquare.ppm").data(), imageOut, h, w);
	std::cout << "Shift using K-mean for split complementary colors\n";
	KMeanShiftRegions(HarmonyType::SPLIT_COMPLEMENTARY, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftSplitComplementary.ppm").data(), imageOut, h, w);
	std::cout << "Shift using K-mean for analogous colors\n";
	KMeanShiftRegions(HarmonyType::ANALOGOUS, h, w, image, imageOut, regions, regionIds);
	ecrire_image_ppm((baseName+"_segmented_kmeanShiftAnalogous.ppm").data(), imageOut, h, w);


	delete [] imageOut;
	delete [] regionsAssociations;
	delete [] regionNeighbours;
	delete [] regions;
	delete [] regionIds;
	delete [] imageYCbCr;
	delete [] imageSeed;
	delete [] image;

    return 0;
}
