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
	uchar * imageSeed = new uchar[size];
	std::memset(imageSeed, 0u, sizeof(uchar)*size);
	computeSeeds(imageYCbCr, imageSeed, h, w, 100, 0.02f);

	// create regions based solely on seeds
	int * regionIds = new int[size];
	int regionCount = createRegionIDFromSeeds(h, w, imageSeed, regionIds);

	// compute region size and average
	Region * regions = new Region[regionCount];
	computeRegionSizeAndAvg(regions, regionCount, regionIds, imageYCbCr, h, w);

	// compute frontiers
	std::multiset<PxDist, CmpPxDist> frontier;
	computeFrontier(h, w, imageYCbCr, regionIds, regions, frontier);

	std::cout << "growing regions" << std::endl;

	// grow regions
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

	std::cout << "\t1. by size" << std::endl;
	const int THRESHOLD_SIZE = (int)(size/1000.0f);
	mergeBySize(regionCount, regions, regionNeighbours, regionsAssociations, regionMin, THRESHOLD_SIZE);


	// update regionsIds
	updateRegionsID(size, regionIds, regionsAssociations);

	int test = *std::max_element(regionIds, regionIds+size);
	ushort * imTest = new ushort[size];
	for (int i=0; i<size; ++i)
		imTest[i] = (regionIds[i]/(float)test)*65535;
	ecrire_image_pgm_2o((baseName+"testmerge.pgm").data(), imTest, h, w);
	delete [] imTest;

	uchar * imTest2 = new uchar[size*3];
	for (int i=0; i<size; ++i)
	{
		ColorRGB rgb{};
		ycbcr_to_rgb(regions[regionIds[i]].avg, rgb);
		imTest2[3*i] = 255.0f*rgb.r;
		imTest2[3*i+1] = 255.0f*rgb.g;
		imTest2[3*i+2] = 255.0f*rgb.b;
	}
//		imTest[i] = (regionIds[i]/(float)test)*65535;
	ecrire_image_ppm((baseName+"testmerge_color.ppm").data(), imTest2, h, w);
	delete [] imTest;
	delete [] imTest2;

	delete [] regionsAssociations;
	delete [] regionNeighbours;
	delete [] regions;
	delete [] regionIds;
	delete [] imageYCbCr;
	delete [] imageSeed;
	delete [] image;

    return 0;
}
