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
	std::vector<PxDist> frontier;
	computeFrontier(h, w, imageYCbCr, regionIds, regions, frontier);

	std::cout << "growing regions" << std::endl;

	// grow regions
	growRegions(h, w, imageYCbCr, regionIds, regions, frontier);

	std::cout << "merging regions" << std::endl;
	std::set<int> * regionNeighbours = new std::set<int>[regionCount];
	for (int row=0; row<h; ++row)
	{
		for (int col=0; col<w-1; ++col)
		{
			int reg1 = regionIds[row*w+col];
			int reg2 = regionIds[row*w+col+1];
			if (reg1==reg2)
				continue;

			regionNeighbours[reg1].insert(reg2);
			regionNeighbours[reg2].insert(reg1);
		}
	}

	for (int row=0; row<h-1; ++row)
	{
		for (int col=0; col<w; ++col)
		{
			int reg1 = regionIds[row*w+col];
			int reg2 = regionIds[(row+1)*w+col];
			if (reg1==reg2)
				continue;

			regionNeighbours[reg1].insert(reg2);
			regionNeighbours[reg2].insert(reg1);
		}
	}

	float * closestRegions = new float[regionCount];
	for (int i=0; i<regionCount; ++i)
		closestRegions[i] = FLT_MAX;
	for (int reg1=0; reg1 < regionCount; ++reg1)
	{
		for (int reg2 : regionNeighbours[reg1])
		{
			float val = distanceRegion(reg1, reg2, regions);
			if (val < closestRegions[reg1])
				closestRegions[reg1] = val;
			if (val < closestRegions[reg2])
				closestRegions[reg2] = val;
		}
	}
	constexpr float THRESHOLD = 0.1;
	int * regionsAssociations = new int[regionCount];
	std::memset(regionsAssociations, -1, sizeof(int)*regionCount);

	// we now have the closest region distance for each region
	int regionMin = std::min_element(closestRegions, closestRegions+regionCount) - closestRegions;
	while (closestRegions[regionMin] < THRESHOLD)
	{
		int minNeighbor = -100000;
		float min = FLT_MAX;
		for (int neigh : regionNeighbours[regionMin])
		{
			float val = distanceRegion(regionMin, neigh, regions);
			if (val < min)
			{
				min = val;
				minNeighbor = neigh;
			}
		}
		// new average
		int newSize = regions[regionMin].size + regions[minNeighbor].size;
		float newAvgY = regions[regionMin].avg.y*regions[regionMin].size + regions[minNeighbor].avg.y*regions[minNeighbor].size;
		float newAvgCb = regions[regionMin].avg.cr*regions[regionMin].size + regions[minNeighbor].avg.cb*regions[minNeighbor].size;
		float newAvgCr = regions[regionMin].avg.cr*regions[regionMin].size + regions[minNeighbor].avg.cr*regions[minNeighbor].size;
		newAvgY/=newSize;
		newAvgCb/=newSize;
		newAvgCr/=newSize;
		for (int reg : regionNeighbours[minNeighbor])
		{
			regionNeighbours[regionMin].insert(reg);
			regionNeighbours[reg].insert(regionMin);
			regionNeighbours[reg].erase(minNeighbor);
		}
		regionNeighbours[regionMin].erase(minNeighbor);
		regionNeighbours[regionMin].erase(regionMin);
		regionNeighbours[minNeighbor].clear();

		regions[regionMin].avg = ColorYCbCr{newAvgY, newAvgCb, newAvgCr};
		regions[regionMin].size = newSize;

		// layer of indirection to avoid iterating over the entire image
		regionsAssociations[minNeighbor] = regionMin;
		//
		// compute new distance to closest neighbour
		closestRegions[minNeighbor] = FLT_MAX;
		closestRegions[regionMin] = FLT_MAX;
		for (int neigh : regionNeighbours[regionMin])
		{
			float val = distanceRegion(regionMin, neigh, regions);
			if (val < closestRegions[regionMin])
				closestRegions[regionMin] = val;
			if (val < closestRegions[neigh])
				closestRegions[neigh] = val;
		}

		regionMin = std::min_element(closestRegions, closestRegions+regionCount) - closestRegions;
	}

	// update regionsIds
	for (int i=0; i<size; ++i)
	{
		int current = regionIds[i];
		while (regionsAssociations[current] != -1)
		{
			current=regionsAssociations[current];
		}
		regionIds[i] = current;
	}




	int test = *std::max_element(regionIds, regionIds+size);
	ushort * imTest = new ushort[size];
	for (int i=0; i<size; ++i)
		imTest[i] = (regionIds[i]/(float)test)*65535;
	ecrire_image_pgm_2o((baseName+"testmerge.pgm").data(), imTest, h, w);
	delete [] imTest;

	delete [] regionNeighbours;
	delete [] regions;
	delete [] regionIds;
	delete [] imageYCbCr;
	delete [] imageSeed;
	delete [] image;

    return 0;
}
