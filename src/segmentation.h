#ifndef SEGMENTATION_H
#define SEGMENTATION_H

#include "ColorSpaces.h"
#include <vector>
#include <set>

enum class OtsuCriterion
{
	BCW_WCV,
	TV_WCV,
	BCV_TV
};

// seeds
float getPxStdDev(const float * img, int h, int w, int row, int col);
void computeSimilarityImage(const float * imageYCbCr, float * imageSimilarity, int h, int w);
double otsu(int k, double *ddp, int grayLevelCount, OtsuCriterion criterion);
void computeSeedsFromSimilarity(const float *imageYCbCr, const float *imageSimilarity, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance);
void filterOtsuSeedsEuclidean(const float *imageYCbCr, const unsigned char *imageSeedOtsu, unsigned char * out, int h, int w, float maxEuclideanDistance);
void computeSeeds(const float *imageYCbCr, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance);

// base regions
struct Px
{
	int row, col;
};
struct PxDist
{
	int row, col;
	double minDist;
};
struct Region
{
	ColorYCbCr avg;
	int size;
};
struct MinReg
{
	int label;
	double dist;
};


inline std::vector<Px> getNeighbours(Px p, int h, int w)
{
	std::vector<Px> res;
	int topRow = p.row-1; int topCol = p.col;
	int bottomRow = p.row + 1; int bottomCol = p.col;
	int leftRow = p.row; int leftCol = p.col-1;
	int rightRow = p.row; int rightCol = p.col + 1;
	if (topRow >= 0)
		res.push_back({topRow, topCol});
	if (bottomRow < h)
		res.push_back({bottomRow, bottomCol});
	if (leftCol >= 0)
		res.push_back({leftRow, leftCol});
	if (rightCol < w)
		res.push_back({rightRow, rightCol});
	return res;
}

inline std::vector<Px> getNeighbours(PxDist p, int h, int w)
{
	return getNeighbours(Px{p.row, p.col}, h, w);
}
struct CmpPxDist {
	bool operator() (const PxDist & a, const PxDist & b) const {
		return a.minDist > b.minDist;
	}
};

int createRegionIDFromSeeds(int h, int w, const unsigned char *imageSeed, int *regionIds);
void computeRegionSizeAndAvg(Region * regions, int regionCount, const int * regionIds, const float * imageYCbCr, int h, int w);
void computeFrontier(int h, int w, const float *imageYCbCr, const int *regionIds, const Region *regions,
					 std::multiset<PxDist, CmpPxDist> & frontier);
MinReg getClosestRegion(const float *imageYCbCr, const int *regionIds, const Region * regions, Px px, int h, int w);
void growRegions(int h, int w, const float *imageYCbCr, int *regionIds, Region *regions,
				 std::multiset<PxDist, CmpPxDist> & frontier);

// merge
float distanceRegion(int region1, int region2, const Region * regions);
void computeRegionNeighbors(int h, int w, const int *regionIds, std::set<int> * regionNeighbours);


void mergeByAverage(int regionCount, Region *regions, std::set<int> *regionNeighbours, const float threshold,
					int *regionsAssociations);
void mergeBySize(int regionCount, Region *regions, std::set<int> *regionNeighbours, int *regionsAssociations, int regionMin,
				 const int threshold);
void updateRegionsID(int size, int *regionIds, const int * regionsAssociations);

#endif // SEGMENTATION_H
