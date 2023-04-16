#include "segmentation.h"
#include "util.h"
#include "ColorSpaces.h"
#include <cmath>
#include <cstring>
#include <cstdint>
#include <algorithm>
#include <queue>
#include <cfloat>
#include <set>

float getPxStdDev(const float * img, int h, int w, int row, int col)
{
	int w3 = 3*w;
	float avgY, avgCb, avgCr;
	avgY = avgCb = avgCr = 0.0f;
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			avgY += img[(row + i) * w3 + 3*(col + j)];
			avgCb += img[(row + i) * w3 + 3*(col + j) + 1];
			avgCr += img[(row + i) * w3 + 3*(col + j) + 2];
		}
	}
	avgY /= 9.0f;
	avgCb /= 9.0f;
	avgCr /= 9.0f;

	float stdY, stdCb, stdCr;
	stdY = stdCb = stdCr = 0.0f;
	for (int i=-1; i<=1; ++i)
	{
		for (int j=-1; j<=1; ++j)
		{
			stdY += (img[(row + i) * w3 + 3*(col + j)] - avgY) * (img[(row + i) * w3 + 3*(col + j)] - avgY);
			stdCb += (img[(row + i) * w3 + 3*(col + j) + 1] - avgCb) * (img[(row + i) * w3 + 3*(col + j) + 1] - avgCb);
			stdCr += (img[(row + i) * w3 + 3*(col + j) + 2] - avgCr) * (img[(row + i) * w3 + 3*(col + j) + 2] - avgCr);
		}
	}
	stdY /= 9.0f;
	stdCb /= 9.0f;
	stdCr /= 9.0f;
	stdY = std::sqrt(stdY);
	stdCb = std::sqrt(stdCb);
	stdCr = std::sqrt(stdCr);

	return stdY + stdCb + stdCr;
}

void computeSimilarityImage(const float * imageYCbCr, float * imageSimilarity, int h, int w)
{
	int size = h*w;
	std::memset(imageSimilarity, 0, size * sizeof(float));
	for (int row=1; row<h-1; ++row)
	{
		for (int col=1; col<w-1; ++col)
		{
			imageSimilarity[row * w + col] = getPxStdDev(imageYCbCr, h, w, row, col);
		}
	}
	float maxStdDev = *std::max_element(imageSimilarity, imageSimilarity + size);
	for (int i=0; i<size; ++i)
		imageSimilarity[i] = 1.0f - imageSimilarity[i] / maxStdDev;
}

double otsu(int k, double *ddp, int grayLevelCount, OtsuCriterion criterion)
{
	double w0, u;
	w0 = u = 0.0;
	for(int i=0; i<=k; ++i)
	{
		w0+=ddp[i];
		u+=i*ddp[i];
	}
	double w1 = 1.0-w0;

	double ut = u;
	for (int i=k+1; i < grayLevelCount; i++)
		ut += i * ddp[i];

	double u0 = u/w0;
	double u1 = (ut - u) / w1;

	double v0 = 0;
	double v1 = 0;
	for(int i=0; i<=k; i++)
		v0 += ((double)i-u0)*((double)i-u0)
			  * (ddp[i] / w0);
	for(int i=k+1; i < grayLevelCount; i++)
		v1 += ((double)i-u1)*((double)i-u1)
			  * (ddp[i] / w1);

	switch (criterion)
	{
		case OtsuCriterion::BCW_WCV:
		{
			double WCV = w0 * v0 + w1 * v1;
			double BCV = w0 * w1 * (u1 - u0) * (u1 - u0);
			return BCV / WCV;
		}
		case OtsuCriterion::TV_WCV:
		{
			double TV = 0;
			for (int i=0; i<grayLevelCount; ++i)
				TV += ((double)i-ut)*((double)i-ut)*ddp[i];
			double WCV = w0 * v0 + w1 * v1;
			return TV/WCV;
		}
		case OtsuCriterion::BCV_TV:
		{
			double TV = 0;
			for (int i=0; i<grayLevelCount; ++i)
				TV += ((double)i-ut)*((double)i-ut)*ddp[i];
			double BCV = w0*w1*(u1-u0)*(u1-u0);
			return BCV/TV;
		}
		default:
			return 0;
	}
}

void computeSeedsFromSimilarity(const float *imageYCbCr, const float *imageSimilarity, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance)
{
	int size = h*w;
	// split similarity image into buckets
	uint * imageOtsu = new uint[size];
	for (int i=0; i<size; ++i)
		imageOtsu[i] = static_cast<float>(otsuBuckets-1)*imageSimilarity[i];

	// compute ddp
	double * ddp = new double[otsuBuckets];
	create_ddp(imageOtsu, h, w, ddp, otsuBuckets);

	// get threshold using Otsu's method
	int threshold = 0;
	double thresholdVal = 0;
	for (int i=0; i<otsuBuckets; ++i)
	{
		double otsuVal = otsu(i, ddp, otsuBuckets, OtsuCriterion::BCW_WCV);
		if (otsuVal > thresholdVal)
		{
			threshold = i;
			thresholdVal = otsuVal;
		}
	}

	// create binary seed image
	unsigned char * imageSeedOtsu = new unsigned char[size];
	for (int i=0; i<size; ++i)
		imageSeedOtsu[i] = imageOtsu[i] >= threshold ? 255 : 0;

	// 2nd condition : euclidean distance
	filterOtsuSeedsEuclidean(imageYCbCr, imageSeedOtsu, imageSeeds, h, w, maxEuclideanDistance);

	delete [] imageSeedOtsu;
	delete [] imageOtsu;
	delete [] ddp;
}

void filterOtsuSeedsEuclidean(const float *imageYCbCr, const unsigned char *imageSeedOtsu, unsigned char * out, int h, int w, float maxEuclideanDistance)
{
	for (int row=1; row<h-1; ++row)
	{
		for (int col=1; col<w-1; ++col)
		{
			if (imageSeedOtsu[row*w+col] == 0)
			{
				out[row*w+col] = 0;
				continue;
			}
			float max = 0.0f;
			for (int i=-1; i<=1; ++i)
			{
				for (int j=-1; j<=1; ++j)
				{
					float dist = ycbcr_distance_relative_euclidean(
						imageYCbCr[row*w*3+col*3],
						imageYCbCr[row*w*3+col*3+1],
						imageYCbCr[row*w*3+col*3+2],
						imageYCbCr[(row+i)*w*3+(col+j)*3],
						imageYCbCr[(row+i)*w*3+(col+j)*3+1],
						imageYCbCr[(row+i)*w*3+(col+j)*3+2]
					);
					if (dist > max)
						max = dist;
				}
			}
			if (max > maxEuclideanDistance)
				out[row*w+col] = 0;
			else
				out[row*w+col] = 255;
		}
	}
}

void computeSeeds(const float *imageYCbCr, unsigned char *imageSeeds, int h, int w, unsigned int otsuBuckets, float maxEuclideanDistance)
{
	int size = h*w;
	float * imageSimilarity = new float[size];
	computeSimilarityImage(imageYCbCr, imageSimilarity, h, w);

	computeSeedsFromSimilarity(imageYCbCr, imageSimilarity, imageSeeds, h, w, otsuBuckets, maxEuclideanDistance);

	delete [] imageSimilarity;
}

int createRegionIDFromSeeds(int h, int w, const unsigned char *imageSeed, int *regionIds)
{
	int size = h*w;
	memset(regionIds, -1, size * sizeof(int));
	int regionCount = 0;
	std::queue<Px> fifo;
	for (int row=0; row<h; ++row)
	{
		for (int col = 0; col < w; ++col)
		{
			if (imageSeed[row*w+col] != 255)
				continue;
			if (regionIds[row*w+col] != -1)
				continue;
			std::vector<Px> neighs = getNeighbours(Px{row, col}, h, w);
			auto foundRegion = std::find_if(neighs.begin(), neighs.end(), [&regionIds, &w](const Px & n){return regionIds[n.row*w+n.col] != -1;});
			if (foundRegion != neighs.end())
			{
				regionIds[row*w+col] = regionIds[(*foundRegion).row*w + (*foundRegion).col];
				continue;
			}

			regionIds[row*w+col] = regionCount++;
			for (const Px & neigh : neighs)
			{
				if (imageSeed[neigh.row*w + neigh.col] == 255)
					fifo.push(neigh);
			}

			while (!fifo.empty())
			{
				Px p = fifo.front();
				fifo.pop();
				std::vector<Px> pNeighs = getNeighbours(p, h, w);
				for (const Px & neigh : pNeighs)
				{
					if (imageSeed[neigh.row*w + neigh.col] == 255 && regionIds[neigh.row*w+neigh.col] == -1)
					{
						regionIds[neigh.row*w + neigh.col] = regionIds[row*w+col];
						fifo.push(neigh);
					}
				}
			}
		}
	}
	return regionCount;
}

void computeRegionSizeAndAvg(Region * regions, int regionCount, const int * regionIds, const float * imageYCbCr, int h, int w)
{
	int size = h*w;
	std::memset(regions, 0, sizeof(Region)*regionCount);
	for (int i=0; i<size; ++i)
	{
		if (regionIds[i] == -1) continue;
		Region & r = regions[regionIds[i]];
		++(r.size);
		r.avg.y+=imageYCbCr[3*i];
		r.avg.cb+=imageYCbCr[3*i+1];
		r.avg.cr+=imageYCbCr[3*i+2];
	}
	for (int i=0; i<regionCount; ++i)
	{
		Region & r = regions[i];
		r.avg.y/=r.size;
		r.avg.cb/=r.size;
		r.avg.cr/=r.size;
	}
}

void computeFrontier(int h, int w, const float *imageYCbCr, const int *regionIds, const Region *regions,
					 std::multiset<PxDist, CmpPxDist> & frontier)
{
	for (int row=0; row < h; ++row)
	{
		for (int col = 0; col < w; ++col)
		{
			if (regionIds[row * w + col] != -1)
				continue;
			std::vector<Px> neighs = getNeighbours(Px{row, col}, h, w);
			MinReg minReg = getClosestRegion(imageYCbCr, regionIds, regions, Px{row, col}, h, w);
			if (minReg.label != -1)
				frontier.emplace(row, col, minReg.dist);
		}
	}
//	std::sort(frontier.begin(), frontier.end(),
//			  [](const PxDist & px1, const PxDist & px2)
//			  {
//				  return px1.minDist > px2.minDist;
//			  }
//	);
}

MinReg getClosestRegion(const float *imageYCbCr, const int *regionIds, const Region * regions, Px px, int h, int w)
{
	double minDist = DBL_MAX;
	int minRegion = -1;
	const float * valPx = imageYCbCr+px.row*w*3+px.col*3;
	for (const Px & neigh : getNeighbours(px, h, w))
	{
		int region = regionIds[neigh.row * w + neigh.col];
		if (region != -1)
		{
			double dist = ycbcr_distance_relative_euclidean_squared(
				valPx[0], valPx[1], valPx[2],
				regions[region].avg.val[0], regions[region].avg.val[1], regions[region].avg.val[2]
			);
			if (dist < minDist)
			{
				minRegion = region;
				minDist = dist;
			}
		}
	}
	return {minRegion, minDist};
}

void growRegions(int h, int w, const float *imageYCbCr, int *regionIds, Region *regions,
				 std::multiset<PxDist, CmpPxDist> & frontier)
{
	while (!frontier.empty())
	{
		PxDist px = *(--frontier.end());
		frontier.erase(--frontier.end());
		std::vector<Px> neighs = getNeighbours(px, h, w);

		// check if all labels are the same
		int label = -1;
		bool sameLabel = true;
		for (int i=0; i<neighs.size(); ++i)
		{
			if (label == -1 && regionIds[neighs[i].row*w+neighs[i].col] != -1)
				label = regionIds[neighs[i].row*w+neighs[i].col];
			else if (label != -1 && regionIds[neighs[i].row*w+neighs[i].col] != label)
			{
				sameLabel = false;
				break;
			}
		}
		if (label != -1 && sameLabel)
		{
			regionIds[px.row*w+px.col] = label;
		}
		else
		{
			// check nearest region
			std::vector<int> labels;
			std::vector<double> distances;
			for (int i=0; i<neighs.size(); ++i)
			{
				int curLabel = regionIds[neighs[i].row*w+neighs[i].col];
				if (curLabel == -1) continue;
				if (std::find(labels.begin(), labels.end(), curLabel) == labels.end())
				{
					labels.push_back(curLabel);
					double dist = ycbcr_distance_relative_euclidean_squared(imageYCbCr+px.row*w*3+px.col*3, regions[curLabel].avg.val);
					distances.push_back(dist);
				}
			}
			label = labels[std::min_element(distances.begin(), distances.end()) - distances.begin()];
			regionIds[px.row*w+px.col] = label;
		}

		// update region's mean
		{
			const float * pxVal = imageYCbCr+px.row*w*3+px.col*3;
			Region & region = regions[label];
			region.avg.y*= region.size;
			region.avg.cb*= region.size;
			region.avg.cr*= region.size;
			region.avg.y+=pxVal[0];
			region.avg.cb+=pxVal[1];
			region.avg.cr+=pxVal[2];
			++region.size;
			region.avg.y/= region.size;
			region.avg.cb/= region.size;
			region.avg.cr/= region.size;
		}

		for (int i=0; i<neighs.size(); ++i)
		{
			Px & neigh = neighs[i];
			if (regionIds[neigh.row * w + neigh.col] == -1)
			{
				MinReg minReg = getClosestRegion(imageYCbCr, regionIds, regions, neigh, h, w);
				if (minReg.label != -1)
				{
					// find if the pixel has a labelled neighbour other than the pixel that we just processed
					std::vector<Px> neighs2 = getNeighbours(neigh, h, w);
					bool alreadyInFrontier = false;
					for (const Px & neigh2 : neighs2)
					{
						if (neigh2.col == px.col && neigh2.row == px.row)
							continue;
						if (regionIds[neigh2.row*w + neigh2.col] != -1)
							alreadyInFrontier = true;
					}
					if (alreadyInFrontier)
						continue;
					PxDist newPx{neigh.row, neigh.col, minReg.dist};
//					frontier.insert(
//						std::upper_bound(frontier.begin(), frontier.end(), newPx,
//										 [](const PxDist & px1, const PxDist & px2){return px1.minDist < px2.minDist;}),
//						newPx);
					frontier.insert(newPx);
				}
			}
		}
	}
}

float distanceRegion(int region1, int region2, const Region * regions)
{
	const ColorYCbCr & cReg1 = regions[region1].avg;
	const ColorYCbCr & cReg2 = regions[region2].avg;
	float denom = std::min(
		std::sqrt(cReg1.y*cReg1.y + cReg1.cb*cReg1.cb + cReg1.cr*cReg1.cr),
		std::sqrt(cReg2.y*cReg2.y + cReg2.cb*cReg2.cb + cReg2.cr*cReg2.cr)
		);
	float num = std::sqrt(
		(cReg1.y-cReg2.y)*(cReg1.y-cReg2.y)
		+ (cReg1.cb-cReg2.cb)*(cReg1.cb-cReg2.cb)
		+ (cReg1.cr-cReg2.cr)*(cReg1.cr-cReg2.cr)
		);
	return num/denom;
}

void computeRegionNeighbors(int h, int w, const int *regionIds, std::set<int> * regionNeighbours)
{
	for (int row=0; row < h; ++row)
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
}

void mergeByAverage(int regionCount, Region *regions, std::set<int> *regionNeighbours, const float threshold,
					int *regionsAssociations)
{
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

	// we now have the closest region distance for each region
	int regionMin = std::min_element(closestRegions, closestRegions+regionCount) - closestRegions;
	while (closestRegions[regionMin] < threshold)
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
		regions[minNeighbor].size = INT32_MAX; // for the size merge
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
	delete [] closestRegions;
}


void mergeBySize(int regionCount, Region *regions, std::set<int> *regionNeighbours, int *regionsAssociations, int regionMin,
			const int threshold)
{
	regionMin = std::min_element(regions, regions + regionCount, [](const Region & r1, const Region & r2){return r1.size < r2.size;}) - regions;
	while (regions[regionMin].size < threshold)
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
		for (int reg : regionNeighbours[regionMin])
		{
			regionNeighbours[minNeighbor].insert(reg);
			regionNeighbours[reg].insert(minNeighbor);
			regionNeighbours[reg].erase(regionMin);
		}
		regionNeighbours[minNeighbor].erase(regionMin);
		regionNeighbours[minNeighbor].erase(minNeighbor);
		regionNeighbours[regionMin].clear();

		regions[minNeighbor].avg = ColorYCbCr{newAvgY, newAvgCb, newAvgCr};
		regions[minNeighbor].size = newSize;

		// layer of indirection to avoid iterating over the entire image
		regionsAssociations[regionMin] = minNeighbor;

		regions[regionMin].size = INT32_MAX;

		regionMin = std::min_element(regions, regions+regionCount, [](const Region & r1, const Region & r2){return r1.size < r2.size;}) - regions;
	}
}

void updateRegionsID(int size, int *regionIds, const int *regionsAssociations)
{
	for (int i=0; i < size; ++i)
	{
		int current = regionIds[i];
		while (regionsAssociations[current] != -1)
		{
			current=regionsAssociations[current];
		}
		regionIds[i] = current;
	}
}
