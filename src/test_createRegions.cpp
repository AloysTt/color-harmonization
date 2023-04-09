#include <iostream>
#include <vector>
#include <cstring>
#include <cfloat>
#include <algorithm>
#include <queue>
#include "image_ppm.h"
#include "ColorSpaces.h"
#include "segmentation.h"
#include "util.h"

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;
/*using namespace std;
typedef unsigned char uchar;


bool pixelSimilar(uchar* ImgIn, int seed, int pt, int seuil)
{
	ColorRGB rgbPixel{(float)ImgIn[3*pt] / 255.0f, (float)ImgIn[3*pt + 1] / 255.0f, (float)ImgIn[3*pt + 2] / 255.0f};
	ColorHSV hsvPixel{};
	rgb_to_hsv(rgbPixel, hsvPixel);
	ColorRGB rgbSeed{(float)ImgIn[3*seed] / 255.0f, (float)ImgIn[3*seed + 1] / 255.0f, (float)ImgIn[3*seed + 2] / 255.0f};
	ColorHSV hsvSeed{};
	rgb_to_hsv(rgbSeed, hsvSeed);

    float hDiff = min(abs(hsvPixel.h - hsvSeed.h), 180 - abs(hsvPixel.h - hsvSeed.h));
    hDiff /= 360;
    float sDiff = abs(hsvPixel.s - hsvSeed.s);
    float vDiff = abs(hsvPixel.v - hsvSeed.v);
    float dist = sqrt(hDiff * hDiff + sDiff * sDiff + vDiff * vDiff);
    dist*=360;

    return (dist < seuil);
}

void regionGrow(uchar* imgIn, uchar* reg, int seed, int seuil, int nH, int nW)
{
    vector<int> points;
    points.push_back(seed);
    int i, j;

    while (!points.empty())
    {
        int pt = points.back();
        points.pop_back();
        i = pt%nW;
        j = pt/nW;

        if (i < 0 || j < 0 || i >= nW || j >= nH || reg[pt] == 255)
            continue;

        if (pixelSimilar(imgIn, seed, pt, seuil))
        {
            reg[pt] = 255;
            points.push_back(j*nW+i+1);
            points.push_back(j*nW+i-1);
            points.push_back((j+1)*nW+i);
            points.push_back((j-1)*nW+i);
        }
    }
}

void rgbToYcbcr(uchar* ImgIn, uchar* ImgOut, int nH, int nW){
    int nW3 = nW*3;
    for (int i=0; i<nH; i++){
        for(int j=0; j<nW; j++){
            int place = i*nW3+j*3;
            float R = ImgIn[place];
            float G = ImgIn[place+1];
            float B = ImgIn[place+2];
            ImgOut[place] = 0.299f*R + 0.587*G + 0.114f*B;
            ImgOut[place+1] = -0.1687f*R - 0.3313f*G + 0.5f*B +128.0f;
            ImgOut[place+2] = 0.5f*R - 0.4187f*G - 0.0813f*B + 128.0f;
        }
    }
}

void ycbcrTorgb(uchar* ImgIn, uchar* ImgOut){
    int nW3 = nW*3;
    for (int i=0; i<nH; i++){
        for(int j=0; j<nW; j++){
            int place = i*nW3+j*3;
            float Y = ImgIn[place];
            float Cb = ImgIn[place+1];
            float Cr = ImgIn[place+2];
            ImgOut[place] = clamp(Y + 1.402f*(Cr - 128.0f));
            ImgOut[place+1] = clamp(Y - 0.34414f*(Cb - 128.0f) - 0.71414f*(Cr -128.0f));
            ImgOut[place+2] = clamp(Y + 1.772f*(Cb-128.0f));
        }
    }
}*/


/*
double autoSeuil(uchar* ImgIn, int nH, int nW){
    //int seuil;
    double max=0;
    double otsu;
    for(int i=0; i<256; i++){
        otsu = otsu(ImgIn, i, nH, nW);
        if(otsu>max){
            max = otsu;
            //seuil = i;
        }
    }
    return max;
}

void seedList(uchar* ImgIn, uchar* ImgOut, int nH, int nW){
    double similiSeuil = autoSeuil(ImgIn, nH, nW);
    int taille = nH*nW;
    int nW3 = nW*3;
    int placeY, placeCb, placeCr;
	uchar * ImgCond1 = new uchar[taille];
    double vY, vCb, vCr, moy;
    double tabV[taille];
    double tabVNorm[taille];
    for (int i=1; i<nH-1; i++){
        for(int j=1; j<nW-1; j++){
            placeY = i*nW3+j*3;
            moy=(ImgIn[placeY]
                +ImgIn[placeY+3]
                +ImgIn[placeY-3]
                +ImgIn[placeY-nW3]
                +ImgIn[placeY-nW3-3]
                +ImgIn[placeY-nW3+3]
                +ImgIn[placeY+nW3]
                +ImgIn[placeY+nW3+3]
                +ImgIn[placeY+nW3-3]
                )/9;
            vY=sqrt((pow(ImgIn[placeY]-moy,2)
                    +pow(ImgIn[placeY-3]-moy,2)
                    +pow(ImgIn[placeY+3]-moy,2)
                    +pow(ImgIn[placeY-nW3]-moy,2)
                    +pow(ImgIn[placeY-3-nW3]-moy,2)
                    +pow(ImgIn[placeY+3-nW3]-moy,2)
                    +pow(ImgIn[placeY+nW3]-moy,2)
                    +pow(ImgIn[placeY-3+nW3]-moy,2)
                    +pow(ImgIn[placeY+3+nW3]-moy,2)
                    )/9
                    );
            placeCb = placeY+1;
            moy=(ImgIn[placeCb]
                +ImgIn[placeCb+3]
                +ImgIn[placeCb-3]
                +ImgIn[placeCb-nW3]
                +ImgIn[placeCb-nW3-3]
                +ImgIn[placeCb-nW3+3]
                +ImgIn[placeCb+nW3]
                +ImgIn[placeCb+nW3+3]
                +ImgIn[placeCb+nW3-3]
                )/9;
            vCb=sqrt((pow(ImgIn[placeCb]-moy,2)
                    +pow(ImgIn[placeCb-3]-moy,2)
                    +pow(ImgIn[placeCb+3]-moy,2)
                    +pow(ImgIn[placeCb-nW3]-moy,2)
                    +pow(ImgIn[placeCb-3-nW3]-moy,2)
                    +pow(ImgIn[placeCb+3-nW3]-moy,2)
                    +pow(ImgIn[placeCb+nW3]-moy,2)
                    +pow(ImgIn[placeCb-3+nW3]-moy,2)
                    +pow(ImgIn[placeCb+3+nW3]-moy,2)
                    )/9
                    );
            placeCr = placeY+2;
            moy=(ImgIn[placeCr]
                +ImgIn[placeCr+3]
                +ImgIn[placeCr-3]
                +ImgIn[placeCr-nW3]
                +ImgIn[placeCr-nW3-3]
                +ImgIn[placeCr-nW3+3]
                +ImgIn[placeY+nW3]
                +ImgIn[placeCr+nW3+3]
                +ImgIn[placeCr+nW3-3]
                )/9;
            vCb=sqrt((pow(ImgIn[placeCr]-moy,2)
                    +pow(ImgIn[placeCr-3]-moy,2)
                    +pow(ImgIn[placeCr+3]-moy,2)
                    +pow(ImgIn[placeCr-nW3]-moy,2)
                    +pow(ImgIn[placeCr-3-nW3]-moy,2)
                    +pow(ImgIn[placeCr+3-nW3]-moy,2)
                    +pow(ImgIn[placeCr+nW3]-moy,2)
                    +pow(ImgIn[placeCr-3+nW3]-moy,2)
                    +pow(ImgIn[placeCr+3+nW3]-moy,2)
                    )/9
                    );
            tabV[i*nW+j]= vY + vCb + vCr;
        }
    }
}*/

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

std::vector<Px> getNeighbours(Px p, int h, int w)
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

std::vector<Px> getNeighbours(PxDist p, int h, int w)
{
	return getNeighbours(Px{p.row, p.col}, h, w);
}

struct MinReg
{
	int label;
	double dist;
};

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

	uchar * imageSeed = new uchar[size];
	std::memset(imageSeed, 0u, sizeof(uchar)*size);
	computeSeeds(imageYCbCr, imageSeed, h, w, 100, 0.02f);

	// create regions based solely on seeds
	int * regionIds = new int[size];
	std::memset(regionIds, -1, size*sizeof(int));
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

	Region * regions = new Region[regionCount];
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

	// compute frontiers
	std::vector<PxDist> frontier;
	for (int row=0; row<h; ++row)
	{
		for (int col = 0; col < w; ++col)
		{
			if (regionIds[row * w + col] != -1)
				continue;
			std::vector<Px> neighs = getNeighbours(Px{row, col}, h, w);
			MinReg minReg = getClosestRegion(imageYCbCr, regionIds, regions, Px{row, col}, h, w);
			if (minReg.label != -1)
				frontier.push_back({row, col, minReg.dist});
		}
	}
	std::sort(frontier.begin(), frontier.end(),
			[](const PxDist & px1, const PxDist & px2)
			{
				return px1.minDist > px2.minDist;
			}
		  );

	// growing regions
	while (!frontier.empty())
	{
		PxDist px = frontier.back();
		frontier.pop_back();
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
			float *pxVal = imageYCbCr+px.row*w*3+px.col*3;
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
					// find if the pixel ha a labelled neighbour other than the pixel that we just processed
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
					frontier.insert(
						std::upper_bound(frontier.begin(), frontier.end(), newPx,
										 [](const PxDist & px1, const PxDist & px2){return px1.minDist < px2.minDist;}),
						newPx);
				}
			}
		}
	}

	// debug
	uchar * tmpTest = new uchar[size];
	std::memset(tmpTest, 0u, sizeof(uchar)*size);
	for (const PxDist & px : frontier)
	{
		tmpTest[px.row*w+px.col] = 255u;
	}
	ecrire_image_pgm((baseName+"testfrontier.pgm").data(), tmpTest, h, w);
	delete [] tmpTest;

	int test = *std::max_element(regionIds, regionIds+size);
	ushort * imTest = new ushort[size];
	for (int i=0; i<size; ++i)
		imTest[i] = (regionIds[i]/(float)test)*65535;
	delete [] imTest;
	ecrire_image_pgm_2o((baseName+"testregions.pgm").data(), imTest, h, w);

	delete [] regions;
	delete [] regionIds;
	delete [] imageYCbCr;
	delete [] imageSeed;
	delete [] image;

    return 0;
}



